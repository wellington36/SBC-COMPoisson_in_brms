# Script para validar métodos de cálculo da constante de normalização da COM-Poisson

# Instalação e carregamento do pacote
if (!require("COMPoissonReg")) {
  install.packages("COMPoissonReg", repos = "https://cloud.r-project.org")
}
library(COMPoissonReg)

# Função para calcular Z(lambda, nu) usando truncamento
calc_Z_truncated <- function(lambda, nu, M = 100) {
  sum_terms <- numeric(M + 1)
  for (y in 0:M) {
    sum_terms[y + 1] <- lambda^y / (factorial(y))^nu
  }
  return(sum(sum_terms))
}

# Função para calcular Z(lambda, nu) usando aproximação assintótica
calc_Z_asymptotic <- function(lambda, nu) {
  term1 <- exp(nu * lambda^(1/nu))
  term2 <- lambda^((nu-1)/2) * (2*pi)^((nu-1)/2) * nu^(1/2)
  return(term1 / term2)
}

# Função híbrida que escolhe o método apropriado
calc_Z_hybrid <- function(lambda, nu, delta = 0.01, M = 100) {
  if (lambda^(-1/nu) < delta) {
    return(calc_Z_asymptotic(lambda, nu))
  } else {
    return(calc_Z_truncated(lambda, nu, M))
  }
}

# Validação com diferentes valores de parâmetros
test_cases <- data.frame(
  lambda = c(0.5, 1.0, 1.5, 2.0, 5.0, 10.0),
  nu = c(0.5, 1.0, 1.2, 1.5, 2.0, 3.0)
)

results <- data.frame(
  lambda = numeric(0),
  nu = numeric(0),
  Z_package = numeric(0),
  Z_trunc = numeric(0),
  Z_asymp = numeric(0),
  Z_hybrid = numeric(0)
)

for (i in 1:nrow(test_cases)) {
  lambda <- test_cases$lambda[i]
  nu <- test_cases$nu[i]
  
  # Cálculo usando o pacote
  Z_package <- exp(ncmp(lambda, nu, log = TRUE))
  
  # Cálculo usando métodos implementados
  Z_trunc <- calc_Z_truncated(lambda, nu, M = 100)
  Z_asymp <- calc_Z_asymptotic(lambda, nu)
  Z_hybrid <- calc_Z_hybrid(lambda, nu)
  
  # Armazenar resultados
  results <- rbind(results, data.frame(
    lambda = lambda,
    nu = nu,
    Z_package = Z_package,
    Z_trunc = Z_trunc,
    Z_asymp = Z_asymp,
    Z_hybrid = Z_hybrid
  ))
}

# Calcular erro relativo
results$err_trunc <- abs(results$Z_trunc - results$Z_package) / results$Z_package
results$err_asymp <- abs(results$Z_asymp - results$Z_package) / results$Z_package
results$err_hybrid <- abs(results$Z_hybrid - results$Z_package) / results$Z_package

# Imprimir resultados
print(results)

# Verificar quando a aproximação assintótica é adequada
results$lambda_inv_nu <- results$lambda^(-1/results$nu)
results$asymp_adequate <- results$lambda_inv_nu < 0.01

# Imprimir recomendações
cat("\nRecomendações baseadas nos resultados:\n")
for (i in 1:nrow(results)) {
  cat(sprintf("Para lambda=%.1f, nu=%.1f:\n", results$lambda[i], results$nu[i]))
  
  if (results$asymp_adequate[i]) {
    cat("  - Aproximação assintótica é adequada (lambda^(-1/nu) < 0.01)\n")
  } else {
    cat("  - Truncamento da série é recomendado (lambda^(-1/nu) >= 0.01)\n")
  }
  
  if (results$err_trunc[i] < 1e-6) {
    cat("  - Truncamento tem alta precisão (erro < 1e-6)\n")
  }
  
  if (results$err_asymp[i] < 1e-2) {
    cat("  - Aproximação assintótica tem boa precisão (erro < 1e-2)\n")
  }
  
  cat("\n")
}

# Exemplo de uso em contexto MCMC
cat("Exemplo de uso em contexto MCMC:\n")

# Função de log-verossimilhança para COM-Poisson
com_poisson_log_likelihood <- function(y, lambda, nu) {
  # Cálculo do log da constante de normalização
  log_Z <- ncmp(lambda, nu, log = TRUE)
  
  # Cálculo da log-verossimilhança para cada observação
  log_lik <- y * log(lambda) - nu * log(factorial(y)) - log_Z
  
  # Retorna a soma das log-verossimilhanças
  return(sum(log_lik))
}

# Dados simulados
set.seed(123)
y <- rpois(100, lambda = 2)

# Calcular log-verossimilhança para diferentes parâmetros
log_liks <- data.frame(
  lambda = seq(1.5, 2.5, by = 0.1),
  nu = rep(1.0, 11)
)

for (i in 1:nrow(log_liks)) {
  log_liks$log_lik[i] <- com_poisson_log_likelihood(y, log_liks$lambda[i], log_liks$nu[i])
}

# Imprimir resultados
print(log_liks)
cat("Valor máximo da log-verossimilhança:", max(log_liks$log_lik), "\n")
cat("Parâmetros correspondentes: lambda =", log_liks$lambda[which.max(log_liks$log_lik)], 
    ", nu =", log_liks$nu[which.max(log_liks$log_lik)], "\n")

# Salvar resultados
write.csv(results, "validation_results.csv", row.names = FALSE)
write.csv(log_liks, "mcmc_example_results.csv", row.names = FALSE)

cat("\nValidação concluída. Resultados salvos em validation_results.csv e mcmc_example_results.csv\n")
