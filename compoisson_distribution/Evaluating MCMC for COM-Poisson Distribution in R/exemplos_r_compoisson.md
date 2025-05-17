# Exemplos em R para Cálculo da Constante de Normalização da COM-Poisson

Este documento contém exemplos práticos em R para calcular a constante de normalização da distribuição Conway-Maxwell-Poisson (COM-Poisson), utilizando principalmente o pacote `COMPoissonReg`.

## 1. Instalação e Carregamento do Pacote

```r
# Instalação do pacote (se necessário)
if (!require("COMPoissonReg")) {
  install.packages("COMPoissonReg")
}

# Carregamento do pacote
library(COMPoissonReg)
```

## 2. Configuração de Controle para Cálculos

O pacote `COMPoissonReg` permite configurar parâmetros de controle para o cálculo da constante de normalização:

```r
# Criação de um objeto de controle personalizado
control <- get.control(
  ymax = 100000,        # Limite superior para truncamento
  hybrid.tol = 1e-2,    # Tolerância para método híbrido
  truncate.tol = 1e-6   # Tolerância para truncamento
)

# Verificação dos valores padrão
control$ymax
# [1] 1e+06
control$hybrid.tol
# [1] 0.01
control$truncate.tol
# [1] 1e-06

# Alteração dos valores padrão na sessão atual
options(COMPoissonReg.control = control)
```

## 3. Cálculo da Constante de Normalização

A função principal para calcular a constante de normalização é `ncmp()`:

```r
# Cálculo básico da constante
ncmp(lambda = 1.5, nu = 1.2)
# [1] 4.01341

# Cálculo na escala logarítmica (recomendado para estabilidade numérica)
ncmp(lambda = 1.5, nu = 1.2, log = TRUE)
# [1] 1.389641

# Especificando controle personalizado
ncmp(lambda = 1.5, nu = 1.2, log = TRUE, control = get.control(hybrid.tol = 1e-10))
# [1] 1.373642

# Usando controle com tolerância diferente
ncmp(lambda = 1.5, nu = 1.2, log = TRUE, control = get.control(hybrid.tol = 1e-10))
# [1] 1.389641
```

## 4. Verificação do Método de Truncamento

A função `tcmp()` retorna o valor de truncamento M obtido do algoritmo:

```r
# Definição de uma função para exibir avisos
print_warning <- function(x) { print(strwrap(x), quote = FALSE) }

# Exemplo com diferentes valores de parâmetros
nu_seq <- c(1, 0.5, 0.2, 0.1, 0.05, 0.03)
tryCatch({ tcmp(lambda = 1.5, nu = nu_seq) }, warning = print_warning)
# [1] simpleWarning in y_trunc(prep$lambda, prep$nu, tol = truncate.tol, ymax = ymax)
```

## 5. Implementação Manual da Constante de Normalização

Para fins educacionais ou quando maior controle é necessário, pode-se implementar o cálculo manualmente:

```r
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

# Exemplo de uso
lambda <- 1.5
nu <- 1.2
Z_trunc <- calc_Z_truncated(lambda, nu, 50)
Z_asymp <- calc_Z_asymptotic(lambda, nu)
Z_hybrid <- calc_Z_hybrid(lambda, nu)

# Comparação com o resultado do pacote
Z_package <- exp(ncmp(lambda, nu, log = TRUE))
```

## 6. Uso em Contexto MCMC

Em análises MCMC, a constante de normalização é frequentemente necessária para calcular a verossimilhança:

```r
# Função de log-verossimilhança para COM-Poisson
com_poisson_log_likelihood <- function(y, lambda, nu) {
  # Cálculo do log da constante de normalização
  log_Z <- ncmp(lambda, nu, log = TRUE)
  
  # Cálculo da log-verossimilhança para cada observação
  log_lik <- y * log(lambda) - nu * log(factorial(y)) - log_Z
  
  # Retorna a soma das log-verossimilhanças
  return(sum(log_lik))
}

# Exemplo com dados simulados
set.seed(123)
y <- rpois(100, lambda = 2)  # Dados simulados de uma Poisson
lambda_test <- 1.8
nu_test <- 1.1
log_lik <- com_poisson_log_likelihood(y, lambda_test, nu_test)
```

## 7. Implementação Eficiente para MCMC

Para MCMC, onde a constante precisa ser calculada repetidamente, pode-se usar uma abordagem de cache:

```r
# Função que implementa cache para constantes já calculadas
cached_ncmp <- function() {
  # Inicializa cache vazio
  cache <- new.env(parent = emptyenv())
  
  # Retorna função com cache
  function(lambda, nu, log = TRUE) {
    # Cria chave para o cache
    key <- paste(lambda, nu, sep = "_")
    
    # Verifica se o valor já está no cache
    if (!exists(key, envir = cache)) {
      # Calcula e armazena no cache
      cache[[key]] <- ncmp(lambda, nu, log = log)
    }
    
    # Retorna valor do cache
    return(cache[[key]])
  }
}

# Cria função com cache
ncmp_cached <- cached_ncmp()

# Uso em MCMC (exemplo simplificado)
mcmc_com_poisson <- function(y, iterations = 1000) {
  # Valores iniciais
  lambda <- mean(y)
  nu <- 1
  
  # Armazena resultados
  results <- matrix(0, nrow = iterations, ncol = 2)
  colnames(results) <- c("lambda", "nu")
  
  # Loop MCMC simplificado
  for (i in 1:iterations) {
    # Propõe novos valores (simplificado)
    lambda_new <- rnorm(1, lambda, 0.1)
    nu_new <- rnorm(1, nu, 0.1)
    
    # Restringe a valores positivos
    if (lambda_new <= 0 || nu_new <= 0) next
    
    # Calcula log-verossimilhança atual e proposta
    log_lik_current <- sum(y * log(lambda) - nu * log(factorial(y)) - ncmp_cached(lambda, nu))
    log_lik_proposed <- sum(y * log(lambda_new) - nu_new * log(factorial(y)) - ncmp_cached(lambda_new, nu_new))
    
    # Razão de aceitação (simplificada, sem prior)
    log_ratio <- log_lik_proposed - log_lik_current
    
    # Aceita ou rejeita
    if (log(runif(1)) < log_ratio) {
      lambda <- lambda_new
      nu <- nu_new
    }
    
    # Armazena resultados
    results[i, ] <- c(lambda, nu)
  }
  
  return(results)
}
```

Estes exemplos demonstram diferentes abordagens para calcular a constante de normalização da distribuição COM-Poisson em R, desde o uso direto do pacote `COMPoissonReg` até implementações manuais e otimizações para MCMC.
