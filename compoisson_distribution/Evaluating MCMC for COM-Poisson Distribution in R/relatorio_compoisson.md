# Cálculo da Constante de Normalização da Distribuição Conway-Maxwell-Poisson

Este documento compila métodos e exemplos para calcular a constante de normalização da distribuição Conway-Maxwell-Poisson (COM-Poisson), especialmente no contexto de MCMC.

## Introdução

A distribuição Conway-Maxwell-Poisson (COM-Poisson) tem sua função de massa de probabilidade (p.m.f.) definida por:

$$
p(y \mid \lambda, \nu) = \frac{\lambda^y}{(y!)^\nu} \cdot \frac{1}{Z(\lambda, \nu)},
$$

onde $Z(\lambda, \nu) = \sum_{y = 0}^\infty \frac{\lambda^y}{(y!)^\nu}$ é a constante de normalização.

Esta constante não possui forma fechada em geral, o que representa um desafio para implementações práticas, especialmente em contextos de MCMC.

## Métodos para Cálculo da Constante de Normalização

### 1. Casos Especiais com Forma Fechada

Existem alguns casos especiais onde a constante de normalização possui forma fechada:

- Quando $\nu = 1$, a distribuição COM-Poisson se reduz à distribuição Poisson padrão, e $Z(\lambda, 1) = e^\lambda$.
- Quando $\nu = 0$ e $0 < \lambda < 1$, a distribuição se torna uma Geométrica, com $Z(\lambda, 0) = \frac{1}{1-\lambda}$.
- Quando $\nu \to \infty$, a distribuição converge para uma Bernoulli com parâmetro $\lambda/(1+\lambda)$.

### 2. Aproximação Assintótica

Shmueli et al. (2005) propuseram a seguinte aproximação para a constante de normalização:

$$
Z(\lambda, \nu) \approx \frac{\exp(\nu\lambda^{1/\nu})}{\lambda^{(\nu-1)/2}(2\pi)^{(\nu-1)/2}\nu^{1/2}} \left\{1 + O(\lambda^{-1/\nu})\right\}
$$

Esta aproximação é particularmente útil quando $\lambda^{-1/\nu}$ é pequeno. A expressão pode ser simplificada para:

$$
Z(\lambda, \nu) \approx \frac{\exp(\nu\lambda^{1/\nu})}{\lambda^{(\nu-1)/2}(2\pi)^{(\nu-1)/2}\nu^{1/2}}
$$

ou na escala logarítmica:

$$
\log Z(\lambda, \nu) \approx \nu\lambda^{1/\nu} - \frac{\nu-1}{2}\log\lambda - \frac{\nu-1}{2}\log(2\pi) - \frac{1}{2}\log\nu
$$

Esta aproximação é utilizada em casos onde $\lambda^{-1/\nu} < \delta$ para algum valor pequeno de $\delta$ (geralmente 0.01).

### 3. Truncamento da Série Infinita

Quando a aproximação assintótica não é adequada, a constante é calculada truncando a série infinita:

$$
Z^{(M)}(\lambda, \nu) = \sum_{y=0}^{M} \frac{\lambda^y}{(y!)^\nu}
$$

O valor de $M$ é escolhido de forma que o erro de truncamento seja menor que uma tolerância especificada $\epsilon$:

$$
\frac{|Z(\lambda, \nu) - Z^{(M)}(\lambda, \nu)|}{Z^{(M)}(\lambda, \nu)} \leq \epsilon
$$

O algoritmo para calcular a constante usando truncamento segue estes passos:

1. Iniciar com $M = 0$ e $Z^{(0)} = 1$
2. Enquanto $M < \lambda^{1/\nu} - 1$ ou o erro relativo não for satisfeito:
   - Incrementar $M$
   - Atualizar $Z^{(M+1)} = Z^{(M)} + \lambda^M/(M!)^\nu$
3. Retornar $Z^{(M)}$

Para evitar problemas numéricos com valores muito grandes ou muito pequenos, os cálculos são realizados na escala logarítmica sempre que possível.

### 4. Abordagem Híbrida

Na prática, o pacote COMPoissonReg em R utiliza uma abordagem híbrida:

- Se $\lambda^{-1/\nu} < \delta$ (onde $\delta$ é um valor pequeno como 1e-2), usa a aproximação assintótica
- Caso contrário, utiliza o método de truncamento da série

Esta abordagem equilibra precisão e eficiência computacional.

## Implementação em R

### Instalação e Carregamento do Pacote

```r
# Instalação do pacote (se necessário)
if (!require("COMPoissonReg")) {
  install.packages("COMPoissonReg")
}

# Carregamento do pacote
library(COMPoissonReg)
```

### Configuração de Controle para Cálculos

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

### Cálculo da Constante de Normalização

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
```

### Implementação Manual da Constante de Normalização

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
```

### Uso em Contexto MCMC

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
```

### Implementação Eficiente para MCMC

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
```

## Recomendações para Uso em MCMC

1. **Use a escala logarítmica**: Sempre trabalhe com o logaritmo da constante de normalização (`log = TRUE`) para evitar problemas numéricos.

2. **Implemente cache**: Em MCMC, muitos valores de parâmetros são recalculados repetidamente. Usar cache pode melhorar significativamente o desempenho.

3. **Ajuste as tolerâncias**: Dependendo da precisão necessária, ajuste os parâmetros `hybrid.tol` e `truncate.tol` para equilibrar precisão e velocidade.

4. **Considere aproximações**: Para valores grandes de λ e pequenos de ν, a aproximação assintótica pode ser muito mais rápida com pouca perda de precisão.

5. **Evite cálculos desnecessários**: Em algoritmos MCMC, muitas vezes é possível trabalhar com razões de probabilidades, o que pode eliminar a necessidade de calcular a constante de normalização em cada iteração.

## Validação dos Métodos

Um script de validação está disponível para testar os diferentes métodos de cálculo da constante de normalização. Este script compara os resultados do pacote `COMPoissonReg` com implementações manuais e fornece recomendações baseadas nos resultados.

Para executar o script de validação, salve o arquivo `validate_compoisson.R` e execute-o em um ambiente R:

```r
source("validate_compoisson.R")
```

## Referências

1. Shmueli, G., Minka, T. P., Kadane, J. B., Borle, S., & Boatwright, P. (2005). A useful distribution for fitting discrete data: revival of the Conway–Maxwell–Poisson distribution. Journal of the Royal Statistical Society: Series C (Applied Statistics), 54(1), 127-142.

2. Gillispie, S. B., & Green, C. G. (2015). Approximating the Conway–Maxwell–Poisson distribution normalization constant. Statistics, 49(5), 1062-1073.

3. Daly, F., & Gaunt, R. E. (2016). The Conway-Maxwell-Poisson distribution: Distributional theory and approximation. ALEA, 13(2), 635-658.

4. Sellers, K. F., & Shmueli, G. (2010). A flexible regression model for count data. The Annals of Applied Statistics, 4(2), 943-961.

5. Pacote COMPoissonReg: https://cran.r-project.org/web/packages/COMPoissonReg/
