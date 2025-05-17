# Métodos para Cálculo da Constante de Normalização da Distribuição Conway-Maxwell-Poisson

A distribuição Conway-Maxwell-Poisson (COM-Poisson) tem sua função de massa de probabilidade (p.m.f.) definida por:

$$
p(y \mid \lambda, \nu) = \frac{\lambda^y}{(y!)^\nu} \cdot \frac{1}{Z(\lambda, \nu)},
$$

onde $Z(\lambda, \nu) = \sum_{y = 0}^\infty \frac{\lambda^y}{(y!)^\nu}$ é a constante de normalização.

Esta constante não possui forma fechada em geral, o que representa um desafio para implementações práticas, especialmente em contextos de MCMC. Abaixo estão os principais métodos para calcular esta constante.

## 1. Casos Especiais com Forma Fechada

Existem alguns casos especiais onde a constante de normalização possui forma fechada:

- Quando $\nu = 1$, a distribuição COM-Poisson se reduz à distribuição Poisson padrão, e $Z(\lambda, 1) = e^\lambda$.
- Quando $\nu = 0$ e $0 < \lambda < 1$, a distribuição se torna uma Geométrica, com $Z(\lambda, 0) = \frac{1}{1-\lambda}$.
- Quando $\nu \to \infty$, a distribuição converge para uma Bernoulli com parâmetro $\lambda/(1+\lambda)$.

## 2. Aproximação Assintótica

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

Esta aproximação é utilizada em casos onde $\lambda^{-1/\nu} < \delta$ para algum valor pequeno de $\delta$.

## 3. Truncamento da Série Infinita

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

## 4. Abordagem Híbrida

Na prática, o pacote COMPoissonReg em R utiliza uma abordagem híbrida:

- Se $\lambda^{-1/\nu} < \delta$ (onde $\delta$ é um valor pequeno como 1e-2), usa a aproximação assintótica
- Caso contrário, utiliza o método de truncamento da série

Esta abordagem equilibra precisão e eficiência computacional.
