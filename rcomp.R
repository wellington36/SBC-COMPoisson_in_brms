log_k_term <- function(log_lambda, nu, k) {
  return((k - 1) * log_lambda - nu * lgamma(k))
}

bound_remainder <- function(k_current_term, k_previous_term) {
  return(k_current_term - log(-expm1(k_current_term - k_previous_term)))
}

stopping_criterio_bucket <- function(k_current_term, k_previous_term, k, leps) {
  if (k %% 50 == 0) {
    return(bound_remainder(k_current_term, k_previous_term) >= leps)
  }
  return(TRUE)
}

log_Z_com_poisson <- function(log_lambda, nu, eps) {
  if (nu == 1) {
    return(log_lambda)
  }
  if (nu <= 0) {
    stop("nu must be positive")
  }
  if (is.infinite(nu)) {
    stop("nu must be finite")
  }

  M <- 1000000
  leps <- log(eps)

  # direct computation of the truncated series
  # check if the Mth term of the series pass in the stopping criteria
  if ((log_k_term(log_lambda, nu, M) > log_k_term(log_lambda, nu, M - 1)) ||
      (bound_remainder(log_k_term(log_lambda, nu, M),
                      log_k_term(log_lambda, nu, M - 1)) >= leps)) {
    stop("nu is too close to zero.")
  }

  log_Z_terms <- numeric(M)
  log_Z_terms[1] <- log_k_term(log_lambda, nu, 1)
  log_Z_terms[2] <- log_k_term(log_lambda, nu, 2)

  k <- 2
  while (((log_Z_terms[k] >= log_Z_terms[k-1]) ||
          (stopping_criterio_bucket(log_Z_terms[k], log_Z_terms[k-1], k, leps))) &&
         k < M) {
    k <- k + 1
    log_Z_terms[k] <- log_k_term(log_lambda, nu, k)
  }
  log_Z <- logsumexp(log_Z_terms[1:k])

  return(log_Z)
}

com_poisson_log_lpmf <- function(y, log_lambda, nu, eps) {
  if (nu == 1) {
    return(dpois(y, exp(log_lambda), log = TRUE))
  }
  return(y * log_lambda - nu * lgamma(y + 1) - log_Z_com_poisson(log_lambda, nu, eps))
}

# Helper function for logsumexp, as it's not native in base R
logsumexp <- function(x) {
  xmax <- max(x)
  xmax + log(sum(exp(x - xmax)))
}

rcomp <- function(n, lam, nu, sumTo = 10000, eps = 10^-6) {
  if (lam <= 0) {
    stop("lambda must be positive")
  }
  if (nu <= 0) {
    stop("nu must be positive")
  }

  log_lambda <- log(lam)

  # Calculate log-PMF for y from 0 to sumTo
  y_values <- 0:sumTo
  log_pmf_values <- sapply(y_values, function(y) com_poisson_log_lpmf(y, log_lambda, nu, eps))

  # Convert log-PMF to PMF (probabilities)
  pmf_values <- exp(log_pmf_values)

  # Normalize probabilities
  pmf_values <- pmf_values / sum(pmf_values)

  # Sample from the probabilities
  samples <- sample(y_values, size = n, replace = TRUE, prob = pmf_values)

  return(samples)
}

# rcomp <- function(n, lambda, nu, max_y = 100000, tol = 1e-10) {
#   # Step 1: Support
#   y <- 0:max_y
# 
#   # Step 2: Compute log-probabilities for stability
#   log_pmf <- y * log(lambda) - nu * lfactorial(y)
#   log_pmf <- log_pmf - max(log_pmf)  # prevent overflow
#   pmf <- exp(log_pmf)
#   pmf <- pmf / sum(pmf)              # normalize
# 
#   # Step 3: Check if PMF sums to 1 (otherwise increase max_y)
#   cum_pmf <- cumsum(pmf)
#   if (tail(cum_pmf, 1) < 1 - tol) {
#     rcomp(n, lambda, nu, 10*max_y, tol)
#   }
# 
#   # Step 4: Inverse transform sampling
#   u <- runif(n)
#   samples <- sapply(u, function(ui) which(cum_pmf >= ui)[1] - 1)
#   return(samples)
# }