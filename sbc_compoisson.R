# --- Load Libraries ---
library(rstan)
library(brms)

# --- Configuration ---
# set.seed(789) # for reproducibility

# SBC parameters
N_simulations <- 5000 # Number of SBC iterations (datasets/fits)
J = 10000 # Number of samples to step 2

# --- Auxiliary Functions ---
rcomp <- function(n, lambda, nu, max_y = 10000, tol = 1e-16) {
  # Support
  y <- 0:max_y
  
  # Compute log-probabilities for stability
  log_pmf <- y * log(lambda) - nu * lfactorial(y)
  log_pmf <- log_pmf - max(log_pmf)  # prevent overflow
  pmf <- exp(log_pmf)
  pmf <- pmf / sum(pmf)              # normalize
  
  # Check if PMF sums to 1 (otherwise increase max_y)
  cum_pmf <- cumsum(pmf)
  if (tail(cum_pmf, 1) < 1 - tol) {
    warning("PMF may be truncated; increase max_y.")
  }
  
  # Inverse transform sampling
  u <- runif(n)
  samples <- sapply(u, function(ui) which(cum_pmf >= ui)[1] - 1)
  return(samples)
}

# --- 1. Draw the Joint Prior ---
lambda_sim <- rlnorm(1, meanlog = 0, sdlog = 1)
nu_sim <- rlnorm(1, meanlog = 0, sdlog = 1)

#Z_sim <- # sla

# --- 2. Draw a dataset y ~ COMPoisson(lambda_sim, nu_sim) ---
dataset <- rcomp(J, lambda_sim, nu_sim)

df <- as.data.frame(table(dataset))
colnames(df) <- c("y", "x")

df$y <- as.numeric(df$y)
df$x <- as.numeric(df$x)


# --- 3. MCMC in Stan do obtain samples ---
fit <- brm(y ~ x,
           data = df,
           prior = c(
             set_prior("lognormal(0, 1)", class = "mu"),
             set_prior("lognormal(0, 1)", class = "nu")
           ),
           cores = 4,
           iter = 10000,
           refresh = 2000,
           backend = "cmdstanr",
           family = "com_poisson")



