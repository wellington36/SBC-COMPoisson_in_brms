library(parallel)
library(brms)
source("rcomp.R")

zinb <- read.csv("https://paul-buerkner.github.io/data/fish.csv")

leps <- 1 * log(2)
N_simulations <- 1   # increase later
stan_chains <- 1
stan_iter <- 1000     # reduce while testing
stan_warmup <- 800

one_test <- function(i) {
  t0 <- Sys.time()
  
  # --- Stan fit ---
  scode_string <- sprintf("real leps_custom() { return %f; }", leps)
  custom_stanvars <- stanvar(scode = scode_string, block = "functions")
  
  fit <- brm(
    count ~ 1,
    data = zinb,
    chains = stan_chains,
    iter = stan_iter,
    warmup = stan_warmup,
    prior = prior(normal(0, 1), class = "Intercept") +
      prior(lognormal(0, 1), class = "shape"),
    stanvars = custom_stanvars,
    backend = "cmdstanr",
    family = "com_poisson",
    control = list()
  )
  
  t1 <- Sys.time()
  
  # --- Extract draws ---
  draws <- as_draws_df(fit)
  lambda_draws <- exp(draws$b_Intercept)
  nu_draws <- draws$shape
  
  # --- Diagnostics with posterior package ---
  coef_tab <- summary(fit)$fixed
  rhat_intercept <- coef_tab["Intercept", "Rhat"]
  estimate_intercept <- coef_tab["Intercept", "Estimate"]
  ess_bulk_intercept <- coef_tab["Intercept", "Bulk_ESS"]
  
  elapsed_sec <- as.numeric(difftime(t1, t0, units = "secs"))
  
  list(
    lambda_mean = log(mean(lambda_draws)),
    nu_mean     = mean(nu_draws),
    rhat        = rhat_intercept,
    ess         = ess_bulk_intercept,
    time_sec    = elapsed_sec,
    ess_per_sec = ess_bulk_intercept / elapsed_sec,
    converged   = as.numeric(ifelse(
      (rhat_intercept <= 1.05) &
        (!(estimate_intercept < -0.26 | estimate_intercept > -0.25)), 
      1, 0
    ))
  )
}

# --- Run in parallel ---
results_list <- mclapply(1:N_simulations,
                        one_test,
                        mc.cores = 8)
results <- do.call(rbind, lapply(results_list, as.data.frame))

cat("\n=== Overall diagnostics ===\n")
cat("Mean ESS:", mean(results$ess), "\n")
cat("Mean time (s):", mean(results$time_sec), "\n")
cat("Mean ESS/sec:", mean(results$ess_per_sec), "\n")
cat("MCMCs converged:", sum(results$converged), "of", N_simulations, "\n")








