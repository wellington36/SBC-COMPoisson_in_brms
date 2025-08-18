library(parallel)
library(brms)
source("rcomp.R")

# Delete old models
files <- list.files("models", full.names = TRUE)
unlink(files)

# --- SBC parameters ---
leps <- -8 * log(2)
N_simulations <- 10   # increase later
J <- 1000
stan_chains <- 1
stan_iter <- 10000     # reduce while testing
stan_warmup <- 8000
M_posterior_draws <- stan_chains * (stan_iter - stan_warmup)

# helper: one SBC run
sbc_one <- function(i) {
  t0 <- Sys.time()
  
  # --- 1. Prior draw ---
  lambda_sim <- rnorm(1)
  nu_sim <- rlnorm(1)
  while (exp(lambda_sim)^(1/nu_sim) >= 100) {
    lambda_sim <- rnorm(1)
    nu_sim <- rlnorm(1)
  }
  
  # --- 2. Simulated data ---
  dataset <- rcomp(J, exp(lambda_sim), nu_sim, sumTo = 10000, eps = 1e-6)
  df <- data.frame(y = dataset)
  
  # --- 3. Stan fit ---
  scode_string <- sprintf("real leps_custom() { return %f; }", leps)
  custom_stanvars <- stanvar(scode = scode_string, block = "functions")
  
  fit <- brm(
    y ~ 1,
    data = df,
    chains = stan_chains,
    iter = stan_iter,
    warmup = stan_warmup,
    prior = prior(normal(0, 1), class = "Intercept") +
      prior(lognormal(0, 1), class = "shape"),
    stanvars = custom_stanvars,
    backend = "cmdstanr",
    family = "com_poisson",
    file = paste0("models/stan_model_", i, ".rds"),
    file_refit = "always"
  )
  
  # --- 4. Extract draws ---
  draws <- as_draws_df(fit)
  lambda_draws <- exp(draws$b_Intercept)
  nu_draws <- draws$shape
  
  # --- 5. Diagnostics with posterior package ---
  rhat_intercept <- coef_tab["Intercept", "Rhat"]
  ess_bulk_intercept <- coef_tab["Intercept", "Bulk_ESS"]
  
  ess_vals <- posterior::ess_basic(fit)
  print(ess_vals)
  ess_mean <- mean(ess_vals, na.rm = TRUE)
  
  rhats <- posterior::rhat(fit)
  print(summary(fit))
  n_converged <- sum(rhats <= 1.01, na.rm = TRUE)
  
  t1 <- Sys.time()
  elapsed_sec <- as.numeric(difftime(t1, t0, units = "secs"))
  ess_per_sec <- ess_mean / elapsed_sec
  
  list(
    rank_lambda = sum(lambda_draws < exp(lambda_sim)),
    rank_nu     = sum(nu_draws < nu_sim),
    lambda_mean = mean(lambda_draws),
    nu_mean     = mean(nu_draws),
    lambda_true = exp(lambda_sim),
    nu_true     = nu_sim,
    ess_mean    = ess_mean,
    time_sec    = elapsed_sec,
    ess_per_sec = ess_per_sec,
    n_converged = n_converged
  )
}


# --- run in parallel ---
results_list <- mclapply(1:N_simulations, 
                         sbc_one, 
                         mc.cores = detectCores(),
                         mc.preschedule = FALSE)
results <- do.call(rbind, lapply(results_list, as.data.frame))

# --- summary diagnostics ---
cat("\n=== Overall diagnostics ===\n")
cat("Mean ESS:", mean(results$ess_mean), "\n")
cat("Mean time (s):", mean(results$time_sec), "\n")
cat("Mean ESS/sec:", mean(results$ess_per_sec), "\n")
cat("MCMCs converged (Rhat <= 1.01):", sum(results$n_converged > 0), "of", N_simulations, "\n")

# --- SBC plots ---
valid_idx <- !is.na(results$rank_lambda) & !is.na(results$rank_nu)
N_valid <- sum(valid_idx)

if (N_valid > 0) {
  cat("Plotting histograms for", N_valid, "valid simulations\n")
  
  breaks <- seq(-0.5, M_posterior_draws + 0.5, by = 1)
  expected_count <- N_valid / (length(breaks) - 1)
  
  par(mfrow = c(1, 2))
  
  hist(results$rank_lambda[valid_idx],
       breaks = breaks,
       main = "SBC Rank Distribution (lambda)",
       xlab = paste("Rank (0 to", M_posterior_draws, ")"),
       ylab = "Frequency",
       col = "black")
  abline(h = expected_count, col = "red", lwd = 2, lty = 2)
  
  hist(results$rank_nu[valid_idx],
       breaks = breaks,
       main = "SBC Rank Distribution (nu)",
       xlab = paste("Rank (0 to", M_posterior_draws, ")"),
       ylab = "Frequency",
       col = "black")
  abline(h = expected_count, col = "red", lwd = 2, lty = 2)
  
  par(mfrow = c(1, 1))
}
