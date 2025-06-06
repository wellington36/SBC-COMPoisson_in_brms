# --- Load Libraries ---
library(rstan)
library(brms)

# --- Configuration ---
# set.seed(789) # for reproducibility

# --- Remove Salved Compiled Stan Model if exists ---
if (file.exists("stan_model.rds")) {
  #Delete file if it exists
  file.remove("stan_model.rds")
}

# --- To have the posterior in brms ---
stanvars <- stanvar(scode = "
  vector[2] lt_sim;
  lt_sim[1] = Intercept;
  lt_sim[2] = shape;
", block = "genquant")

# SBC parameters
N_simulations <- 10 # Number of SBC iterations (datasets/fits)
J = 200 # Number of samples to step 2
stan_chains <- 1
stan_iter <- 10000 # Iterations per chain
stan_warmup <- 9000 # Warmup iterations per chain
M_posterior_draws <- stan_chains * (stan_iter - stan_warmup)

sbc_ranks <- matrix(NA, nrow = N_simulations, ncol = 2, dimnames = list(NULL, c("lambda", "nu")))

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

# --- SBC for COMPoisson ---
for (i in 1:N_simulations) {
  # --- 1. Draw the Joint Prior ---
  lambda_sim <- rlnorm(1, meanlog = 0, sdlog = 1)
  nu_sim <- rlnorm(1, meanlog = 0, sdlog = 1)

  #Z_sim <- # sla

  # --- 2. Draw a dataset y ~ COMPoisson(lambda_sim, nu_sim) ---
  dataset <- rcomp(J, lambda_sim, nu_sim)

  df <- as.data.frame(table(dataset))
  colnames(df) <- c("x", "y")

  df$x <- as.numeric(df$x)
  df$y <- as.numeric(df$y)


  # --- 3. MCMC in Stan do obtain samples ---
  fit <- brm(y ~ 1,
             data = df,
             chains = stan_chains,
             iter = stan_iter,
             cores = 8,
             warmup = stan_warmup,
             file = "stan_model.rds",     # saves compiled Stan
             file_refit = "always",
             #save_model = "saved_model",    # optional: save the .stan file
             #recompile = FALSE,             # critical: don't recompile unless needed
             backend = "cmdstanr",
             family = "com_poisson",
             control = list(adapt_delta = 0.99),
             stanvars = stanvars
             )

  # -- 4. Obtain ranks
  draws <- as_draws_df(fit)
  
  # Check for lt_sim columns
  lt_cols <- grep("^lt_sim", names(draws), value = TRUE)
  
  if (length(lt_cols) == 0) {
    stop("No lt_sim variables found. Did you define them in 'brms' with `bf(..., nl = TRUE)` or via custom Stan code?")
  }
  
  lt_sim_draws <- draws[, lt_cols]
  
  lambda_draws <- lt_sim_draws[[1]]
  nu_draws     <- lt_sim_draws[[2]]
  
  sbc_ranks[i, "lambda"] <- sum(lambda_draws < log(lambda_sim))
  sbc_ranks[i, "nu"]     <- sum(nu_draws < nu_sim)
  
  if (i %% 1 == 0 || i == 1) {
    print(paste("Simulation", i, "of", N_simulations, "completed."))
  }
}

print("SBC simulations finished.")

# --- Check Uniformity ---

# Remove any rows with NA ranks (e.g., if extraction failed)
valid_ranks <- !is.na(sbc_ranks[, 1]) & !is.na(sbc_ranks[, 2])
sbc_ranks_valid <- sbc_ranks[valid_ranks, , drop = FALSE]
N_valid_simulations <- nrow(sbc_ranks_valid)

if (N_valid_simulations < N_simulations) {
  print(paste("Warning:", N_simulations - N_valid_simulations, "simulations had issues and were excluded."))
}

if (N_valid_simulations > 0) {
  print(paste("Plotting histograms for", N_valid_simulations, "valid simulations."))
  
  # Define breaks for the histogram (M+1 bins from -0.5 to M+0.5)
  # Ranks range from 0 to M_posterior_draws
  breaks <- seq(-0.5, M_posterior_draws + 0.5, by = 1)
  num_bins <- length(breaks) - 1
  
  # Calculate expected count per bin for uniform distribution
  expected_count <- N_valid_simulations / num_bins
  
  # Set up plot area for two histograms
  par(mfrow = c(1, 2)) # Arrange plots in 1 row, 2 columns
  
  # Plot histogram for lambda ranks
  hist(sbc_ranks_valid[, "lambda"],
       breaks = breaks,
       main = "SBC Rank Distribution (lambda)",
       xlab = paste("Rank (0 to", M_posterior_draws, ")"),
       ylab = "Frequency",
       col = "black")
  abline(h = expected_count, col = "red", lwd = 2, lty = 2)
  legend("topright", legend = "Expected Uniform", lty = 2, col = "red", bty = "n")
  
  # Plot histogram for nu ranks
  hist(sbc_ranks_valid[, "nu"],
       breaks = breaks,
       main = "SBC Rank Distribution (nu)",
       xlab = paste("Rank (0 to", M_posterior_draws, ")"),
       ylab = "Frequency",
       col = "black")
  abline(h = expected_count, col = "red", lwd = 2, lty = 2)
  legend("topright", legend = "Expected Uniform", lty = 2, col = "red", bty = "n")
  
  # Reset plot layout
  par(mfrow = c(1, 1))
  
  print("Histograms plotted. If the bars are roughly flat and close to the red dashed line for both lambda and nu, the Stan sampler is well-calibrated for this model.")
  
} else {
  print("No valid simulations completed. Cannot plot histograms.")
}


