# --- Load Libraries ---
library(brms)
source("rcomp.R")


# --- Remove Salved Compiled Stan Model if exists ---
if (file.exists("stan_model.rds")) {
  #Delete file if it exists
  file.remove("stan_model.rds")
}

# SBC parameters
N_simulations <- 1000 # Number of SBC iterations (datasets/fits)
J = 200 # Number of samples to step 2
stan_chains <- 1
stan_iter <- 3000 # Iterations per chain
stan_warmup <- 2980 # Warmup iterations per chain
M_posterior_draws <- stan_chains * (stan_iter - stan_warmup)

sbc_ranks <- matrix(NA, nrow = N_simulations, ncol = 2, dimnames = list(NULL, c("lambda", "nu")))

# --- Auxiliary Functions ---


# --- SBC for COMPoisson ---
for (i in 1:N_simulations) {
  # --- 1. Draw the Joint Prior ---
  lambda_sim <- rnorm(1)
  nu_sim <- rlnorm(1)
  
  while (nu_sim < 0.5) {  
    nu_sim <- rlnorm(1)
  }

  #Z_sim <- # Todo

  # --- 2. Draw a dataset y ~ COMPoisson(lambda_sim, nu_sim) ---
  dataset <- rcomp(J, exp(lambda_sim), nu_sim, sumTo = 100000, eps = 10^-10)

  df <- data.frame(y = dataset)


  # --- 3. MCMC in Stan do obtain samples ---
  leps <- - 3 * log(2)
  
  scode_string <- sprintf("real leps_custom() { return %f; }", leps)
  
  custom_stanvars <- stanvar(scode = scode_string, block = "functions")
  
  fit <- brm(y ~ 1,
             data = df,
             chains = stan_chains,
             iter = stan_iter,
             cores = 4,
             warmup = stan_warmup,
             prior = prior(normal(0, 1), class="Intercept") +
               prior(lognormal(0, 1), class="shape"),
             stanvars = custom_stanvars,
             file = "stan_model.rds",     # saves compiled Stan
             file_refit = "always",
             backend = "cmdstanr",
             family = "com_poisson"
             )
  
  # -- 4. Obtain ranks
  draws <- as_draws_df(fit)
  
  # Extract draws for mu and shape
  lambda_draws <- exp(draws$b_Intercept)
  nu_draws     <- draws$shape
  
  sbc_ranks[i, "lambda"] <- sum(lambda_draws < exp(lambda_sim))
  sbc_ranks[i, "nu"]     <- sum(nu_draws < nu_sim)
  
  print(summary(fit))
  print(sbc_ranks[i, "lambda"])
  print(sbc_ranks[i, "nu"])
  print(mean(lambda_draws))
  print(mean(nu_draws))
  print(exp(lambda_sim))
  print(nu_sim)
  
  if (i %% 1 == 0 || i == 1) {
    print(paste("Simulation", i, "of", N_simulations, "completed."))
  }
  
  
  ##### Plot #####
  valid_ranks <- !is.na(sbc_ranks[, 1]) & !is.na(sbc_ranks[, 2])
  sbc_ranks_valid <- sbc_ranks[valid_ranks, , drop = FALSE]
  N_valid_simulations <- nrow(sbc_ranks_valid)
  
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
}

print("SBC simulations finished.")

if (N_valid_simulations < N_simulations) {
  print(paste("Warning:", N_simulations - N_valid_simulations, "simulations had issues and were excluded."))
}
