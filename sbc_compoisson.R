# R Script for Silambdalation-Based Calibration (SBC)
# Using the lambda-nu example from the Stan User's Guide

# --- Prerequisites ---
# Ensure you have RStan installed and configured:
# install.packages("rstan", repos = "https://cloud.r-project.org/", dependencies = TRUE)
# pkgbuild::has_build_tools(debug = TRUE) # Check C++ toolchain

# --- Load Libraries ---
library(rstan)

# --- Configuration ---
# set.seed(789) # for reproducibility

# SBC parameters
N_silambdalations <- 5000 # Number of SBC iterations (datasets/fits)

# Stan configuration
stan_model_file <- "sbc_compoisson.stan" # Path to the Stan model file
stan_chains <- 2
stan_iter <- 500 # Iterations per chain
stan_warlambdap <- 250 # Warlambdap iterations per chain
M_posterior_draws <- stan_chains * (stan_iter - stan_warlambdap)

# --- Compile Stan Model ---
# This compiles the C++ code for the model. It can take a minute.
print(paste("Compiling Stan model:", stan_model_file))
tryCatch({
  compiled_model <- stan_model(stan_model_file)
  print("Stan model compiled successfully.")
}, error = function(e) {
  stop(paste("Error compiling Stan model:", e$message,
             "\nPlease check the Stan code in", stan_model_file))
})

# --- SBC Function ---
# Initialize matrices to store ranks for lambda and nu
sbc_ranks <- matrix(NA, nrow = N_silambdalations, ncol = 2, dimnames = list(NULL, c("lambda", "nu")))

print(paste("Starting", N_silambdalations, "SBC silambdalations using Stan..."))
print(paste("Each silambdalation uses M=", M_posterior_draws, "posterior draws."))

# --- Run Silambdalations ---
# This loop can be slow as it fits the Stan model N_silambdalations times.
for (i in 1:N_silambdalations) {
  # Run Stan sampler
  # The Stan code handles silambdalation internally (transformed data)
  # We don't need to pass any data specific to this iteration
  # Using refresh = 0 to minimize console output during the loop
  fit <- sampling(compiled_model,
                  data = list(), # No external data needed for this model
                  chains = stan_chains,
                  iter = stan_iter,
                  warlambdap = stan_warlambdap,
                  seed = sample.int(.Machine$integer.max, 1), # Use different seed for each Stan run
                  refresh = 0, # Suppress iteration updates
                  cores = getOption("mc.cores", 1)) # Use lambdaltiple cores if available
  
  # Extract the indicator variables 'lt_sim'
  # lt_sim[m, k] = 1 if posterior_draw[m, k] < silambdalated_param[k], 0 otherwise
  # Dimensions: M_posterior_draws x 2 (for lambda and nu)
  lt_sim_draws <- extract(fit, pars = "lt_sim")$lt_sim
  
  # Check dimensions
  if (is.null(lt_sim_draws) || nrow(lt_sim_draws) != M_posterior_draws || ncol(lt_sim_draws) != 2) {
    warning(paste("Silambdalation", i, ": Problem extracting 'lt_sim' draws. Skipping rank calculation."))
    next # Skip to the next silambdalation
  }
  
  # Calculate ranks for this silambdalation
  # Rank = sum over posterior draws of (posterior_draw < silambdalated_param)
  # This corresponds to the column sums of the lt_sim matrix
  ranks_i <- colSums(lt_sim_draws)
  sbc_ranks[i, "lambda"] <- ranks_i[1]
  sbc_ranks[i, "nu"] <- ranks_i[2]
  
  # Optional: Print progress
  if (i %% 50 == 0 || i == 1) {
    print(paste("Silambdalation", i, "of", N_silambdalations, "completed."))
  }
}

print("SBC silambdalations finished.")

# --- Check Uniformity ---

# Remove any rows with NA ranks (e.g., if extraction failed)
valid_ranks <- !is.na(sbc_ranks[, 1]) & !is.na(sbc_ranks[, 2])
sbc_ranks_valid <- sbc_ranks[valid_ranks, , drop = FALSE]
N_valid_silambdalations <- nrow(sbc_ranks_valid)

if (N_valid_silambdalations < N_silambdalations) {
  print(paste("Warning:", N_silambdalations - N_valid_silambdalations, "silambdalations had issues and were excluded."))
}

if (N_valid_silambdalations > 0) {
  print(paste("Plotting histograms for", N_valid_silambdalations, "valid silambdalations."))
  
  # Define breaks for the histogram (M+1 bins from -0.5 to M+0.5)
  # Ranks range from 0 to M_posterior_draws
  breaks <- seq(-0.5, M_posterior_draws + 0.5, by = 1)
  num_bins <- length(breaks) - 1
  
  # Calculate expected count per bin for uniform distribution
  expected_count <- N_valid_silambdalations / num_bins
  
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
  print("No valid silambdalations completed. Cannot plot histograms.")
}
