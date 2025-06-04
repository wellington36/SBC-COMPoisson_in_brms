// Model: y ~ Normal(mu, sigma)
// Priors: mu ~ Normal(0, 1), sigma ~ Lognormal(0, 1)

transformed data {
  // --- Simulation Step ---
  // Simulate parameters from their priors
  real mu_sim = normal_rng(0, 1);
  real<lower=0> sigma_sim = lognormal_rng(0, 1);

  // Simulate data given the simulated parameters
  int<lower=0> J = 10; // Number of data points per simulation
  array[J] real y_sim;
  for (j in 1:J) {
    y_sim[j] = normal_rng(mu_sim, sigma_sim);
  }

  // Pass simulated values to generated quantities through transformed data
  array[2] real sim_params;
  sim_params[1] = mu_sim;
  sim_params[2] = sigma_sim;
}

parameters {
  // --- Parameters for Inference Step ---
  real mu;
  real<lower=0> sigma;
}

model {
  // --- Model Definition for Inference Step ---
  // Priors
  mu ~ normal(0, 1);
  sigma ~ lognormal(0, 1);

  // Likelihood
  y_sim ~ normal(mu, sigma);
}

generated quantities {
  // --- Rank Calculation Step ---
  // Compare posterior draws to the original simulated values
  array[2] int<lower=0, upper=1> lt_sim;
  lt_sim[1] = mu < sim_params[1];    // Compare mu with mu_sim
  lt_sim[2] = sigma < sim_params[2]; // Compare sigma with sigma_sim
}