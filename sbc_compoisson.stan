// Model: y ~ CMP(lambda, nu)
// Priors: lambda ~ Normal(0, 1), nu ~ Lognormal(0, 1)

functions {
  // log of kth term of the normalizing series of the COM Poisson distribution
  // Args:
  //   log_lambda: log location parameter
  //   shape: positive shape parameter
  //   k: k-th term
  real log_k_term(real log_lambda, real nu, int k) {
    return (k - 1) * log_lambda - nu * lgamma(k);
  }

  // bound for the remainder of the normalizing series of the COM Poisson
  // distribution given the last two terms in log-scale
  // Args:
  //   k_current_term: the log of a_k term
  //   k_previous_term: the log of a_(k-1) term
  real bound_remainder(real k_current_term, real k_previous_term) {
    return k_current_term - log(- expm1(k_current_term - k_previous_term));
  }

  // stopping criterio with bucket
  // Args:
  //   k_current_term: the log of a_k term
  //   k_previous_term: the log of a_(k-1) term
  //   k: k term of the series
  //   leps: log(eps)
  int stopping_criterio_bucket(real k_current_term, real k_previous_term, int k, real leps) {
    if (k % 50 == 0) {
      return (bound_remainder(k_current_term, k_previous_term) >= leps);
    }
    return (1e300 >= leps); // Int > leps
  }

  // log normalizing constant of the COM Poisson distribution
  // implementation inspired by code of Ben Goodrich
  // improved following suggestions of Sebastian Weber (#892)
  // Args:
  //   log_lambda: log location parameter
  //   shape: positive shape parameter
  real log_Z_com_poisson(real log_lambda, real nu) {
    real log_Z;
    int k = 2;
    int M = 10000;
    real leps = -23 * log2();
    vector[M] log_Z_terms;

    if (nu == 1) {
      return exp(log_lambda);
    }
    // nu == 0 or Inf will fail in this parameterization
    if (nu <= 0) {
      reject("nu lambdast be positive");
    }
    if (nu == positive_infinity()) {
      reject("nu lambdast be finite");
    }

    // direct computation of the truncated series
    // check if the Mth term of the series pass in the stopping criteria
    if (bound_remainder(log_k_term(log_lambda, nu, M),
                        log_k_term(log_lambda, nu, M - 1)) >= leps) {
      reject("nu is too close to zero.");
    }

    // first 2 terms of the series
    log_Z_terms[1] = log_k_term(log_lambda, nu, 1);
    log_Z_terms[2] = log_k_term(log_lambda, nu, 2);

    while (((log_Z_terms[k] >= log_Z_terms[k-1]) ||
      (stopping_criterio_bucket(log_Z_terms[k], log_Z_terms[k-1], k, leps))) &&
      k < M) {
      k += 1;
      log_Z_terms[k] = log_k_term(log_lambda, nu, k);
    }
    log_Z = log_sum_exp(log_Z_terms[1:k]);

    return log_Z;
  }
  // COM Poisson log-PMF for a single response (log parameterization)
  // Args:
  //   y: the response value
  //   log_lambda: log location parameter
  //   shape: positive shape parameter
  real com_poisson_log_lpmf(int y, real log_lambda, real nu) {
    if (nu == 1) return poisson_log_lpmf(y | log_lambda);
    return y * log_lambda - nu * lgamma(y + 1) - log_Z_com_poisson(log_lambda, nu);
  }
  // COM Poisson log-PMF for a single response
  real com_poisson_lpmf(int y, real lambda, real nu) {
    if (nu == 1) return poisson_lpmf(y | lambda);
    return com_poisson_log_lpmf(y | log(lambda), nu);
  }
  // COM Poisson log-CDF for a single response
  real com_poisson_lcdf(int y, real lambda, real nu) {
    real log_lambda;
    real log_Z;  // log denominator
    vector[y] log_num_terms;  // terms of the log numerator
    if (nu == 1) {
      return poisson_lcdf(y | lambda);
    }
    // nu == 0 or Inf will fail in this parameterization
    if (nu <= 0) {
      reject("nu lambdast be positive");
    }
    if (nu == positive_infinity()) {
      reject("nu lambdast be finite");
    }
    if (y > 10000) {
      reject("cannot handle y > 10000");
    }
    log_lambda = log(lambda);
    if (y * log_lambda - nu * lgamma(y + 1) <= -36.0) {
      // y is large enough for the CDF to be very close to 1;
      return 0;
    }
    log_Z = log_Z_com_poisson(log_lambda, nu);
    if (y == 0) {
      return -log_Z;
    }
    // first 2 terms of the series
    log_num_terms[1] = log1p_exp(nu * log_lambda);
    // remaining terms of the series until y
    for (k in 2:y) {
      log_num_terms[k] = k * log_lambda - nu * lgamma(k + 1);
    }
    return log_sum_exp(log_num_terms) - log_Z;
  }
  // COM Poisson log-CCDF for a single response
  real com_poisson_lccdf(int y, real lambda, real nu) {
    return log1m_exp(com_poisson_lcdf(y | lambda, nu));
  }
}

transformed data {
  // --- Silambdalation Step ---
  // Silambdalate parameters from their priors
  real<lower=0> lambda_sim = lognormal_rng(0, 1);
  real<lower=0> nu_sim = lognormal_rng(0, 1);
  
  // Silambdalate data given the silambdalated parameters
  int<lower=0> J = 500;
  array[J] real y_sim;
  for (j in 1:J) {
    y_sim[j] = exp(com_poisson_lpmf(lambda_sim, nu_sim));
  }
  
  // Pass simulated values to generated quantities through transformed data
  array[2] real sim_params;
  sim_params[1] = lambda_sim;
  sim_params[2] = nu_sim;
}

parameters {
  real<lower=0> lambda;
  real<lower=0> nu;
}

model {
  // --- Model Definition for Inference Step ---
  // Priors
  lambda ~ lognormal(0, 1);
  nu ~ lognormal(0, 1);
  
  // Likelihood
  y_sim ~ com_poisson_lcdf(lambda, nu);
}

generated quantities {
  // --- Rank Calculation Step ---
  // Compare posterior draws to the original simulated values
  array[2] int<lower=0, upper=1> lt_sim;
  lt_sim[1] = lambda < sim_params[1];    // Compare lambda with lambda_sim
  lt_sim[2] = nu < sim_params[2]; // Compare nu with nu_sim
}