data {
  int<lower=1> N;             // Number of observations
  int<lower=1> J;             // Number of inputs
  int<lower=1> M;             // Number of z variables
  vector[N] y;                // Output (log)
  matrix[N, J] X;             // Inputs (log)
  matrix[N, M] Z;             // Exogenous variables for inefficiency
}

parameters {
  real alpha0;                // Intercept
  vector[J] alpha;            // First order coefficients
  matrix[J, J] B_raw;         // Second order coefficients (symmetric)
  vector[M] delta;            // Coefficients for z in production function
  matrix[J, M] Delta;         // Coefficients for z in input slacks
  
  cov_matrix[J] Omega;        // Covariance of input slacks
  real<lower=0> sigma_sq;     // Noise variance
  real<lower=0> sigma_u_sq;   // Output inefficiency variance
  
  matrix<upper=0>[N, J] theta; // Input slacks (latent, negative)
  vector<upper=0>[N] u0;       // Output inefficiency (latent, negative)
}

transformed parameters {
  matrix[J, J] B;
  // Enforce symmetry for B
  B = 0.5 * (B_raw + B_raw');
}

model {
  // Priors
  alpha0 ~ normal(0, 10);
  alpha ~ normal(0, 10);
  to_vector(B_raw) ~ normal(0, 10);
  delta ~ normal(0, 10);
  to_vector(Delta) ~ normal(0, 10);
  
  sigma_sq ~ inv_gamma(1, 1);
  sigma_u_sq ~ inv_gamma(1, 1);
  // Omega prior is implicit (uniform on cov matrices) or can be Wishart
  
  // Likelihood for latent variables
  for (i in 1:N) {
    theta[i] ~ multi_normal(Delta * Z[i]', Omega);
    u0[i] ~ normal(0, sqrt(sigma_u_sq));
  }
  
  // Likelihood for data
  for (i in 1:N) {
    vector[J] X_eff = X[i]' + theta[i]';
    real mu = alpha0 + dot_product(delta, Z[i]) + 
              dot_product(alpha, X_eff) + 
              0.5 * quad_form(B, X_eff) + 
              u0[i];
    y[i] ~ normal(mu, sqrt(sigma_sq));
  }
}
