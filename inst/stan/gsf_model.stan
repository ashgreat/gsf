data {
  int<lower=1> N;                    // number of observations
  int<lower=1> J;                    // number of inputs
  int<lower=1> M_slack;              // z vars driving slacks (includes intercept)
  int<lower=0> M_neutral;            // neutral shifters

  vector[N] y;                       // log output
  matrix[N, J] X;                    // log inputs
  matrix[N, M_slack] Z_slack;        // determinants of inefficiency (with intercept)
  matrix[N, M_neutral] Z_neutral;    // neutral technology shifters

  // Prior hyperparameters
  real<lower=0> prior_sigma_scale;
  real<lower=0> prior_sigma_u_scale;
  real<lower=0> prior_delta_scale;
  real<lower=0> prior_omega_scale;
}

transformed data {
  // Pre-compute for efficiency
  vector[J*(J+1)/2] zero_vech = rep_vector(0.0, J*(J+1)/2);
}

parameters {
  // Production function parameters
  real alpha0;                       // intercept
  vector[J] alpha;                   // first-order coefficients
  vector[J*(J+1)/2] B_vech;         // second-order coefficients (lower triangle + diagonal)

  // Neutral technology shifters
  vector[M_neutral] delta0;          // coefficients for neutral shifters

  // Input slack parameters
  matrix[J, M_slack] Delta_raw;     // coefficients for z in slack equations
  cholesky_factor_corr[J] L_Omega;  // correlation matrix for slacks
  vector<lower=0>[J] sigma_slack;   // scale of slack distributions

  // Latent variables
  matrix<upper=0>[N, J] theta;      // input slacks (all negative)
  vector<upper=0>[N] u0;            // output technical inefficiency (negative)

  // Variance parameters
  real<lower=0> sigma;              // noise standard deviation
  real<lower=0> sigma_u;            // output inefficiency scale
}

transformed parameters {
  // Construct symmetric B matrix from vech
  matrix[J, J] B;
  matrix[J, M_slack] Delta = Delta_raw;
  matrix[J, J] L_Sigma;

  // Fill in B matrix (symmetric)
  {
    int idx = 1;
    for (i in 1:J) {
      for (j in i:J) {
        B[i, j] = B_vech[idx];
        B[j, i] = B_vech[idx];
        idx = idx + 1;
      }
    }
  }

  // Covariance factor for slacks
  L_Sigma = diag_pre_multiply(sigma_slack, L_Omega);
}

model {
  // Priors for production function parameters
  alpha0 ~ normal(0, 10);
  alpha ~ normal(0, 2);
  B_vech ~ normal(0, 1);

  // Priors for technology shifters
  delta0 ~ normal(0, prior_delta_scale);

  // Priors for slack equation parameters
  to_vector(Delta_raw) ~ normal(0, prior_delta_scale);
  L_Omega ~ lkj_corr_cholesky(2);
  sigma_slack ~ exponential(1.0 / prior_omega_scale);

  // Priors for variance parameters
  sigma ~ exponential(1.0 / prior_sigma_scale);
  sigma_u ~ exponential(1.0 / prior_sigma_u_scale);

  // Priors for latent effects
  for (n in 1:N) {
    vector[J] mu_theta = Delta * Z_slack[n]';
    theta[n]' ~ multi_normal_cholesky(mu_theta, L_Sigma);
  }
  u0 ~ normal(0, sigma_u);

  // Likelihood
  for (n in 1:N) {
    vector[J] x_n = X[n]';
    vector[J] theta_n = theta[n]';
    vector[J] x_eff = x_n + theta_n;
    real mu_y;

    // Translog production function
    mu_y = alpha0;
    if (M_neutral > 0) {
      mu_y += dot_product(delta0, Z_neutral[n]');
    }
    mu_y += dot_product(alpha, x_eff) +
           0.5 * quad_form(B, x_eff);

    // Add output inefficiency
    mu_y = mu_y + u0[n];

    // Likelihood contribution
    y[n] ~ normal(mu_y, sigma);
  }
}

generated quantities {
  // Efficiency measures
  vector[N] technical_inefficiency;  // exp(-u0) - 1, percentage
  matrix[N, J] input_slacks;         // -theta, percentage over-use
  vector[N] output_loss_from_slacks; // loss due to input slacks
  vector[N] total_output_loss;       // total loss
  matrix[N, J] elasticities;         // input elasticities
  vector[N] returns_to_scale;        // sum of elasticities
  vector[N] y_pred;                  // predicted output
  vector[N] log_lik;                 // log likelihood for model comparison

  for (n in 1:N) {
    vector[J] x_n = X[n]';
    vector[J] theta_n = theta[n]';
    vector[J] x_eff = x_n + theta_n;
    vector[J] q_n;  // elasticity vector
    real mu_y;

    // Technical inefficiency (percentage loss)
    technical_inefficiency[n] = (exp(-u0[n]) - 1) * 100;

    // Input slacks (percentage over-use)
    for (j in 1:J) {
      input_slacks[n, j] = (exp(-theta[n, j]) - 1) * 100;
    }

    // Output elasticities: q = alpha + B*x_eff
    q_n = alpha + B * x_eff;
    elasticities[n] = q_n';

    // Returns to scale
    returns_to_scale[n] = sum(q_n);

    // Output loss from slacks: theta'q + 0.5*theta'B*theta
    output_loss_from_slacks[n] = -(dot_product(theta_n, q_n) +
                                    0.5 * quad_form(B, theta_n)) * 100;

    // Total output loss
    total_output_loss[n] = technical_inefficiency[n] + output_loss_from_slacks[n];

    // Predicted value
    mu_y = alpha0;
    if (M_neutral > 0) {
      mu_y += dot_product(delta0, Z_neutral[n]');
    }
    mu_y += dot_product(alpha, x_eff) +
           0.5 * quad_form(B, x_eff) + u0[n];
    y_pred[n] = mu_y;

    // Log likelihood
    log_lik[n] = normal_lpdf(y[n] | mu_y, sigma);
  }
}
