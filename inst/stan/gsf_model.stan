/**
 * Generalized Stochastic Frontier Model
 *
 * Implementation of Kumbhakar & Tsionas (2021) "Dissections of input and
 * output efficiency: A generalized stochastic frontier model"
 * International Journal of Production Economics, 232, 107940.
 *
 * This model estimates:
 * - Output technical inefficiency (u0 <= 0, half-normal)
 * - Input-specific slacks (theta <= 0, truncated multivariate normal)
 *
 * Key equations from the paper:
 * - Production function (eq. 2): y = alpha0 + delta'z + alpha'(x+theta) + 0.5*(x+theta)'B(x+theta) + u0 + epsilon
 * - Input slacks (eq. 3b): theta = Delta*z + u, theta <= 0
 * - Output loss from slacks (eq. 5): theta'q + 0.5*theta'B*theta, where q = alpha + B*x
 */

functions {
  /**
   * Compute the bivariate normal CDF using Drezner & Wesolowsky (1990) approximation
   * This is needed for the normalizing constant of the truncated bivariate normal
   *
   * @param x Upper limit for first variable
   * @param y Upper limit for second variable
   * @param rho Correlation coefficient
   * @return P(X <= x, Y <= y) where (X,Y) ~ BVN(0, 0, 1, 1, rho)
   */
  real binormal_cdf(real x, real y, real rho) {
    // Handle edge cases
    if (rho == 0) {
      return Phi(x) * Phi(y);
    }
    if (rho == 1) {
      return Phi(fmin(x, y));
    }
    if (rho == -1) {
      if (x + y >= 0) {
        return fmax(0.0, Phi(x) - Phi(-y));
      } else {
        return 0.0;
      }
    }

    // Drezner & Wesolowsky (1990) approximation
    // Using Gauss-Legendre quadrature with 6 points
    real result;
    real abs_rho = fabs(rho);

    if (abs_rho < 0.925) {
      // Use direct formula for moderate correlations
      real h = -x;
      real k = -y;
      real hk = h * k;
      real hs = (h * h + k * k) / 2;
      real asr = asin(rho);

      // Gauss-Legendre weights and abscissas for 6 points
      array[3] real w = {0.1713244923791705, 0.3607615730481384, 0.4679139345726904};
      array[3] real xi = {0.9324695142031522, 0.6612093864662647, 0.2386191860831970};

      real bvn = 0;
      for (i in 1:3) {
        for (j in 1:2) {
          real sn = sin(asr * (j == 1 ? 1 - xi[i] : 1 + xi[i]) / 2);
          bvn += w[i] * exp((sn * hk - hs) / (1 - sn * sn));
        }
      }
      bvn = bvn * asr / (4 * pi());
      result = bvn + Phi(-h) * Phi(-k);
    } else {
      // Use alternative formula for high correlations
      real h = -x;
      real k = -y;

      if (rho < 0) {
        k = -k;
      }

      real hk = h * k;
      real bvn = 0;

      if (abs_rho < 1) {
        real ass = (1 - rho) * (1 + rho);
        real a = sqrt(ass);
        real bs = (h - k) * (h - k);
        real c = (4 - hk) / 8;
        real d = (12 - hk) / 16;
        real asr = -(bs / ass + hk) / 2;

        if (asr > -100) {
          bvn = a * exp(asr) * (1 - c * (bs - ass) * (1 - d * bs / 5) / 3 + c * d * ass * ass / 5);
        }

        if (-hk < 100) {
          real b = sqrt(bs);
          bvn = bvn - exp(-hk / 2) * sqrt(2 * pi()) * Phi(-b / a) * b * (1 - c * bs * (1 - d * bs / 5) / 3);
        }

        a = a / 2;

        // Gauss-Legendre weights and abscissas
        array[5] real w = {0.04717533638651177, 0.1069393259953183, 0.1600783285433464,
                           0.2031674267230659, 0.2334925365383547};
        array[5] real xi = {0.9815606342467191, 0.9041172563704750, 0.7699026741943050,
                            0.5873179542866171, 0.3678314989981802};

        for (i in 1:5) {
          for (j in 1:2) {
            real xs = a * (j == 1 ? 1 - xi[i] : 1 + xi[i]);
            xs = xs * xs;
            real rs = sqrt(1 - xs);
            asr = -(bs / xs + hk) / 2;
            if (asr > -100) {
              bvn += a * w[i] * exp(asr) * (exp(-hk * (1 - rs) / (2 * (1 + rs))) / rs
                     - (1 + c * xs * (1 + d * xs)));
            }
          }
        }
        bvn = -bvn / (2 * pi());
      }

      if (rho > 0) {
        result = bvn + Phi(-fmax(h, k));
      } else {
        result = -bvn;
        if (k > h) {
          result += Phi(k) - Phi(h);
        }
      }
    }

    // Ensure result is in valid range
    return fmax(1e-10, fmin(1 - 1e-10, result));
  }

  /**
   * Log of bivariate normal CDF for truncated MVN normalizing constant
   * Computes log P(X <= a, Y <= b) where (X,Y) ~ N(mu, Sigma)
   *
   * @param a Upper limits (length 2)
   * @param mu Mean vector (length 2)
   * @param L_Sigma Cholesky factor of covariance matrix
   * @return log P(X <= a)
   */
  real log_binormal_cdf(vector a, vector mu, matrix L_Sigma) {
    // Standardize
    vector[2] sigma = sqrt(diagonal(L_Sigma * L_Sigma'));
    real rho = (L_Sigma * L_Sigma')[1, 2] / (sigma[1] * sigma[2]);

    real z1 = (a[1] - mu[1]) / sigma[1];
    real z2 = (a[2] - mu[2]) / sigma[2];

    return log(binormal_cdf(z1, z2, rho));
  }
}

data {
  int<lower=1> N;                    // number of observations
  int<lower=1> J;                    // number of inputs (must be 2 for bivariate normal CDF)
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
  vector[J] zero_vec = rep_vector(0.0, J);  // Upper bound for truncation
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

  // Prior for output technical inefficiency u0
  // u0 ~ N(0, sigma_u) truncated at u0 <= 0 (half-normal)
  // The truncation probability is 0.5 (by symmetry), so we add log(2) per observation
  u0 ~ normal(0, sigma_u);
  target += N * log(2);  // Normalizing constant for half-normal (truncated at 0)

  // Prior for input slacks theta
  // theta_i ~ N_J(Delta * z_i, Sigma) truncated at theta <= 0
  // Must account for the normalizing constant C_i(Omega, Delta)
  for (n in 1:N) {
    vector[J] mu_theta = Delta * Z_slack[n]';

    // Log density of multivariate normal
    theta[n]' ~ multi_normal_cholesky(mu_theta, L_Sigma);

    // Subtract log of truncation probability (normalizing constant)
    // C_i = P(theta <= 0 | theta ~ N(mu_theta, Sigma))
    // For J=2, use bivariate normal CDF
    target += -log_binormal_cdf(zero_vec, mu_theta, L_Sigma);
  }

  // Likelihood
  for (n in 1:N) {
    vector[J] x_n = X[n]';
    vector[J] theta_n = theta[n]';
    vector[J] x_eff = x_n + theta_n;
    real mu_y;

    // Translog production function
    // y = alpha0 + delta0'z_neutral + alpha'(x + theta) + 0.5*(x + theta)'B(x + theta) + u0 + epsilon
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
  matrix[N, J] elasticities;         // input elasticities (at effective inputs)
  vector[N] returns_to_scale;        // sum of elasticities
  vector[N] y_pred;                  // predicted output
  vector[N] log_lik;                 // log likelihood for model comparison

  for (n in 1:N) {
    vector[J] x_n = X[n]';
    vector[J] theta_n = theta[n]';
    vector[J] x_eff = x_n + theta_n;
    vector[J] q_observed;  // elasticity at observed inputs (for output loss)
    vector[J] q_effective; // elasticity at effective inputs (for reporting)
    real mu_y;

    // Technical inefficiency (percentage loss)
    // For small u0, exp(-u0) - 1 ≈ -u0
    technical_inefficiency[n] = -u0[n] * 100;

    // Input slacks (percentage over-use)
    // For small theta, exp(-theta) - 1 ≈ -theta
    for (j in 1:J) {
      input_slacks[n, j] = -theta[n, j] * 100;
    }

    // Output elasticities at EFFECTIVE inputs (for reporting, as in Table 2)
    // Elas = alpha + B*(x + theta)
    q_effective = alpha + B * x_eff;
    elasticities[n] = q_effective';

    // Returns to scale (at effective inputs)
    returns_to_scale[n] = sum(q_effective);

    // Output loss from slacks (Equation 5 in paper)
    // Loss = theta'q + 0.5*theta'B*theta, where q = alpha + B*x (OBSERVED inputs)
    // This is the composed error term contribution from slacks
    q_observed = alpha + B * x_n;
    output_loss_from_slacks[n] = -(dot_product(theta_n, q_observed) +
                                    0.5 * quad_form(B, theta_n)) * 100;

    // Total output loss (technical inefficiency + loss from slacks)
    total_output_loss[n] = technical_inefficiency[n] + output_loss_from_slacks[n];

    // Predicted value (at effective inputs with inefficiency)
    mu_y = alpha0;
    if (M_neutral > 0) {
      mu_y += dot_product(delta0, Z_neutral[n]');
    }
    mu_y += dot_product(alpha, x_eff) +
           0.5 * quad_form(B, x_eff) + u0[n];
    y_pred[n] = mu_y;

    // Log likelihood (for model comparison via WAIC/LOO)
    log_lik[n] = normal_lpdf(y[n] | mu_y, sigma);
  }
}
