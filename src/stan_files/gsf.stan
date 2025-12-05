/**
 * Generalized Stochastic Frontier Model (Legacy/Alternative version)
 *
 * Implementation of Kumbhakar & Tsionas (2021) "Dissections of input and
 * output efficiency: A generalized stochastic frontier model"
 *
 * NOTE: This is a simplified version. For full implementation with proper
 * truncation corrections, use inst/stan/gsf_model.stan
 */

functions {
  /**
   * Compute the bivariate normal CDF using Drezner & Wesolowsky (1990) approximation
   */
  real binormal_cdf(real x, real y, real rho) {
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

    real result;
    real abs_rho = fabs(rho);

    if (abs_rho < 0.925) {
      real h = -x;
      real k = -y;
      real hk = h * k;
      real hs = (h * h + k * k) / 2;
      real asr = asin(rho);

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

    return fmax(1e-10, fmin(1 - 1e-10, result));
  }

  real log_binormal_cdf(vector a, vector mu, matrix Omega) {
    vector[2] sigma = sqrt(diagonal(Omega));
    real rho = Omega[1, 2] / (sigma[1] * sigma[2]);

    real z1 = (a[1] - mu[1]) / sigma[1];
    real z2 = (a[2] - mu[2]) / sigma[2];

    return log(binormal_cdf(z1, z2, rho));
  }
}

data {
  int<lower=1> N;             // Number of observations
  int<lower=1> J;             // Number of inputs
  int<lower=1> M;             // Number of z variables
  vector[N] y;                // Output (log)
  matrix[N, J] X;             // Inputs (log)
  matrix[N, M] Z;             // Exogenous variables for inefficiency
}

transformed data {
  vector[J] zero_vec = rep_vector(0.0, J);
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
  alpha ~ normal(0, 2);
  to_vector(B_raw) ~ normal(0, 1);
  delta ~ normal(0, 10);
  to_vector(Delta) ~ normal(0, 10);

  sigma_sq ~ inv_gamma(1, 1);
  sigma_u_sq ~ inv_gamma(1, 1);
  Omega ~ wishart(J + 1, diag_matrix(rep_vector(0.1, J)));

  // Output technical inefficiency (half-normal, truncated at 0)
  u0 ~ normal(0, sqrt(sigma_u_sq));
  target += N * log(2);  // Normalizing constant for half-normal

  // Input slacks (truncated multivariate normal)
  for (i in 1:N) {
    vector[J] mu_theta = Delta * Z[i]';
    theta[i]' ~ multi_normal(mu_theta, Omega);
    // Add truncation normalizing constant
    target += -log_binormal_cdf(zero_vec, mu_theta, Omega);
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

generated quantities {
  vector[N] technical_inefficiency;
  matrix[N, J] input_slacks;
  vector[N] output_loss_from_slacks;
  matrix[N, J] elasticities;
  vector[N] returns_to_scale;

  for (i in 1:N) {
    vector[J] x_n = X[i]';
    vector[J] theta_n = theta[i]';
    vector[J] x_eff = x_n + theta_n;
    vector[J] q_observed;
    vector[J] q_effective;

    technical_inefficiency[i] = -u0[i] * 100;

    for (j in 1:J) {
      input_slacks[i, j] = -theta[i, j] * 100;
    }

    // Elasticities at effective inputs
    q_effective = alpha + B * x_eff;
    elasticities[i] = q_effective';
    returns_to_scale[i] = sum(q_effective);

    // Output loss using observed inputs (per equation 5)
    q_observed = alpha + B * x_n;
    output_loss_from_slacks[i] = -(dot_product(theta_n, q_observed) +
                                    0.5 * quad_form(B, theta_n)) * 100;
  }
}
