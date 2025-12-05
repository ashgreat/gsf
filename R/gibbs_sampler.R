#' Gibbs Sampler for GSF Model
#'
#' @description
#' Implements the Gibbs sampler with data augmentation for the Generalized
#' Stochastic Frontier model following Kumbhakar and Tsionas (2021), Appendix B.
#'
#' @keywords internal

#' Sample from Truncated Univariate Normal Distribution
#'
#' Samples from N(mu, sigma^2) truncated to (-Inf, upper]
#'
#' @param n Number of samples
#' @param mu Mean of the untruncated normal
#' @param sigma Standard deviation
#' @param upper Upper bound (default 0)
#' @return Vector of samples
#' @keywords internal
rtruncnorm_upper <- function(n, mu, sigma, upper = 0) {
  # Using inverse CDF method
  u <- runif(n)
  # P(X <= upper) = Phi((upper - mu) / sigma)
  p_upper <- pnorm((upper - mu) / sigma)
  # Sample from U(0, p_upper) and invert
  qnorm(u * p_upper) * sigma + mu
}

#' Sample from Truncated Multivariate Normal (Geweke 1991)
#'
#' Samples from N_J(mu, Sigma) truncated to x <= upper using Gibbs sampling.
#' Based on Geweke (1991) "Efficient simulation from the multivariate normal
#' and student-t distributions subject to linear constraints."
#'
#' @param mu Mean vector (J x 1)
#' @param Sigma Covariance matrix (J x J)
#' @param upper Upper bounds (J x 1), default all zeros
#' @param n_gibbs Number of Gibbs iterations for sampling
#' @param init Initial values (if NULL, starts at conditional means)
#' @return A vector of length J
#' @keywords internal
rtmvnorm_geweke <- function(mu, Sigma, upper = rep(0, length(mu)),
                             n_gibbs = 10, init = NULL) {
  J <- length(mu)

  # Precompute conditional parameters
  # For each j, we need: mu_j|x_{-j} and sigma^2_j|x_{-j}
  # These are linear functions of x_{-j}

  # Cholesky decomposition
  L <- chol(Sigma)
  Sigma_inv <- chol2inv(L)

  # Initialize
  if (is.null(init)) {
    # Start at truncated conditional means
    x <- pmin(mu, upper - 0.01)
  } else {
    x <- init
  }

  # Gibbs sampling
  for (iter in 1:n_gibbs) {
    for (j in 1:J) {
      # Conditional distribution of x_j | x_{-j}
      # mu_j|{-j} = mu_j - Sigma_j{-j} * Sigma_{-j}{-j}^{-1} * (x_{-j} - mu_{-j})
      # sigma^2_j|{-j} = Sigma_jj - Sigma_j{-j} * Sigma_{-j}{-j}^{-1} * Sigma_{-j}j

      if (J == 1) {
        cond_mean <- mu[1]
        cond_var <- Sigma[1, 1]
      } else {
        idx_j <- j
        idx_notj <- setdiff(1:J, j)

        # Sigma_j{-j}: row j, columns not j
        Sigma_j_notj <- Sigma[idx_j, idx_notj, drop = FALSE]
        # Sigma_{-j}{-j}: submatrix excluding row and col j
        Sigma_notj_notj <- Sigma[idx_notj, idx_notj, drop = FALSE]
        Sigma_notj_notj_inv <- solve(Sigma_notj_notj)

        # Conditional mean and variance
        cond_mean <- mu[j] + Sigma_j_notj %*% Sigma_notj_notj_inv %*%
          (x[idx_notj] - mu[idx_notj])
        cond_var <- Sigma[j, j] - Sigma_j_notj %*% Sigma_notj_notj_inv %*%
          t(Sigma_j_notj)
      }

      cond_sd <- sqrt(as.numeric(cond_var))

      # Sample from truncated normal
      x[j] <- rtruncnorm_upper(1, as.numeric(cond_mean), cond_sd, upper[j])
    }
  }

  return(x)
}

#' Compute Log of Multivariate Normal CDF
#'
#' Computes log P(X <= a) where X ~ N(mu, Sigma) using Monte Carlo.
#'
#' @param a Upper bounds (J x 1)
#' @param mu Mean vector (J x 1)
#' @param Sigma Covariance matrix (J x J)
#' @param n_mc Number of Monte Carlo samples
#' @return Log probability
#' @keywords internal
log_mvnorm_cdf <- function(a, mu, Sigma, n_mc = 1000) {
  J <- length(mu)

  if (J == 1) {
    return(pnorm((a - mu) / sqrt(Sigma), log.p = TRUE))
  }

  if (J == 2) {
    # Use bivariate normal CDF (more accurate)
    return(log_binormal_cdf_r(a, mu, Sigma))
  }

  # Monte Carlo for J > 2
  L <- t(chol(Sigma))
  count <- 0
  for (i in 1:n_mc) {
    z <- rnorm(J)
    x <- mu + L %*% z
    if (all(x <= a)) {
      count <- count + 1
    }
  }

  # Avoid log(0)
  p <- max(count / n_mc, 1e-10)
  return(log(p))
}

#' Bivariate Normal CDF (R implementation)
#'
#' Computes log P(X <= a) for bivariate normal using Drezner & Wesolowsky (1990).
#'
#' @param a Upper bounds (2 x 1)
#' @param mu Mean vector (2 x 1)
#' @param Sigma Covariance matrix (2 x 2)
#' @return Log probability
#' @keywords internal
log_binormal_cdf_r <- function(a, mu, Sigma) {
  # Standardize
  sigma1 <- sqrt(Sigma[1, 1])
  sigma2 <- sqrt(Sigma[2, 2])
  rho <- Sigma[1, 2] / (sigma1 * sigma2)

  z1 <- (a[1] - mu[1]) / sigma1
  z2 <- (a[2] - mu[2]) / sigma2

  # Use pbivnorm if available, otherwise approximate
  if (requireNamespace("pbivnorm", quietly = TRUE)) {
    p <- pbivnorm::pbivnorm(z1, z2, rho)
  } else {
    # Approximate using pmvnorm from mvtnorm
    if (requireNamespace("mvtnorm", quietly = TRUE)) {
      corr <- matrix(c(1, rho, rho, 1), 2, 2)
      p <- mvtnorm::pmvnorm(upper = c(z1, z2), corr = corr)[1]
    } else {
      # Very rough approximation using independence
      p <- pnorm(z1) * pnorm(z2)
    }
  }

  return(log(max(p, 1e-10)))
}

#' Main Gibbs Sampler for GSF Model
#'
#' @param y Log output vector (N x 1)
#' @param X Log input matrix (N x J)
#' @param Z_slack Slack determinants with intercept (N x M)
#' @param Z_neutral Neutral shifters (N x M_neutral)
#' @param n_iter Number of MCMC iterations
#' @param n_burnin Number of burn-in iterations
#' @param n_thin Thinning rate
#' @param prior Prior hyperparameters list
#' @param seed Random seed
#' @param verbose Print progress
#' @return List with posterior samples
#' @keywords internal
gsf_gibbs_sampler <- function(y, X, Z_slack, Z_neutral = NULL,
                               n_iter = 10000,
                               n_burnin = 5000,
                               n_thin = 1,
                               prior = NULL,
                               seed = NULL,
                               verbose = TRUE) {

  if (!is.null(seed)) {
    set.seed(seed)
  }

  N <- length(y)
  J <- ncol(X)
  M <- ncol(Z_slack)
  M_neutral <- if (is.null(Z_neutral)) 0 else ncol(Z_neutral)

  # Number of translog parameters: intercept + J first-order + J(J+1)/2 second-order
  n_B_vech <- J * (J + 1) / 2

  # Set default priors
  if (is.null(prior)) {
    prior <- list(
      # For delta (production function parameters)
      delta_mean = rep(0, 1 + J + n_B_vech + M_neutral),
      delta_var = diag(100, 1 + J + n_B_vech + M_neutral),
      # For sigma^2 (noise variance): IG(a_sigma, b_sigma)
      a_sigma = 1,
      b_sigma = 1,
      # For sigma_u^2 (output inefficiency): IG(a_u, b_u)
      a_u = 1,
      b_u = 1,
      # For Delta (slack coefficients): each row ~ N(0, prior_Delta_var * I)
      Delta_var = 100,
      # For Omega (slack covariance): Wishart(nu_Omega, S_Omega)
      nu_Omega = J + 1,
      S_Omega = diag(0.1, J)
    )
  }

  # Initialize parameters
  # Production function parameters (delta = [alpha0, alpha, B_vech, delta0])
  n_delta <- 1 + J + n_B_vech + M_neutral
  delta <- rep(0, n_delta)
  delta[1] <- mean(y)  # intercept
  delta[2:(J + 1)] <- 1 / J  # equal elasticities

  # Variance parameters - start with reasonable values
  sigma_sq <- var(y) * 0.05  # Assume 5% of variance is noise
  sigma_u_sq <- 0.02^2  # Start with 2% SD for inefficiency

  # Slack parameters
  Delta <- matrix(0, J, M)
  Delta[, 1] <- -0.05  # Small negative intercept for slacks (mean slack ~5%)

  # Slack covariance - start with small values
  Omega <- diag(0.03^2, J)  # 3% SD for slacks
  # Add small correlation
  if (J == 2) {
    Omega[1, 2] <- Omega[2, 1] <- 0.5 * 0.03^2
  }
  L_Omega <- chol(Omega)

  # Latent variables - initialize to prior means
  theta <- matrix(0, N, J)
  for (i in 1:N) {
    theta[i, ] <- pmin(as.vector(Delta %*% Z_slack[i, ]), -0.01)  # Ensure negative
  }
  u0 <- rep(-0.02, N)  # Output inefficiency (negative), ~2% mean

  # Storage for samples
  n_samples <- floor((n_iter - n_burnin) / n_thin)
  samples <- list(
    alpha0 = numeric(n_samples),
    alpha = matrix(0, n_samples, J),
    B_vech = matrix(0, n_samples, n_B_vech),
    delta0 = if (M_neutral > 0) matrix(0, n_samples, M_neutral) else NULL,
    sigma = numeric(n_samples),
    sigma_u = numeric(n_samples),
    Delta = array(0, c(n_samples, J, M)),
    Omega = array(0, c(n_samples, J, J)),
    theta = array(0, c(n_samples, N, J)),
    u0 = matrix(0, n_samples, N)
  )

  # Acceptance rates for Metropolis steps
  accept_theta <- rep(0, N)
  accept_Omega <- 0
  accept_Delta <- 0
  n_mh_theta <- 0
  n_mh_Omega <- 0
  n_mh_Delta <- 0

  # Adaptive proposal scales
  theta_prop_scale <- rep(1.0, N)
  Delta_prop_sd <- 0.02  # Start smaller
  Omega_prop_scale <- 0.1  # Scale for Omega proposals

  # Helper function: construct B matrix from vech
  vech_to_B <- function(B_vech, J) {
    B <- matrix(0, J, J)
    idx <- 1
    for (i in 1:J) {
      for (j in i:J) {
        B[i, j] <- B_vech[idx]
        B[j, i] <- B_vech[idx]
        idx <- idx + 1
      }
    }
    return(B)
  }

  # Helper: compute effective inputs
  compute_x_eff <- function(X, theta) {
    X + theta
  }

  # Helper: compute production function value
  compute_mu_y <- function(X_eff, delta, Z_neutral, J, n_B_vech, M_neutral) {
    alpha0 <- delta[1]
    alpha <- delta[2:(J + 1)]
    B_vech <- delta[(J + 2):(J + 1 + n_B_vech)]
    B <- vech_to_B(B_vech, J)

    N <- nrow(X_eff)
    mu <- numeric(N)

    for (i in 1:N) {
      x_eff <- X_eff[i, ]
      mu[i] <- alpha0 + sum(alpha * x_eff) + 0.5 * t(x_eff) %*% B %*% x_eff
    }

    if (M_neutral > 0) {
      delta0 <- delta[(J + 2 + n_B_vech):length(delta)]
      mu <- mu + Z_neutral %*% delta0
    }

    return(mu)
  }

  # Main Gibbs loop
  sample_idx <- 0

  if (verbose) {
    message("Starting Gibbs sampler...")
    message("Iterations: ", n_iter, " (burn-in: ", n_burnin, ", thin: ", n_thin, ")")
    pb <- txtProgressBar(min = 0, max = n_iter, style = 3)
  }

  for (iter in 1:n_iter) {
    # Current effective inputs
    X_eff <- compute_x_eff(X, theta)

    # ==================================================================
    # Step 1: Sample delta (production function parameters)
    # ==================================================================
    # y_i = mu_i + u0_i + epsilon_i
    # where mu_i = alpha0 + alpha'x_eff_i + 0.5*x_eff'B*x_eff + delta0'z_neutral
    # This is linear in delta if we construct the design matrix appropriately

    # Construct design matrix for delta
    W <- matrix(0, N, n_delta)
    W[, 1] <- 1  # intercept
    W[, 2:(J + 1)] <- X_eff  # first-order terms

    # Second-order terms (vech of x_eff %*% t(x_eff))
    for (i in 1:N) {
      x_eff <- X_eff[i, ]
      idx <- J + 2
      for (j1 in 1:J) {
        for (j2 in j1:J) {
          if (j1 == j2) {
            W[i, idx] <- 0.5 * x_eff[j1]^2
          } else {
            W[i, idx] <- x_eff[j1] * x_eff[j2]  # off-diagonal: x1*x2 (not 0.5)
          }
          idx <- idx + 1
        }
      }
    }

    # Neutral shifters
    if (M_neutral > 0) {
      W[, (n_delta - M_neutral + 1):n_delta] <- Z_neutral
    }

    # Response for delta regression: y - u0
    y_adj <- y - u0

    # Posterior for delta: normal conjugate
    # Prior: delta ~ N(delta_mean, delta_var)
    # Likelihood: y_adj ~ N(W*delta, sigma_sq*I)
    # Posterior: delta | ... ~ N(delta_post_mean, delta_post_var)

    delta_prior_prec <- solve(prior$delta_var)
    delta_lik_prec <- (1 / sigma_sq) * crossprod(W)

    delta_post_var <- solve(delta_prior_prec + delta_lik_prec)
    delta_post_mean <- delta_post_var %*% (delta_prior_prec %*% prior$delta_mean +
                                             (1 / sigma_sq) * crossprod(W, y_adj))

    delta <- as.vector(mvtnorm::rmvnorm(1, delta_post_mean, delta_post_var))

    # ==================================================================
    # Step 2: Sample sigma^2 (noise variance)
    # ==================================================================
    # epsilon_i = y_i - mu_i - u0_i
    mu_y <- compute_mu_y(X_eff, delta, Z_neutral, J, n_B_vech, M_neutral)
    epsilon <- y - mu_y - u0

    # Posterior: IG(a_post, b_post)
    a_post <- prior$a_sigma + N / 2
    b_post <- prior$b_sigma + sum(epsilon^2) / 2

    sigma_sq <- 1 / rgamma(1, a_post, b_post)

    # ==================================================================
    # Step 3: Sample sigma_u^2 (output inefficiency variance)
    # ==================================================================
    # u0_i ~ N(0, sigma_u^2) truncated at u0 <= 0
    # For half-normal, the sufficient statistic is sum(u0^2)
    # Note: E[|u0|] = sigma_u * sqrt(2/pi) for half-normal

    a_u_post <- prior$a_u + N / 2
    b_u_post <- prior$b_u + sum(u0^2) / 2

    sigma_u_sq <- 1 / rgamma(1, a_u_post, b_u_post)

    # Bound sigma_u to prevent it from becoming too large
    # (helps with identification when signal is weak)
    sigma_u_sq <- min(sigma_u_sq, 0.5^2)  # Cap at 50% inefficiency SD

    # ==================================================================
    # Step 4: Sample u0 (output technical inefficiency)
    # ==================================================================
    # u0_i | ... ~ N(u0_hat_i, sigma_star^2) truncated at u0 <= 0
    # where sigma_star^2 = sigma^2 * sigma_u^2 / (sigma^2 + sigma_u^2)
    # and u0_hat_i = sigma_star^2 * epsilon_i / sigma^2

    sigma_star_sq <- sigma_sq * sigma_u_sq / (sigma_sq + sigma_u_sq)
    sigma_star <- sqrt(sigma_star_sq)

    for (i in 1:N) {
      # Compute residual without u0
      e_i <- y[i] - mu_y[i]
      u0_hat <- sigma_star_sq * e_i / sigma_sq
      u0[i] <- rtruncnorm_upper(1, u0_hat, sigma_star, upper = 0)
    }

    # ==================================================================
    # Step 5: Sample theta (input slacks) - Metropolis step
    # ==================================================================
    # theta_i | ... is complex because theta appears in production function
    # Use Metropolis-Hastings with random walk proposal

    n_mh_theta <- n_mh_theta + 1  # Count iterations, not N*iterations

    # Extract current production function parameters
    alpha_pf <- delta[2:(J + 1)]
    B_vech_curr <- delta[(J + 2):(J + 1 + n_B_vech)]
    B_curr <- vech_to_B(B_vech_curr, J)

    for (i in 1:N) {
      # Current values
      theta_curr <- theta[i, ]
      x_eff_curr <- X_eff[i, ]

      # Mean of theta prior
      mu_theta <- as.vector(Delta %*% Z_slack[i, ])

      # Random walk proposal with adaptive scale
      prop_sd <- sqrt(diag(Omega)) * theta_prop_scale[i] * 0.3
      theta_prop <- theta_curr + rnorm(J, 0, prop_sd)

      # Reject immediately if any theta > 0 (truncation constraint)
      if (any(theta_prop > 0)) {
        next
      }

      # Compute log acceptance ratio
      # Log likelihood contribution
      x_eff_prop <- X[i, ] + theta_prop
      mu_y_curr <- delta[1] + sum(alpha_pf * x_eff_curr) + 0.5 * t(x_eff_curr) %*% B_curr %*% x_eff_curr
      mu_y_prop <- delta[1] + sum(alpha_pf * x_eff_prop) + 0.5 * t(x_eff_prop) %*% B_curr %*% x_eff_prop

      if (M_neutral > 0) {
        delta0 <- delta[(J + 2 + n_B_vech):length(delta)]
        neutral_contrib <- sum(delta0 * Z_neutral[i, ])
        mu_y_curr <- mu_y_curr + neutral_contrib
        mu_y_prop <- mu_y_prop + neutral_contrib
      }

      ll_curr <- dnorm(y[i], mu_y_curr + u0[i], sqrt(sigma_sq), log = TRUE)
      ll_prop <- dnorm(y[i], mu_y_prop + u0[i], sqrt(sigma_sq), log = TRUE)

      # Log prior contribution (MVN, truncation handled by rejection above)
      lp_curr <- mvtnorm::dmvnorm(theta_curr, mu_theta, Omega, log = TRUE)
      lp_prop <- mvtnorm::dmvnorm(theta_prop, mu_theta, Omega, log = TRUE)

      # Symmetric random walk proposal, so proposal ratio = 1 (log = 0)
      log_alpha <- (ll_prop + lp_prop) - (ll_curr + lp_curr)

      if (log(runif(1)) < log_alpha) {
        theta[i, ] <- theta_prop
        accept_theta[i] <- accept_theta[i] + 1
      }
    }

    # Update effective inputs after theta update
    X_eff <- compute_x_eff(X, theta)

    # ==================================================================
    # Step 6: Sample Delta (slack determinant coefficients)
    # ==================================================================
    # theta_i = Delta * z_i + u_i, where u_i ~ N(0, Omega), theta_i <= 0
    # For each row of Delta (each input j), this is a multivariate regression
    # Use Gibbs sampling for each row treating truncation as data augmented

    n_mh_Delta <- n_mh_Delta + 1

    # For the truncated MVN model, we can use a conjugate update if we
    # ignore the truncation correction (approximate). For more accuracy,
    # use Metropolis with small proposal.

    # Proposal: random walk with smaller step size
    Delta_prop <- Delta + matrix(rnorm(J * M, 0, Delta_prop_sd), J, M)

    # Log prior for Delta
    lp_Delta_curr <- -sum(Delta^2) / (2 * prior$Delta_var)
    lp_Delta_prop <- -sum(Delta_prop^2) / (2 * prior$Delta_var)

    # Likelihood contribution from theta (truncated MVN)
    # Use vectorized computation for speed
    ll_Delta_curr <- 0
    ll_Delta_prop <- 0

    Omega_inv <- solve(Omega)
    log_det_Omega <- log(det(Omega))

    for (i in 1:N) {
      mu_curr <- as.vector(Delta %*% Z_slack[i, ])
      mu_prop <- as.vector(Delta_prop %*% Z_slack[i, ])

      # MVN log density (without constant)
      resid_curr <- theta[i, ] - mu_curr
      resid_prop <- theta[i, ] - mu_prop

      ll_Delta_curr <- ll_Delta_curr - 0.5 * (t(resid_curr) %*% Omega_inv %*% resid_curr)
      ll_Delta_prop <- ll_Delta_prop - 0.5 * (t(resid_prop) %*% Omega_inv %*% resid_prop)

      # Truncation correction (compute less frequently for speed)
      # The normalizing constant P(theta <= 0 | mu, Omega) changes with mu
      if (iter %% 10 == 0 || iter <= n_burnin) {
        ll_Delta_curr <- ll_Delta_curr - log_mvnorm_cdf(rep(0, J), mu_curr, Omega)
        ll_Delta_prop <- ll_Delta_prop - log_mvnorm_cdf(rep(0, J), mu_prop, Omega)
      }
    }

    log_alpha_Delta <- (ll_Delta_prop + lp_Delta_prop) - (ll_Delta_curr + lp_Delta_curr)

    if (is.finite(log_alpha_Delta) && log(runif(1)) < log_alpha_Delta) {
      Delta <- Delta_prop
      accept_Delta <- accept_Delta + 1
    }

    # Adapt proposal SD during burn-in
    if (iter <= n_burnin && iter %% 100 == 0) {
      accept_rate_Delta <- accept_Delta / n_mh_Delta
      if (accept_rate_Delta < 0.2) {
        Delta_prop_sd <- Delta_prop_sd * 0.8
      } else if (accept_rate_Delta > 0.5) {
        Delta_prop_sd <- Delta_prop_sd * 1.2
      }
    }

    # ==================================================================
    # Step 7: Sample Omega (slack covariance) - Random walk on log-variances
    # ==================================================================
    n_mh_Omega <- n_mh_Omega + 1

    # Use a simpler approach: random walk on log standard deviations
    # and correlation (for J=2 case)

    # Current parameters
    sigma_omega <- sqrt(diag(Omega))
    if (J == 2) {
      rho_omega <- Omega[1, 2] / (sigma_omega[1] * sigma_omega[2])
    }

    # Proposal on log scale for variances
    log_sigma_curr <- log(sigma_omega)
    log_sigma_prop <- log_sigma_curr + rnorm(J, 0, Omega_prop_scale)
    sigma_prop <- exp(log_sigma_prop)

    # Proposal for correlation (if J=2)
    if (J == 2) {
      # Use Fisher z-transform for correlation
      z_curr <- atanh(rho_omega)
      z_prop <- z_curr + rnorm(1, 0, Omega_prop_scale)
      rho_prop <- tanh(z_prop)  # Ensures -1 < rho < 1

      # Construct proposed Omega
      Omega_prop <- diag(sigma_prop) %*% matrix(c(1, rho_prop, rho_prop, 1), 2, 2) %*% diag(sigma_prop)
    } else {
      # For J > 2, just update diagonal (simplification)
      Omega_prop <- diag(sigma_prop^2)
      # Add small off-diagonal if needed
      for (j1 in 1:(J-1)) {
        for (j2 in (j1+1):J) {
          Omega_prop[j1, j2] <- Omega[j1, j2] * sigma_prop[j1] * sigma_prop[j2] / (sigma_omega[j1] * sigma_omega[j2])
          Omega_prop[j2, j1] <- Omega_prop[j1, j2]
        }
      }
    }

    # Check that Omega_prop is positive definite
    eig <- eigen(Omega_prop, symmetric = TRUE, only.values = TRUE)
    if (all(eig$values > 1e-8)) {
      # Compute acceptance ratio

      # Log likelihood for theta under current and proposed Omega
      ll_Omega_curr <- 0
      ll_Omega_prop <- 0

      Omega_inv_curr <- solve(Omega)
      Omega_inv_prop <- solve(Omega_prop)
      log_det_curr <- log(det(Omega))
      log_det_prop <- log(det(Omega_prop))

      for (i in 1:N) {
        mu_theta <- as.vector(Delta %*% Z_slack[i, ])
        resid <- theta[i, ] - mu_theta

        ll_Omega_curr <- ll_Omega_curr - 0.5 * log_det_curr - 0.5 * t(resid) %*% Omega_inv_curr %*% resid
        ll_Omega_prop <- ll_Omega_prop - 0.5 * log_det_prop - 0.5 * t(resid) %*% Omega_inv_prop %*% resid

        # Truncation correction (subsample for speed)
        if (i %% 20 == 1) {
          ll_Omega_curr <- ll_Omega_curr - log_mvnorm_cdf(rep(0, J), mu_theta, Omega) * min(20, N - i + 1)
          ll_Omega_prop <- ll_Omega_prop - log_mvnorm_cdf(rep(0, J), mu_theta, Omega_prop) * min(20, N - i + 1)
        }
      }

      # Prior: Inverse Wishart (approximately flat on reasonable scale)
      # Use weak prior favoring small variances
      lp_Omega_curr <- -sum(log(sigma_omega)) - sum(sigma_omega^2) / (2 * 0.1^2)
      lp_Omega_prop <- -sum(log(sigma_prop)) - sum(sigma_prop^2) / (2 * 0.1^2)

      # Jacobian for log transform
      log_jacobian <- sum(log_sigma_prop) - sum(log_sigma_curr)

      log_alpha_Omega <- (ll_Omega_prop + lp_Omega_prop) - (ll_Omega_curr + lp_Omega_curr) + log_jacobian

      if (is.finite(log_alpha_Omega) && log(runif(1)) < log_alpha_Omega) {
        Omega <- Omega_prop
        accept_Omega <- accept_Omega + 1
      }
    }

    # Adapt proposal scale during burn-in
    if (iter <= n_burnin && iter %% 100 == 0) {
      accept_rate_Omega <- accept_Omega / n_mh_Omega
      if (accept_rate_Omega < 0.15) {
        Omega_prop_scale <- Omega_prop_scale * 0.8
      } else if (accept_rate_Omega > 0.4) {
        Omega_prop_scale <- Omega_prop_scale * 1.2
      }
    }

    # ==================================================================
    # Store samples
    # ==================================================================
    if (iter > n_burnin && (iter - n_burnin) %% n_thin == 0) {
      sample_idx <- sample_idx + 1

      samples$alpha0[sample_idx] <- delta[1]
      samples$alpha[sample_idx, ] <- delta[2:(J + 1)]
      samples$B_vech[sample_idx, ] <- delta[(J + 2):(J + 1 + n_B_vech)]
      if (M_neutral > 0) {
        samples$delta0[sample_idx, ] <- delta[(J + 2 + n_B_vech):length(delta)]
      }
      samples$sigma[sample_idx] <- sqrt(sigma_sq)
      samples$sigma_u[sample_idx] <- sqrt(sigma_u_sq)
      samples$Delta[sample_idx, , ] <- Delta
      samples$Omega[sample_idx, , ] <- Omega
      samples$theta[sample_idx, , ] <- theta
      samples$u0[sample_idx, ] <- u0
    }

    if (verbose && iter %% 100 == 0) {
      setTxtProgressBar(pb, iter)
    }
  }

  if (verbose) {
    close(pb)
    message("\nGibbs sampler completed.")
    message("Acceptance rates:")
    message("  theta (avg): ", round(mean(accept_theta) / n_mh_theta * 100, 1), "%")
    message("  Delta: ", round(accept_Delta / n_mh_Delta * 100, 1), "%")
    message("  Omega: ", round(accept_Omega / n_mh_Omega * 100, 1), "%")
  }

  # Add acceptance rates (per observation for theta)
  samples$accept_rates <- list(
    theta = accept_theta / n_mh_theta,
    Delta = accept_Delta / n_mh_Delta,
    Omega = accept_Omega / n_mh_Omega
  )

  return(samples)
}

#' Density of Inverse Wishart Distribution
#'
#' @param X The matrix
#' @param nu Degrees of freedom
#' @param S Scale matrix
#' @return Density value
#' @keywords internal
dWishart_inv <- function(X, nu, S) {
  J <- nrow(X)
  X_inv <- solve(X)

  # Log normalizing constant
  log_const <- (nu * J / 2) * log(2) + (J * (J - 1) / 4) * log(pi) +
    sum(lgamma((nu + 1 - 1:J) / 2)) - (nu / 2) * log(det(S))

  # Log kernel
  log_kernel <- -((nu + J + 1) / 2) * log(det(X)) - 0.5 * sum(diag(S %*% X_inv))

  exp(log_kernel - log_const)
}

#' Density of Wishart Distribution
#'
#' @param X The matrix
#' @param nu Degrees of freedom
#' @param S Scale matrix
#' @return Density value
#' @keywords internal
dWishart <- function(X, nu, S) {
  J <- nrow(X)
  S_inv <- solve(S)

  # Log normalizing constant
  log_const <- (nu * J / 2) * log(2) + (J * (J - 1) / 4) * log(pi) +
    sum(lgamma((nu + 1 - 1:J) / 2)) + (nu / 2) * log(det(S))

  # Log kernel
  log_kernel <- ((nu - J - 1) / 2) * log(det(X)) - 0.5 * sum(diag(S_inv %*% X))

  exp(log_kernel - log_const)
}

#' Convert Gibbs Samples to Stan-like Format
#'
#' Converts the output of gsf_gibbs_sampler to a format compatible with
#' the existing extraction functions.
#'
#' @param samples Output from gsf_gibbs_sampler
#' @param N Number of observations
#' @param J Number of inputs
#' @return A list mimicking rstan extract() output
#' @keywords internal
gibbs_to_stan_format <- function(samples, N, J) {
  n_samples <- length(samples$sigma)
  n_B_vech <- J * (J + 1) / 2

  # Construct B matrix for each sample
  B <- array(0, c(n_samples, J, J))
  for (s in 1:n_samples) {
    idx <- 1
    for (i in 1:J) {
      for (j in i:J) {
        B[s, i, j] <- samples$B_vech[s, idx]
        B[s, j, i] <- samples$B_vech[s, idx]
        idx <- idx + 1
      }
    }
  }

  # Compute efficiency measures
  technical_inefficiency <- -samples$u0 * 100  # Convert to percentage
  input_slacks <- -samples$theta * 100

  # Compute elasticities and returns to scale
  elasticities <- array(0, c(n_samples, N, J))
  returns_to_scale <- matrix(0, n_samples, N)
  output_loss_from_slacks <- matrix(0, n_samples, N)

  for (s in 1:n_samples) {
    alpha <- samples$alpha[s, ]
    B_s <- B[s, , ]

    for (i in 1:N) {
      x_eff <- samples$theta[s, i, ]  # This should be X + theta but we need X
      # For now, skip detailed elasticity computation
      # This will be handled properly when we have access to X
    }
  }

  list(
    alpha0 = samples$alpha0,
    alpha = samples$alpha,
    B_vech = samples$B_vech,
    B = B,
    sigma = samples$sigma,
    sigma_u = samples$sigma_u,
    sigma_slack = apply(samples$Omega, 1, function(O) sqrt(diag(O))),
    Delta_raw = samples$Delta,
    L_Omega = apply(samples$Omega, 1, function(O) t(chol(O))),
    theta = samples$theta,
    u0 = samples$u0,
    technical_inefficiency = technical_inefficiency,
    input_slacks = input_slacks
  )
}
