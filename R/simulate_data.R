#' Simulate Data from the GSF Model
#'
#' @description
#' Simulates data from the Generalized Stochastic Frontier model for testing
#' and demonstration purposes. This function generates data similar to the
#' UK manufacturing data used in Kumbhakar and Tsionas (2021).
#'
#' The key innovation is proper sampling from the truncated multivariate normal
#' distribution for input slacks theta, following equation (3b) in the paper:
#' theta_i = Delta * z_i + u_i, where u_i ~ N_J(0, Omega) and theta_i <= 0
#'
#' @param n_firms Number of firms.
#' @param n_periods Number of time periods per firm (for panel data).
#' @param n_inputs Number of inputs (default: 2 for labor and capital).
#' @param n_z Number of z variables (slack determinants).
#' @param alpha0 Intercept of production function.
#' @param alpha Vector of first-order coefficients.
#' @param B Symmetric matrix of second-order coefficients.
#' @param Delta Matrix of slack determinant coefficients (J x (n_z + 1)).
#' @param Omega Covariance matrix for input slacks (J x J).
#' @param sigma Noise standard deviation.
#' @param sigma_u Output inefficiency scale.
#' @param seed Random seed for reproducibility.
#'
#' @return A data frame with simulated data.
#'
#' @examples
#' # Simulate data similar to UK manufacturing
#' sim_data <- simulate_gsf_data(
#'   n_firms = 100,
#'   n_periods = 10,
#'   n_inputs = 2,
#'   n_z = 3,
#'   seed = 123
#' )
#'
#' # Check structure
#' str(sim_data)
#'
#' @export
simulate_gsf_data <- function(n_firms = 100,
                               n_periods = 10,
                               n_inputs = 2,
                               n_z = 3,
                               alpha0 = 2.0,
                               alpha = NULL,
                               B = NULL,
                               Delta = NULL,
                               Omega = NULL,
                               sigma = 0.1,
                               sigma_u = 0.05,
                               seed = NULL) {

  if (!is.null(seed)) {
    set.seed(seed)
  }

  N <- n_firms * n_periods

  # Default alpha (labor and capital elasticities)
  if (is.null(alpha)) {
    if (n_inputs == 2) {
      alpha <- c(0.6, 0.3)  # Labor, Capital
    } else {
      alpha <- rep(1/n_inputs, n_inputs)
    }
  }

  # Default B matrix (symmetric, with curvature for identification)
  if (is.null(B)) {
    B <- matrix(0, n_inputs, n_inputs)
    diag(B) <- -0.1  # Decreasing marginal products
    if (n_inputs > 1) {
      # Small cross-effects
      for (i in 1:(n_inputs-1)) {
        for (j in (i+1):n_inputs) {
          B[i, j] <- 0.05
          B[j, i] <- 0.05
        }
      }
    }
  }

  # Default Delta matrix (J x (n_z + 1), includes intercept)
  # Design so that mean theta is moderately negative
  if (is.null(Delta)) {
    Delta <- matrix(0, n_inputs, n_z + 1)
    # Intercept: ensures mean slack is negative even when z = 0
    Delta[, 1] <- -0.15
    # Z coefficients: moderate effects
    for (j in 1:n_z) {
      Delta[, j + 1] <- runif(n_inputs, -0.3, 0.3)
    }
  }

  # Default Omega (covariance matrix for slack errors)
  # Per paper, slacks can be correlated across inputs

  if (is.null(Omega)) {
    sigma_slack <- rep(0.08, n_inputs)  # SD for each input slack
    rho <- 0.3  # Moderate positive correlation
    # Construct correlation matrix
    R <- diag(n_inputs)
    if (n_inputs > 1) {
      for (i in 1:(n_inputs-1)) {
        for (j in (i+1):n_inputs) {
          R[i, j] <- rho
          R[j, i] <- rho
        }
      }
    }
    Omega <- diag(sigma_slack) %*% R %*% diag(sigma_slack)
  }

  # Generate firm and time identifiers
  firm_id <- rep(1:n_firms, each = n_periods)
  time_id <- rep(1:n_periods, times = n_firms)

  # Generate log inputs
  X <- matrix(0, N, n_inputs)
  colnames(X) <- paste0("X", 1:n_inputs)
  if (n_inputs == 2) {
    colnames(X) <- c("L", "K")  # Labor, Capital
  }

  # Simulate inputs with firm-specific heterogeneity
  for (i in 1:n_firms) {
    idx <- ((i-1)*n_periods + 1):(i*n_periods)
    # Firm-specific means
    firm_mean <- rnorm(n_inputs, mean = c(3, 4)[1:n_inputs], sd = 0.5)
    # Time variation
    for (j in 1:n_inputs) {
      X[idx, j] <- firm_mean[j] + cumsum(rnorm(n_periods, 0, 0.1))
    }
  }

  # Generate Z variables (slack determinants)
  Z <- matrix(0, N, n_z)
  z_names <- c("CR", "MKSH", "FP", "IMP", "TREND")[1:n_z]
  if (n_z > 5) {
    z_names <- c(z_names, paste0("Z", 6:n_z))
  }
  colnames(Z) <- z_names

  # Simulate z variables similar to UK manufacturing data
  for (col_idx in 1:n_z) {
    if (z_names[col_idx] == "CR") {
      # Concentration ratio (0-1)
      for (i in 1:n_firms) {
        idx <- ((i-1)*n_periods + 1):(i*n_periods)
        Z[idx, col_idx] <- runif(1, 0.3, 0.8) + rnorm(n_periods, 0, 0.05)
      }
    } else if (z_names[col_idx] == "MKSH") {
      # Market share (0-1, typically small)
      for (i in 1:n_firms) {
        idx <- ((i-1)*n_periods + 1):(i*n_periods)
        Z[idx, col_idx] <- runif(1, 0.01, 0.15) + rnorm(n_periods, 0, 0.01)
      }
    } else if (z_names[col_idx] == "FP") {
      # Financial pressure (0-1)
      for (i in 1:n_firms) {
        idx <- ((i-1)*n_periods + 1):(i*n_periods)
        Z[idx, col_idx] <- rbeta(n_periods, 2, 5)
      }
    } else if (z_names[col_idx] == "IMP") {
      # Import penetration (0-1)
      for (i in 1:n_firms) {
        idx <- ((i-1)*n_periods + 1):(i*n_periods)
        Z[idx, col_idx] <- runif(1, 0.1, 0.4) + rnorm(n_periods, 0, 0.02)
      }
    } else if (z_names[col_idx] == "TREND") {
      # Time trend
      Z[, col_idx] <- time_id / n_periods
    } else {
      # Generic z variable
      Z[, col_idx] <- rnorm(N, 0, 1)
    }
  }

  # Clip Z to reasonable ranges
  Z <- pmax(pmin(Z, 1), 0)

  # Generate input slacks theta (N x J) using PROPER truncated MVN
  # theta_i ~ N_J(Delta * z_i, Omega) truncated at theta <= 0
  Z_with_intercept <- cbind(1, Z)
  theta <- matrix(0, N, n_inputs)

  for (i in 1:N) {
    mu_theta <- as.vector(Delta %*% Z_with_intercept[i, ])
    # Sample from truncated multivariate normal
    theta[i, ] <- sample_truncated_mvn(mu_theta, Omega, upper = rep(0, n_inputs))
  }

  # Generate output technical inefficiency u0 (half-normal, u0 <= 0)
  u0 <- -abs(rnorm(N, 0, sigma_u))

  # Generate noise
  epsilon <- rnorm(N, 0, sigma)

  # Compute output using translog production function (Equation 2 in paper)
  # y_i = alpha0 + alpha'(x_i + theta_i) + 0.5*(x_i + theta_i)'B(x_i + theta_i) + u0_i + epsilon_i
  y <- numeric(N)
  for (i in 1:N) {
    x_eff <- X[i, ] + theta[i, ]
    y[i] <- alpha0 + sum(alpha * x_eff) + 0.5 * t(x_eff) %*% B %*% x_eff +
            u0[i] + epsilon[i]
  }

  # Create output data frame
  data <- data.frame(
    firm_id = firm_id,
    time_id = time_id,
    Y = exp(y),  # Convert to levels
    stringsAsFactors = FALSE
  )

  # Add inputs (in levels)
  for (j in 1:n_inputs) {
    data[[colnames(X)[j]]] <- exp(X[, j])
  }

  # Add Z variables
  for (j in 1:n_z) {
    data[[colnames(Z)[j]]] <- Z[, j]
  }

  # Store true parameters as attributes (including Omega for recovery testing)
  attr(data, "true_params") <- list(
    alpha0 = alpha0,
    alpha = alpha,
    B = B,
    Delta = Delta,
    Omega = Omega,
    sigma = sigma,
    sigma_u = sigma_u,
    sigma_slack = sqrt(diag(Omega))  # For backward compatibility
  )

  # Store true efficiency measures
  attr(data, "true_efficiency") <- list(
    theta = theta,
    u0 = u0,
    technical_inefficiency = -u0 * 100,
    input_slacks = -theta * 100
  )

  return(data)
}


#' Sample from Truncated Multivariate Normal Distribution
#'
#' @description
#' Samples from N_J(mu, Sigma) truncated to x <= upper using Gibbs sampling.
#' Based on Geweke (1991) "Efficient simulation from the multivariate normal
#' and student-t distributions subject to linear constraints."
#'
#' This is essential for proper simulation of input slacks in the GSF model,
#' where theta ~ N_J(Delta * z, Omega) truncated at theta <= 0.
#'
#' @param mu Mean vector (J x 1)
#' @param Sigma Covariance matrix (J x J)
#' @param upper Upper bounds (J x 1), default all zeros
#' @param n_gibbs Number of Gibbs iterations for sampling (default 25)
#' @return A vector of length J sampled from the truncated distribution
#'
#' @keywords internal
sample_truncated_mvn <- function(mu, Sigma, upper = rep(0, length(mu)),
                                  n_gibbs = 25) {
  J <- length(mu)

  # For J = 1, use simple univariate truncated normal

  if (J == 1) {
    return(sample_truncated_normal(mu, sqrt(Sigma[1,1]), upper[1]))
  }

  # Initialize at a feasible point
  # Start at conditional means truncated to satisfy constraint
  x <- pmin(mu, upper - abs(mu) * 0.1 - 0.01)

  # Ensure we start in the feasible region
  if (any(x > upper)) {
    x <- upper - 0.1
  }

  # Precompute for conditional distributions
  # For Gibbs sampling, we cycle through each component

  # Gibbs sampling
  for (iter in 1:n_gibbs) {
    for (j in 1:J) {
      # Conditional distribution of x_j | x_{-j}
      # x_j | x_{-j} ~ N(mu_cond, sigma_cond^2) truncated at upper[j]

      if (J == 2) {
        # Explicit formulas for bivariate case (more efficient)
        k <- 3 - j  # The other index
        sigma_j <- sqrt(Sigma[j, j])
        sigma_k <- sqrt(Sigma[k, k])
        rho <- Sigma[j, k] / (sigma_j * sigma_k)

        cond_mean <- mu[j] + rho * (sigma_j / sigma_k) * (x[k] - mu[k])
        cond_sd <- sigma_j * sqrt(1 - rho^2)
      } else {
        # General case for J > 2
        idx_j <- j
        idx_notj <- setdiff(1:J, j)

        Sigma_j_notj <- Sigma[idx_j, idx_notj, drop = FALSE]
        Sigma_notj_notj <- Sigma[idx_notj, idx_notj, drop = FALSE]
        Sigma_notj_notj_inv <- solve(Sigma_notj_notj)

        cond_mean <- mu[j] + Sigma_j_notj %*% Sigma_notj_notj_inv %*%
          (x[idx_notj] - mu[idx_notj])
        cond_var <- Sigma[j, j] - Sigma_j_notj %*% Sigma_notj_notj_inv %*%
          t(Sigma_j_notj)
        cond_sd <- sqrt(as.numeric(cond_var))
        cond_mean <- as.numeric(cond_mean)
      }

      # Sample from truncated univariate normal
      x[j] <- sample_truncated_normal(cond_mean, cond_sd, upper[j])
    }
  }

  return(x)
}


#' Sample from Truncated Univariate Normal Distribution
#'
#' @description
#' Samples from N(mu, sigma^2) truncated to (-Inf, upper] using inverse CDF method.
#'
#' @param mu Mean of the untruncated normal
#' @param sigma Standard deviation
#' @param upper Upper bound
#' @return A single sample from the truncated distribution
#'
#' @keywords internal
sample_truncated_normal <- function(mu, sigma, upper) {
  # P(X <= upper) for X ~ N(mu, sigma^2)
  p_upper <- pnorm((upper - mu) / sigma)

  # Handle edge case where truncation probability is very small
  if (p_upper < 1e-10) {
    # Almost all mass is above upper, return upper - small value
    return(upper - abs(sigma) * 0.01)
  }

  # Sample u ~ Uniform(0, p_upper) and invert
  u <- runif(1, 0, p_upper)
  return(qnorm(u) * sigma + mu)
}


#' Simulate UK Manufacturing-like Data
#'
#' @description
#' Simulates data that mimics the structure of the UK manufacturing data
#' used in Kumbhakar and Tsionas (2021). Uses proper truncated multivariate
#' normal sampling for input slacks.
#'
#' @param n_firms Number of firms (default: 582 as in the paper).
#' @param avg_periods Average number of periods per firm (default: 9).
#' @param seed Random seed.
#'
#' @return A data frame similar to the UK manufacturing data.
#'
#' @export
simulate_uk_manufacturing <- function(n_firms = 582,
                                       avg_periods = 9,
                                       seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Create unbalanced panel
  n_periods_per_firm <- sample(5:13, n_firms, replace = TRUE)

  # Adjust to get approximately the target average
  while (abs(mean(n_periods_per_firm) - avg_periods) > 0.5) {
    if (mean(n_periods_per_firm) < avg_periods) {
      idx <- sample(1:n_firms, 10)
      n_periods_per_firm[idx] <- pmin(n_periods_per_firm[idx] + 1, 13)
    } else {
      idx <- sample(1:n_firms, 10)
      n_periods_per_firm[idx] <- pmax(n_periods_per_firm[idx] - 1, 5)
    }
  }

  total_N <- sum(n_periods_per_firm)

  # Initialize data
  firm_id <- integer(total_N)
  time_id <- integer(total_N)
  idx <- 1

  for (i in 1:n_firms) {
    n_t <- n_periods_per_firm[i]
    firm_id[idx:(idx + n_t - 1)] <- i
    time_id[idx:(idx + n_t - 1)] <- 1:n_t
    idx <- idx + n_t
  }

  # Generate inputs (Labor and Capital)
  L <- numeric(total_N)
  K <- numeric(total_N)

  for (i in 1:n_firms) {
    firm_idx <- which(firm_id == i)
    n_t <- length(firm_idx)

    L_base <- exp(rnorm(1, 5, 1))
    K_base <- exp(rnorm(1, 8, 1.5))

    L[firm_idx] <- L_base * exp(cumsum(rnorm(n_t, 0.01, 0.05)))
    K[firm_idx] <- K_base * exp(cumsum(rnorm(n_t, 0.02, 0.08)))
  }

  # Z variables
  CR <- numeric(total_N)
  MKSH <- numeric(total_N)
  FP <- numeric(total_N)
  IMP <- numeric(total_N)
  TREND <- time_id / 13

  for (i in 1:n_firms) {
    firm_idx <- which(firm_id == i)
    n_t <- length(firm_idx)

    CR[firm_idx] <- runif(1, 0.3, 0.7) + rnorm(n_t, 0, 0.02)
    MKSH[firm_idx] <- runif(1, 0.01, 0.1) * (1 + rnorm(n_t, 0, 0.1))
    FP[firm_idx] <- rbeta(n_t, 2, 8)
    IMP[firm_idx] <- runif(1, 0.15, 0.35) + cumsum(rnorm(n_t, 0.005, 0.01))
  }

  CR <- pmax(pmin(CR, 1), 0)
  MKSH <- pmax(pmin(MKSH, 1), 0)
  FP <- pmax(pmin(FP, 1), 0)
  IMP <- pmax(pmin(IMP, 1), 0)

  # True parameters (similar to paper's estimates, Table 2 and 3)
  alpha0 <- 2.5
  alpha <- c(0.52, 0.42)
  B <- matrix(c(-0.08, 0.03, 0.03, -0.05), 2, 2)

  # Delta coefficients (Table 3 in paper)
  # Columns: intercept, CR, MKSH, FP, IMP, TREND
  Delta <- matrix(c(
    -0.10, -0.15, 0.12, 0.10, 0.08, -0.02,  # Labor
    -0.12, -0.08, -0.05, 0.06, 0.05, -0.01  # Capital
  ), 2, 6, byrow = TRUE)

  sigma <- 0.08
  sigma_u <- 0.025  # ~2.5% mean technical inefficiency

  # Covariance matrix for slacks (with correlation)
  sigma_slack <- c(0.04, 0.05)  # SDs for labor and capital slacks
  rho_slack <- 0.4  # Positive correlation between slacks
  Omega <- matrix(c(
    sigma_slack[1]^2, rho_slack * sigma_slack[1] * sigma_slack[2],
    rho_slack * sigma_slack[1] * sigma_slack[2], sigma_slack[2]^2
  ), 2, 2)

  # Generate slacks using proper truncated MVN
  theta <- matrix(0, total_N, 2)
  Z_with_intercept <- cbind(1, CR, MKSH, FP, IMP, TREND)

  for (i in 1:total_N) {
    mu_theta <- as.vector(Delta %*% Z_with_intercept[i, ])
    theta[i, ] <- sample_truncated_mvn(mu_theta, Omega, upper = c(0, 0))
  }

  # Generate output inefficiency (half-normal)
  u0 <- -abs(rnorm(total_N, 0, sigma_u))

  # Generate noise
  epsilon <- rnorm(total_N, 0, sigma)

  # Compute output (Sales in logs)
  X <- cbind(log(L), log(K))
  y <- numeric(total_N)

  for (i in 1:total_N) {
    x_eff <- X[i, ] + theta[i, ]
    y[i] <- alpha0 + sum(alpha * x_eff) + 0.5 * t(x_eff) %*% B %*% x_eff +
            u0[i] + epsilon[i]
  }

  # Create output data frame
  data <- data.frame(
    firm_id = firm_id,
    time_id = time_id,
    Y = exp(y),
    L = L,
    K = K,
    CR = CR,
    MKSH = MKSH,
    FP = FP,
    IMP = IMP,
    TREND = TREND
  )

  # Store true parameters (including Omega)
  attr(data, "true_params") <- list(
    alpha0 = alpha0,
    alpha = alpha,
    B = B,
    Delta = Delta,
    Omega = Omega,
    sigma = sigma,
    sigma_u = sigma_u,
    sigma_slack = sigma_slack
  )

  attr(data, "true_efficiency") <- list(
    theta = theta,
    u0 = u0,
    technical_inefficiency = -u0 * 100,
    input_slacks = -theta * 100
  )

  message("Generated ", total_N, " observations from ", n_firms, " firms")
  message("Average periods per firm: ", round(mean(n_periods_per_firm), 2))

  return(data)
}


#' Simulate Data with Stronger Identification Signals
#'
#' @description
#' Generates a synthetic data set designed for parameter recovery testing.
#' Identification between output inefficiency (u0) and input slacks (theta)
#' is strengthened by:
#' \itemize{
#'   \item Stronger curvature in the production function (larger B matrix effects)
#'   \item High-variance, independent Z drivers that strongly predict slacks
#'   \item Proper truncated multivariate normal sampling for theta
#'   \item Smaller output-inefficiency scale relative to slack scale
#' }
#'
#' @param n_firms Number of firms (default: 100).
#' @param n_periods Number of periods per firm (default: 5).
#' @param seed Random seed for reproducibility.
#'
#' @return A data frame with levels of output and inputs plus the Z variables.
#'   True parameters and efficiencies are stored as attributes `true_params`
#'   and `true_efficiency`.
#'
#' @details
#' The key identification insight from Kumbhakar & Tsionas (2021) is that
#' theta can be separated from u0 because:
#' 1. theta appears both linearly (alpha'*theta) and quadratically (theta'*B*theta)
#' 2. theta interacts with x through the term theta'*B*x
#' 3. The Z variables predict theta but not u0 directly
#'
#' This function designs the DGP to maximize these identification signals.
#'
#' @examples
#' sim <- simulate_gsf_data_strong_id(n_firms = 100, n_periods = 5, seed = 123)
#' str(sim)
#' attr(sim, "true_params")
#'
#' @export
simulate_gsf_data_strong_id <- function(n_firms = 100,
                                         n_periods = 5,
                                         seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  N <- n_firms * n_periods
  J <- 2
  n_z <- 2

  # Production function parameters with STRONG curvature for identification
  # The B matrix creates the interaction theta'*B*x that helps identify theta
  alpha0 <- 2.0
  alpha <- c(0.55, 0.40)  # Reasonable elasticities summing to ~0.95
  B <- matrix(c(-0.15, 0.06,
                 0.06, -0.12), 2, 2, byrow = TRUE)

 # Delta: slack determinants with STRONG, DISTINCT effects
 # Column order: intercept, z1, z2
 # Design: z1 strongly affects labor slack, z2 strongly affects capital slack
 # This differential response helps identify input-specific slacks
  Delta <- matrix(
    c(-0.08,  0.12, -0.03,   # Labor: intercept, strong z1 effect, weak z2
      -0.10, -0.02,  0.15),  # Capital: intercept, weak z1, strong z2 effect
    nrow = J,
    byrow = TRUE
  )

  # Variance parameters designed for identification
  # Key: sigma_u should be SMALLER than slack variation so u0 doesn't dominate
  sigma <- 0.06       # Noise SD
  sigma_u <- 0.02     # Output inefficiency SD (small: ~2% mean inefficiency)

  # Covariance matrix for slacks - moderate correlation
  sigma_slack <- c(0.06, 0.07)  # ~6-7% SD for slacks
  rho_slack <- 0.35
  Omega <- matrix(c(
    sigma_slack[1]^2, rho_slack * sigma_slack[1] * sigma_slack[2],
    rho_slack * sigma_slack[1] * sigma_slack[2], sigma_slack[2]^2
  ), 2, 2)

  # Firm/time identifiers
  firm_id <- rep(1:n_firms, each = n_periods)
  time_id <- rep(1:n_periods, times = n_firms)

  # Inputs in logs with firm heterogeneity
  X <- matrix(0, N, J)
  colnames(X) <- c("X1", "X2")
  for (i in 1:n_firms) {
    idx <- ((i - 1) * n_periods + 1):(i * n_periods)
    base <- rnorm(J, mean = c(4.0, 4.5), sd = 0.4)
    drift <- matrix(rnorm(n_periods * J, 0, 0.06), n_periods, J)
    X[idx, ] <- matrix(base, n_periods, J, byrow = TRUE) + apply(drift, 2, cumsum)
  }

  # Z drivers: HIGH VARIANCE and INDEPENDENT for strong identification
  # z1: varies at observation level (time-varying)
  # z2: varies at firm level (firm-specific) with small time noise
  z1 <- rnorm(N, 0, 1.0)  # Strong variation
  z2 <- rep(rnorm(n_firms, 0, 0.8), each = n_periods) + rnorm(N, 0, 0.15)
  Z <- cbind(z1, z2)
  colnames(Z) <- c("z1", "z2")

  # Generate slacks using PROPER truncated MVN
  Z_with_intercept <- cbind(1, Z)
  theta <- matrix(0, N, J)

  for (i in 1:N) {
    mu_theta <- as.vector(Delta %*% Z_with_intercept[i, ])
    theta[i, ] <- sample_truncated_mvn(mu_theta, Omega, upper = c(0, 0))
  }

  # Output inefficiency (half-normal, small relative to slacks)
  u0 <- -abs(rnorm(N, 0, sigma_u))

  # Noise
  epsilon <- rnorm(N, 0, sigma)

  # Compute output using translog (Equation 2)
  y_log <- numeric(N)
  for (i in 1:N) {
    x_eff <- X[i, ] + theta[i, ]
    y_log[i] <- alpha0 +
      sum(alpha * x_eff) +
      0.5 * t(x_eff) %*% B %*% x_eff +
      u0[i] + epsilon[i]
  }

  data <- data.frame(
    firm_id = firm_id,
    time_id = time_id,
    Y = exp(y_log),
    X1 = exp(X[, 1]),
    X2 = exp(X[, 2]),
    z1 = z1,
    z2 = z2
  )

  # Store ALL true parameters for recovery testing
  attr(data, "true_params") <- list(
    alpha0 = alpha0,
    alpha = alpha,
    B = B,
    Delta = Delta,
    Omega = Omega,
    sigma = sigma,
    sigma_u = sigma_u,
    sigma_slack = sigma_slack,
    rho_slack = rho_slack
  )

  attr(data, "true_efficiency") <- list(
    theta = theta,
    u0 = u0,
    technical_inefficiency = -u0 * 100,
    input_slacks = -theta * 100,
    # Also store log inputs for computing elasticities
    X_log = X
  )

  message("Generated ", N, " observations from ", n_firms, " firms")
  message("Mean |theta|: ", round(mean(abs(theta)) * 100, 2), "% (input slack)")
  message("Mean |u0|: ", round(mean(abs(u0)) * 100, 2), "% (output inefficiency)")

  return(data)
}


#' Validate Simulated Data for Parameter Recovery
#'
#' @description
#' Checks that simulated data has properties conducive to parameter recovery.
#' Reports diagnostics about the data generating process.
#'
#' @param data Simulated data from simulate_gsf_data or related functions
#' @return Invisibly returns a list of diagnostics
#'
#' @export
validate_simulation <- function(data) {
  true_params <- attr(data, "true_params")
  true_eff <- attr(data, "true_efficiency")

  if (is.null(true_params) || is.null(true_eff)) {
    stop("Data does not have true_params or true_efficiency attributes")
  }

  cat("=== Simulation Diagnostics ===\n\n")

  cat("Sample size:", nrow(data), "\n\n")

  cat("Production function parameters:\n")
  cat("  alpha0:", true_params$alpha0, "\n")
  cat("  alpha:", paste(round(true_params$alpha, 3), collapse = ", "), "\n")
  cat("  B diagonal:", paste(round(diag(true_params$B), 3), collapse = ", "), "\n")
  cat("  Returns to scale:", round(sum(true_params$alpha), 3), "\n\n")

  cat("Variance parameters:\n")
  cat("  sigma (noise):", round(true_params$sigma, 4), "\n")
  cat("  sigma_u (output inefficiency):", round(true_params$sigma_u, 4), "\n")

  if (!is.null(true_params$Omega)) {
    cat("  Omega (slack covariance):\n")
    print(round(true_params$Omega, 6))
    cat("  Slack correlation:", round(true_params$Omega[1,2] /
        sqrt(true_params$Omega[1,1] * true_params$Omega[2,2]), 3), "\n")
  }
  cat("\n")

  cat("Efficiency measures (true values):\n")
  cat("  Mean |u0|:", round(mean(abs(true_eff$u0)) * 100, 2), "%\n")
  cat("  SD |u0|:", round(sd(abs(true_eff$u0)) * 100, 2), "%\n")

  theta <- true_eff$theta
  cat("  Mean |theta_1| (input 1 slack):", round(mean(abs(theta[,1])) * 100, 2), "%\n")
  cat("  Mean |theta_2| (input 2 slack):", round(mean(abs(theta[,2])) * 100, 2), "%\n")
  cat("  Correlation(theta_1, theta_2):", round(cor(theta[,1], theta[,2]), 3), "\n\n")

  # Check identification conditions
  cat("Identification checks:\n")

  # 1. Slack variation vs u0 variation
  slack_var <- var(rowSums(theta))
  u0_var <- var(true_eff$u0)
  cat("  Var(sum theta) / Var(u0):", round(slack_var / u0_var, 2),
      "(should be > 1 for good identification)\n")

  # 2. B matrix effects
  B_norm <- sqrt(sum(true_params$B^2))
  cat("  ||B|| (Frobenius norm):", round(B_norm, 3),
      "(larger helps identification)\n")

  # 3. Z variation
  if ("z1" %in% names(data)) {
    cat("  SD(z1):", round(sd(data$z1), 3), "\n")
    cat("  SD(z2):", round(sd(data$z2), 3), "\n")

    # Check if Z predicts theta
    z_mat <- cbind(1, data$z1, data$z2)
    r2_theta1 <- summary(lm(theta[,1] ~ data$z1 + data$z2))$r.squared
    r2_theta2 <- summary(lm(theta[,2] ~ data$z1 + data$z2))$r.squared
    cat("  R^2 of theta_1 ~ Z:", round(r2_theta1, 3), "\n")
    cat("  R^2 of theta_2 ~ Z:", round(r2_theta2, 3),
        "(higher = stronger identification)\n")
  }

  cat("\n")

  invisible(list(
    slack_var = slack_var,
    u0_var = u0_var,
    B_norm = B_norm
  ))
}
