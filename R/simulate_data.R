#' Simulate Data from the GSF Model
#'
#' @description
#' Simulates data from the Generalized Stochastic Frontier model for testing
#' and demonstration purposes. This function generates data similar to the
#' UK manufacturing data used in Kumbhakar and Tsionas (2021).
#'
#' @param n_firms Number of firms.
#' @param n_periods Number of time periods per firm (for panel data).
#' @param n_inputs Number of inputs (default: 2 for labor and capital).
#' @param n_z Number of z variables (slack determinants).
#' @param alpha0 Intercept of production function.
#' @param alpha Vector of first-order coefficients.
#' @param B Symmetric matrix of second-order coefficients.
#' @param delta_mean Mean of Delta coefficients for slack determinants.
#' @param sigma Noise standard deviation.
#' @param sigma_u Output inefficiency scale.
#' @param sigma_slack Scale of input slacks.
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
                               delta_mean = -0.3,
                               sigma = 0.1,
                               sigma_u = 0.05,
                               sigma_slack = 0.1,
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

  # Default B matrix (symmetric, small second-order effects)
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

  # Generate Delta matrix (J x M, where M = n_z + 1 for intercept)
  Delta <- matrix(rnorm(n_inputs * (n_z + 1), delta_mean, 0.3),
                  n_inputs, n_z + 1)

  # Generate input slacks theta (N x J)
  # theta_i = Delta * z_i + u_i, theta_i <= 0
  theta <- matrix(0, N, n_inputs)
  Z_with_intercept <- cbind(1, Z)

  for (i in 1:N) {
    mu_theta <- Delta %*% Z_with_intercept[i, ]
    # Add noise
    theta_i <- mu_theta + rnorm(n_inputs, 0, sigma_slack)
    # Truncate at 0 from above
    theta[i, ] <- pmin(theta_i, -0.001)
  }

  # Generate output technical inefficiency u0 (half-normal, u0 <= 0)
  u0 <- -abs(rnorm(N, 0, sigma_u))

  # Generate noise
  epsilon <- rnorm(N, 0, sigma)

  # Compute output using translog production function
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

  # Store true parameters as attributes
  attr(data, "true_params") <- list(
    alpha0 = alpha0,
    alpha = alpha,
    B = B,
    Delta = Delta,
    sigma = sigma,
    sigma_u = sigma_u,
    sigma_slack = sigma_slack
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


#' Simulate UK Manufacturing-like Data
#'
#' @description
#' Simulates data that mimics the structure of the UK manufacturing data
#' used in Kumbhakar and Tsionas (2021).
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
      # Increase some firms' periods
      idx <- sample(1:n_firms, 10)
      n_periods_per_firm[idx] <- pmin(n_periods_per_firm[idx] + 1, 13)
    } else {
      # Decrease some firms' periods
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
  L <- numeric(total_N)  # Number of employees
  K <- numeric(total_N)  # Capital

  for (i in 1:n_firms) {
    firm_idx <- which(firm_id == i)
    n_t <- length(firm_idx)

    # Firm-specific baseline
    L_base <- exp(rnorm(1, 5, 1))  # 50-1000 employees
    K_base <- exp(rnorm(1, 8, 1.5))  # Capital

    # Time series with growth
    L[firm_idx] <- L_base * exp(cumsum(rnorm(n_t, 0.01, 0.05)))
    K[firm_idx] <- K_base * exp(cumsum(rnorm(n_t, 0.02, 0.08)))
  }

  # Z variables
  CR <- numeric(total_N)    # Concentration ratio
  MKSH <- numeric(total_N)  # Market share
  FP <- numeric(total_N)    # Financial pressure
  IMP <- numeric(total_N)   # Import penetration
  TREND <- time_id / 13     # Normalized time trend

  for (i in 1:n_firms) {
    firm_idx <- which(firm_id == i)
    n_t <- length(firm_idx)

    # Concentration ratio (industry-level, varies slowly)
    CR[firm_idx] <- runif(1, 0.3, 0.7) + rnorm(n_t, 0, 0.02)

    # Market share
    MKSH[firm_idx] <- runif(1, 0.01, 0.1) * (1 + rnorm(n_t, 0, 0.1))

    # Financial pressure
    FP[firm_idx] <- rbeta(n_t, 2, 8)

    # Import penetration
    IMP[firm_idx] <- runif(1, 0.15, 0.35) + cumsum(rnorm(n_t, 0.005, 0.01))
  }

  # Clip to valid ranges
  CR <- pmax(pmin(CR, 1), 0)
  MKSH <- pmax(pmin(MKSH, 1), 0)
  FP <- pmax(pmin(FP, 1), 0)
  IMP <- pmax(pmin(IMP, 1), 0)

  # True parameters (similar to paper's estimates)
  alpha0 <- 2.5
  alpha <- c(0.52, 0.42)  # Labor and capital elasticities
  B <- matrix(c(-0.08, 0.03, 0.03, -0.05), 2, 2)

  # Delta coefficients (Table 3 in paper)
  # Note: positive coefficients mean less slack (more efficiency)
  Delta <- matrix(c(
    -0.8, -1.2, 1.1, 1.2, 1.4, -0.02,  # Labor
    0.2, -0.2, -0.3, 0.15, 0.24, -0.01  # Capital
  ), 2, 6, byrow = TRUE)

  sigma <- 0.08
  sigma_u <- 0.03
  sigma_slack <- c(0.04, 0.03)  # Different for each input

  # Generate slacks
  theta <- matrix(0, total_N, 2)
  Z_with_intercept <- cbind(1, CR, MKSH, FP, IMP, TREND)

  for (i in 1:total_N) {
    mu_theta <- Delta %*% Z_with_intercept[i, ]
    theta_i <- mu_theta + rnorm(2, 0, sigma_slack)
    theta[i, ] <- pmin(theta_i, -0.001)
  }

  # Generate output inefficiency
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
    Y = exp(y),  # Sales
    L = L,       # Labor
    K = K,       # Capital
    CR = CR,
    MKSH = MKSH,
    FP = FP,
    IMP = IMP,
    TREND = TREND
  )

  # Store true parameters
  attr(data, "true_params") <- list(
    alpha0 = alpha0,
    alpha = alpha,
    B = B,
    Delta = Delta,
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
