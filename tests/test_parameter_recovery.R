#!/usr/bin/env Rscript
#' Parameter Recovery Test for GSF Model
#'
#' This script tests whether the GSF model can accurately recover
#' the true parameters from simulated data.
#'
#' Based on Kumbhakar & Tsionas (2021) methodology.

library(rstan)

# Set options for Stan
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores() - 1)

# Set seed for reproducibility
set.seed(42)

cat("=======================================================\n")
cat("GSF Model Parameter Recovery Test\n")
cat("=======================================================\n\n")

#' Simulate data from the GSF model
#'
#' @param N Number of observations
#' @param true_params List of true parameter values
#' @return List with simulated data and true values
simulate_gsf_data <- function(N = 200, true_params = NULL) {

  if (is.null(true_params)) {
    true_params <- list(
      # Production function parameters
      alpha0 = 2.0,
      alpha = c(0.5, 0.4),  # Labor, Capital elasticities
      B = matrix(c(-0.08, 0.02,
                   0.02, -0.05), 2, 2, byrow = TRUE),

      # Variance parameters
      sigma = 0.1,          # Noise SD
      sigma_u = 0.04,       # Output inefficiency SD

      # Input slack parameters
      sigma_slack = c(0.05, 0.04),  # SD for each input slack
      rho_slack = 0.3,              # Correlation between slacks

      # Delta coefficients (intercept only for simplicity)
      Delta = matrix(c(-0.15, -0.12), 2, 1)  # Mean slacks (negative)
    )
  }

  J <- 2  # Number of inputs
  M_slack <- 1  # Just intercept

  # Generate log inputs
  X <- matrix(0, N, J)
  X[, 1] <- rnorm(N, 3, 0.5)  # log(Labor)
  X[, 2] <- rnorm(N, 4, 0.6)  # log(Capital)

  # Z matrix (just intercept)
  Z_slack <- matrix(1, N, M_slack)

  # Construct covariance matrix for slacks
  Sigma_slack <- diag(true_params$sigma_slack) %*%
    matrix(c(1, true_params$rho_slack,
             true_params$rho_slack, 1), 2, 2) %*%
    diag(true_params$sigma_slack)

  # Generate input slacks (truncated multivariate normal, theta <= 0)
  theta <- matrix(0, N, J)
  mu_theta <- as.vector(true_params$Delta)  # Same for all obs (intercept only)

  for (i in 1:N) {
    # Sample until we get theta <= 0
    accepted <- FALSE
    while (!accepted) {
      theta_i <- MASS::mvrnorm(1, mu_theta, Sigma_slack)
      if (all(theta_i <= 0)) {
        theta[i, ] <- theta_i
        accepted <- TRUE
      }
    }
  }

  # Generate output inefficiency (half-normal, u0 <= 0)
  u0 <- -abs(rnorm(N, 0, true_params$sigma_u))

  # Generate noise
  epsilon <- rnorm(N, 0, true_params$sigma)

  # Compute output using translog production function
  y <- numeric(N)
  for (i in 1:N) {
    x_eff <- X[i, ] + theta[i, ]
    y[i] <- true_params$alpha0 +
      sum(true_params$alpha * x_eff) +
      0.5 * t(x_eff) %*% true_params$B %*% x_eff +
      u0[i] + epsilon[i]
  }

  return(list(
    y = y,
    X = X,
    Z_slack = Z_slack,
    N = N,
    J = J,
    M_slack = M_slack,
    true_params = true_params,
    true_theta = theta,
    true_u0 = u0
  ))
}

#' Run parameter recovery test
#'
#' @param sim_data Simulated data from simulate_gsf_data
#' @param iter Number of MCMC iterations
#' @param warmup Number of warmup iterations
#' @param chains Number of chains
#' @return Stan fit object and comparison results
run_recovery_test <- function(sim_data,
                               iter = 2000,
                               warmup = 1000,
                               chains = 4) {

  # Prepare Stan data
  stan_data <- list(
    N = sim_data$N,
    J = sim_data$J,
    M_slack = sim_data$M_slack,
    M_neutral = 0,
    y = sim_data$y,
    X = sim_data$X,
    Z_slack = sim_data$Z_slack,
    Z_neutral = matrix(0, sim_data$N, 0),
    prior_sigma_scale = 1,
    prior_sigma_u_scale = 0.5,
    prior_delta_scale = 10,
    prior_omega_scale = 0.1
  )

  # Compile Stan model
  stan_file <- file.path(getwd(), "inst", "stan", "gsf_model.stan")
  if (!file.exists(stan_file)) {
    stan_file <- system.file("stan", "gsf_model.stan", package = "gsf")
  }

  cat("Compiling Stan model from:", stan_file, "\n")
  model <- stan_model(stan_file)

  # Run MCMC
  cat("\nRunning MCMC with", chains, "chains,", iter, "iterations each...\n")
  cat("This may take several minutes...\n\n")

  fit <- sampling(
    model,
    data = stan_data,
    chains = chains,
    iter = iter,
    warmup = warmup,
    cores = chains,
    seed = 123,
    control = list(adapt_delta = 0.95, max_treedepth = 12)
  )

  return(fit)
}

#' Compare estimated and true parameters
#'
#' @param fit Stan fit object
#' @param sim_data Simulated data with true parameters
#' @return Data frame with comparison
compare_parameters <- function(fit, sim_data) {

  true <- sim_data$true_params
  summ <- summary(fit)$summary

  # Extract key parameters
  results <- data.frame(
    parameter = character(),
    true_value = numeric(),
    est_mean = numeric(),
    est_sd = numeric(),
    ci_lower = numeric(),
    ci_upper = numeric(),
    covered = logical(),
    stringsAsFactors = FALSE
  )

  # alpha0
  results <- rbind(results, data.frame(
    parameter = "alpha0",
    true_value = true$alpha0,
    est_mean = summ["alpha0", "mean"],
    est_sd = summ["alpha0", "sd"],
    ci_lower = summ["alpha0", "2.5%"],
    ci_upper = summ["alpha0", "97.5%"],
    covered = true$alpha0 >= summ["alpha0", "2.5%"] & true$alpha0 <= summ["alpha0", "97.5%"]
  ))

  # alpha
  for (j in 1:2) {
    param_name <- paste0("alpha[", j, "]")
    results <- rbind(results, data.frame(
      parameter = param_name,
      true_value = true$alpha[j],
      est_mean = summ[param_name, "mean"],
      est_sd = summ[param_name, "sd"],
      ci_lower = summ[param_name, "2.5%"],
      ci_upper = summ[param_name, "97.5%"],
      covered = true$alpha[j] >= summ[param_name, "2.5%"] & true$alpha[j] <= summ[param_name, "97.5%"]
    ))
  }

  # B matrix (vech form)
  B_vech_true <- c(true$B[1,1], true$B[1,2], true$B[2,2])
  for (idx in 1:3) {
    param_name <- paste0("B_vech[", idx, "]")
    results <- rbind(results, data.frame(
      parameter = param_name,
      true_value = B_vech_true[idx],
      est_mean = summ[param_name, "mean"],
      est_sd = summ[param_name, "sd"],
      ci_lower = summ[param_name, "2.5%"],
      ci_upper = summ[param_name, "97.5%"],
      covered = B_vech_true[idx] >= summ[param_name, "2.5%"] & B_vech_true[idx] <= summ[param_name, "97.5%"]
    ))
  }

  # sigma
  results <- rbind(results, data.frame(
    parameter = "sigma",
    true_value = true$sigma,
    est_mean = summ["sigma", "mean"],
    est_sd = summ["sigma", "sd"],
    ci_lower = summ["sigma", "2.5%"],
    ci_upper = summ["sigma", "97.5%"],
    covered = true$sigma >= summ["sigma", "2.5%"] & true$sigma <= summ["sigma", "97.5%"]
  ))

  # sigma_u
  results <- rbind(results, data.frame(
    parameter = "sigma_u",
    true_value = true$sigma_u,
    est_mean = summ["sigma_u", "mean"],
    est_sd = summ["sigma_u", "sd"],
    ci_lower = summ["sigma_u", "2.5%"],
    ci_upper = summ["sigma_u", "97.5%"],
    covered = true$sigma_u >= summ["sigma_u", "2.5%"] & true$sigma_u <= summ["sigma_u", "97.5%"]
  ))

  # sigma_slack
  for (j in 1:2) {
    param_name <- paste0("sigma_slack[", j, "]")
    results <- rbind(results, data.frame(
      parameter = param_name,
      true_value = true$sigma_slack[j],
      est_mean = summ[param_name, "mean"],
      est_sd = summ[param_name, "sd"],
      ci_lower = summ[param_name, "2.5%"],
      ci_upper = summ[param_name, "97.5%"],
      covered = true$sigma_slack[j] >= summ[param_name, "2.5%"] & true$sigma_slack[j] <= summ[param_name, "97.5%"]
    ))
  }

  # Delta coefficients
  for (j in 1:2) {
    param_name <- paste0("Delta_raw[", j, ",1]")
    results <- rbind(results, data.frame(
      parameter = param_name,
      true_value = true$Delta[j, 1],
      est_mean = summ[param_name, "mean"],
      est_sd = summ[param_name, "sd"],
      ci_lower = summ[param_name, "2.5%"],
      ci_upper = summ[param_name, "97.5%"],
      covered = true$Delta[j, 1] >= summ[param_name, "2.5%"] & true$Delta[j, 1] <= summ[param_name, "97.5%"]
    ))
  }

  rownames(results) <- NULL
  return(results)
}

#' Compare efficiency estimates
#'
#' @param fit Stan fit object
#' @param sim_data Simulated data with true efficiencies
#' @return List with correlation and RMSE
compare_efficiencies <- function(fit, sim_data) {

  # Extract posterior means
  samples <- extract(fit)

  # Technical inefficiency
  est_u0 <- colMeans(samples$u0)
  true_u0 <- sim_data$true_u0
  u0_cor <- cor(est_u0, true_u0)
  u0_rmse <- sqrt(mean((est_u0 - true_u0)^2))

  # Input slacks
  est_theta <- apply(samples$theta, c(2, 3), mean)
  true_theta <- sim_data$true_theta

  theta1_cor <- cor(est_theta[, 1], true_theta[, 1])
  theta1_rmse <- sqrt(mean((est_theta[, 1] - true_theta[, 1])^2))

  theta2_cor <- cor(est_theta[, 2], true_theta[, 2])
  theta2_rmse <- sqrt(mean((est_theta[, 2] - true_theta[, 2])^2))

  return(list(
    u0 = list(correlation = u0_cor, rmse = u0_rmse),
    theta1 = list(correlation = theta1_cor, rmse = theta1_rmse),
    theta2 = list(correlation = theta2_cor, rmse = theta2_rmse)
  ))
}

# ============================================
# Main execution
# ============================================

cat("Step 1: Simulating data from GSF model...\n")
sim_data <- simulate_gsf_data(N = 200)

cat("\nTrue parameter values:\n")
cat("  alpha0:", sim_data$true_params$alpha0, "\n")
cat("  alpha:", sim_data$true_params$alpha, "\n")
cat("  B[1,1]:", sim_data$true_params$B[1,1], "\n")
cat("  B[1,2]:", sim_data$true_params$B[1,2], "\n")
cat("  B[2,2]:", sim_data$true_params$B[2,2], "\n")
cat("  sigma:", sim_data$true_params$sigma, "\n")
cat("  sigma_u:", sim_data$true_params$sigma_u, "\n")
cat("  sigma_slack:", sim_data$true_params$sigma_slack, "\n")
cat("  Delta:", sim_data$true_params$Delta, "\n")

cat("\nTrue efficiency summary:\n")
cat("  Mean technical inefficiency (u0):", mean(-sim_data$true_u0) * 100, "%\n")
cat("  Mean input slack 1 (theta_L):", mean(-sim_data$true_theta[,1]) * 100, "%\n")
cat("  Mean input slack 2 (theta_K):", mean(-sim_data$true_theta[,2]) * 100, "%\n")

cat("\nStep 2: Fitting GSF model...\n")
fit <- run_recovery_test(sim_data, iter = 2000, warmup = 1000, chains = 4)

cat("\nStep 3: Checking convergence...\n")
summ <- summary(fit)$summary
rhats <- summ[, "Rhat"]
rhats <- rhats[!is.na(rhats) & is.finite(rhats)]
max_rhat <- max(rhats)
cat("  Max Rhat:", round(max_rhat, 3), "\n")

if (max_rhat > 1.1) {
  cat("  WARNING: Some parameters have Rhat > 1.1\n")
  bad_params <- names(rhats)[rhats > 1.1]
  cat("  Parameters:", paste(head(bad_params, 10), collapse = ", "), "\n")
} else {
  cat("  All parameters converged (Rhat < 1.1)\n")
}

n_eff <- summ[, "n_eff"]
n_eff <- n_eff[!is.na(n_eff) & is.finite(n_eff)]
min_neff <- min(n_eff)
cat("  Min effective sample size:", round(min_neff), "\n")

cat("\nStep 4: Comparing parameter estimates...\n")
comparison <- compare_parameters(fit, sim_data)
print(comparison, digits = 3)

coverage_rate <- mean(comparison$covered)
cat("\n95% CI Coverage rate:", round(coverage_rate * 100, 1), "%\n")
cat("(Expected: ~95% if model is correct)\n")

cat("\nStep 5: Comparing efficiency estimates...\n")
eff_comparison <- compare_efficiencies(fit, sim_data)

cat("\nLatent variable recovery:\n")
cat("  Technical inefficiency (u0):\n")
cat("    Correlation with true:", round(eff_comparison$u0$correlation, 3), "\n")
cat("    RMSE:", round(eff_comparison$u0$rmse, 4), "\n")

cat("  Input slack 1 (theta_L):\n")
cat("    Correlation with true:", round(eff_comparison$theta1$correlation, 3), "\n")
cat("    RMSE:", round(eff_comparison$theta1$rmse, 4), "\n")

cat("  Input slack 2 (theta_K):\n")
cat("    Correlation with true:", round(eff_comparison$theta2$correlation, 3), "\n")
cat("    RMSE:", round(eff_comparison$theta2$rmse, 4), "\n")

cat("\n=======================================================\n")
cat("Parameter Recovery Test Complete\n")
cat("=======================================================\n")

# Save results
results <- list(
  simulation_data = sim_data,
  stanfit = fit,
  parameter_comparison = comparison,
  efficiency_comparison = eff_comparison,
  convergence = list(max_rhat = max_rhat, min_neff = min_neff),
  coverage_rate = coverage_rate
)

saveRDS(results, "parameter_recovery_results.rds")
cat("\nResults saved to: parameter_recovery_results.rds\n")

# Return pass/fail status
if (coverage_rate >= 0.80 && max_rhat < 1.1 &&
    eff_comparison$u0$correlation > 0.5 &&
    eff_comparison$theta1$correlation > 0.5 &&
    eff_comparison$theta2$correlation > 0.5) {
  cat("\n*** TEST PASSED ***\n")
  quit(status = 0)
} else {
  cat("\n*** TEST FAILED ***\n")
  cat("Possible issues:\n")
  if (coverage_rate < 0.80) cat("  - Low coverage rate\n")
  if (max_rhat >= 1.1) cat("  - Convergence issues\n")
  if (eff_comparison$u0$correlation <= 0.5) cat("  - Poor u0 recovery\n")
  if (eff_comparison$theta1$correlation <= 0.5) cat("  - Poor theta1 recovery\n")
  if (eff_comparison$theta2$correlation <= 0.5) cat("  - Poor theta2 recovery\n")
  quit(status = 1)
}
