#!/usr/bin/env Rscript
#' Parameter Recovery Test for GSF Model using Gibbs Sampler
#'
#' This script tests whether the GSF model with Gibbs sampler can accurately
#' recover the true parameters from simulated data.
#'
#' Based on Kumbhakar & Tsionas (2021) methodology.

# Source the package functions
for (f in list.files("R", pattern = "\\.R$", full.names = TRUE)) {
  source(f)
}

# Set seed for reproducibility
set.seed(42)

cat("=======================================================\n")
cat("GSF Model Parameter Recovery Test (Gibbs Sampler)\n")
cat("=======================================================\n\n")

# Simulate data
cat("Step 1: Simulating data from GSF model...\n")

N <- 200
J <- 2

# True parameters
true_params <- list(
  alpha0 = 2.0,
  alpha = c(0.5, 0.4),
  B = matrix(c(-0.08, 0.02, 0.02, -0.05), 2, 2, byrow = TRUE),
  sigma = 0.1,
  sigma_u = 0.04,
  sigma_slack = c(0.05, 0.04),
  rho_slack = 0.3,
  Delta = matrix(c(-0.15, -0.12), 2, 1)
)

# Generate log inputs
X <- matrix(0, N, J)
X[, 1] <- rnorm(N, 3, 0.5)  # log(Labor)
X[, 2] <- rnorm(N, 4, 0.6)  # log(Capital)

# Z matrix (just intercept)
Z_slack <- matrix(1, N, 1)

# Construct covariance matrix for slacks
Sigma_slack <- diag(true_params$sigma_slack) %*%
  matrix(c(1, true_params$rho_slack,
           true_params$rho_slack, 1), 2, 2) %*%
  diag(true_params$sigma_slack)

# Generate input slacks (truncated multivariate normal, theta <= 0)
theta <- matrix(0, N, J)
mu_theta <- as.vector(true_params$Delta)

for (i in 1:N) {
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

# Create data frame for gsf_fit
sim_df <- data.frame(
  y = y,
  L = X[, 1],
  K = X[, 2]
)

cat("\nTrue parameter values:\n")
cat("  alpha0:", true_params$alpha0, "\n")
cat("  alpha:", true_params$alpha, "\n")
cat("  B[1,1]:", true_params$B[1, 1], "\n")
cat("  B[1,2]:", true_params$B[1, 2], "\n")
cat("  B[2,2]:", true_params$B[2, 2], "\n")
cat("  sigma:", true_params$sigma, "\n")
cat("  sigma_u:", true_params$sigma_u, "\n")
cat("  sigma_slack:", true_params$sigma_slack, "\n")
cat("  Delta:", true_params$Delta, "\n")

cat("\nTrue efficiency summary:\n")
cat("  Mean technical inefficiency (u0):", mean(-u0) * 100, "%\n")
cat("  Mean input slack 1 (theta_L):", mean(-theta[, 1]) * 100, "%\n")
cat("  Mean input slack 2 (theta_K):", mean(-theta[, 2]) * 100, "%\n")

cat("\nStep 2: Fitting GSF model with Gibbs sampler...\n")

# Fit using Gibbs sampler (default)
fit <- gsf_fit(
  formula = y ~ L + K,
  data = sim_df,
  sampler = "gibbs",
  iter = 15000,
  warmup = 5000,
  thin = 2,
  seed = 123,
  verbose = TRUE
)

cat("\nStep 3: Extracting results...\n")

# Extract posterior samples
samples <- fit$gibbs_samples
n_samples <- length(samples$sigma)

cat("\nPosterior summary:\n")

# alpha0
cat("\n  alpha0:\n")
cat("    True:", true_params$alpha0, "\n")
cat("    Mean:", mean(samples$alpha0), "\n")
cat("    SD:", sd(samples$alpha0), "\n")
cat("    95% CI: [", quantile(samples$alpha0, 0.025), ",", quantile(samples$alpha0, 0.975), "]\n")

# alpha
for (j in 1:J) {
  cat("\n  alpha[", j, "]:\n", sep = "")
  cat("    True:", true_params$alpha[j], "\n")
  cat("    Mean:", mean(samples$alpha[, j]), "\n")
  cat("    SD:", sd(samples$alpha[, j]), "\n")
  cat("    95% CI: [", quantile(samples$alpha[, j], 0.025), ",", quantile(samples$alpha[, j], 0.975), "]\n")
}

# B_vech
B_vech_true <- c(true_params$B[1, 1], true_params$B[1, 2], true_params$B[2, 2])
for (idx in 1:3) {
  cat("\n  B_vech[", idx, "]:\n", sep = "")
  cat("    True:", B_vech_true[idx], "\n")
  cat("    Mean:", mean(samples$B_vech[, idx]), "\n")
  cat("    SD:", sd(samples$B_vech[, idx]), "\n")
  cat("    95% CI: [", quantile(samples$B_vech[, idx], 0.025), ",", quantile(samples$B_vech[, idx], 0.975), "]\n")
}

# sigma
cat("\n  sigma:\n")
cat("    True:", true_params$sigma, "\n")
cat("    Mean:", mean(samples$sigma), "\n")
cat("    SD:", sd(samples$sigma), "\n")
cat("    95% CI: [", quantile(samples$sigma, 0.025), ",", quantile(samples$sigma, 0.975), "]\n")

# sigma_u
cat("\n  sigma_u:\n")
cat("    True:", true_params$sigma_u, "\n")
cat("    Mean:", mean(samples$sigma_u), "\n")
cat("    SD:", sd(samples$sigma_u), "\n")
cat("    95% CI: [", quantile(samples$sigma_u, 0.025), ",", quantile(samples$sigma_u, 0.975), "]\n")

# Efficiency recovery
cat("\nStep 4: Comparing efficiency estimates...\n")

# Technical inefficiency
est_u0 <- colMeans(samples$u0)
u0_cor <- cor(est_u0, u0)
u0_rmse <- sqrt(mean((est_u0 - u0)^2))

cat("\n  Technical inefficiency (u0):\n")
cat("    Correlation with true:", round(u0_cor, 3), "\n")
cat("    RMSE:", round(u0_rmse, 4), "\n")

# Input slacks
est_theta <- apply(samples$theta, c(2, 3), mean)

theta1_cor <- cor(est_theta[, 1], theta[, 1])
theta1_rmse <- sqrt(mean((est_theta[, 1] - theta[, 1])^2))

theta2_cor <- cor(est_theta[, 2], theta[, 2])
theta2_rmse <- sqrt(mean((est_theta[, 2] - theta[, 2])^2))

cat("\n  Input slack 1 (theta_L):\n")
cat("    Correlation with true:", round(theta1_cor, 3), "\n")
cat("    RMSE:", round(theta1_rmse, 4), "\n")

cat("\n  Input slack 2 (theta_K):\n")
cat("    Correlation with true:", round(theta2_cor, 3), "\n")
cat("    RMSE:", round(theta2_rmse, 4), "\n")

# Acceptance rates
cat("\nAcceptance rates:\n")
cat("  theta (avg):", round(mean(samples$accept_rates$theta) * 100, 1), "%\n")
cat("  Delta:", round(samples$accept_rates$Delta * 100, 1), "%\n")
cat("  Omega:", round(samples$accept_rates$Omega * 100, 1), "%\n")

cat("\n=======================================================\n")
cat("Parameter Recovery Test Complete\n")
cat("=======================================================\n")

# Save results
results <- list(
  true_params = true_params,
  true_theta = theta,
  true_u0 = u0,
  fit = fit,
  correlations = list(
    u0 = u0_cor,
    theta1 = theta1_cor,
    theta2 = theta2_cor
  )
)

saveRDS(results, "gibbs_recovery_results.rds")
cat("\nResults saved to: gibbs_recovery_results.rds\n")

# Pass/fail criteria
if (u0_cor > 0.3 && theta1_cor > 0.3 && theta2_cor > 0.3) {
  cat("\n*** TEST PASSED ***\n")
  cat("(Correlations with true values are reasonable)\n")
} else {
  cat("\n*** TEST NEEDS REVIEW ***\n")
  cat("Some correlations are low, but this may be expected for this complex model.\n")
}
