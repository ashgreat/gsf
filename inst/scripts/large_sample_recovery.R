#!/usr/bin/env Rscript

# Large sample parameter recovery test using Gibbs sampler
# Uses N=10,000 observations for more precise estimation of quadratic terms (B matrix)

# Source R files
r_files <- list.files("R", pattern = "\\.R$", full.names = TRUE)
invisible(lapply(r_files, source))

set.seed(42)

# Simulate large dataset
cat("=== Simulating Large Dataset ===\n")
cat("This may take a moment...\n\n")

# Use more firms and periods for N=10,000
sim <- simulate_gsf_data_strong_id(n_firms = 1000, n_periods = 10, seed = 42)
true_params <- attr(sim, "true_params")
true_eff <- attr(sim, "true_efficiency")

cat("=== Simulated Data ===\n")
cat("Observations:", nrow(sim), "\n")
cat("Firms:", length(unique(sim$firm_id)), "\n")
cat("Periods:", length(unique(sim$time_id)), "\n\n")

# Show true parameters
cat("=== True Parameters ===\n")
cat("alpha0:", true_params$alpha0, "\n")
cat("alpha:", true_params$alpha, "\n")
cat("B diagonal:", diag(true_params$B), "\n")
cat("B off-diagonal:", true_params$B[1,2], "\n")
cat("sigma:", true_params$sigma, "\n")
cat("sigma_u:", true_params$sigma_u, "\n")
cat("sigma_slack:", sqrt(diag(true_params$Omega)), "\n\n")

# Fit using Gibbs sampler
cat("=== Fitting with Gibbs Sampler ===\n")
cat("This will take several minutes for N=10,000...\n\n")

fit <- gsf_fit(
  formula = log(Y) ~ log(X1) + log(X2),
  data = sim,
  z_vars = c("z1", "z2"),
  sampler = "gibbs",
  iter = 15000,
  warmup = 5000,
  seed = 42,
  verbose = TRUE
)

# Extract posterior means
samples <- fit$gibbs_samples
alpha0_est <- mean(samples$alpha0)
alpha_est <- colMeans(samples$alpha)
B_vech_est <- colMeans(samples$B_vech)
sigma_est <- mean(samples$sigma)
sigma_u_est <- mean(samples$sigma_u)
sigma_slack_est <- sqrt(colMeans(apply(samples$Omega, 1, diag)))

# Reconstruct B matrix for comparison
B_est <- matrix(0, 2, 2)
B_est[1,1] <- B_vech_est[1]
B_est[1,2] <- B_est[2,1] <- B_vech_est[2]
B_est[2,2] <- B_vech_est[3]

cat("\n=== Parameter Recovery ===\n\n")
cat(sprintf("%-15s %10s %10s %10s\n", "Parameter", "True", "Estimated", "Error%"))
cat(sprintf("%-15s %10.4f %10.4f %10.1f%%\n", "alpha0", true_params$alpha0, alpha0_est, abs(alpha0_est-true_params$alpha0)/abs(true_params$alpha0)*100))
cat(sprintf("%-15s %10.4f %10.4f %10.1f%%\n", "alpha[1]", true_params$alpha[1], alpha_est[1], abs(alpha_est[1]-true_params$alpha[1])/abs(true_params$alpha[1])*100))
cat(sprintf("%-15s %10.4f %10.4f %10.1f%%\n", "alpha[2]", true_params$alpha[2], alpha_est[2], abs(alpha_est[2]-true_params$alpha[2])/abs(true_params$alpha[2])*100))
cat(sprintf("%-15s %10.4f %10.4f %10.1f%%\n", "B[1,1]", true_params$B[1,1], B_est[1,1], abs(B_est[1,1]-true_params$B[1,1])/abs(true_params$B[1,1])*100))
cat(sprintf("%-15s %10.4f %10.4f %10.1f%%\n", "B[1,2]", true_params$B[1,2], B_est[1,2], abs(B_est[1,2]-true_params$B[1,2])/abs(true_params$B[1,2])*100))
cat(sprintf("%-15s %10.4f %10.4f %10.1f%%\n", "B[2,2]", true_params$B[2,2], B_est[2,2], abs(B_est[2,2]-true_params$B[2,2])/abs(true_params$B[2,2])*100))
cat(sprintf("%-15s %10.4f %10.4f %10.1f%%\n", "sigma", true_params$sigma, sigma_est, abs(sigma_est-true_params$sigma)/abs(true_params$sigma)*100))
cat(sprintf("%-15s %10.4f %10.4f %10.1f%%\n", "sigma_u", true_params$sigma_u, sigma_u_est, abs(sigma_u_est-true_params$sigma_u)/abs(true_params$sigma_u)*100))
cat(sprintf("%-15s %10.4f %10.4f %10.1f%%\n", "sigma_slack[1]", sqrt(true_params$Omega[1,1]), sigma_slack_est[1], abs(sigma_slack_est[1]-sqrt(true_params$Omega[1,1]))/sqrt(true_params$Omega[1,1])*100))
cat(sprintf("%-15s %10.4f %10.4f %10.1f%%\n", "sigma_slack[2]", sqrt(true_params$Omega[2,2]), sigma_slack_est[2], abs(sigma_slack_est[2]-sqrt(true_params$Omega[2,2]))/sqrt(true_params$Omega[2,2])*100))

# Efficiency recovery
u0_hat <- colMeans(samples$u0)
theta_hat <- apply(samples$theta, c(2, 3), mean)

cat("\n=== Efficiency Recovery ===\n")
cat("u0 correlation:", round(cor(u0_hat, true_eff$u0), 4), "\n")
cat("theta[1] correlation:", round(cor(theta_hat[,1], true_eff$theta[,1]), 4), "\n")
cat("theta[2] correlation:", round(cor(theta_hat[,2], true_eff$theta[,2]), 4), "\n")

cat("\n=== Summary Statistics ===\n")
cat("True mean |u0|:", round(mean(abs(true_eff$u0))*100, 2), "%\n")
cat("Estimated mean |u0|:", round(mean(abs(u0_hat))*100, 2), "%\n")
cat("True mean |theta|:", round(mean(abs(true_eff$theta))*100, 2), "%\n")
cat("Estimated mean |theta|:", round(mean(abs(theta_hat))*100, 2), "%\n")

# Save results
results <- list(
  true_params = true_params,
  estimated = list(
    alpha0 = alpha0_est,
    alpha = alpha_est,
    B = B_est,
    sigma = sigma_est,
    sigma_u = sigma_u_est,
    sigma_slack = sigma_slack_est
  ),
  efficiency_correlations = list(
    u0 = cor(u0_hat, true_eff$u0),
    theta1 = cor(theta_hat[,1], true_eff$theta[,1]),
    theta2 = cor(theta_hat[,2], true_eff$theta[,2])
  ),
  samples = samples
)

saveRDS(results, "large_sample_recovery_results.rds")
cat("\nResults saved to large_sample_recovery_results.rds\n")

cat("\nDone.\n")
