#!/usr/bin/env Rscript

# Diagnostic script: Can we recover production function parameters
# if we KNOW the true theta and u0?

# Source R files
r_files <- list.files("R", pattern = "\\.R$", full.names = TRUE)
invisible(lapply(r_files, source))

set.seed(42)

# Simulate data
sim <- simulate_gsf_data_strong_id(n_firms = 100, n_periods = 5, seed = 42)
true_params <- attr(sim, "true_params")
true_eff <- attr(sim, "true_efficiency")

N <- nrow(sim)
J <- 2

cat("=== Checking Identification with Known Latents ===\n\n")

# Get logged values
y <- log(sim$Y)
X <- cbind(log(sim$X1), log(sim$X2))

# True values
true_theta <- true_eff$theta
true_u0 <- true_eff$u0

cat("True parameters:\n")
cat("  alpha0:", true_params$alpha0, "\n")
cat("  alpha:", true_params$alpha, "\n")
cat("  B diagonal:", diag(true_params$B), "\n")
cat("  sigma:", true_params$sigma, "\n")
cat("  sigma_u:", true_params$sigma_u, "\n\n")

# ========================================================================
# Test 1: Can we recover alpha, B from y, X, theta, u0?
# ========================================================================
cat("Test 1: OLS on y - u0 with effective inputs (true theta known)\n")

X_eff <- X + true_theta

# Build design matrix for translog
W <- cbind(
  1,  # intercept (alpha0)
  X_eff,  # first-order terms (alpha)
  0.5 * X_eff[,1]^2,  # B[1,1]
  X_eff[,1] * X_eff[,2],  # B[1,2]
  0.5 * X_eff[,2]^2   # B[2,2]
)

# Dependent variable: y - u0
y_adj <- y - true_u0

# OLS regression
fit_ols <- lm(y_adj ~ W - 1)
coefs <- coef(fit_ols)
names(coefs) <- c("alpha0", "alpha1", "alpha2", "B11", "B12", "B22")

cat("  OLS estimates:\n")
cat("    alpha0:", round(coefs["alpha0"], 4), "(true:", true_params$alpha0, ")\n")
cat("    alpha1:", round(coefs["alpha1"], 4), "(true:", true_params$alpha[1], ")\n")
cat("    alpha2:", round(coefs["alpha2"], 4), "(true:", true_params$alpha[2], ")\n")
cat("    B[1,1]:", round(coefs["B11"], 4), "(true:", true_params$B[1,1], ")\n")
cat("    B[1,2]:", round(coefs["B12"], 4), "(true:", true_params$B[1,2], ")\n")
cat("    B[2,2]:", round(coefs["B22"], 4), "(true:", true_params$B[2,2], ")\n")
cat("    sigma (residual SD):", round(sd(residuals(fit_ols)), 4), "(true:", true_params$sigma, ")\n")
cat("\n")

# ========================================================================
# Test 2: Can we recover Delta from theta, Z?
# ========================================================================
cat("Test 2: Regression of theta on Z (true theta known)\n")

Z <- cbind(1, sim$z1, sim$z2)

# For each input j, regress theta_j on Z
for (j in 1:J) {
  fit_j <- lm(true_theta[, j] ~ Z - 1)
  cat("  Delta[", j, ",]:", round(coef(fit_j), 4), "\n")
  cat("         True:", round(true_params$Delta[j,], 4), "\n")
}
cat("\n")

# ========================================================================
# Test 3: Can we estimate Omega from theta residuals?
# ========================================================================
cat("Test 3: Covariance of theta residuals\n")

theta_fitted <- Z %*% t(true_params$Delta)
theta_resid <- true_theta - theta_fitted
Omega_hat <- cov(theta_resid)

cat("  Estimated Omega:\n")
print(round(Omega_hat, 6))
cat("  True Omega:\n")
print(round(true_params$Omega, 6))
cat("\n")

# ========================================================================
# Test 4: What happens if theta = 0? How big is u0?
# ========================================================================
cat("Test 4: If we assume theta = 0, what do residuals look like?\n")

# With theta = 0, effective inputs = raw inputs
W_raw <- cbind(
  1,
  X,
  0.5 * X[,1]^2,
  X[,1] * X[,2],
  0.5 * X[,2]^2
)

fit_raw <- lm(y ~ W_raw - 1)
resid_raw <- residuals(fit_raw)

cat("  Residual SD (theta=0 assumed):", round(sd(resid_raw), 4), "\n")
cat("  True sigma + sigma_u combined:", round(sqrt(true_params$sigma^2 + true_params$sigma_u^2), 4), "\n")
cat("  True mean(|theta|):", round(mean(abs(true_theta)) * 100, 2), "%\n\n")

# The key insight: if theta is truly ~14% on average, but we assume theta=0,
# the residuals will be much larger because the production function is misspecified

cat("=== Key Insight ===\n")
cat("If the sampler assumes theta ~ 0, it will attribute the production shortfall\n")
cat("to u0 (output inefficiency) instead of theta (input slacks).\n\n")

cat("The identification of theta vs u0 depends on:\n")
cat("1. Z variables that predict theta but not u0\n")
cat("2. The quadratic B terms that create theta*x interactions\n")
cat("3. The truncation constraint (theta <= 0)\n\n")

# ========================================================================
# Test 5: Information content of Z
# ========================================================================
cat("Test 5: How much does Z explain theta?\n")

for (j in 1:J) {
  fit_j <- lm(true_theta[, j] ~ sim$z1 + sim$z2)
  cat("  R^2 for theta[", j, "]:", round(summary(fit_j)$r.squared, 4), "\n")
}
cat("\nHigher R^2 = stronger identification of theta through Z\n")
