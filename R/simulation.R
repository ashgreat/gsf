#' Simulate Data for GSF Model
#'
#' @param N Number of observations
#' @param J Number of inputs
#' @param M Number of z variables
#' @return A list containing simulated data and parameters
#' @export
simulate_gsf_data <- function(N = 100, J = 2, M = 2) {
    set.seed(123)

    # Exogenous variables
    Z <- matrix(rnorm(N * M), N, M)
    X <- matrix(rnorm(N * J), N, J)

    # Parameters
    alpha0 <- 1.0
    alpha <- rep(0.5, J)
    B <- matrix(0.1, J, J)
    diag(B) <- 0.2
    delta <- rep(0.1, M)
    Delta <- matrix(0.1, J, M)
    Omega <- diag(0.05, J)
    sigma_sq <- 0.01
    sigma_u_sq <- 0.04

    # Latent variables
    theta <- matrix(0, N, J)
    u0 <- numeric(N)

    for (i in 1:N) {
        mu_theta <- as.vector(Delta %*% Z[i, ])
        theta[i, ] <- -abs(MASS::mvrnorm(1, mu_theta, Omega)) # Ensure negative
        u0[i] <- -abs(rnorm(1, 0, sqrt(sigma_u_sq))) # Ensure negative
    }

    # Output
    y <- numeric(N)
    for (i in 1:N) {
        X_eff <- X[i, ] + theta[i, ]
        mu <- alpha0 + sum(delta * Z[i, ]) +
            sum(alpha * X_eff) +
            0.5 * t(X_eff) %*% B %*% X_eff +
            u0[i]
        y[i] <- rnorm(1, mu, sqrt(sigma_sq))
    }

    list(
        y = y, X = X, Z = Z,
        params = list(
            alpha0 = alpha0, alpha = alpha, B = B, delta = delta,
            Delta = Delta, Omega = Omega, sigma_sq = sigma_sq, sigma_u_sq = sigma_u_sq
        ),
        latent = list(theta = theta, u0 = u0)
    )
}
