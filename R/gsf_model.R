#' Fit Generalized Stochastic Frontier Model
#'
#' @param y Vector of log outputs
#' @param X Matrix of log inputs
#' @param Z Matrix of exogenous variables
#' @param ... Additional arguments passed to rstan::sampling
#' @return An object of class `stanfit`
#' @export
gsf <- function(y, X, Z, ...) {
    if (!requireNamespace("rstan", quietly = TRUE)) {
        stop("Package 'rstan' is required for this function.")
    }

    data_list <- list(
        N = length(y),
        J = ncol(X),
        M = ncol(Z),
        y = as.vector(y),
        X = as.matrix(X),
        Z = as.matrix(Z)
    )

    # Find the Stan model file
    stan_file <- system.file("stan_files", "gsf.stan", package = "GSF")
    if (stan_file == "") {
        # Fallback for development (when not installed)
        stan_file <- "src/stan_files/gsf.stan"
    }

    fit <- rstan::stan(file = stan_file, data = data_list, ...)
    return(fit)
}
