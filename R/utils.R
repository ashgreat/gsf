#' Extract Inefficiencies
#'
#' @param fit An object of class `stanfit`
#' @return A list containing posterior means of input slacks (theta) and output inefficiency (u0)
#' @export
extract_inefficiencies <- function(fit) {
    if (!inherits(fit, "stanfit")) {
        stop("Object must be of class 'stanfit'")
    }

    samples <- rstan::extract(fit)

    theta_mean <- apply(samples$theta, c(2, 3), mean)
    u0_mean <- apply(samples$u0, 2, mean)

    list(
        input_slacks = theta_mean,
        output_inefficiency = u0_mean
    )
}
