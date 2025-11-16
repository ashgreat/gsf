#' Simulated Data for Parameter Recovery Checks
#'
#' `gsf_recovery` is a compact panel generated from the same data-generating
#' process used in the Stan model. It is meant for unit tests or vignettes that
#' demonstrate how the package can recover production-function parameters,
#' slack determinants, and inefficiency profiles when the model is correctly
#' specified.
#'
#' @format A data frame with 50 observations and the following columns:
#' \describe{
#'   \item{firm_id}{Integer panel identifier (10 firms).}
#'   \item{time_id}{Integer period indicator (5 periods per firm).}
#'   \item{Y}{Output in levels.}
#'   \item{L, K}{Input levels (labor and capital).}
#'   \item{CR, MKSH, FP, IMP}{Exogenous determinants of input slacks.}
#' }
#'
#' @details
#' The object carries two attributes:
#' \itemize{
#'   \item `attr(., "true_params")` — list containing the production-function
#'         coefficients, covariance parameters, and slack determinants used to
#'         generate the data.
#'   \item `attr(., "true_efficiency")` — list containing the exact input
#'         slacks and output technical inefficiency for every observation,
#'         expressed in percentage terms.
#' }
#'
#' Analysts can fit the model via:
#' \preformatted{
#'   fit <- gsf_fit(
#'     log(Y) ~ log(L) + log(K),
#'     data = gsf_recovery,
#'     z_vars = c("CR", "MKSH", "FP", "IMP")
#'   )
#' }
#' and then compare `efficiency_summary(fit)` or `coef(fit)` against the stored
#' ground truth.
#'
#' @source Generated via `simulate_gsf_data(seed = 2025)`; see
#'   `inst/scripts/simulate_data.R` for implementation details.
"gsf_recovery"
