#' gsf: Generalized Stochastic Frontier Model with Input-Specific Inefficiencies
#'
#' @description
#' Implements the Generalized Stochastic Frontier (GSF) model from
#' Kumbhakar and Tsionas (2021) "Dissections of input and output efficiency:
#' A generalized stochastic frontier model". The model jointly estimates output
#' technical inefficiency and input-specific inefficiencies (slacks) using
#' Bayesian MCMC methods via Stan.
#'
#' @section Main Features:
#' \itemize{
#'   \item Joint estimation of output technical inefficiency and input slacks
#'   \item Translog production function with flexible substitution patterns
#'   \item Exogenous determinants of input slacks (z variables)
#'   \item Full Bayesian inference via Markov Chain Monte Carlo (Stan)
#'   \item Comprehensive diagnostics and visualization tools
#' }
#'
#' @section Key Functions:
#' \describe{
#'   \item{\code{\link{gsf_fit}}}{Main function to fit the GSF model}
#'   \item{\code{\link{extract_technical_inefficiency}}}{Extract output technical inefficiency}
#'   \item{\code{\link{extract_input_slacks}}}{Extract input-specific slacks}
#'   \item{\code{\link{extract_elasticities}}}{Extract input elasticities}
#'   \item{\code{\link{efficiency_summary}}}{Summary table of efficiency measures}
#'   \item{\code{\link{plot.gsf_fit}}}{Various diagnostic plots}
#' }
#'
#' @section Model Specification:
#' The GSF model uses a translog production function:
#' \deqn{y_i = \alpha_0 + \delta'z_i + \alpha'(x_i + \theta_i) + \frac{1}{2}(x_i + \theta_i)'B(x_i + \theta_i) + u_{0i} + \varepsilon_{0i}}
#'
#' Where:
#' \itemize{
#'   \item \eqn{y_i} is log output
#'   \item \eqn{x_i} is the vector of log inputs
#'   \item \eqn{\theta_i \leq 0} are input-specific slacks (negative values indicate over-use)
#'   \item \eqn{u_{0i} \leq 0} is output technical inefficiency
#'   \item \eqn{\varepsilon_{0i}} is random noise
#'   \item \eqn{z_i} are exogenous variables that determine slacks
#' }
#'
#' The input slacks follow:
#' \deqn{\theta_i = \Delta z_i + u_i, \quad u_i \sim N_J(0, \Omega), \quad \theta_i \leq 0}
#'
#' @section References:
#' Kumbhakar, S.C. and Tsionas, M.G. (2021). "Dissections of input and output
#' efficiency: A generalized stochastic frontier model." International Journal
#' of Production Economics, 232, 107940.
#'
#' @docType package
#' @name gsf-package
#' @aliases gsf
#'
#' @import methods
#' @importFrom rstan sampling extract stan_model traceplot get_divergent_iterations
#' @importFrom stats model.frame model.matrix model.response complete.cases
#' @importFrom stats na.pass median var sd quantile logLik rnorm runif rbeta
#' @importFrom ggplot2 ggplot aes geom_density geom_vline geom_pointrange
#' @importFrom ggplot2 facet_wrap labs theme_minimal theme coord_flip
#' @importFrom ggplot2 position_dodge annotate
#' @importFrom tidyr pivot_wider
#' @importFrom parallel detectCores
#'
"_PACKAGE"

# The following block is used by usethis to automatically manage
# temporary build directory when running R CMD check
# See: https://r-pkgs.org/src.html#c-cpp-best-practices
.onLoad <- function(libname, pkgname) {
  # Compile Stan model on load (optional, can be slow)
  # This is handled dynamically in gsf_fit instead
}
