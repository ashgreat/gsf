#' Fit Generalized Stochastic Frontier Model
#'
#' @description
#' Fits the Generalized Stochastic Frontier (GSF) model from Kumbhakar and
#' Tsionas (2021) that jointly estimates output technical inefficiency and
#' input-specific inefficiencies (slacks).
#'
#' @param formula A formula specifying the production function. The left-hand
#'   side should be the log output, and the right-hand side should contain
#'   the log inputs (e.g., log(Y) ~ log(L) + log(K)).
#' @param data A data frame containing the variables.
#' @param z_vars Character vector of variable names that determine input
#'   slacks. These are the exogenous variables that affect input effectiveness.
#'   If NULL, only stochastic slacks are estimated.
#' @param neutral_vars Character vector of variable names that act as neutral
#'   technology shifters (affect output directly). If NULL, no neutral shifters.
#' @param panel_id Character string naming the panel identifier variable.
#'   If NULL, cross-sectional data is assumed.
#' @param time_id Character string naming the time identifier variable.
#'   Required if panel_id is specified.
#' @param chains Number of MCMC chains (default: 4).
#' @param iter Number of iterations per chain (default: 2000).
#' @param warmup Number of warmup iterations (default: 1000).
#' @param thin Thinning rate (default: 1).
#' @param cores Number of cores for parallel computation (default: parallel::detectCores() - 1).
#' @param prior_sigma_scale Prior scale for noise standard deviation (default: 1).
#' @param prior_sigma_u_scale Prior scale for output inefficiency (default: 0.5).
#' @param prior_delta_scale Prior scale for regression coefficients (default: 10).
#' @param prior_omega_scale Prior scale for slack covariance (default: 0.1).
#' @param seed Random seed for reproducibility.
#' @param ... Additional arguments passed to rstan::sampling.
#'
#' @return An object of class "gsf_fit" containing:
#'   \item{stanfit}{The rstan fit object}
#'   \item{data}{The processed data}
#'   \item{model_info}{Information about the model specification}
#'   \item{call}{The original function call}
#'
#' @examples
#' \dontrun{
#' # Simulate data
#' set.seed(123)
#' n <- 100
#' data <- data.frame(
#'   Y = exp(rnorm(n, 5, 0.5)),
#'   L = exp(rnorm(n, 3, 0.3)),
#'   K = exp(rnorm(n, 4, 0.4)),
#'   competition = runif(n, 0, 1)
#' )
#'
#' # Fit model
#' fit <- gsf_fit(
#'   formula = log(Y) ~ log(L) + log(K),
#'   data = data,
#'   z_vars = "competition",
#'   chains = 2,
#'   iter = 1000
#' )
#'
#' # Summary
#' summary(fit)
#' }
#'
#' @export
gsf_fit <- function(formula,
                    data,
                    z_vars = NULL,
                    neutral_vars = NULL,
                    panel_id = NULL,
                    time_id = NULL,
                    chains = 4,
                    iter = 2000,
                    warmup = 1000,
                    thin = 1,
                    cores = parallel::detectCores() - 1,
                    prior_sigma_scale = 1,
                    prior_sigma_u_scale = 0.5,
                    prior_delta_scale = 10,
                    prior_omega_scale = 0.1,
                    seed = NULL,
                    ...) {

  call <- match.call()

  # Parse formula
  parsed <- parse_gsf_formula(formula, data)
  y <- parsed$y
  X <- parsed$X
  input_names <- parsed$input_names
  y_name <- parsed$y_name
  clean_data <- parsed$data

  # Number of observations and inputs
  N <- nrow(X)
  J <- ncol(X)

  # Prepare Z matrix (determinants of slacks)
  if (is.null(z_vars) || length(z_vars) == 0) {
    Z_slack <- matrix(1, nrow = N, ncol = 1)
    z_slack_names <- "intercept"
  } else {
    missing_z <- setdiff(z_vars, names(clean_data))
    if (length(missing_z) > 0) {
      stop("Variables not found in data: ", paste(missing_z, collapse = ", "))
    }
    Z_slack <- cbind(1, as.matrix(clean_data[, z_vars, drop = FALSE]))
    z_slack_names <- c("intercept", z_vars)
  }
  M_slack <- ncol(Z_slack)

  # Prepare neutral shifters
  if (is.null(neutral_vars)) {
    Z_neutral <- matrix(0, nrow = N, ncol = 0)
    neutral_names <- character(0)
  } else {
    missing_n <- setdiff(neutral_vars, names(clean_data))
    if (length(missing_n) > 0) {
      stop("Neutral variables not found in data: ", paste(missing_n, collapse = ", "))
    }
    Z_neutral <- as.matrix(clean_data[, neutral_vars, drop = FALSE])
    neutral_names <- neutral_vars
  }
  M_neutral <- ncol(Z_neutral)

  # Panel identifiers (optional)
  panel_info <- NULL
  if (!is.null(panel_id)) {
    if (is.null(time_id)) {
      stop("time_id must be supplied when panel_id is provided.")
    }
    if (!(panel_id %in% names(clean_data))) {
      stop("panel_id '", panel_id, "' not found in data.")
    }
    if (!(time_id %in% names(clean_data))) {
      stop("time_id '", time_id, "' not found in data.")
    }
    panel_info <- list(
      panel_id = panel_id,
      time_id = time_id,
      ids = clean_data[[panel_id]],
      times = clean_data[[time_id]]
    )
  }

  # Prepare Stan data
  stan_data <- list(
    N = N,
    J = J,
    M_slack = M_slack,
    M_neutral = M_neutral,
    y = as.vector(y),
    X = X,
    Z_slack = Z_slack,
    Z_neutral = Z_neutral,
    prior_sigma_scale = prior_sigma_scale,
    prior_sigma_u_scale = prior_sigma_u_scale,
    prior_delta_scale = prior_delta_scale,
    prior_omega_scale = prior_omega_scale
  )

  # Get Stan model
  stan_file <- system.file("stan", "gsf_model.stan", package = "gsf")
  if (stan_file == "") {
    # For development, use local path
    stan_file <- file.path(getwd(), "inst", "stan", "gsf_model.stan")
    if (!file.exists(stan_file)) {
      stop("Stan model file not found. Please ensure gsf package is properly installed.")
    }
  }

  # Compile and fit model
  message("Compiling Stan model...")
  stan_model <- rstan::stan_model(stan_file)

  message("Fitting GSF model with ", chains, " chains, ", iter, " iterations each...")
  message("This may take a while for large datasets...")

  if (!is.null(seed)) {
    set.seed(seed)
  }

  stanfit <- rstan::sampling(
    stan_model,
    data = stan_data,
    chains = chains,
    iter = iter,
    warmup = warmup,
    thin = thin,
    cores = cores,
    seed = seed,
    control = list(adapt_delta = 0.95, max_treedepth = 12),
    ...
  )

  # Build brmsfit shell for downstream interoperability
  brms_shell <- build_brms_shell(
    stanfit = stanfit,
    formula = formula,
    data = clean_data,
    stan_file = stan_file,
    stan_data = stan_data
  )

  # Store model information
  model_info <- list(
    formula = formula,
    terms = stats::terms(formula),
    y_name = y_name,
    input_names = input_names,
    z_slack_names = z_slack_names,
    z_slack_vars = if (is.null(z_vars)) character(0) else z_vars,
    neutral_names = neutral_names,
    neutral_vars = if (is.null(neutral_vars)) character(0) else neutral_vars,
    N = N,
    J = J,
    M_slack = M_slack,
    panel = panel_info
  )

  # Store processed data
  processed_data <- list(
    y = y,
    X = X,
    Z_slack = Z_slack,
    Z_neutral = Z_neutral,
    original_data = clean_data
  )

  # Create output object
  result <- list(
    stanfit = stanfit,
    brmsfit = brms_shell,
    data = processed_data,
    model_info = model_info,
    call = call
  )

  class(result) <- "gsf_fit"
  return(result)
}


#' Parse GSF Formula
#'
#' @param formula The model formula
#' @param data The data frame
#' @return A list with parsed components
#' @keywords internal
parse_gsf_formula <- function(formula, data) {
  # Get model frame
  mf <- model.frame(formula, data = data, na.action = na.pass)
  row_ids <- row.names(mf)
  if (is.null(row_ids)) {
    row_ids <- seq_len(nrow(mf))
  }
  data_in_model <- data[row_ids, , drop = FALSE]

  # Extract response
  y <- model.response(mf)
  y_name <- as.character(formula[[2]])

  # Extract predictors
  X <- model.matrix(formula, data = mf)

  # Remove intercept if present (we handle it separately)
  if ("(Intercept)" %in% colnames(X)) {
    X <- X[, -1, drop = FALSE]
  }

  input_names <- colnames(X)
  if (is.null(input_names) || length(input_names) == 0) {
    input_names <- paste0("X", 1:ncol(X))
  }

  complete_cases <- complete.cases(y, X)
  if (any(!complete_cases)) {
    warning("Missing values detected. Rows with NAs were removed.")
    y <- y[complete_cases]
    X <- X[complete_cases, , drop = FALSE]
    mf <- mf[complete_cases, , drop = FALSE]
    data_in_model <- data_in_model[complete_cases, , drop = FALSE]
  }

  list(
    y = as.vector(y),
    X = as.matrix(X),
    input_names = input_names,
    y_name = y_name,
    data = data_in_model
  )
}


#' Print method for gsf_fit
#'
#' @param x A gsf_fit object
#' @param ... Additional arguments (ignored)
#' @export
print.gsf_fit <- function(x, ...) {
  cat("Generalized Stochastic Frontier Model\n")
  cat("======================================\n\n")
  cat("Call:\n")
  print(x$call)
  cat("\n")
  cat("Number of observations:", x$model_info$N, "\n")
  cat("Number of inputs:", x$model_info$J, "\n")
  cat("Input names:", paste(x$model_info$input_names, collapse = ", "), "\n")

  if (length(x$model_info$z_slack_names) > 1) {
    cat("Slack determinants:", paste(x$model_info$z_slack_names[-1], collapse = ", "), "\n")
  } else {
    cat("Slack determinants: None (stochastic only)\n")
  }

  if (length(x$model_info$neutral_names) > 0) {
    cat("Neutral shifters:", paste(x$model_info$neutral_names, collapse = ", "), "\n")
  }

  cat("\nMCMC Sampling:\n")
  cat("Chains:", length(x$stanfit@stan_args), "\n")
  cat("Iterations:", x$stanfit@stan_args[[1]]$iter, "\n")
  cat("Warmup:", x$stanfit@stan_args[[1]]$warmup, "\n")

  # Quick convergence check
  rhats <- rstan::summary(x$stanfit)$summary[, "Rhat"]
  rhats <- rhats[!is.na(rhats) & is.finite(rhats)]
  max_rhat <- max(rhats)
  cat("\nMax Rhat:", round(max_rhat, 3), "\n")
  if (max_rhat > 1.1) {
    cat("WARNING: Some parameters have Rhat > 1.1. Consider more iterations.\n")
  }

  invisible(x)
}


#' Construct a minimal brmsfit object for interoperability
#' @keywords internal
build_brms_shell <- function(stanfit, formula, data, stan_file, stan_data) {
  if (!requireNamespace("brms", quietly = TRUE)) {
    return(NULL)
  }

  stan_code <- paste(readLines(stan_file), collapse = "\n")
  threads_stub <- structure(list(threads = 1L), class = "brmsthreads")
  opencl_stub <- structure(list(opencl = FALSE), class = "brmsopencl")

  structure(
    list(
      formula = brms::brmsformula(formula),
      data = data,
      prior = brms::prior(NULL),
      data2 = list(),
      stanvars = structure(list(), class = "stanvars"),
      model = stan_code,
      save_pars = brms::save_pars(all = TRUE),
      algorithm = "sampling",
      backend = "rstan",
      threads = threads_stub,
      opencl = opencl_stub,
      stan_args = stanfit@stan_args,
      fit = stanfit,
      basis = list(),
      criteria = list(),
      file = NULL,
      version = utils::packageVersion("brms"),
      family = brms::brmsfamily("gaussian"),
      autocor = NULL,
      ranef = structure(data.frame(), class = c("reframe", "data.frame")),
      cov_ranef = NULL,
      stan_funs = NULL,
      data.name = "",
      standata = stan_data
    ),
    class = "brmsfit"
  )
}
