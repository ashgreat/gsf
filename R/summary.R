#' Summary Method for gsf_fit Objects
#'
#' @description
#' Provides a comprehensive summary of the fitted GSF model including
#' parameter estimates, efficiency measures, and model diagnostics.
#'
#' @param object A gsf_fit object.
#' @param ... Additional arguments (ignored).
#'
#' @return Invisibly returns a list with summary components.
#'
#' @export
summary.gsf_fit <- function(object, ...) {
  cat("Generalized Stochastic Frontier Model\n")
  cat("=====================================\n\n")

  # Model info
  cat("Model Specification:\n")
  cat("  Formula:", deparse(object$model_info$formula), "\n")
  cat("  Observations:", object$model_info$N, "\n")
  cat("  Inputs:", paste(object$model_info$input_names, collapse = ", "), "\n")
  if (length(object$model_info$z_slack_names) > 1) {
    cat("  Slack determinants:", paste(object$model_info$z_slack_names[-1], collapse = ", "), "\n")
  }
  cat("\n")

  # MCMC diagnostics
  cat("MCMC Diagnostics:\n")
  # Check convergence
  summ <- rstan::summary(object$stanfit)$summary
  key_params <- c("alpha0", "sigma", "sigma_u")

  # Add alpha parameters
  alpha_params <- paste0("alpha[", 1:object$model_info$J, "]")
  key_params <- c(key_params, alpha_params)

  # Filter to available parameters
  available_params <- rownames(summ)
  key_params <- intersect(key_params, available_params)

  if (length(key_params) > 0) {
    param_summ <- summ[key_params, c("mean", "sd", "2.5%", "97.5%", "n_eff", "Rhat"), drop = FALSE]
    cat("\nKey Parameter Estimates:\n")
    print(round(param_summ, 4))
  }

  # Efficiency summary
  cat("\nEfficiency Measures (Sample Averages):\n")
  eff_summ <- efficiency_summary(object)
  print(eff_summ, row.names = FALSE)

  # Delta coefficients (slack determinants)
  if (length(object$model_info$z_slack_names) > 1) {
    cat("\nSlack Determinants (Delta Coefficients):\n")
    delta_coef <- extract_delta_coefficients(object)
    print(delta_coef, row.names = FALSE)
  }

  # Convergence warnings
  rhats <- summ[, "Rhat"]
  rhats <- rhats[!is.na(rhats) & is.finite(rhats)]
  if (any(rhats > 1.1)) {
    cat("\nWARNING: Some parameters have Rhat > 1.1. Consider more iterations.\n")
    bad_params <- names(rhats)[rhats > 1.1]
    if (length(bad_params) <= 10) {
      cat("Parameters with poor convergence:", paste(bad_params, collapse = ", "), "\n")
    } else {
      cat("Number of parameters with poor convergence:", length(bad_params), "\n")
    }
  }

  # Check for divergences
  divergences <- sum(rstan::get_divergent_iterations(object$stanfit))
  if (divergences > 0) {
    cat("\nWARNING:", divergences, "divergent transitions after warmup.\n")
    cat("Consider increasing adapt_delta or reparameterizing the model.\n")
  }

  invisible(list(
    model_info = object$model_info,
    efficiency_summary = eff_summ,
    param_summary = if (exists("param_summ")) param_summ else NULL
  ))
}


#' Extract Coefficients from gsf_fit
#'
#' @param object A gsf_fit object.
#' @param type Character string specifying which coefficients to extract:
#'   "production" for production function parameters, "delta" for slack
#'   determinants, or "all" for everything.
#' @param ... Additional arguments (ignored).
#'
#' @return A data frame with coefficient estimates.
#'
#' @export
coef.gsf_fit <- function(object, type = "all", ...) {
  if (!inherits(object, "gsf_fit")) {
    stop("object must be of class 'gsf_fit'")
  }

  summ <- rstan::summary(object$stanfit)$summary

  results <- list()

  if (type %in% c("production", "all")) {
    # Production function parameters
    prod_params <- c("alpha0")
    alpha_params <- paste0("alpha[", 1:object$model_info$J, "]")
    prod_params <- c(prod_params, alpha_params)

    # B matrix parameters
    idx <- 1
    for (i in 1:object$model_info$J) {
      for (j in i:object$model_info$J) {
        prod_params <- c(prod_params, paste0("B_vech[", idx, "]"))
        idx <- idx + 1
      }
    }

    available <- intersect(prod_params, rownames(summ))
    if (length(available) > 0) {
      prod_df <- data.frame(
        parameter = available,
        mean = summ[available, "mean"],
        sd = summ[available, "sd"],
        lower = summ[available, "2.5%"],
        upper = summ[available, "97.5%"],
        stringsAsFactors = FALSE
      )
      results$production <- prod_df
    }
  }

  if (type %in% c("delta", "all")) {
    # Delta coefficients
    if (length(object$model_info$z_slack_names) > 1) {
      delta_coef <- extract_delta_coefficients(object)
      results$delta <- delta_coef
    }
  }

  if (type %in% c("variance", "all")) {
    # Variance parameters
    var_params <- c("sigma", "sigma_u")
    available <- intersect(var_params, rownames(summ))
    if (length(available) > 0) {
      var_df <- data.frame(
        parameter = available,
        mean = summ[available, "mean"],
        sd = summ[available, "sd"],
        lower = summ[available, "2.5%"],
        upper = summ[available, "97.5%"],
        stringsAsFactors = FALSE
      )
      results$variance <- var_df
    }
  }

  if (type == "all") {
    return(results)
  } else if (type == "production") {
    return(results$production)
  } else if (type == "delta") {
    return(results$delta)
  } else if (type == "variance") {
    return(results$variance)
  }
}


#' Predict Method for gsf_fit
#'
#' @param object A gsf_fit object.
#' @param newdata Optional new data for prediction. If NULL, returns fitted values.
#' @param type Character string: "response" for posterior means of observed output,
#'   "frontier" for frontier values (without inefficiency or noise), or "draws"
#'   for posterior draws.
#' @param draws Optional number of posterior draws to use for prediction when
#'   `newdata` is provided. Defaults to all draws.
#' @param seed Optional seed for reproducibility of simulated predictions when
#'   `newdata` is provided.
#' @param ... Additional arguments (ignored).
#'
#' @return A vector or data frame of predictions.
#'
#' @export
predict.gsf_fit <- function(object,
                            newdata = NULL,
                            type = c("response", "frontier", "draws"),
                            draws = NULL,
                            seed = NULL,
                            ...) {
  if (!inherits(object, "gsf_fit")) {
    stop("object must be of class 'gsf_fit'")
  }

  type <- match.arg(type)

  if (is.null(newdata)) {
    if (type == "frontier") {
      frontier_draws <- compute_frontier_draws(object)
      if (!is.null(draws) && draws < nrow(frontier_draws)) {
        keep <- sample(nrow(frontier_draws), draws)
        frontier_draws <- frontier_draws[keep, , drop = FALSE]
      }
      return(colMeans(frontier_draws))
    }

    y_pred <- rstan::extract(object$stanfit, pars = "y_pred")[[1]]
    if (!is.null(draws) && draws < nrow(y_pred)) {
      keep <- sample(seq_len(nrow(y_pred)), draws)
      y_pred <- y_pred[keep, , drop = FALSE]
    }
    if (type == "draws") {
      return(y_pred)
    }
    return(colMeans(y_pred))
  }

  prediction_data <- prepare_prediction_data(object, newdata)
  preds <- simulate_new_predictions(
    object = object,
    pred_data = prediction_data,
    draws = draws,
    type = type,
    seed = seed
  )

  if (type == "draws") {
    return(preds)
  }
  colMeans(preds)
}


#' Log-Likelihood for Model Comparison
#'
#' @param object A gsf_fit object.
#' @param ... Additional arguments (ignored).
#'
#' @return A matrix of log-likelihood values (iterations x observations).
#'
#' @export
logLik.gsf_fit <- function(object, ...) {
  if (!inherits(object, "gsf_fit")) {
    stop("object must be of class 'gsf_fit'")
  }

  log_lik <- rstan::extract(object$stanfit, pars = "log_lik")[[1]]
  return(log_lik)
}


#' WAIC for Model Comparison
#'
#' @description
#' Compute the Widely Applicable Information Criterion (WAIC) for model
#' comparison.
#'
#' @param object A gsf_fit object.
#'
#' @return A list with WAIC value and effective number of parameters.
#'
#' @export
waic <- function(object) {
  if (!inherits(object, "gsf_fit")) {
    stop("object must be of class 'gsf_fit'")
  }

  log_lik <- logLik(object)

  # Compute WAIC components
  # lppd: log pointwise predictive density
  lppd <- sum(log(colMeans(exp(log_lik))))

  # p_waic: effective number of parameters
  p_waic <- sum(apply(log_lik, 2, var))

  # WAIC
  waic_val <- -2 * (lppd - p_waic)

  return(list(
    waic = waic_val,
    lppd = lppd,
    p_waic = p_waic
  ))
}


#' LOO-CV for Model Comparison
#'
#' @description
#' Compute Leave-One-Out Cross-Validation using the loo package.
#'
#' @param object A gsf_fit object.
#'
#' @return A loo object.
#'
#' @export
loo_cv <- function(object) {
  if (!inherits(object, "gsf_fit")) {
    stop("object must be of class 'gsf_fit'")
  }

  if (!requireNamespace("loo", quietly = TRUE)) {
    stop("Package 'loo' is required for LOO-CV. Please install it.")
  }

  log_lik <- logLik(object)
  loo_result <- loo::loo(log_lik)

  return(loo_result)
}
