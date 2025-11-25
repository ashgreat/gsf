#' Extract Technical Inefficiency Estimates
#'
#' @description
#' Extract posterior estimates of output technical inefficiency from a fitted
#' GSF model.
#'
#' @param object A gsf_fit object.
#' @param type Character string: "mean" for posterior means, "median" for
#'   posterior medians, or "samples" for all posterior samples.
#' @param interval Numeric between 0 and 1, credible interval level (default: 0.95).
#'
#' @return A data frame with inefficiency estimates for each observation.
#'
#' @export
extract_technical_inefficiency <- function(object, type = "mean", interval = 0.95) {
  if (!inherits(object, "gsf_fit")) {
    stop("object must be of class 'gsf_fit'")
  }

  # Extract samples
  samples <- rstan::extract(object$stanfit, pars = "technical_inefficiency")[[1]]

  if (type == "samples") {
    return(samples)
  }

  # Compute summary statistics
  N <- ncol(samples)
  alpha <- (1 - interval) / 2

  result <- data.frame(
    observation = 1:N,
    mean = colMeans(samples),
    median = apply(samples, 2, median),
    sd = apply(samples, 2, sd),
    lower = apply(samples, 2, quantile, probs = alpha),
    upper = apply(samples, 2, quantile, probs = 1 - alpha)
  )

  if (type == "mean") {
    result$estimate <- result$mean
  } else if (type == "median") {
    result$estimate <- result$median
  }

  return(result)
}


#' Extract Input-Specific Slack Estimates
#'
#' @description
#' Extract posterior estimates of input-specific slacks (inefficiencies) from
#' a fitted GSF model.
#'
#' @param object A gsf_fit object.
#' @param input Character string or integer specifying which input's slack to
#'   extract. If NULL (default), returns slacks for all inputs.
#' @param type Character string: "mean" for posterior means, "median" for
#'   posterior medians, or "samples" for all posterior samples.
#' @param interval Numeric between 0 and 1, credible interval level (default: 0.95).
#'
#' @return A data frame with slack estimates, or a 3D array if type = "samples".
#'
#' @export
extract_input_slacks <- function(object, input = NULL, type = "mean", interval = 0.95) {
  if (!inherits(object, "gsf_fit")) {
    stop("object must be of class 'gsf_fit'")
  }

  # Extract samples (iterations x N x J)
  samples <- rstan::extract(object$stanfit, pars = "input_slacks")[[1]]

  if (type == "samples") {
    if (!is.null(input)) {
      if (is.character(input)) {
        input_idx <- which(object$model_info$input_names == input)
        if (length(input_idx) == 0) {
          stop("Input '", input, "' not found. Available: ",
               paste(object$model_info$input_names, collapse = ", "))
        }
      } else {
        input_idx <- input
      }
      return(samples[, , input_idx])
    }
    return(samples)
  }

  # Compute summary statistics
  N <- dim(samples)[2]
  J <- dim(samples)[3]
  alpha <- (1 - interval) / 2

  if (!is.null(input)) {
    if (is.character(input)) {
      input_idx <- which(object$model_info$input_names == input)
      if (length(input_idx) == 0) {
        stop("Input '", input, "' not found. Available: ",
             paste(object$model_info$input_names, collapse = ", "))
      }
    } else {
      input_idx <- input
    }

    slack_samples <- samples[, , input_idx]
    result <- data.frame(
      observation = 1:N,
      input = object$model_info$input_names[input_idx],
      mean = colMeans(slack_samples),
      median = apply(slack_samples, 2, median),
      sd = apply(slack_samples, 2, sd),
      lower = apply(slack_samples, 2, quantile, probs = alpha),
      upper = apply(slack_samples, 2, quantile, probs = 1 - alpha)
    )
  } else {
    # All inputs
    results_list <- list()
    for (j in 1:J) {
      slack_samples <- samples[, , j]
      results_list[[j]] <- data.frame(
        observation = 1:N,
        input = object$model_info$input_names[j],
        mean = colMeans(slack_samples),
        median = apply(slack_samples, 2, median),
        sd = apply(slack_samples, 2, sd),
        lower = apply(slack_samples, 2, quantile, probs = alpha),
        upper = apply(slack_samples, 2, quantile, probs = 1 - alpha)
      )
    }
    result <- do.call(rbind, results_list)
  }

  if (type == "mean") {
    result$estimate <- result$mean
  } else if (type == "median") {
    result$estimate <- result$median
  }

  return(result)
}


#' Extract Output Loss Estimates
#'
#' @description
#' Extract estimates of output loss due to technical inefficiency and/or
#' input slacks.
#'
#' @param object A gsf_fit object.
#' @param component Character string: "total" for total loss, "slacks" for
#'   loss from input slacks only, or "technical" for output technical
#'   inefficiency only.
#' @param type Character string: "mean" for posterior means, or "samples"
#'   for all posterior samples.
#'
#' @return A data frame or matrix with output loss estimates.
#'
#' @export
extract_output_loss <- function(object, component = "total", type = "mean") {
  if (!inherits(object, "gsf_fit")) {
    stop("object must be of class 'gsf_fit'")
  }

  if (component == "total") {
    samples <- rstan::extract(object$stanfit, pars = "total_output_loss")[[1]]
  } else if (component == "slacks") {
    samples <- rstan::extract(object$stanfit, pars = "output_loss_from_slacks")[[1]]
  } else if (component == "technical") {
    samples <- rstan::extract(object$stanfit, pars = "technical_inefficiency")[[1]]
  } else {
    stop("component must be 'total', 'slacks', or 'technical'")
  }

  if (type == "samples") {
    return(samples)
  }

  N <- ncol(samples)
  result <- data.frame(
    observation = 1:N,
    mean = colMeans(samples),
    median = apply(samples, 2, median),
    sd = apply(samples, 2, sd),
    lower = apply(samples, 2, quantile, probs = 0.025),
    upper = apply(samples, 2, quantile, probs = 0.975)
  )

  return(result)
}


#' Extract Input Elasticities
#'
#' @description
#' Extract posterior estimates of input elasticities from the translog
#' production function.
#'
#' @param object A gsf_fit object.
#' @param type Character string: "mean" for posterior means, or "samples"
#'   for all posterior samples.
#'
#' @return A data frame with elasticity estimates for each observation and input.
#'
#' @export
extract_elasticities <- function(object, type = "mean") {
  if (!inherits(object, "gsf_fit")) {
    stop("object must be of class 'gsf_fit'")
  }

  # Extract samples (iterations x N x J)
  samples <- rstan::extract(object$stanfit, pars = "elasticities")[[1]]

  if (type == "samples") {
    return(samples)
  }

  N <- dim(samples)[2]
  J <- dim(samples)[3]

  results_list <- list()
  for (j in 1:J) {
    elas_samples <- samples[, , j]
    results_list[[j]] <- data.frame(
      observation = 1:N,
      input = object$model_info$input_names[j],
      mean = colMeans(elas_samples),
      median = apply(elas_samples, 2, median),
      sd = apply(elas_samples, 2, sd),
      lower = apply(elas_samples, 2, quantile, probs = 0.025),
      upper = apply(elas_samples, 2, quantile, probs = 0.975)
    )
  }

  result <- do.call(rbind, results_list)
  return(result)
}


#' Extract Returns to Scale
#'
#' @description
#' Extract posterior estimates of returns to scale from the translog
#' production function.
#'
#' @param object A gsf_fit object.
#' @param type Character string: "mean" for posterior means, or "samples"
#'   for all posterior samples.
#'
#' @return A data frame with RTS estimates for each observation.
#'
#' @export
extract_returns_to_scale <- function(object, type = "mean") {
  if (!inherits(object, "gsf_fit")) {
    stop("object must be of class 'gsf_fit'")
  }

  samples <- rstan::extract(object$stanfit, pars = "returns_to_scale")[[1]]

  if (type == "samples") {
    return(samples)
  }

  N <- ncol(samples)
  result <- data.frame(
    observation = 1:N,
    mean = colMeans(samples),
    median = apply(samples, 2, median),
    sd = apply(samples, 2, sd),
    lower = apply(samples, 2, quantile, probs = 0.025),
    upper = apply(samples, 2, quantile, probs = 0.975)
  )

  return(result)
}


#' Get Efficiency Summary Table
#'
#' @description
#' Creates a summary table of all efficiency measures similar to Table 2
#' in Kumbhakar and Tsionas (2021).
#'
#' @param object A gsf_fit object.
#'
#' @return A data frame with summary statistics for efficiency measures.
#'
#' @export
efficiency_summary <- function(object) {
  if (!inherits(object, "gsf_fit")) {
    stop("object must be of class 'gsf_fit'")
  }

  # Technical inefficiency
  ti <- extract_technical_inefficiency(object)

  # Input slacks
  slacks <- extract_input_slacks(object)

  # Elasticities
  elas <- extract_elasticities(object)

  # Returns to scale
  rts <- extract_returns_to_scale(object)

  # Create summary table
  summary_list <- list(
    data.frame(
      measure = "Technical Inefficiency (u0)",
      mean = mean(ti$mean),
      sd = sd(ti$mean)
    )
  )

  # Add input slacks
  for (input_name in unique(slacks$input)) {
    slack_input <- slacks[slacks$input == input_name, ]
    summary_list[[length(summary_list) + 1]] <- data.frame(
      measure = paste0("Slack (", input_name, ")"),
      mean = mean(slack_input$mean),
      sd = sd(slack_input$mean)
    )
  }

  # Add elasticities
  for (input_name in unique(elas$input)) {
    elas_input <- elas[elas$input == input_name, ]
    summary_list[[length(summary_list) + 1]] <- data.frame(
      measure = paste0("Elasticity (", input_name, ")"),
      mean = mean(elas_input$mean),
      sd = sd(elas_input$mean)
    )
  }

  # Add returns to scale
  summary_list[[length(summary_list) + 1]] <- data.frame(
    measure = "Returns to Scale",
    mean = mean(rts$mean),
    sd = sd(rts$mean)
  )

  result <- do.call(rbind, summary_list)
  rownames(result) <- NULL

  return(result)
}


#' Extract Delta Coefficients (Slack Determinants)
#'
#' @description
#' Extract posterior estimates of the Delta coefficients that determine
#' input-specific slacks.
#'
#' @param object A gsf_fit object.
#'
#' @return A data frame with coefficient estimates.
#'
#' @export
extract_delta_coefficients <- function(object) {
  if (!inherits(object, "gsf_fit")) {
    stop("object must be of class 'gsf_fit'")
  }

  # Extract Delta matrix samples
  samples <- rstan::extract(object$stanfit, pars = "Delta_raw")[[1]]

  # samples is (iterations x J x M)
  J <- dim(samples)[2]
  M <- dim(samples)[3]

  results_list <- list()
  for (j in 1:J) {
    for (m in 1:M) {
      delta_samples <- samples[, j, m]
      results_list[[length(results_list) + 1]] <- data.frame(
        input = object$model_info$input_names[j],
        variable = object$model_info$z_slack_names[m],
        mean = mean(delta_samples),
        sd = sd(delta_samples),
        lower = quantile(delta_samples, 0.025),
        upper = quantile(delta_samples, 0.975)
      )
    }
  }

  result <- do.call(rbind, results_list)
  rownames(result) <- NULL
  return(result)
}
