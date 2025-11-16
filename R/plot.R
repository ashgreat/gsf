#' Plot Method for gsf_fit Objects
#'
#' @description
#' Creates various diagnostic and summary plots for the fitted GSF model.
#'
#' @param x A gsf_fit object.
#' @param type Character string specifying the plot type:
#'   \describe{
#'     \item{"efficiency"}{Distribution of efficiency measures (default)}
#'     \item{"slacks"}{Input-specific slack distributions}
#'     \item{"output_loss"}{Output loss distributions}
#'     \item{"elasticities"}{Input elasticity distributions}
#'     \item{"rts"}{Returns to scale distribution}
#'     \item{"trace"}{MCMC trace plots for key parameters}
#'     \item{"pairs"}{Pairs plot for efficiency measures}
#'   }
#' @param ... Additional arguments passed to plotting functions.
#'
#' @return A ggplot object.
#'
#' @export
plot.gsf_fit <- function(x, type = "efficiency", ...) {
  if (!inherits(x, "gsf_fit")) {
    stop("x must be of class 'gsf_fit'")
  }

  if (type == "efficiency") {
    plot_efficiency_distributions(x, ...)
  } else if (type == "slacks") {
    plot_input_slacks(x, ...)
  } else if (type == "output_loss") {
    plot_output_loss(x, ...)
  } else if (type == "elasticities") {
    plot_elasticities(x, ...)
  } else if (type == "rts") {
    plot_returns_to_scale(x, ...)
  } else if (type == "trace") {
    plot_trace(x, ...)
  } else if (type == "pairs") {
    plot_efficiency_pairs(x, ...)
  } else {
    stop("Unknown plot type: ", type)
  }
}


#' Plot Efficiency Distributions
#'
#' @param object A gsf_fit object.
#' @param ... Additional arguments.
#'
#' @return A ggplot object.
#'
#' @keywords internal
plot_efficiency_distributions <- function(object, ...) {
  # Extract efficiency measures
  ti <- extract_technical_inefficiency(object)
  slacks <- extract_input_slacks(object)

  # Prepare data for plotting
  plot_data <- data.frame(
    value = c(ti$mean, slacks$mean),
    type = c(
      rep("Technical Inefficiency", nrow(ti)),
      rep(paste0("Slack (", slacks$input, ")"), nrow(slacks))
    )
  )

  # Create plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = value, fill = type)) +
    ggplot2::geom_density(alpha = 0.6) +
    ggplot2::labs(
      title = "Distribution of Efficiency Measures",
      subtitle = "Technical inefficiency and input-specific slacks (%)",
      x = "Inefficiency/Slack (%)",
      y = "Density",
      fill = "Measure"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "bottom")

  return(p)
}


#' Plot Input Slacks
#'
#' @param object A gsf_fit object.
#' @param ... Additional arguments.
#'
#' @return A ggplot object.
#'
#' @keywords internal
plot_input_slacks <- function(object, ...) {
  slacks <- extract_input_slacks(object)

  p <- ggplot2::ggplot(slacks, ggplot2::aes(x = mean, fill = input)) +
    ggplot2::geom_density(alpha = 0.6) +
    ggplot2::facet_wrap(~ input, scales = "free_y") +
    ggplot2::labs(
      title = "Input-Specific Slack Distributions",
      subtitle = "Percentage over-use of each input",
      x = "Slack (%)",
      y = "Density"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "none")

  return(p)
}


#' Plot Output Loss
#'
#' @param object A gsf_fit object.
#' @param ... Additional arguments.
#'
#' @return A ggplot object.
#'
#' @keywords internal
plot_output_loss <- function(object, ...) {
  # Extract output loss components
  ti <- extract_technical_inefficiency(object)
  slack_loss <- extract_output_loss(object, component = "slacks")

  plot_data <- data.frame(
    value = c(ti$mean, slack_loss$mean),
    type = c(
      rep("Technical Inefficiency", nrow(ti)),
      rep("Output Loss from Slacks", nrow(slack_loss))
    )
  )

  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = value, fill = type)) +
    ggplot2::geom_density(alpha = 0.6) +
    ggplot2::labs(
      title = "Output Loss Distributions",
      subtitle = "Comparison of technical inefficiency and loss from input slacks",
      x = "Output Loss (%)",
      y = "Density",
      fill = "Component"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "bottom")

  return(p)
}


#' Plot Input Elasticities
#'
#' @param object A gsf_fit object.
#' @param ... Additional arguments.
#'
#' @return A ggplot object.
#'
#' @keywords internal
plot_elasticities <- function(object, ...) {
  elas <- extract_elasticities(object)

  p <- ggplot2::ggplot(elas, ggplot2::aes(x = mean, fill = input)) +
    ggplot2::geom_density(alpha = 0.6) +
    ggplot2::facet_wrap(~ input, scales = "free") +
    ggplot2::labs(
      title = "Input Elasticity Distributions",
      subtitle = "Output elasticities for each input",
      x = "Elasticity",
      y = "Density"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "none")

  return(p)
}


#' Plot Returns to Scale
#'
#' @param object A gsf_fit object.
#' @param ... Additional arguments.
#'
#' @return A ggplot object.
#'
#' @keywords internal
plot_returns_to_scale <- function(object, ...) {
  rts <- extract_returns_to_scale(object)

  p <- ggplot2::ggplot(rts, ggplot2::aes(x = mean)) +
    ggplot2::geom_density(fill = "steelblue", alpha = 0.6) +
    ggplot2::geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
    ggplot2::labs(
      title = "Returns to Scale Distribution",
      subtitle = "Dashed line indicates constant returns to scale",
      x = "Returns to Scale",
      y = "Density"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::annotate(
      "text",
      x = 1.05,
      y = 0,
      label = "CRS",
      color = "red",
      hjust = 0
    )

  return(p)
}


#' Plot MCMC Trace Plots
#'
#' @param object A gsf_fit object.
#' @param pars Character vector of parameter names to plot.
#' @param ... Additional arguments.
#'
#' @return A plot.
#'
#' @keywords internal
plot_trace <- function(object, pars = NULL, ...) {
  if (is.null(pars)) {
    # Default parameters
    pars <- c("alpha0", "sigma", "sigma_u")
    # Add first two alpha parameters
    pars <- c(pars, paste0("alpha[", 1:min(2, object$model_info$J), "]"))
  }

  rstan::traceplot(object$stanfit, pars = pars, inc_warmup = FALSE)
}


#' Plot Efficiency Pairs
#'
#' @param object A gsf_fit object.
#' @param ... Additional arguments.
#'
#' @return A ggplot object.
#'
#' @keywords internal
plot_efficiency_pairs <- function(object, ...) {
  # Extract efficiency measures
  ti <- extract_technical_inefficiency(object)
  slacks <- extract_input_slacks(object)

  # Reshape slacks to wide format
  slacks_wide <- tidyr::pivot_wider(
    slacks[, c("observation", "input", "mean")],
    names_from = input,
    values_from = mean,
    names_prefix = "Slack_"
  )

  # Combine
  plot_data <- merge(
    data.frame(observation = ti$observation, TI = ti$mean),
    slacks_wide,
    by = "observation"
  )

  # Create pairs plot
  if (ncol(plot_data) > 2) {
    if ("GGally" %in% loadedNamespaces()) {
      ggpairs <- getExportedValue("GGally", "ggpairs")
      p <- ggpairs(
        plot_data[, -1],
        title = "Relationships Between Efficiency Measures",
        lower = list(continuous = "points"),
        diag = list(continuous = "densityDiag"),
        upper = list(continuous = "cor")
      )
    } else {
      warning("Package 'GGally' not loaded; using base pairs plot instead.")
      graphics::pairs(plot_data[, -1, drop = FALSE])
      return(invisible(NULL))
    }
  } else {
    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = TI)) +
      ggplot2::geom_density(fill = "steelblue", alpha = 0.6) +
      ggplot2::theme_minimal()
  }

  return(p)
}


#' Plot Marginal Effects of Z Variables
#'
#' @description
#' Plots the marginal effects of z variables on input slacks.
#'
#' @param object A gsf_fit object.
#' @param ... Additional arguments.
#'
#' @return A ggplot object.
#'
#' @export
plot_marginal_effects <- function(object, ...) {
  if (!inherits(object, "gsf_fit")) {
    stop("object must be of class 'gsf_fit'")
  }

  delta_coef <- extract_delta_coefficients(object)

  # Remove intercept
  delta_coef <- delta_coef[delta_coef$variable != "intercept", ]

  if (nrow(delta_coef) == 0) {
    message("No z variables in the model.")
    return(invisible(NULL))
  }

  p <- ggplot2::ggplot(
    delta_coef,
    ggplot2::aes(
      x = variable,
      y = mean,
      ymin = lower,
      ymax = upper,
      color = input
    )
  ) +
    ggplot2::geom_pointrange(position = ggplot2::position_dodge(width = 0.5)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title = "Effects of Z Variables on Input Slacks",
      subtitle = "Posterior means with 95% credible intervals",
      x = "Variable",
      y = "Coefficient",
      color = "Input"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "bottom")

  return(p)
}


#' Compare Best and Worst Firms
#'
#' @description
#' Creates comparison plots for the most and least efficient firms,
#' similar to Figure 3 in Kumbhakar and Tsionas (2021).
#'
#' @param object A gsf_fit object.
#' @param n_firms Number of best/worst firms to compare (default: 5).
#' @param ... Additional arguments.
#'
#' @return A list of ggplot objects.
#'
#' @export
plot_best_worst_firms <- function(object, n_firms = 5, ...) {
  if (!inherits(object, "gsf_fit")) {
    stop("object must be of class 'gsf_fit'")
  }

  # Get total output loss
  total_loss <- extract_output_loss(object, component = "total")

  # Identify best and worst firms
  worst_idx <- order(total_loss$mean, decreasing = TRUE)[1:n_firms]
  best_idx <- order(total_loss$mean, decreasing = FALSE)[1:n_firms]

  # Extract samples for these firms
  ti_samples <- rstan::extract(object$stanfit, pars = "technical_inefficiency")[[1]]
  slack_samples <- rstan::extract(object$stanfit, pars = "input_slacks")[[1]]

  # Prepare data for worst firms
  worst_ti <- data.frame(
    iteration = rep(1:nrow(ti_samples), n_firms),
    firm = rep(paste("Firm", worst_idx), each = nrow(ti_samples)),
    value = as.vector(ti_samples[, worst_idx])
  )

  # Plot worst firms TI
  p_worst <- ggplot2::ggplot(worst_ti, ggplot2::aes(x = value, fill = firm)) +
    ggplot2::geom_density(alpha = 0.5) +
    ggplot2::labs(
      title = paste("Technical Inefficiency:", n_firms, "Worst Firms"),
      x = "Technical Inefficiency (%)",
      y = "Posterior Density"
    ) +
    ggplot2::theme_minimal()

  # Prepare data for best firms
  best_ti <- data.frame(
    iteration = rep(1:nrow(ti_samples), n_firms),
    firm = rep(paste("Firm", best_idx), each = nrow(ti_samples)),
    value = as.vector(ti_samples[, best_idx])
  )

  # Plot best firms TI
  p_best <- ggplot2::ggplot(best_ti, ggplot2::aes(x = value, fill = firm)) +
    ggplot2::geom_density(alpha = 0.5) +
    ggplot2::labs(
      title = paste("Technical Inefficiency:", n_firms, "Best Firms"),
      x = "Technical Inefficiency (%)",
      y = "Posterior Density"
    ) +
    ggplot2::theme_minimal()

  return(list(worst = p_worst, best = p_best))
}
