#' Convert vech vector to symmetric matrix
#' @keywords internal
#' @noRd
vech_to_matrix <- function(vech, J) {
  mat <- matrix(0, nrow = J, ncol = J)
  idx <- 1
  for (i in seq_len(J)) {
    for (j in i:J) {
      mat[i, j] <- vech[idx]
      mat[j, i] <- vech[idx]
      idx <- idx + 1
    }
  }
  mat
}

#' Draw from a multivariate normal truncated at zero from above
#' @keywords internal
#' @noRd
sample_negative_mvn <- function(mu, L, max_attempts = 2000) {
  J <- length(mu)
  for (attempt in seq_len(max_attempts)) {
    z <- stats::rnorm(J)
    candidate <- as.vector(mu + L %*% z)
    if (all(candidate <= 0)) {
      return(candidate)
    }
  }
  warning("Failed to draw truncated slacks after ", max_attempts,
          " attempts. Using deterministic clipping.")
  pmin(mu, -1e-6)
}

#' Prepare design matrices for prediction
#' @keywords internal
#' @noRd
prepare_prediction_data <- function(object, newdata) {
  rhs <- stats::delete.response(object$model_info$terms)
  mf <- stats::model.frame(rhs, data = newdata, na.action = stats::na.fail)
  X <- stats::model.matrix(rhs, mf)
  if ("(Intercept)" %in% colnames(X)) {
    X <- X[, -1, drop = FALSE]
  }

  # Reorder to match training columns
  X <- X[, object$model_info$input_names, drop = FALSE]

  # Slack determinants
  slack_vars <- object$model_info$z_slack_vars
  if (length(slack_vars) == 0) {
    Z_slack <- matrix(1, nrow = nrow(X), ncol = 1)
  } else {
    missing <- setdiff(slack_vars, names(newdata))
    if (length(missing) > 0) {
      stop("Variables not found in newdata: ", paste(missing, collapse = ", "))
    }
    Z_slack <- cbind(1, as.matrix(newdata[, slack_vars, drop = FALSE]))
  }

  # Neutral shifters
  neutral_vars <- object$model_info$neutral_vars
  if (length(neutral_vars) == 0) {
    Z_neutral <- matrix(0, nrow = nrow(X), ncol = 0)
  } else {
    missing_n <- setdiff(neutral_vars, names(newdata))
    if (length(missing_n) > 0) {
      stop("Neutral variables missing from newdata: ",
           paste(missing_n, collapse = ", "))
    }
    Z_neutral <- as.matrix(newdata[, neutral_vars, drop = FALSE])
  }

  list(
    X = X,
    Z_slack = Z_slack,
    Z_neutral = Z_neutral,
    N = nrow(X)
  )
}

#' Simulate posterior predictions for new data
#' @keywords internal
#' @noRd
simulate_new_predictions <- function(object,
                                     pred_data,
                                     draws = NULL,
                                     type = c("response", "frontier", "draws"),
                                     seed = NULL) {
  type <- match.arg(type)
  if (!is.null(seed)) {
    set.seed(seed)
  }

  samples <- rstan::extract(object$stanfit, pars = c(
    "alpha0", "alpha", "B_vech", "delta0",
    "Delta_raw", "sigma", "sigma_u",
    "sigma_slack", "L_Omega"
  ))

  total_draws <- length(samples$alpha0)
  draw_ids <- if (is.null(draws) || draws >= total_draws) {
    seq_len(total_draws)
  } else {
    sample(total_draws, draws)
  }

  n_sims <- length(draw_ids)
  N_new <- pred_data$N
  preds <- matrix(NA_real_, nrow = n_sims, ncol = N_new)

  for (s in seq_len(n_sims)) {
    idx <- draw_ids[s]
    alpha0 <- samples$alpha0[idx]
    alpha <- samples$alpha[idx, ]
    B <- vech_to_matrix(samples$B_vech[idx, ], object$model_info$J)

    delta_draw <- samples$delta0
    if (length(delta_draw)) {
      delta_vec <- if (is.matrix(delta_draw)) {
        delta_draw[idx, ]
      } else {
        delta_draw[idx]
      }
    } else {
      delta_vec <- numeric(0)
    }

    Delta_mat <- samples$Delta_raw[idx, , ]
    sigma <- samples$sigma[idx]
    sigma_u <- samples$sigma_u[idx]
    sigma_slack <- samples$sigma_slack[idx, ]
    L_Omega <- samples$L_Omega[idx, , ]
    L_Sigma <- diag(sigma_slack) %*% L_Omega

    for (n in seq_len(N_new)) {
      mu_theta <- as.vector(Delta_mat %*% pred_data$Z_slack[n, ])
      theta_draw <- sample_negative_mvn(mu_theta, L_Sigma)
      x_eff <- pred_data$X[n, ] + theta_draw

      mu_y <- alpha0 +
        sum(alpha * x_eff) +
        0.5 * as.numeric(t(x_eff) %*% B %*% x_eff)

      if (ncol(pred_data$Z_neutral) > 0 && length(delta_vec) > 0) {
        mu_y <- mu_y + sum(delta_vec * pred_data$Z_neutral[n, ])
      }

      if (type == "frontier") {
        preds[s, n] <- mu_y
      } else {
        u0 <- -abs(stats::rnorm(1, mean = 0, sd = sigma_u))
        eps <- stats::rnorm(1, mean = 0, sd = sigma)
        preds[s, n] <- mu_y + u0 + eps
      }
    }
  }

  preds
}

#' Frontier draws for observed data
#' @keywords internal
#' @noRd
compute_frontier_draws <- function(object) {
  samples <- rstan::extract(object$stanfit, pars = c(
    "alpha0", "alpha", "B_vech", "delta0",
    "theta"
  ))

  X <- object$data$X
  Z_neutral <- object$data$Z_neutral
  n_draws <- length(samples$alpha0)
  N <- nrow(X)
  preds <- matrix(NA_real_, nrow = n_draws, ncol = N)

  for (s in seq_len(n_draws)) {
    alpha0 <- samples$alpha0[s]
    alpha <- samples$alpha[s, ]
    B <- vech_to_matrix(samples$B_vech[s, ], object$model_info$J)
    theta <- samples$theta[s, , ]

    delta_draw <- samples$delta0
    if (length(delta_draw)) {
      delta_vec <- if (is.matrix(delta_draw)) {
        delta_draw[s, ]
      } else {
        delta_draw[s]
      }
    } else {
      delta_vec <- numeric(0)
    }

    for (n in seq_len(N)) {
      x_eff <- X[n, ] + theta[n, ]
      mu <- alpha0 +
        sum(alpha * x_eff) +
        0.5 * as.numeric(t(x_eff) %*% B %*% x_eff)
      if (ncol(Z_neutral) > 0 && length(delta_vec) > 0) {
        mu <- mu + sum(delta_vec * Z_neutral[n, ])
      }
      preds[s, n] <- mu
    }
  }

  preds
}
