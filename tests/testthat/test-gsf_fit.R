test_that("gsf_fit runs on small simulated data", {
  skip_on_cran()
  skip_if_not_installed("rstan")

  set.seed(123)
  toy <- simulate_gsf_data(
    n_firms = 3,
    n_periods = 2,
    seed = 123
  )

  fit <- gsf_fit(
    formula = log(Y) ~ log(L) + log(K),
    data = toy,
    z_vars = c("CR", "MKSH"),
    neutral_vars = NULL,
    chains = 1,
    iter = 400,
    warmup = 200,
    thin = 1,
    cores = 1,
    seed = 2024
  )

  expect_s3_class(fit, "gsf_fit")
  expect_true(is.list(efficiency_summary(fit)))

  preds_train <- predict(fit)
  expect_length(preds_train, nrow(fit$data$X))

  new_obs <- toy[1:2, ]
  preds_new <- predict(fit, newdata = new_obs, draws = 10, seed = 99)
  expect_length(preds_new, nrow(new_obs))
})
