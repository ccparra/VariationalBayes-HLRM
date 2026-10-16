test_that("as_yrep returns the draw by observation orientation expected by bayesplot", {
  X <- cbind(intercept = 1, slope = seq(-1, 1, length.out = 4))
  fit <- fit_model(
    lrm_model(2), list(X = X, y = drop(X %*% c(1, 2))),
    lrm_prior(c(0, 0), 10 * diag(2), 4, 1), "cavi",
    initialization_ols(), cavi_control(3, evaluation_draws = 3),
    seed = 750, verbose = FALSE)
  prediction <- posterior_predict(
    fit, X, n_draws = 3, seed = 751, return_draws = TRUE)
  yrep <- as_yrep(prediction)
  expect_equal(dim(yrep), c(3, 4))
  expect_equal(yrep, t(do.call(rbind, prediction$draws)), ignore_attr = TRUE)
  expect_named(as.data.frame(yrep), c(
    "group_1.observation_1", "group_1.observation_2",
    "group_1.observation_3", "group_1.observation_4"))

  mean_prediction <- prediction
  mean_prediction$target <- "mean"
  expect_error(as_yrep(mean_prediction), "response predictions")
  prediction$draws <- NULL
  expect_error(as_yrep(prediction), "not retained")
})

test_that("PPC yrep and log likelihood adapters transpose exactly once", {
  ppc <- structure(list(
    model = "LRM", method = "CAVI", predictive_target = "observed_data",
    y_rep = matrix(1:12, nrow = 4, ncol = 3)),
    class = c("vb_posterior_predictive_check", "list"))
  expect_equal(as_yrep(ppc), t(ppc$y_rep), ignore_attr = TRUE)
  ppc$y_rep <- NULL
  expect_error(as_yrep(ppc), "not retained")

  likelihood <- log_likelihood(
    "LRM", y = 1:4, X = matrix(1, 4, 1),
    draws = list(Beta = matrix(c(0, 1, 2), ncol = 1),
                 sigma2 = rep(1, 3)),
    unit = "observation", conditioning = "parameter_conditional")
  external <- as_log_lik_matrix(likelihood)
  expect_equal(dim(external), c(3, 4))
  expect_equal(external, t(likelihood$values), ignore_attr = TRUE)
  expect_identical(attr(external, "predictive_units"),
                   likelihood$predictive_units)
  expect_identical(attr(external, "estimand")$unit, "observation")
})
