.idiomatic_method_fits <- function() {
  X <- cbind(intercept = 1, slope = seq(-1, 1, length.out = 12))
  lrm_data <- list(X = X, y = drop(X %*% c(1, 2)))
  group_rows <- split(seq_len(12), rep(seq_len(4), each = 3))
  group_coefficients <- list(c(0, 1), c(.2, .8), c(2, -1), c(1.8, -.8))
  grouped_data <- list(
    X = lapply(group_rows, function(rows) X[rows, , drop = FALSE]),
    y = Map(function(rows, coefficients)
      drop(X[rows, , drop = FALSE] %*% coefficients),
      group_rows, group_coefficients))
  control <- cavi_control(
    max_iterations = 4, evaluation_draws = 4,
    stopping = cavi_stopping(
      objective_tolerance = 1e20, parameter_tolerance = 1e20,
      min_iterations = 2, patience = 1))
  list(
    LRM = fit_model(
      lrm_model(2), lrm_data, lrm_prior(c(0, 0), 10 * diag(2), 4, 1),
      "cavi", initialization_ols(), control, seed = 740, verbose = FALSE),
    HLRM = fit_model(
      hlrm_model(2), grouped_data,
      hlrm_prior(c(0, 0), 10 * diag(2), 4, diag(2), 4, 2, 2),
      "cavi", initialization_ols(), control, seed = 741, verbose = FALSE),
    CHLRM = fit_model(
      chlrm_model(2, 2), grouped_data,
      chlrm_prior(c(0, 0), 10 * diag(2), 4, diag(2), 4, 2, 2, c(1, 1)),
      "cavi", initialization_kmeans(), control, seed = 742, verbose = FALSE))
}

test_that("standard fitted-model methods expose only applicable coefficients", {
  fits <- .idiomatic_method_fits()
  expect_named(coef(fits$LRM), c("intercept", "slope"))
  expect_named(coef(fits$HLRM),
               c("population_coefficients", "group_coefficients"))
  expect_named(coef(fits$CHLRM),
               c("population_coefficients", "component_coefficients"))
  expect_equal(dim(coef(fits$HLRM)$group_coefficients), c(4, 2))
  expect_equal(dim(coef(fits$CHLRM)$component_coefficients), c(2, 2))
  expect_false("mixture_weights" %in% names(coef(fits$LRM)))
  expect_s3_class(summary(fits$LRM), "vb_posterior_summary")
})

test_that("fitted values, residuals and nobs follow input observation order", {
  fits <- .idiomatic_method_fits()
  for (fit in fits) {
    expect_length(fitted(fit), 12)
    expect_length(residuals(fit), 12)
    expect_identical(nobs(fit), 12L)
    expect_true(all(is.finite(fitted(fit))))
    expect_true(all(is.finite(residuals(fit))))
  }
  expect_match(capture.output(print(fits$LRM))[1L], "Bayesian LRM fit")
})
