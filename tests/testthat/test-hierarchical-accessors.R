.accessor_hlrm_fit <- function() {
  X <- replicate(3, cbind(intercept = 1, slope = c(-1, 0, 1)),
                 simplify = FALSE)
  coefficients <- list(c(-1, 1), c(0, .5), c(2, -1))
  y <- Map(function(Xj, beta) drop(Xj %*% beta), X, coefficients)
  fit_model(
    hlrm_model(2), list(X = X, y = y),
    hlrm_prior(c(0, 0), 10 * diag(2), 4, diag(2), 4, 2, 2),
    "cavi", initialization_ols(),
    cavi_control(4, evaluation_draws = 4), seed = 770, verbose = FALSE)
}

test_that("hierarchical accessors preserve population and group responsibilities", {
  fit <- .accessor_hlrm_fit()
  population <- population_parameters(fit)
  groups <- group_parameters(fit)
  expect_named(population, c(
    "population_coefficients", "population_covariance", "variance_scale"))
  expect_named(groups, c("group_coefficients", "group_variances"))
  deviations <- group_deviations(fit)
  expect_equal(dim(deviations), c(3, 2))
  expect_named(as.data.frame(deviations), c("intercept", "slope"))
  expected <- sweep(
    groups$group_coefficients$mean, 2L,
    population$population_coefficients$mean, `-`)
  expect_equal(deviations, expected, ignore_attr = TRUE)
  expect_match(attr(deviations, "estimand"), "posterior mean")
})

test_that("hierarchical accessors reject nonapplicable model families", {
  X <- cbind(1, seq(-1, 1, length.out = 8))
  fit <- fit_model(
    lrm_model(2), list(X = X, y = drop(X %*% c(1, 2))),
    lrm_prior(c(0, 0), diag(2), 4, 1), "cavi",
    initialization_ols(), cavi_control(3, evaluation_draws = 3),
    seed = 771, verbose = FALSE)
  expect_error(population_parameters(fit), "only to HLRM or CHLRM")
  expect_error(group_parameters(fit), "only to HLRM")
  expect_error(group_deviations(fit), "only to HLRM")
})
