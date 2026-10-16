.accessor_chlrm_fit <- function() {
  X <- replicate(4, cbind(intercept = 1, slope = c(-1, 0, 1)),
                 simplify = FALSE)
  coefficients <- list(c(-1, 1), c(-.8, .9), c(2, -1), c(2.2, -.8))
  y <- Map(function(Xj, beta) drop(Xj %*% beta), X, coefficients)
  fit_model(
    chlrm_model(2, 2), list(X = X, y = y),
    chlrm_prior(c(0, 0), 10 * diag(2), 4, diag(2), 4, 2, 2, c(1, 1)),
    "cavi", initialization_kmeans(),
    cavi_control(4, evaluation_draws = 4), seed = 760, verbose = FALSE)
}

test_that("CHLRM accessors expose clustering without private-list navigation", {
  fit <- .accessor_chlrm_fit()
  component_parameters <- components(fit)
  expect_named(component_parameters, c(
    "component_coefficients", "component_variances", "mixture_weights"))
  probability <- responsibilities(fit)
  expect_equal(dim(probability), c(4, 2))
  expect_equal(unname(rowSums(probability)), rep(1, 4), tolerance = 1e-10)
  classification <- classify_groups(fit)
  expect_named(classification, c(
    "group", "component", "assignment_probability", "uncertainty",
    "entropy"))
  expect_equal(classification$uncertainty,
               unname(1 - apply(probability, 1, max)))
  expect_equal(classification_uncertainty(fit), classification)
  similarity <- co_clustering(fit)
  expect_equal(unname(diag(similarity)), rep(1, 4))
  expect_equal(similarity, t(similarity))
})

test_that("CHLRM-only accessors reject models without mixture latents", {
  X <- cbind(1, seq(-1, 1, length.out = 8))
  fit <- fit_model(
    lrm_model(2), list(X = X, y = drop(X %*% c(1, 2))),
    lrm_prior(c(0, 0), diag(2), 4, 1), "cavi",
    initialization_ols(), cavi_control(3, evaluation_draws = 3),
    seed = 761, verbose = FALSE)
  expect_error(responsibilities(fit), "only to CHLRM")
  expect_error(components(fit), "only to CHLRM")
  expect_error(co_clustering(fit), "only to CHLRM")
})
