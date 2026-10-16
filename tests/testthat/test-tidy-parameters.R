test_that("tidy_parameters represents every applicable parameter exactly once", {
  X <- cbind(intercept = 1, slope = seq(-1, 1, length.out = 10))
  fit <- fit_model(
    lrm_model(2), list(X = X, y = drop(X %*% c(1, 2))),
    lrm_prior(c(0, 0), diag(2), 4, 1), "cavi",
    initialization_ols(), cavi_control(4, evaluation_draws = 4),
    seed = 780, verbose = FALSE)
  from_fit <- tidy_parameters(fit)
  from_summary <- tidy_parameters(summary(fit))
  expect_s3_class(from_fit, "vb_tidy_parameters")
  expect_identical(from_fit, from_summary)
  expect_named(from_fit, c(
    "model", "method", "parameter", "index_1", "index_2", "estimate",
    "sd", "lower", "upper", "mean_exists", "contains_zero",
    "distribution", "interval_level"))
  counts <- table(from_fit$parameter)
  expect_equal(unname(as.integer(counts)), c(2, 1))
  expect_equal(names(counts), c("coefficients", "residual_variance"))
  expect_true(all(from_fit$model == "LRM"))
  expect_true(all(from_fit$method == "CAVI"))
})

test_that("tidy_parameters keeps both indices of matrix-valued summaries", {
  X <- replicate(3, cbind(1, c(-1, 0, 1)), simplify = FALSE)
  y <- lapply(seq_along(X), function(group)
    drop(X[[group]] %*% c(group, (-1)^group)))
  fit <- fit_model(
    hlrm_model(2), list(X = X, y = y),
    hlrm_prior(c(0, 0), diag(2), 4, diag(2), 4, 2, 2), "cavi",
    initialization_ols(), cavi_control(4, evaluation_draws = 4),
    seed = 781, verbose = FALSE)
  tidy <- tidy_parameters(fit)
  covariance <- tidy[tidy$parameter == "population_covariance", ]
  groups <- tidy[tidy$parameter == "group_coefficients", ]
  expect_equal(nrow(covariance), 4)
  expect_true(all(!is.na(covariance$index_2)))
  expect_equal(nrow(groups), 6)
  expect_equal(sort(unique(groups$index_1)), 1:3)
  expect_equal(sort(unique(groups$index_2)), 1:2)
})
