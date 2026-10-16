test_that("LRM scaled minibatch sufficient statistics are unbiased", {
  X <- cbind(1, c(-2, -1, .5, 3))
  y <- c(-1, 0, 2, 5)
  stats <- list(X = X, y = y, n = length(y),
                XtX = crossprod(X), Xty = crossprod(X, y),
                yty = drop(crossprod(y)))
  batches <- combn(seq_along(y), 2, simplify = FALSE)
  estimates <- lapply(batches, function(B)
    .lrm_scaled_sufficient_statistics(stats, B))
  expect_equal(Reduce(`+`, lapply(estimates, `[[`, "XtX")) / length(batches),
               stats$XtX)
  expect_equal(Reduce(`+`, lapply(estimates, `[[`, "Xty")) / length(batches),
               stats$Xty)
  expect_equal(mean(vapply(estimates, `[[`, numeric(1), "yty")), stats$yty)
})

test_that("LRM full-batch natural update with rate one equals CAVI", {
  X <- cbind(1, seq(-1, 1, length.out = 10))
  y <- drop(X %*% c(1, 2)) + cos(seq_len(10)) / 10
  hyper <- hyperparams_lrm(list(X = X, y = y), "OLS")
  cavi <- vi_lrm$new(y, X, hyper,
    list(iterations = 1, stochastic = FALSE, epsilon = NULL, seed = 720))
  svi <- vi_lrm$new(y, X, hyper,
    list(iterations = 1, stochastic = TRUE, minibatch_m = nrow(X),
         tau = 0, kappa = .7, epsilon = NULL, elbo_every = 1,
         min_iter = 1, seed = 720))
  expect_identical(cavi$lrm_vi_object$params, svi$lrm_vi_object$params)
  cavi$fit(FALSE)
  svi$fit(FALSE)
  expect_equal(cavi$results, svi$results, tolerance = 1e-12)
})

test_that("public LRM SVI is reproducible and exposes the uniform operations", {
  X <- cbind(intercept = 1, predictor = seq(-2, 2, length.out = 30))
  y <- drop(X %*% c(.5, 1.5)) + sin(seq_len(30)) / 5
  prior <- lrm_prior(c(0, 0), 10 * diag(2), 4, 1)
  control <- svi_control(
    max_iterations = 20, minibatch_size = 8, evaluation_draws = 6,
    stopping = svi_stopping(
      objective_tolerance = 0, parameter_tolerance = NULL,
      min_iterations = 1, patience = 2, check_every = 5,
      smoothing_window = 2))
  before <- rng_snapshot()
  a <- fit_model(lrm_model(2), list(X = X, y = y), prior, "svi",
                 initialization_ols(), control, seed = 721, verbose = FALSE)
  b <- fit_model(lrm_model(2), list(X = X, y = y), prior, "svi",
                 initialization_ols(), control, seed = 721, verbose = FALSE)
  expect_identical(rng_snapshot(), before)
  expect_equal(a$parameters, b$parameters)
  expect_identical(a$method, "SVI")
  expect_identical(a$metadata$sampling_unit, "observation")
  expect_equal(a$metadata$minibatch_size, 8)
  expect_identical(summarize_posterior(a)$method, "SVI")
  expect_identical(predict(a, X[1:2, ], n_draws = 4, seed = 722)$method,
                   "SVI")
  expect_identical(posterior_predictive_check(
    a, n_draws = 4, seed = 723)$method, "SVI")
  expect_identical(predictive_criteria(
    a, n_draws = 4, seed = 724)$method, "SVI")
})
