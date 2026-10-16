test_that("public results share a core schema without nonapplicable fields", {
  dl <- .test_lrm_simulation(12, c(1, .5), 1, c(-1, 1), seed = 1001)
  hl <- hyperparams_lrm(dl$train, "OLS")
  lrm <- vi_lrm$new(dl$train$y, dl$train$X, hl,
    list(iterations = 3, epsilon = NULL, stochastic = FALSE))
  fit_lrm <- lrm$fit(FALSE)

  dh <- .test_hlrm_simulation(rep(6, 3),
    matrix(rep(c(1, .5), 3), 3, 2, byrow = TRUE), rep(1, 3), seed = 1002)
  hh <- hyperparams_hlrm(dh$train, "OLS")
  hlrm <- vi_hlrm$new(dh$train$y, dh$train$X, hh,
    list(iterations = 3, epsilon = NULL, stochastic = FALSE))
  fit_hlrm <- hlrm$fit(FALSE)

  dc <- .test_chlrm_simulation(4, 6, 2, 2,
    component_coefficients = rbind(c(-1, .5), c(2, -.2)),
    component_variances = c(1, 1), allocations = c(1, 1, 2, 2),
    seed = 1003)
  hc <- hyperparams_chlrm(dc$train, 2, "OLS")
  chlrm <- vi_chlrm$new(dc$train$y, dc$train$X, 2, hc,
    list(iterations = 3, epsilon = NULL, stochastic = FALSE))
  fit_chlrm <- chlrm$fit(FALSE)

  fits <- list(fit_lrm, fit_hlrm, fit_chlrm)
  expect_true(all(vapply(fits, inherits, logical(1), "vb_fit_result")))
  fit_core <- c("schema_version", "model", "method", "status",
                "representation", "parameters", "objective", "diagnostics",
                "timing", "metadata")
  expect_true(all(vapply(fits, function(x)
    all(fit_core %in% names(x)), logical(1))))
  expect_true(all(vapply(seq_along(fits), function(i)
    identical(names(fits[[i]]$parameters),
              parameter_names_for_model(c("LRM", "HLRM", "CHLRM")[i])),
    logical(1))))
  expect_false("latent_variables" %in% names(fit_lrm))
  expect_false("latent_variables" %in% names(fit_hlrm))
  expect_identical(names(fit_chlrm$latent_variables),
                   "allocation_probabilities")
  expect_true(all(vapply(fits, function(x)
    identical(names(x$objective), c("elbo", "evaluation_iterations")),
    logical(1))))
  expect_identical(vapply(fits, `[[`, character(1), "model"),
                   c("LRM", "HLRM", "CHLRM"))

  summaries <- list(lrm$summarize_posterior(), hlrm$summarize_posterior(),
                    chlrm$summarize_posterior(include_incidence = FALSE))
  expect_true(all(vapply(summaries, inherits, logical(1), "vb_posterior_summary")))
  summary_core <- c("schema_version", "model", "method", "interval", "parameters")
  expect_true(all(vapply(summaries, function(x)
    all(summary_core %in% names(x)), logical(1))))
  expect_true(all(vapply(seq_along(summaries), function(i)
    identical(names(summaries[[i]]$parameters),
              parameter_names_for_model(c("LRM", "HLRM", "CHLRM")[i])),
    logical(1))))
  expect_false("allocations" %in% names(summaries[[1]]))
  expect_false("relabeling" %in% names(summaries[[2]]))
  expect_true(all(c("allocations", "relabeling") %in% names(summaries[[3]])))
  expect_false(any(vapply(summaries[[3]]$allocations, is.null, logical(1L))))
  expect_false(any(vapply(summaries[[3]]$relabeling, is.null, logical(1L))))
  expect_false(any(c("dahl_labels", "co_clustering", "component_count") %in%
                   names(summaries[[3]]$allocations)))
  expect_identical(vapply(summaries, `[[`, character(1), "model"),
                   c("LRM", "HLRM", "CHLRM"))

  fits_in_sample <- list(lrm$evaluate_in_sample_fit(),
                         hlrm$evaluate_in_sample_fit(),
                         chlrm$evaluate_in_sample_fit())
  expect_true(all(vapply(fits_in_sample, inherits, logical(1),
                         "vb_in_sample_fit")))
  ppcs <- list(lrm$posterior_predictive_check(3, seed = 1004),
               hlrm$posterior_predictive_check(3, seed = 1004),
               chlrm$posterior_predictive_check(3, seed = 1004))
  expect_true(all(vapply(ppcs, inherits, logical(1),
                         "vb_posterior_predictive_check")))
  criteria <- list(lrm$predictive_criteria(3, seed = 1004),
                   hlrm$predictive_criteria(3, seed = 1004),
                   chlrm$predictive_criteria(3, seed = 1004))
  expect_true(all(vapply(criteria, inherits, logical(1),
                         "vb_predictive_criteria")))
  expect_identical(names(criteria[[1]]$criteria),
                   "waic_observation")
  expect_identical(names(criteria[[2]]$criteria),
                   c("waic_observation_conditional", "waic_group_conditional"))
  expect_identical(names(criteria[[3]]$criteria),
                   "waic_group_marginal")
  expect_identical(vapply(criteria, `[[`, character(1), "model"),
                   c("LRM", "HLRM", "CHLRM"))
  retained <- list(
    lrm$predictive_criteria(3, seed = 1004,
                            return_log_likelihood = TRUE),
    hlrm$predictive_criteria(3, seed = 1004,
                             return_log_likelihood = TRUE),
    chlrm$predictive_criteria(3, seed = 1004,
                              return_log_likelihood = TRUE))
  expect_true(all(vapply(retained, function(x)
    length(x$log_likelihood) >= 1L &&
      all(vapply(x$log_likelihood, inherits, logical(1L),
                 "vb_log_likelihood")), logical(1L))))
  expect_true(all(vapply(criteria, function(x)
    !"log_likelihood" %in% names(x), logical(1L))))

  predictions <- list(lrm$predict(dl$train$X, n_draws = 3, seed = 1005),
                      hlrm$predict(dh$train$X, n_draws = 3, seed = 1005),
                      chlrm$predict(dc$train$X, n_draws = 3, seed = 1005))
  expect_true(all(vapply(predictions, inherits, logical(1),
                         "vb_prediction_result")))
  prediction_core <- c("schema_version", "model", "method", "target",
                       "interval", "metadata", "summary")
  expect_true(all(vapply(predictions, function(x)
    all(prediction_core %in% names(x)), logical(1))))
  expect_identical(vapply(predictions, `[[`, character(1), "model"),
                   c("LRM", "HLRM", "CHLRM"))
  expect_true(all(vapply(predictions, `[[`, character(1), "target") == "mean"))
  expect_false("allocation" %in% names(predictions[[1]]))
  expect_false("draws" %in% names(predictions[[1]]))
  expect_true("allocation" %in% names(predictions[[3]]))
})

test_that("parameter summaries have one uniform interval vocabulary", {
  draws <- matrix(1:12, 4, 3)
  x <- summarize_draws_parameter(draws, 2)
  expect_identical(names(x), c("mean", "sd", "lower", "upper", "level",
                               "distribution", "mean_exists", "contains_zero"))
  expect_equal(dim(x$mean), NULL)
  expect_length(x$mean, 3)

  ig <- summarize_inverse_gamma_parameter(c(.8, 3), c(1, 2))
  expect_identical(ig$mean_exists, c(FALSE, TRUE))
  expect_true(is.na(ig$mean[1]))
  expect_true(all(is.finite(c(ig$lower, ig$upper))))
})

test_that("model outputs omit nonapplicable quantities", {
  parameters <- new_parameter_slots("LRM")
  parameters$coefficients <- summarize_draws_parameter(matrix(rnorm(20), 10, 2), 2)
  parameters$residual_variance <- summarize_draws_parameter(rgamma(10, 2))
  out <- new_posterior_summary_result("LRM", "MCMC", parameters)
  expect_false("component_coefficients" %in% names(out$parameters))
  expect_false("allocations" %in% names(out))
  expect_error(new_posterior_summary_result("LRM", "MCMC", list()), "canonical")
})

test_that("MCMC, CAVI and SVI classes use identical public method signatures", {
  methods <- c("fit", "summarize_posterior", "evaluate_in_sample_fit",
               "posterior_predictive_check", "predictive_criteria",
               "predict", "posterior_predict")
  pairs <- list(list(mcmc_lrm, vi_lrm), list(mcmc_hlrm, vi_hlrm),
                list(mcmc_chlrm, vi_chlrm))
  for (pair in pairs) for (method in methods) {
    expect_identical(names(formals(pair[[1]]$public_methods[[method]])),
                     names(formals(pair[[2]]$public_methods[[method]])))
  }
})
