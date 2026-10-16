.typed_test_data <- function(seed = 710) {
  scope <- inference_rng_scope(seed)
  on.exit(scope$restore())
  X <- cbind(intercept = 1, predictor = stats::rnorm(30))
  lrm <- list(X = X, y = drop(X %*% c(1, 2) + stats::rnorm(30, sd = .5)))
  groups <- lapply(seq_len(6), function(j) {
    Xj <- cbind(intercept = 1, predictor = stats::rnorm(8))
    coefficient <- if (j <= 3) c(0, 1) else c(2, -1)
    list(X = Xj, y = drop(Xj %*% coefficient + stats::rnorm(8, sd = .5)))
  })
  list(lrm = lrm,
       grouped = list(X = lapply(groups, `[[`, "X"),
                      y = lapply(groups, `[[`, "y")))
}

.typed_test_priors <- function() {
  list(
    lrm = lrm_prior(c(0, 0), diag(10, 2), 4, 1),
    hlrm = hlrm_prior(c(0, 0), diag(10, 2), 4, diag(2), 4, 2, 2),
    chlrm = chlrm_prior(c(0, 0), diag(10, 2), 4, diag(2), 4, 2, 2,
                        c(1, 1)))
}

test_that("typed specifications reject incompatible or incomplete requests", {
  expect_error(chlrm_model(0), "component_count")
  expect_error(mcmc_control(post_warmup_draws = 1), "post_warmup_draws")
  expect_error(mcmc_control(chains = 1, stopping = mcmc_stopping()),
               "R-hat is not defined")
  expect_error(mcmc_control(chains = 2), "require stopping")
  expect_s3_class(mcmc_control(chains = 2, stopping = mcmc_stopping()),
                  "vb_mcmc_control")
  expect_equal(mcmc_control(chains = 1, max_memory_mb = 10)$max_memory_mb, 10)
  expect_equal(cavi_control(max_memory_mb = 10)$max_memory_mb, 10)
  expect_equal(svi_control(minibatch_size = 1, max_memory_mb = 10)$max_memory_mb,
               10)
  expect_error(cavi_control(starts = 1, require_converged = TRUE),
               "only to multiple")
  expect_error(cavi_control(2, stopping = cavi_stopping(min_iterations = 3)),
               "min_iterations")
  expect_error(svi_control(2, 1, stopping = svi_stopping(min_iterations = 3)),
               "min_iterations")
  data <- .typed_test_data()
  priors <- .typed_test_priors()
  expect_error(fit_model(lrm_model(), data$lrm, priors$lrm, "svi",
                         control = svi_control(2, 31,
                           stopping = svi_stopping(min_iterations = 1))),
               "number of observations")
  expect_error(fit_model(lrm_model(), data$lrm, priors$hlrm, "cavi",
                         control = cavi_control(2)), "same model family")
})

test_that("prior_moments reports existence and the package parameterization", {
  lrm <- prior_moments(lrm_prior(c(0, 1), diag(2), 6, 2))
  expect_identical(lrm$schema_version, "prior-moments-1")
  expect_equal(lrm$parameters$residual_variance$mean, 3)
  expect_equal(lrm$parameters$residual_variance$variance, 9)
  heavy_tail <- prior_moments(lrm_prior(0, matrix(1), 2, 1))
  expect_false(heavy_tail$parameters$residual_variance$mean_exists)
  expect_true(is.na(heavy_tail$parameters$residual_variance$mean))

  hierarchical <- prior_moments(hlrm_prior(
    c(0, 0), diag(2), 5, 2 * diag(2), 6, 2, 4))
  expect_equal(hierarchical$parameters$population_covariance$mean, diag(2))
  expect_equal(hierarchical$parameters$variance_scale$mean, .5)
  expect_equal(hierarchical$parameters$group_coefficients$covariance,
               2 * diag(2))
  expect_equal(hierarchical$parameters$group_variance$mean, .75)
  expect_equal(hierarchical$parameters$group_variance$variance, 1.125)

  clustered <- prior_moments(chlrm_prior(
    c(0, 0), diag(2), 5, diag(2), 6, 2, 4, c(2, 3)))
  expect_equal(clustered$parameters$mixture_weights$mean, c(.4, .6))
  expect_equal(rowSums(clustered$parameters$mixture_weights$covariance),
               c(0, 0), tolerance = 1e-15)
})

test_that("explicit initialization uses semantic public names", {
  data <- .typed_test_data()
  prior <- .typed_test_priors()$lrm
  values <- list(
    coefficient_mean = c(0, 0),
    coefficient_covariance = diag(2),
    residual_variance_shape = 2,
    residual_variance_rate = 1)
  fit <- fit_model(
    lrm_model(2), data$lrm, prior, "cavi", initialization_values(values),
    cavi_control(2, 2), seed = 902, verbose = FALSE)
  expect_s3_class(fit, "vb_model_fit")
  expect_error(fit_model(
    lrm_model(2), data$lrm, prior, "cavi",
    initialization_values(list(
      MuBeta = c(0, 0), SigmaBeta = diag(2), asigma2 = 2, bsigma2 = 1)),
    cavi_control(2, 2), seed = 902, verbose = FALSE),
    "public fields")
  invalid <- values
  invalid$residual_variance_rate <- -1
  expect_error(fit_model(
    lrm_model(2), data$lrm, prior, "cavi",
    initialization_values(invalid), cavi_control(2, 2),
    seed = 902, verbose = FALSE), "must contain positive")
})

test_that("fit_model gives one public contract across implemented engines", {
  data <- .typed_test_data()
  priors <- .typed_test_priors()
  fast_cavi <- cavi_control(
    4, 2, cavi_stopping(1e20, 1e20, min_iterations = 2, patience = 1))
  fast_svi <- svi_control(
    4, 3, evaluation_draws = 2,
    stopping = svi_stopping(1e20, 1e20, min_iterations = 2,
                            patience = 1, check_every = 1,
                            smoothing_window = 1))
  fast_mcmc <- mcmc_control(warmup = 1, post_warmup_draws = 3,
                            evaluation_draws = 2)
  requests <- list(
    list(lrm_model(2), data$lrm, priors$lrm, "cavi",
         initialization_ols(), fast_cavi),
    list(lrm_model(2), data$lrm, priors$lrm, "svi",
         initialization_ols(), fast_svi),
    list(lrm_model(2), data$lrm, priors$lrm, "mcmc",
         initialization_ols(), fast_mcmc),
    list(hlrm_model(2), data$grouped, priors$hlrm, "cavi",
         initialization_ols(), fast_cavi),
    list(hlrm_model(2), data$grouped, priors$hlrm, "svi",
         initialization_ols(), fast_svi),
    list(hlrm_model(2), data$grouped, priors$hlrm, "mcmc",
         initialization_ols(), fast_mcmc),
    list(chlrm_model(2, 2), data$grouped, priors$chlrm, "cavi",
         initialization_kmeans(), fast_cavi),
    list(chlrm_model(2, 2), data$grouped, priors$chlrm, "svi",
         initialization_kmeans(), fast_svi),
    list(chlrm_model(2, 2), data$grouped, priors$chlrm, "mcmc",
         initialization_kmeans(), fast_mcmc))
  methods <- c("CAVI", "SVI", "MCMC", "CAVI", "SVI", "MCMC", "CAVI", "SVI", "MCMC")
  for (i in seq_along(requests)) {
    request <- requests[[i]]
    fit <- fit_model(request[[1]], request[[2]], request[[3]], request[[4]],
                     request[[5]], request[[6]], seed = 800 + i,
                     verbose = FALSE)
    expect_s3_class(fit, "vb_model_fit")
    expect_identical(fit$schema_version, "fit-result-2")
    expect_identical(fit$method, methods[i])
    expect_s3_class(fit$diagnostics, "vb_convergence")
    expect_identical(fit$diagnostics$schema_version, "convergence-1")
    expect_identical(fit$specifications$prior, request[[3]])
    expect_true(inherits(attr(fit, "backend"), "R6"))
    draws <- posterior_draws(
      fit, n_draws = 2, seed = 880 + i,
      include_latent = identical(fit$model, "CHLRM"))
    expect_identical(names(draws$parameters),
                     parameter_names_for_model(fit$model))
    expect_true(all(vapply(draws$parameters, function(value)
      all(is.finite(value)), logical(1L))))
  }
})

test_that("typed parameter stopping is visible in fitted diagnostics", {
  data <- .typed_test_data()
  prior <- .typed_test_priors()$lrm
  fit <- fit_model(
    lrm_model(2), data$lrm, prior, "cavi", initialization_ols(),
    cavi_control(5, 2, cavi_stopping(
      objective_tolerance = 1e20, parameter_tolerance = 1e20,
      min_iterations = 2, patience = 1)),
    seed = 901, verbose = FALSE)
  expect_identical(fit$diagnostics$status, "converged")
  expect_identical(fit$diagnostics$reason,
                   "objective_and_parameters_stable")
  expect_true(is.finite(
    fit$diagnostics$metrics$parameter_relative_change))
})

test_that("typed multistart is uniform across all variational engines", {
  data <- .typed_test_data()
  priors <- .typed_test_priors()
  cavi <- cavi_control(
    4, 2, cavi_stopping(1e20, 1e20, 2, 1), starts = 2,
    require_converged = FALSE, max_memory_mb = 50)
  svi <- svi_control(
    4, 3, evaluation_draws = 2,
    stopping = svi_stopping(1e20, 1e20, 2, 1, 1, 1), starts = 2,
    require_converged = FALSE, max_memory_mb = 50)
  requests <- list(
    list(lrm_model(2), data$lrm, priors$lrm, "cavi",
         initialization_ols(), cavi),
    list(lrm_model(2), data$lrm, priors$lrm, "svi",
         initialization_ols(), svi),
    list(hlrm_model(2), data$grouped, priors$hlrm, "cavi",
         initialization_ols(), cavi),
    list(hlrm_model(2), data$grouped, priors$hlrm, "svi",
         initialization_ols(), svi),
    list(chlrm_model(2, 2), data$grouped, priors$chlrm, "cavi",
         initialization_kmeans(), cavi),
    list(chlrm_model(2, 2), data$grouped, priors$chlrm, "svi",
         initialization_kmeans(), svi))
  for (index in seq_along(requests)) {
    request <- requests[[index]]
    fit <- fit_model(
      request[[1]], request[[2]], request[[3]], request[[4]],
      request[[5]], request[[6]], seed = 920 + index, verbose = FALSE)
    expect_identical(fit$metadata$starts, 2L)
    expect_true(fit$metadata$selected_start %in% 1:2)
    expect_equal(nrow(fit$metadata$start_summary), 2)
    orchestration <- fit_orchestration(fit)
    expect_s3_class(orchestration, "vb_fit_orchestration")
    expect_identical(orchestration$schema_version, "fit-orchestration-1")
    expect_false(any(c("attempts", "selected_model", "run_settings") %in%
                     names(orchestration)))
    if (request[[1]]$model == "CHLRM")
      expect_true(is.data.frame(orchestration$start_classifications))
    expect_identical(attr(fit, "backend")$fit_status, "fitted")
    expect_identical(summarize_posterior(fit)$method,
                     toupper(request[[4]]))
  }
})

test_that("typed multichain is uniform across LRM, HLRM and CHLRM", {
  data <- .typed_test_data()
  priors <- .typed_test_priors()
  control <- mcmc_control(
    warmup = 1, post_warmup_draws = 10, chains = 2,
    evaluation_draws = 2,
    stopping = mcmc_stopping(
      min_post_warmup_draws = 4, check_every = 5,
      rhat_max = 100, ess_bulk_min = .1, ess_tail_min = .1),
    max_memory_mb = 50)
  requests <- list(
    list(lrm_model(2), data$lrm, priors$lrm, initialization_ols()),
    list(hlrm_model(2), data$grouped, priors$hlrm, initialization_ols()),
    list(chlrm_model(2, 2), data$grouped, priors$chlrm,
         initialization_kmeans()))
  for (index in seq_along(requests)) {
    request <- requests[[index]]
    fit <- fit_model(
      request[[1]], request[[2]], request[[3]], "mcmc", request[[4]],
      control, seed = 930 + index, verbose = FALSE)
    expect_identical(fit$metadata$chains, 2L)
    orchestration <- fit_orchestration(fit)
    expect_s3_class(orchestration, "vb_fit_orchestration")
    expect_identical(orchestration$type, "multichain")
    expect_false(any(c("chains", "chain_models", "run_settings") %in%
                     names(orchestration)))
    expect_equal(nrow(attr(fit, "backend")$chains$Beta),
                 2 * fit$diagnostics$work_completed)
    expect_identical(
      summarize_posterior(fit)$schema_version,
      "posterior-summary-1")
    expect_identical(summarize_posterior(fit)$method, "MCMC")
  }
})

test_that("public fit operations dispatch without exposing the R6 backend", {
  data <- .typed_test_data()
  fit <- fit_model(
    lrm_model(2), data$lrm, .typed_test_priors()$lrm, "cavi",
    initialization_ols(), cavi_control(4, 2), seed = 940,
    verbose = FALSE)
  expect_error(fit_orchestration(fit), "one chain or one variational start")
  results <- list(
    summarize_posterior(fit),
    evaluate_in_sample_fit(fit),
    posterior_predictive_check(fit, n_draws = 2, seed = 941),
    predictive_criteria(
      fit, n_draws = 2, seed = 942,
      return_log_likelihood = TRUE),
    predict(fit, data$lrm$X[1:3, , drop = FALSE], n_draws = 2,
            seed = 943),
    posterior_predict(
      fit, data$lrm$X[1:3, , drop = FALSE], n_draws = 2,
      seed = 944, return_draws = FALSE))
  expect_true(all(vapply(results, function(result)
    identical(result$method, "CAVI"), logical(1L))))
  expect_identical(results[[4]]$schema_version, "predictive-criteria-2")
  expect_identical(results[[5]]$target, "mean")
  expect_identical(results[[6]]$target, "response")
})

test_that("posterior_draws exposes complete semantic draws for CAVI and SVI", {
  data <- .typed_test_data()
  priors <- .typed_test_priors()
  requests <- list(
    list(lrm_model(2), data$lrm, priors$lrm, "cavi",
         initialization_ols(), cavi_control(3, evaluation_draws = 4)),
    list(hlrm_model(2), data$grouped, priors$hlrm, "svi",
         initialization_ols(), svi_control(
           3, 3, evaluation_draws = 4,
           stopping = svi_stopping(min_iterations = 1))),
    list(chlrm_model(2, 2), data$grouped, priors$chlrm, "cavi",
         initialization_kmeans(), cavi_control(3, evaluation_draws = 4))
  )
  for (index in seq_along(requests)) {
    request <- requests[[index]]
    fit <- fit_model(
      request[[1]], request[[2]], request[[3]], request[[4]],
      request[[5]], request[[6]], seed = 950 + index, verbose = FALSE)
    draws <- posterior_draws(
      fit, seed = 960 + index,
      include_latent = identical(fit$model, "CHLRM"))
    repeated <- posterior_draws(
      fit, seed = 960 + index,
      include_latent = identical(fit$model, "CHLRM"))
    expect_s3_class(draws, "vb_posterior_draws")
    expect_identical(draws$schema_version, "posterior-draws-1")
    expect_identical(draws$distribution, "variational_approximation")
    expect_identical(names(draws$parameters),
                     parameter_names_for_model(fit$model))
    expect_identical(draws$parameters, repeated$parameters)
    expect_equal(draws$metadata$draw_count, 4L)
    expect_true(all(vapply(draws$parameters, function(value)
      if (is.null(dim(value))) length(value) == 4L else dim(value)[1L] == 4L,
      logical(1L))))
    if (fit$model == "CHLRM") {
      expect_equal(dim(draws$latent_variables$group_allocations), c(4, 6))
    } else {
      expect_false("latent_variables" %in% names(draws))
    }
  }
})

test_that("posterior_draws retains MCMC provenance and applies memory limits", {
  data <- .typed_test_data()
  control <- mcmc_control(
    warmup = 1, post_warmup_draws = 6, chains = 2,
    evaluation_draws = 2,
    stopping = mcmc_stopping(
      min_post_warmup_draws = 4, check_every = 3,
      rhat_max = 100, ess_bulk_min = .1, ess_tail_min = .1),
    max_memory_mb = 50)
  fit <- fit_model(
    chlrm_model(2, 2), data$grouped, .typed_test_priors()$chlrm, "mcmc",
    initialization_kmeans(), control, seed = 970, verbose = FALSE)
  draws <- posterior_draws(
    fit, n_draws = 4, seed = 971, include_latent = TRUE)
  expect_identical(draws$distribution, "posterior")
  expect_named(draws$draw_index,
               c("draw", "stored_draw", "chain", "iteration"))
  expect_true(all(draws$draw_index$chain %in% 1:2))
  expect_equal(dim(draws$parameters$component_coefficients), c(4, 2, 2))
  expect_equal(dim(draws$latent_variables$group_allocations), c(4, 6))
  expect_error(posterior_draws(fit, n_draws = 100), "cannot exceed")
  expect_error(posterior_draws(fit, max_output_mb = 1e-8),
               "exceed max_output_mb")
})
