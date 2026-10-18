.end_to_end_simulation_fixtures <- function() {
  list(
    LRM = .test_lrm_simulation(
      sample_size = 30,
      coefficients = c(1, 2),
      residual_variance = 0.5,
      covariate_range = c(-1, 1),
      seed = 1701),
    HLRM = .test_hlrm_simulation(
      group_sizes = rep(8, 6),
      group_coefficients = rbind(
        matrix(rep(c(0, 1), 3), nrow = 3, byrow = TRUE),
        matrix(rep(c(2, -1), 3), nrow = 3, byrow = TRUE)),
      group_variances = rep(0.5, 6),
      seed = 1702),
    CHLRM = .test_chlrm_simulation(
      group_count = 6,
      group_sizes = rep(8, 6),
      component_count = 2,
      coefficient_count = 2,
      component_coefficients = rbind(c(0, 1), c(2, -1)),
      component_variances = c(0.5, 0.5),
      mixture_weights = c(0.5, 0.5),
      allocations = c(1, 1, 1, 2, 2, 2),
      seed = 1703))
}

.end_to_end_priors <- function() {
  list(
    LRM = lrm_prior(c(0, 0), diag(10, 2), 4, 1),
    HLRM = hlrm_prior(c(0, 0), diag(10, 2), 4, diag(2), 4, 2, 2),
    CHLRM = chlrm_prior(
      c(0, 0), diag(10, 2), 4, diag(2), 4, 2, 2, c(1, 1)))
}

.end_to_end_controls <- function() {
  list(
    cavi = cavi_control(
      max_iterations = 4,
      evaluation_draws = 3,
      stopping = cavi_stopping(
        objective_tolerance = 1e20,
        parameter_tolerance = 1e20,
        min_iterations = 2,
        patience = 1)),
    svi = svi_control(
      max_iterations = 4,
      minibatch_size = 3,
      evaluation_draws = 3,
      stopping = svi_stopping(
        objective_tolerance = 1e20,
        parameter_tolerance = 1e20,
        min_iterations = 2,
        patience = 1,
        check_every = 1,
        smoothing_window = 1)),
    mcmc = mcmc_control(
      warmup = 1,
      post_warmup_draws = 3,
      evaluation_draws = 3))
}

test_that("public simulation fixtures exercise every model and engine", {
  simulations <- .end_to_end_simulation_fixtures()
  priors <- .end_to_end_priors()
  controls <- .end_to_end_controls()
  specifications <- list(
    LRM = lrm_model(2),
    HLRM = hlrm_model(2),
    CHLRM = chlrm_model(2, 2))
  initializations <- list(
    LRM = initialization_ols(),
    HLRM = initialization_ols(),
    CHLRM = initialization_kmeans())

  request_index <- 0L
  for (model in names(simulations)) {
    expect_s3_class(simulations[[model]], "vb_simulation_result")
    expect_identical(simulations[[model]]$model, model)
    for (engine in names(controls)) {
      request_index <- request_index + 1L
      fit <- fit_model(
        model = specifications[[model]],
        data = simulations[[model]]$train,
        prior = priors[[model]],
        method = engine,
        initialization = initializations[[model]],
        control = controls[[engine]],
        seed = 1710 + request_index,
        verbose = FALSE)

      expect_s3_class(fit, "vb_model_fit")
      expect_identical(fit$schema_version, "fit-result-2")
      expect_identical(fit$model, model)
      expect_identical(fit$method, toupper(engine))

      draws <- posterior_draws(
        fit,
        n_draws = 3,
        seed = 1730 + request_index,
        include_latent = identical(model, "CHLRM"))
      expect_identical(names(draws$parameters),
                       parameter_names_for_model(model))
      expect_true(all(vapply(draws$parameters, function(parameter)
        all(is.finite(parameter)), logical(1L))))
      expect_identical("latent_variables" %in% names(draws),
                       identical(model, "CHLRM"))

      in_sample <- evaluate_in_sample_fit(fit)
      check <- posterior_predictive_check(
        fit, n_draws = 3, seed = 1750 + request_index)
      criteria <- predictive_criteria(
        fit, n_draws = 3, seed = 1770 + request_index,
        return_log_likelihood = TRUE)
      expect_s3_class(in_sample, "vb_in_sample_fit")
      expect_s3_class(check, "vb_posterior_predictive_check")
      expect_s3_class(criteria, "vb_predictive_criteria")
      expect_true(length(criteria$log_likelihood) >= 1L)
      expect_true(all(vapply(criteria$log_likelihood, inherits, logical(1L),
                             "vb_log_likelihood")))
    }
  }
})
