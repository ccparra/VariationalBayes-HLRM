test_that("scaling is learned only from supplied training data and is reversible", {
  X <- cbind(intercept = 1, continuous = c(2, 4, 6, 8), category = c(0, 1, 0, 1))
  data <- list(X = X, y = c(10, 14, 18, 22))
  specification <- data_scaling(
    predictors = "continuous", response = TRUE,
    intercept_column = "intercept")
  map <- fit_scaling(specification, data)
  transformed <- transform_data(data, map)
  expect_equal(mean(transformed$X[, "continuous"]), 0)
  expect_equal(sd(transformed$X[, "continuous"]), 1)
  expect_identical(transformed$X[, "category"], X[, "category"])
  expect_equal(mean(transformed$y), 0)
  expect_equal(sd(transformed$y), 1)
  expect_equal(map$training$observation_count, 4)

  beta <- c(3, 2, -1)
  beta_model <- transform_coefficients(beta, map)
  expect_equal(restore_coefficients(beta_model, map), beta)
  expect_equal(drop(transformed$X %*% beta_model),
               (drop(X %*% beta) - map$response_center) /
                 map$response_scale)
})

test_that("grouped scaling pools training rows and validates new column order", {
  X <- list(
    cbind(intercept = 1, x = c(1, 2)),
    cbind(intercept = 1, x = c(5, 8, 9)))
  y <- list(c(2, 4), c(6, 10, 12))
  map <- fit_scaling(data_scaling(), list(X = X, y = y))
  expect_equal(map$predictor_center[2], mean(c(1, 2, 5, 8, 9)))
  expect_equal(map$response_center, mean(unlist(y)))
  transformed <- transform_data(list(X = X, y = y), map)
  expect_equal(mean(unlist(lapply(transformed$X, function(value) value[, 2]))), 0)
  bad <- list(X = lapply(X, function(value) value[, 2:1]), y = y)
  expect_error(transform_data(bad, map), "column names or order")
})

test_that("prior transformations follow affine and variance-scale identities", {
  X <- cbind(1, c(-2, -1, 1, 4), c(0, 1, 0, 1))
  colnames(X) <- c("intercept", "x", "binary")
  data <- list(X = X, y = c(2, 3, 7, 12))
  map <- fit_scaling(data_scaling(predictors = "x"), data)

  lrm <- lrm_prior(c(1, 2, -1), diag(c(4, 2, 3)), 6, 5)
  scaled_lrm <- transform_prior(lrm, map)
  expect_equal(scaled_lrm$coefficient_mean,
               transform_coefficients(lrm$coefficient_mean, map))
  expect_equal(scaled_lrm$residual_variance_scale,
               lrm$residual_variance_scale / map$response_scale^2)
  restored_mean <- restore_coefficients(scaled_lrm$coefficient_mean, map)
  restored_covariance <- restore_covariance(
    scaled_lrm$coefficient_covariance, map)
  expect_equal(restored_mean, lrm$coefficient_mean)
  expect_equal(restored_covariance, lrm$coefficient_covariance)

  hierarchical <- chlrm_prior(
    c(1, 2, -1), diag(3), 6, 2 * diag(3), 8, 3, 4, c(2, 1))
  scaled <- transform_prior(hierarchical, map)
  expect_equal(scaled$residual_variance_scale_rate,
               hierarchical$residual_variance_scale_rate *
                 map$response_scale^2)
  expect_equal(scaled$mixture_concentration,
               hierarchical$mixture_concentration)
  expect_equal(restore_covariance(
    scaled$coefficient_dispersion_scale, map),
    hierarchical$coefficient_dispersion_scale)
})

test_that("scaling rejects leakage-prone or noninvertible declarations", {
  data <- list(X = cbind(1, c(1, 1, 1)), y = 1:3)
  expect_error(fit_scaling(data_scaling(), data), "nonzero finite")
  expect_error(data_scaling(intercept_column = NULL), "Centering requires")
  expect_error(fit_scaling(
    data_scaling(predictors = "x", intercept_column = 1), data),
    "no column names")
})

test_that("held-out data cannot alter a training scaling map", {
  training <- list(
    X = cbind(intercept = 1, x = 1:5),
    y = 2 * (1:5))
  held_out <- list(
    X = cbind(intercept = 1, x = c(100, 200)),
    y = c(1000, 2000))
  map <- fit_scaling(data_scaling(), training)
  before <- serialize(map, NULL)
  transformed <- transform_data(held_out, map)
  expect_identical(serialize(map, NULL), before)
  expect_gt(mean(transformed$X[, "x"]), 50)
  expect_gt(mean(transformed$y), 100)
  expect_equal(map$training$observation_count, 5)
})

test_that("fit_model keeps numerical scaling internal and restores LRM outputs", {
  scope <- inference_rng_scope(1201)
  on.exit(scope$restore())
  X <- cbind(intercept = 1, x = stats::rnorm(30))
  data <- list(X = X, y = drop(X %*% c(8, 2.5) + stats::rnorm(30, sd = 1.5)))
  fit <- fit_model(
    lrm_model(2), data, lrm_prior(c(0, 0), diag(100, 2), 4, 4),
    "cavi", initialization_ols(), cavi_control(4, evaluation_draws = 3),
    scaling = data_scaling(), seed = 1202, verbose = FALSE)
  map <- scaling_map(fit)
  backend <- attr(fit, "backend")

  expect_identical(fit$metadata$parameter_scale, "input")
  expect_identical(fit$metadata$objective_scale, "standardized_backend")
  expect_equal(
    as.numeric(fit$parameters$coefficients$mean),
    restore_coefficients(as.numeric(backend$results$MuBeta), map))
  expect_equal(
    fit$parameters$coefficients$covariance,
    restore_covariance(backend$results$SigmaBeta, map))
  expect_equal(
    fit$parameters$residual_variance$scale,
    restore_residual_variance(backend$results$bsigma2, map))
  expect_equal(
    summarize_posterior(fit)$parameters$coefficients$mean,
    as.numeric(fit$parameters$coefficients$mean))

  X_new <- X[1:4, , drop = FALSE]
  internal_prediction <- backend$predict(
    .validate_scaling_X(X_new, map), n_draws = 3, seed = 1203,
    return_draws = TRUE)
  public_prediction <- predict(
    fit, X_new, n_draws = 3, seed = 1203, return_draws = TRUE)
  expect_equal(public_prediction$draws[[1]],
               map$response_center + map$response_scale *
                 internal_prediction$draws[[1]])
  expect_identical(public_prediction$metadata$output_scale, "input")

  internal_fit <- backend$evaluate_in_sample_fit(return_details = TRUE)
  public_fit <- evaluate_in_sample_fit(fit, return_details = TRUE)
  expect_equal(public_fit$statistics$mse,
               internal_fit$statistics$mse * map$response_scale^2)
  expect_equal(public_fit$fitted_values[[1]],
               map$response_center + map$response_scale *
                 internal_fit$fitted_values[[1]])

  internal_criteria <- backend$predictive_criteria(
    n_draws = 3, seed = 1204, return_log_likelihood = TRUE)
  public_criteria <- predictive_criteria(
    fit, n_draws = 3, seed = 1204, return_log_likelihood = TRUE)
  expected_log_likelihood <- internal_criteria$log_likelihood$observation$values -
    log(map$response_scale)
  expect_equal(public_criteria$log_likelihood$observation$values,
               expected_log_likelihood)
  expect_equal(public_criteria$criteria$waic_observation$WAIC,
               waic_from_loglik(expected_log_likelihood)$WAIC)
  expect_false("log_likelihood" %in% names(predictive_criteria(
    fit, n_draws = 3, seed = 1204)))
})

test_that("scaling restores hierarchical and mixture variational parameters", {
  scope <- inference_rng_scope(1210)
  on.exit(scope$restore())
  groups <- lapply(seq_len(5), function(group) {
    X <- cbind(intercept = 1, x = stats::rnorm(6))
    list(X = X, y = drop(X %*% c(12, -3) + stats::rnorm(6)))
  })
  data <- list(X = lapply(groups, `[[`, "X"),
               y = lapply(groups, `[[`, "y"))
  controls <- list(
    cavi = cavi_control(3, evaluation_draws = 2),
    svi = svi_control(
      3, minibatch_size = 3, evaluation_draws = 2,
      stopping = svi_stopping(min_iterations = 1)))
  requests <- list(
    list(hlrm_model(2), hlrm_prior(
      c(0, 0), diag(10, 2), 4, diag(2), 4, 2, 2),
      "svi", initialization_ols(), controls$svi),
    list(chlrm_model(2, 2), chlrm_prior(
      c(0, 0), diag(10, 2), 4, diag(2), 4, 2, 2, c(1, 1)),
      "cavi", initialization_kmeans(), controls$cavi))

  for (index in seq_along(requests)) {
    request <- requests[[index]]
    fit <- fit_model(
      request[[1]], data, request[[2]], request[[3]], request[[4]],
      request[[5]], scaling = data_scaling(), seed = 1210 + index,
      verbose = FALSE)
    map <- scaling_map(fit)
    backend <- attr(fit, "backend")
    expect_equal(
      as.numeric(fit$parameters$population_coefficients$mean),
      restore_coefficients(as.numeric(backend$results$MuBeta), map))
    expect_equal(
      fit$parameters$population_covariance$scale,
      restore_covariance(backend$results$SS, map))
    expect_equal(
      fit$parameters$variance_scale$rate,
      backend$results$bxi2 / map$response_scale^2)
    local_name <- if (fit$model == "HLRM")
      "group_coefficients" else "component_coefficients"
    expect_equal(dim(fit$parameters[[local_name]]$mean),
                 dim(backend$results[[if (fit$model == "HLRM")
                   "MuBetaj" else "MuBetak"]]))
    numeric_values <- unlist(lapply(fit$parameters, function(parameter)
      parameter[setdiff(names(parameter), "distribution")]))
    expect_true(all(is.finite(numeric_values)))
  }
})

test_that("explicit initialization values are transformed without mutation", {
  data <- list(
    X = cbind(intercept = 1, x = c(-2, -1, 1, 4)),
    y = c(5, 6, 10, 15))
  map <- fit_scaling(data_scaling(), data)
  initialization <- initialization_values(list(
    coefficient_mean = c(5, 2),
    coefficient_covariance = diag(c(4, 9)),
    residual_variance_shape = 3,
    residual_variance_rate = 8))
  before <- serialize(initialization, NULL)
  transformed <- transform_initialization(
    initialization, map, "LRM", "CAVI")
  expect_identical(serialize(initialization, NULL), before)
  expect_equal(transformed$values$coefficient_mean,
               transform_coefficients(c(5, 2), map))
  expect_equal(transformed$values$coefficient_covariance,
               .transform_covariance(diag(c(4, 9)), map$to_model_matrix))
  expect_equal(transformed$values$residual_variance_rate,
               8 / map$response_scale^2)
})

test_that("scaled MCMC restoration covers LRM HLRM and CHLRM", {
  scope <- inference_rng_scope(1220)
  on.exit(scope$restore())
  X <- cbind(intercept = 1, x = stats::rnorm(20))
  lrm_data <- list(X = X, y = drop(X %*% c(7, -2) + stats::rnorm(20)))
  groups <- lapply(seq_len(4), function(group) {
    Xj <- cbind(intercept = 1, x = stats::rnorm(6))
    list(X = Xj, y = drop(Xj %*% c(7, -2) + stats::rnorm(6)))
  })
  grouped_data <- list(
    X = lapply(groups, `[[`, "X"), y = lapply(groups, `[[`, "y"))
  requests <- list(
    list(lrm_model(2), lrm_data,
         lrm_prior(c(0, 0), diag(10, 2), 4, 2), initialization_ols()),
    list(hlrm_model(2), grouped_data,
         hlrm_prior(c(0, 0), diag(10, 2), 4, diag(2), 4, 2, 2),
         initialization_ols()),
    list(chlrm_model(2, 2), grouped_data,
         chlrm_prior(c(0, 0), diag(10, 2), 4, diag(2), 4, 2, 2,
                     c(1, 1)), initialization_kmeans()))

  for (index in seq_along(requests)) {
    request <- requests[[index]]
    fit <- fit_model(
      request[[1]], request[[2]], request[[3]], "mcmc", request[[4]],
      mcmc_control(warmup = 1, post_warmup_draws = 3,
                   evaluation_draws = 2),
      scaling = data_scaling(), seed = 1220 + index, verbose = FALSE)
    map <- scaling_map(fit)
    backend <- attr(fit, "backend")
    coefficient_name <- if (fit$model == "LRM")
      "coefficients" else "population_coefficients"
    expect_equal(
      fit$parameters[[coefficient_name]]$draws,
      .transform_coefficient_array(backend$chains$Beta, map, "to_input"))
    if (fit$model == "LRM") {
      expect_equal(fit$parameters$residual_variance$draws,
                   backend$chains$sigma2 * map$response_scale^2)
    } else {
      expect_equal(
        fit$parameters$population_covariance$draws,
        .transform_covariance_collection(
          backend$chains$Sigma, map$restore_matrix))
      local_name <- if (fit$model == "HLRM")
        "group_coefficients" else "component_coefficients"
      native_name <- if (fit$model == "HLRM") "Betaj" else "Betak"
      expect_equal(
        fit$parameters[[local_name]]$draws,
        .transform_coefficient_array(
          backend$chains[[native_name]], map, "to_input"))
    }
    draws <- posterior_draws(fit, n_draws = 2, seed = 1230 + index)
    expect_true(all(is.finite(unlist(draws$parameters))))
  }
})
