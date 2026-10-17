new_in_sample_fit_result <- function(model, method, mse, r2, metadata,
                                     fitted_values = NULL,
                                     residuals = NULL) {
  if (!model %in% c("LRM", "HLRM", "CHLRM")) stop("unknown model")
  if (!method %in% c("MCMC", "CAVI", "SVI")) stop("unknown inference method")
  if (!is.list(metadata)) stop("metadata must be a list")
  if (!is.numeric(mse) || length(mse) != 1L || !is.finite(mse) || mse < 0)
    stop("mse must be a nonnegative finite scalar")
  if (!is.numeric(r2) || length(r2) != 1L || (!is.na(r2) && !is.finite(r2)))
    stop("r2 must be a finite scalar or NA when response variance is zero")
  out <- list(schema_version = "in-sample-fit-1", model = model,
              method = method, target = "observed_data",
              statistics = list(mse = mse, r2 = r2), metadata = metadata)
  if (!is.null(fitted_values)) out$fitted_values <- fitted_values
  if (!is.null(residuals)) out$residuals <- residuals
  structure(out, class = c("vb_in_sample_fit", "list"))
}

new_predictive_criteria_result <- function(model, method, predictive_target,
                                           criteria, metadata,
                                           log_likelihood = NULL) {
  if (!model %in% c("LRM", "HLRM", "CHLRM")) stop("unknown model")
  if (!method %in% c("MCMC", "CAVI", "SVI")) stop("unknown inference method")
  predictive_target <- match.arg(predictive_target,
    c("observed_data", "observed_groups", "new_groups"))
  allowed <- c("waic_observation", "waic_group",
               "waic_observation_conditional",
               "waic_group_conditional", "waic_group_marginal")
  if (!is.list(criteria) || !length(criteria) || is.null(names(criteria)) ||
      anyDuplicated(names(criteria)) || any(!names(criteria) %in% allowed) ||
      any(vapply(criteria, is.null, logical(1L))))
    stop("criteria contains invalid, duplicated or NULL entries")
  if (!is.list(metadata)) stop("metadata must be a list")
  if (!is.null(log_likelihood) &&
      (!is.list(log_likelihood) || !length(log_likelihood) ||
       is.null(names(log_likelihood)) || anyDuplicated(names(log_likelihood)) ||
       any(!vapply(log_likelihood, inherits, logical(1L),
                   "vb_log_likelihood"))))
    stop("log_likelihood must be NULL or a named list of log-likelihood-1 objects")
  out <- list(schema_version = "predictive-criteria-2", model = model,
    method = method, predictive_target = predictive_target,
    criteria = criteria, metadata = metadata)
  if (!is.null(log_likelihood)) out$log_likelihood <- log_likelihood
  structure(out, class = c("vb_predictive_criteria", "list"))
}

new_predictive_score_result <- function(fit, log_likelihood,
                                        return_log_likelihood = FALSE) {
  if (!inherits(fit, "vb_model_fit"))
    stop("fit must be produced by fit_model()")
  if (!inherits(log_likelihood, "vb_log_likelihood"))
    stop("log_likelihood must be a log-likelihood-1 object")
  if (!is.logical(return_log_likelihood) ||
      length(return_log_likelihood) != 1L || is.na(return_log_likelihood))
    stop("return_log_likelihood must be TRUE or FALSE")
  pointwise_log_density <- apply(
    log_likelihood$values, 1L, logmeanexp)
  pointwise <- cbind(
    log_likelihood$predictive_units,
    log_predictive_density = pointwise_log_density)
  out <- list(
    schema_version = "predictive-score-1",
    model = fit$model,
    method = fit$method,
    predictive_target = log_likelihood$metadata$predictive_target,
    unit = log_likelihood$unit,
    conditioning = log_likelihood$conditioning,
    statistics = list(
      elpd = sum(pointwise_log_density),
      mean_log_predictive_density = mean(pointwise_log_density),
      predictive_units = length(pointwise_log_density)),
    pointwise = pointwise,
    metadata = list(
      draws = ncol(log_likelihood$values),
      draw_source = if (fit$method == "MCMC")
        "posterior_mcmc" else "variational_q",
      density_scale = "input",
      interpretation = "out-of-sample log predictive density under the declared target and unit"))
  if (return_log_likelihood) out$log_likelihood <- log_likelihood
  structure(out, class = c("vb_predictive_score", "list"))
}

new_cross_validation_result <- function(model, method, plan, fold_statistics,
                                        pointwise, seeds, fits = NULL,
                                        prediction_pointwise = NULL) {
  if (!model %in% c("LRM", "HLRM", "CHLRM")) stop("unknown model")
  if (!method %in% c("MCMC", "CAVI", "SVI")) stop("unknown method")
  if (!inherits(plan, "vb_cross_validation_plan"))
    stop("plan must be created by cross_validation_plan()")
  if (!is.data.frame(fold_statistics) ||
      nrow(fold_statistics) != plan$fold_count)
    stop("fold_statistics must contain one row per fold")
  if (!is.data.frame(pointwise) || !nrow(pointwise))
    stop("pointwise must be a nonempty data frame")
  total_units <- sum(fold_statistics$predictive_units)
  statistics <- list(
    elpd = sum(fold_statistics$elpd),
    mean_log_predictive_density =
      sum(fold_statistics$elpd) / total_units,
    predictive_units = total_units,
    folds = plan$fold_count)
  if (!is.null(prediction_pointwise)) {
    if (!is.data.frame(prediction_pointwise) ||
        !all(c("squared_error", "absolute_error", "covered",
               "interval_width") %in% names(prediction_pointwise)))
      stop("prediction_pointwise has an invalid schema")
    statistics$rmse <- sqrt(mean(prediction_pointwise$squared_error))
    statistics$mae <- mean(prediction_pointwise$absolute_error)
    statistics$empirical_coverage <- mean(prediction_pointwise$covered)
    statistics$mean_interval_width <- mean(
      prediction_pointwise$interval_width)
    statistics$observations <- nrow(prediction_pointwise)
  }
  out <- list(
    schema_version = "cross-validation-result-1",
    model = model,
    method = method,
    scheme = plan$scheme,
    predictive_target = plan$predictive_target,
    unit = plan$predictive_unit,
    statistics = statistics,
    folds = fold_statistics,
    pointwise = pointwise,
    metadata = list(
      seeds = seeds,
      execution = "sequential",
      scaling = "learned independently inside each training fold"))
  if (!is.null(prediction_pointwise))
    out$prediction_pointwise <- prediction_pointwise
  if (!is.null(fits)) out$fits <- fits
  structure(out, class = c("vb_cross_validation_result", "list"))
}

new_predictive_performance_result <- function(
    fit, predictive_target, level, observed, mean_prediction,
    response_prediction, score) {
  if (!inherits(fit, "vb_model_fit"))
    stop("fit must be produced by fit_model()")
  if (!inherits(score, "vb_predictive_score"))
    stop("score must be produced by predictive_score()")
  observed <- as.numeric(observed)
  mean_summary <- mean_prediction$summary
  response_summary <- response_prediction$summary
  if (length(observed) != nrow(mean_summary) ||
      length(observed) != nrow(response_summary))
    stop("observed responses and prediction rows differ")
  error <- observed - mean_summary$mean
  covered <- observed >= response_summary$lower &
    observed <= response_summary$upper
  width <- response_summary$upper - response_summary$lower
  pointwise <- data.frame(
    mean_summary[c("group", "group_id", "group_status", "row")],
    observed = observed,
    predicted_mean = mean_summary$mean,
    residual = error,
    absolute_error = abs(error),
    squared_error = error^2,
    predictive_lower = response_summary$lower,
    predictive_upper = response_summary$upper,
    covered = covered,
    interval_width = width,
    stringsAsFactors = FALSE)
  structure(list(
    schema_version = "predictive-performance-1",
    model = fit$model,
    method = fit$method,
    predictive_target = predictive_target,
    interval = response_prediction$interval,
    statistics = list(
      rmse = sqrt(mean(error^2)),
      mae = mean(abs(error)),
      r2 = coefficient_of_determination(observed, mean_summary$mean),
      empirical_coverage = mean(covered),
      mean_interval_width = mean(width),
      elpd = score$statistics$elpd,
      mean_log_predictive_density =
        score$statistics$mean_log_predictive_density,
      observations = length(observed),
      log_score_units = score$statistics$predictive_units),
    pointwise = pointwise,
    log_score_pointwise = score$pointwise,
    metadata = list(
      prediction_scale = "input",
      point_prediction = "posterior mean of conditional response mean",
      interval_type = "pointwise posterior predictive",
      interpretation = "held-out predictive performance; not in-sample fit or posterior predictive checking")),
    class = c("vb_predictive_performance", "list"))
}

new_ppc_result <- function(model, method, predictive_target, check,
                           metadata = list()) {
  if (!model %in% c("LRM", "HLRM", "CHLRM")) stop("unknown model")
  if (!method %in% c("MCMC", "CAVI", "SVI")) stop("unknown inference method")
  predictive_target <- match.arg(predictive_target,
    c("observed_data", "observed_groups"))
  required <- c("repdata", "t0", "ppp", "mcse_iid", "target",
                "statistic_definition")
  if (!is.list(check) || !all(required %in% names(check)))
    stop("check must be a posterior predictive check result")
  if (!is.list(metadata)) stop("metadata must be a list")
  check$effective_sample_size <- rep.int(nrow(check$repdata), length(check$ppp))
  names(check$effective_sample_size) <- names(check$ppp)
  check$mcse <- check$mcse_iid
  check$mcse_method <- if (method == "MCMC")
    "nominal_iid_pending_chain_adjustment" else "independent_q_draws"
  retained <- c(required, "effective_sample_size", "mcse", "mcse_method")
  out <- list(schema_version = "posterior-predictive-check-2",
              model = model, method = method,
              predictive_target = predictive_target,
              discrepancies = check[retained], metadata = metadata)
  if (!is.null(check$y_rep)) out$y_rep <- check$y_rep
  structure(out, class = c("vb_posterior_predictive_check", "list"))
}

new_prediction_result <- function(model, method, target, level, metadata,
                                  summary, allocation = NULL, draws = NULL) {
  if (!model %in% c("LRM", "HLRM", "CHLRM")) stop("unknown model")
  if (!method %in% c("MCMC", "CAVI", "SVI")) stop("unknown inference method")
  if (!target %in% c("mean", "response")) stop("unknown prediction target")
  level <- validate_prediction_level(level)
  if (!is.list(metadata)) stop("metadata must be a list")
  if (!is.data.frame(summary)) stop("summary must be a data frame")
  if (!is.null(allocation) && !is.data.frame(allocation))
    stop("allocation must be NULL or a data frame")
  if (!is.null(draws) && !is.list(draws)) stop("draws must be NULL or a list")
  out <- list(
    schema_version = "posterior-prediction-2",
    model = model,
    method = method,
    target = target,
    interval = list(
      level = level,
      type = if (target == "mean") "pointwise posterior credible" else
        "pointwise posterior predictive"
    ),
    metadata = metadata,
    summary = summary)
  if (!is.null(allocation)) out$allocation <- allocation
  if (!is.null(draws)) out$draws <- draws
  structure(out, class = c("vb_prediction_result", "list"))
}
