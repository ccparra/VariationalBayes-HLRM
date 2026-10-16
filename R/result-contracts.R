# Uniform result schemas shared by LRM, HLRM and CHLRM. Each model exposes only
# mathematically applicable quantities; nonapplicable fields are absent rather
# than represented by NULL or fabricated estimates.

inference_method <- function(run_settings, family = c("MCMC", "CAVI")) {
  family <- match.arg(family)
  if (family == "MCMC") return("MCMC")
  if (isTRUE(run_settings$stochastic)) "SVI" else "CAVI"
}

new_parameter_summary <- function(mean, sd, lower, upper, level = 0.95,
                                  distribution, mean_exists = NULL) {
  if (!is.numeric(level) || length(level) != 1L || !is.finite(level) ||
      level <= 0 || level >= 1) stop("level must be strictly between zero and one")
  if (!is.character(distribution) || length(distribution) != 1L ||
      !nzchar(distribution)) stop("distribution must be a nonempty string")
  template <- mean
  fill_like <- function(value) {
    if (length(value) == length(template)) {
      dim(value) <- dim(template)
      dimnames(value) <- dimnames(template)
      return(value)
    }
    stop("summary fields must have compatible dimensions")
  }
  mean <- as.numeric(mean); dim(mean) <- dim(template); dimnames(mean) <- dimnames(template)
  sd <- fill_like(as.numeric(sd))
  lower <- fill_like(as.numeric(lower))
  upper <- fill_like(as.numeric(upper))
  if (is.null(mean_exists)) mean_exists <- is.finite(mean)
  mean_exists <- fill_like(as.logical(mean_exists))
  list(mean = mean, sd = sd, lower = lower, upper = upper,
       level = level, distribution = distribution,
       mean_exists = mean_exists,
       contains_zero = lower <= 0 & 0 <= upper)
}

summarize_draws_parameter <- function(draws, margin = NULL, level = 0.95,
                                      distribution = "posterior draws") {
  if (!is.numeric(draws) || !length(draws) || any(!is.finite(draws)))
    stop("draws must be a nonempty finite numeric object")
  alpha <- (1 - level) / 2
  summarize <- function(x) c(
    mean = mean(x), sd = stats::sd(x),
    lower = stats::quantile(x, alpha, names = FALSE),
    upper = stats::quantile(x, 1 - alpha, names = FALSE))
  values <- if (is.null(margin)) summarize(draws) else apply(draws, margin, summarize)
  take <- function(name) {
    if (is.null(margin)) return(unname(values[name]))
    subs <- c(list(values, name), rep(list(TRUE), length(dim(values)) - 1L),
              list(drop = TRUE))
    do.call(`[`, subs)
  }
  new_parameter_summary(take("mean"), take("sd"), take("lower"), take("upper"),
                        level, distribution)
}

summarize_normal_parameter <- function(mean, covariance, level = 0.95) {
  mean <- as.numeric(mean)
  covariance <- as.matrix(covariance)
  if (!identical(dim(covariance), as.integer(c(length(mean), length(mean)))))
    stop("covariance dimension must match mean")
  chol_spd(covariance, "normal covariance")
  standard_deviation <- sqrt(diag(covariance))
  critical <- stats::qnorm(1 - (1 - level) / 2)
  new_parameter_summary(mean, standard_deviation,
                        mean - critical * standard_deviation,
                        mean + critical * standard_deviation,
                        level, "normal variational factor")
}

summarize_independent_normals <- function(mean, covariances, level = 0.95) {
  mean <- as.matrix(mean)
  if (!is.list(covariances) || length(covariances) != nrow(mean))
    stop("covariances must contain one matrix per row of mean")
  standard_deviation <- t(vapply(seq_len(nrow(mean)), function(i) {
    covariance <- as.matrix(covariances[[i]])
    chol_spd(covariance, paste0("normal covariance ", i))
    sqrt(diag(covariance))
  }, numeric(ncol(mean))))
  critical <- stats::qnorm(1 - (1 - level) / 2)
  new_parameter_summary(mean, standard_deviation,
                        mean - critical * standard_deviation,
                        mean + critical * standard_deviation,
                        level, "independent normal variational factors")
}

summarize_inverse_gamma_parameter <- function(shape, scale, level = 0.95) {
  shape <- as.numeric(shape); scale <- as.numeric(scale)
  if (length(shape) != length(scale) || any(!is.finite(shape)) ||
      any(!is.finite(scale)) || any(shape <= 0) || any(scale <= 0))
    stop("shape and scale must be positive finite vectors of equal length")
  mean_exists <- shape > 1
  mean <- rep(NA_real_, length(shape)); mean[mean_exists] <- scale[mean_exists] / (shape[mean_exists] - 1)
  sd <- rep(NA_real_, length(shape)); variance_exists <- shape > 2
  sd[variance_exists] <- scale[variance_exists] /
    ((shape[variance_exists] - 1) * sqrt(shape[variance_exists] - 2))
  alpha <- (1 - level) / 2
  lower <- 1 / stats::qgamma(1 - alpha, shape = shape, rate = scale)
  upper <- 1 / stats::qgamma(alpha, shape = shape, rate = scale)
  new_parameter_summary(mean, sd, lower, upper, level,
                        "inverse-gamma variational factor", mean_exists)
}

summarize_gamma_parameter <- function(shape, rate, level = 0.95) {
  shape <- as.numeric(shape); rate <- as.numeric(rate)
  if (length(shape) != length(rate) || any(!is.finite(shape)) ||
      any(!is.finite(rate)) || any(shape <= 0) || any(rate <= 0))
    stop("shape and rate must be positive finite vectors of equal length")
  alpha <- (1 - level) / 2
  new_parameter_summary(shape / rate, sqrt(shape) / rate,
                        stats::qgamma(alpha, shape, rate = rate),
                        stats::qgamma(1 - alpha, shape, rate = rate),
                        level, "gamma variational factor")
}

summarize_dirichlet_parameter <- function(alpha, level = 0.95) {
  alpha <- as.numeric(alpha)
  if (!length(alpha) || any(!is.finite(alpha)) || any(alpha <= 0))
    stop("alpha must contain positive finite values")
  total <- sum(alpha)
  if (length(alpha) == 1L) {
    return(new_parameter_summary(1, 0, 1, 1, level,
                                 "degenerate one-component simplex"))
  }
  mean <- alpha / total
  sd <- sqrt(alpha * (total - alpha) / (total^2 * (total + 1)))
  tail <- (1 - level) / 2
  lower <- vapply(seq_along(alpha), function(k)
    stats::qbeta(tail, alpha[k], total - alpha[k]), numeric(1))
  upper <- vapply(seq_along(alpha), function(k)
    stats::qbeta(1 - tail, alpha[k], total - alpha[k]), numeric(1))
  new_parameter_summary(mean, sd, lower, upper, level,
                        "Dirichlet variational factor")
}

summarize_inverse_wishart_mean <- function(df, scale, level = 0.95) {
  scale <- as.matrix(scale); p <- nrow(scale)
  chol_spd(scale, "inverse-Wishart scale")
  mean_exists <- df > p + 1
  mean <- if (mean_exists) scale / (df - p - 1) else matrix(NA_real_, p, p)
  missing_moment <- matrix(NA_real_, p, p)
  new_parameter_summary(mean, missing_moment, missing_moment, missing_moment,
                        level, "inverse-Wishart variational factor",
                        matrix(mean_exists, p, p))
}

parameter_names_for_model <- function(model) {
  switch(model,
    LRM = c("coefficients", "residual_variance"),
    HLRM = c("population_coefficients", "population_covariance",
             "variance_scale", "group_coefficients", "group_variances"),
    CHLRM = c("population_coefficients", "population_covariance",
              "variance_scale", "component_coefficients",
              "component_variances", "mixture_weights"),
    stop("unknown model"))
}

new_parameter_slots <- function(model) {
  stats::setNames(vector("list", length(parameter_names_for_model(model))),
                  parameter_names_for_model(model))
}

posterior_draw_parameter <- function(draws) {
  list(distribution = "posterior draws", draws = draws)
}

normal_variational_parameter <- function(mean, covariance,
                                         independent = FALSE) {
  list(distribution = if (independent) "independent normal variational factors" else
         "normal variational factor",
       mean = mean, covariance = covariance)
}

inverse_gamma_variational_parameter <- function(shape, scale) {
  list(distribution = "inverse-gamma variational factor",
       shape = shape, scale = scale)
}

canonical_fit_parameters <- function(model, method, native_parameters) {
  parameters <- new_parameter_slots(model)
  if (method == "MCMC") {
    if (model == "LRM") {
      parameters$coefficients <- posterior_draw_parameter(native_parameters$Beta)
      parameters$residual_variance <- posterior_draw_parameter(native_parameters$sigma2)
    } else {
      parameters$population_coefficients <- posterior_draw_parameter(native_parameters$Beta)
      parameters$population_covariance <- posterior_draw_parameter(native_parameters$Sigma)
      parameters$variance_scale <- posterior_draw_parameter(native_parameters$xi2)
      if (model == "HLRM") {
        parameters$group_coefficients <- posterior_draw_parameter(native_parameters$Betaj)
        parameters$group_variances <- posterior_draw_parameter(native_parameters$sigma2j)
      } else {
        parameters$component_coefficients <- posterior_draw_parameter(native_parameters$Betak)
        parameters$component_variances <- posterior_draw_parameter(native_parameters$sigma2k)
        parameters$mixture_weights <- posterior_draw_parameter(native_parameters$omega)
      }
    }
    return(parameters)
  }
  if (model == "LRM") {
    parameters$coefficients <- normal_variational_parameter(
      native_parameters$MuBeta, native_parameters$SigmaBeta)
    parameters$residual_variance <- inverse_gamma_variational_parameter(
      native_parameters$asigma2, native_parameters$bsigma2)
  } else {
    parameters$population_coefficients <- normal_variational_parameter(
      native_parameters$MuBeta, native_parameters$SigmaBeta)
    parameters$population_covariance <- list(
      distribution = "inverse-Wishart variational factor",
      df = native_parameters$nuS, scale = native_parameters$SS)
    parameters$variance_scale <- list(
      distribution = "gamma variational factor",
      shape = native_parameters$axi2, rate = native_parameters$bxi2)
    if (model == "HLRM") {
      parameters$group_coefficients <- normal_variational_parameter(
        native_parameters$MuBetaj, native_parameters$SigmaBetaj,
        independent = TRUE)
      parameters$group_variances <- inverse_gamma_variational_parameter(
        native_parameters$aj, native_parameters$bj)
    } else {
      parameters$component_coefficients <- normal_variational_parameter(
        native_parameters$MuBetak, native_parameters$SigmaBetak,
        independent = TRUE)
      parameters$component_variances <- inverse_gamma_variational_parameter(
        native_parameters$ak, native_parameters$bk)
      parameters$mixture_weights <- list(
        distribution = "Dirichlet variational factor",
        concentration = native_parameters$alpha)
    }
  }
  parameters
}

canonical_fit_latent_variables <- function(model, method, native_parameters) {
  if (model != "CHLRM") return(NULL)
  if (method == "MCMC") list(group_allocations = native_parameters$z) else
    list(allocation_probabilities = native_parameters$rho)
}

canonical_fit_objective <- function(method, native_parameters) {
  if (method == "MCMC") {
    return(list(log_likelihood = native_parameters$loglik,
                log_joint = native_parameters$logjoint))
  }
  list(elbo = native_parameters$elbo,
       evaluation_iterations = native_parameters$elbo_iter)
}

new_posterior_summary_result <- function(model, method, parameters,
                                         allocations = NULL,
                                         relabeling = NULL, level = 0.95) {
  models <- c("LRM", "HLRM", "CHLRM")
  methods <- c("MCMC", "CAVI", "SVI")
  if (!model %in% models) stop("unknown model")
  if (!method %in% methods) stop("unknown inference method")
  required <- parameter_names_for_model(model)
  if (!is.list(parameters) || !identical(names(parameters), required))
    stop("parameters must use the canonical slots for the selected model")
  if (any(vapply(parameters, is.null, logical(1L))))
    stop("applicable parameters cannot be NULL")
  out <- list(
    schema_version = "posterior-summary-1",
    model = model,
    method = method,
    interval = list(level = level, type = "pointwise credible"),
    parameters = parameters)
  if (!is.null(allocations)) out$allocations <- allocations
  if (!is.null(relabeling)) out$relabeling <- relabeling
  structure(out, class = c("vb_posterior_summary", "list"))
}

new_posterior_draws_result <- function(fit, sampled, draw_index, seed,
                                       estimated_storage_mb) {
  if (!inherits(fit, "vb_model_fit"))
    stop("fit must be produced by fit_model()")
  if (!is.list(sampled) || !is.list(sampled$parameters) ||
      !identical(names(sampled$parameters), parameter_names_for_model(fit$model)))
    stop("posterior draws must match the parameters of the selected model")
  if (!is.data.frame(draw_index) || !nrow(draw_index))
    stop("draw_index must be a nonempty data frame")
  out <- list(
    schema_version = "posterior-draws-1",
    model = fit$model,
    method = fit$method,
    distribution = if (fit$method == "MCMC") "posterior" else
      "variational_approximation",
    parameters = sampled$parameters,
    draw_index = draw_index,
    metadata = list(
      draw_count = nrow(draw_index), draw_axis = 1L, seed = seed,
      estimated_storage_mb = estimated_storage_mb))
  if (!is.null(sampled$latent_variables))
    out$latent_variables <- sampled$latent_variables
  structure(out, class = c("vb_posterior_draws", "list"))
}

new_fit_result <- function(model, method, representation, parameters,
                           diagnostics, total_seconds, cumulative_seconds,
                           evaluation_iterations = NULL,
                           update_seconds = NULL, objective_seconds = NULL,
                           metadata = list()) {
  if (!model %in% c("LRM", "HLRM", "CHLRM")) stop("unknown model")
  if (!method %in% c("MCMC", "CAVI", "SVI")) stop("unknown inference method")
  if (!representation %in% c("posterior draws", "variational parameters"))
    stop("unknown posterior representation")
  if (!is.list(parameters) || !length(parameters))
    stop("native parameters must be a nonempty list")
  if (!is.list(diagnostics)) stop("diagnostics must be a list")
  if (!is.numeric(total_seconds) || length(total_seconds) != 1L ||
      !is.finite(total_seconds) || total_seconds < 0)
    stop("total_seconds must be a nonnegative finite scalar")
  native_parameters <- parameters
  parameters <- canonical_fit_parameters(model, method, native_parameters)
  timing <- list(total_seconds = total_seconds,
                 cumulative_seconds = cumulative_seconds)
  if (!is.null(evaluation_iterations))
    timing$evaluation_iterations <- evaluation_iterations
  if (!is.null(update_seconds)) timing$update_seconds <- update_seconds
  if (!is.null(objective_seconds)) timing$objective_seconds <- objective_seconds
  out <- list(
    schema_version = "fit-result-1",
    model = model,
    method = method,
    status = "fitted",
    representation = representation,
    parameters = parameters,
    objective = canonical_fit_objective(method, native_parameters),
    diagnostics = diagnostics,
    timing = timing,
    metadata = metadata
  )
  latent_variables <- canonical_fit_latent_variables(
    model, method, native_parameters)
  if (!is.null(latent_variables)) out$latent_variables <- latent_variables
  structure(out, class = c("vb_fit_result", "list"))
}

new_convergence_result <- function(method, status, reason, work_completed,
                                   checks_completed, metrics = list(),
                                   thresholds = list(), history = data.frame(),
                                   work_unit) {
  if (!method %in% c("MCMC", "CAVI", "SVI"))
    stop("unknown inference method")
  if (!status %in% c("converged", "max_work", "failed", "interrupted"))
    stop("unknown convergence status")
  if (!is.character(reason) || length(reason) != 1L || !nzchar(reason))
    stop("reason must be a nonempty string")
  work_completed <- validate_integer_scalar(
    work_completed, "work_completed", 0L)
  checks_completed <- validate_integer_scalar(
    checks_completed, "checks_completed", 0L)
  if (!is.list(metrics) || !is.list(thresholds) || !is.data.frame(history))
    stop("metrics and thresholds must be lists and history a data frame")
  if (!is.character(work_unit) || length(work_unit) != 1L || !nzchar(work_unit))
    stop("work_unit must be a nonempty string")
  structure(list(
    schema_version = "convergence-1", method = method, status = status,
    reason = reason, work_completed = work_completed,
    work_unit = work_unit, checks_completed = checks_completed,
    metrics = metrics, thresholds = thresholds, history = history
  ), class = c("vb_convergence", "list"))
}

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
      r2 = r2(observed, mean_summary$mean),
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
