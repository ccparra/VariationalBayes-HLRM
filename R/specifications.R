# Typed, model-level specifications. These objects contain no data-dependent
# estimates and no backend notation; translation to the current kernels lives
# in fit-model.R.

.new_vb_spec <- function(fields, classes) {
  structure(fields, class = c(classes, "vb_spec", "list"))
}

.optional_count <- function(x, name) {
  if (is.null(x)) return(NULL)
  validate_integer_scalar(x, name, 1L)
}

.logical_scalar <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x))
    stop(name, " must be TRUE or FALSE")
  x
}

.positive_finite_scalar <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x <= 0)
    stop(name, " must be a positive finite scalar")
  as.numeric(x)
}

lrm_model <- function(coefficient_count = NULL) {
  .new_vb_spec(list(
    schema_version = "model-spec-1",
    model = "LRM",
    coefficient_count = .optional_count(coefficient_count, "coefficient_count"),
    likelihood = "gaussian_linear_regression"
  ), c("vb_lrm_model_spec", "vb_model_spec"))
}

hlrm_model <- function(coefficient_count = NULL) {
  .new_vb_spec(list(
    schema_version = "model-spec-1",
    model = "HLRM",
    coefficient_count = .optional_count(coefficient_count, "coefficient_count"),
    likelihood = "gaussian_hierarchical_linear_regression"
  ), c("vb_hlrm_model_spec", "vb_model_spec"))
}

chlrm_model <- function(component_count, coefficient_count = NULL) {
  .new_vb_spec(list(
    schema_version = "model-spec-1",
    model = "CHLRM",
    component_count = validate_integer_scalar(
      component_count, "component_count", 1L),
    coefficient_count = .optional_count(coefficient_count, "coefficient_count"),
    likelihood = "clustered_gaussian_hierarchical_linear_regression"
  ), c("vb_chlrm_model_spec", "vb_model_spec"))
}

.validate_prior_vector <- function(x, name) {
  x <- as.numeric(x)
  if (!length(x) || any(!is.finite(x)))
    stop(name, " must be a nonempty finite numeric vector")
  x
}

lrm_prior <- function(coefficient_mean, coefficient_covariance,
                      residual_degrees_of_freedom,
                      residual_variance_scale) {
  coefficient_mean <- .validate_prior_vector(
    coefficient_mean, "coefficient_mean")
  p <- length(coefficient_mean)
  coefficient_covariance <- .spd_matrix(
    as.matrix(coefficient_covariance), p, "coefficient_covariance")
  residual_degrees_of_freedom <- .positive_scalar(
    residual_degrees_of_freedom, "residual_degrees_of_freedom")
  residual_variance_scale <- .positive_scalar(
    residual_variance_scale, "residual_variance_scale")
  .new_vb_spec(list(
    schema_version = "prior-spec-1",
    model = "LRM",
    coefficient_mean = coefficient_mean,
    coefficient_covariance = coefficient_covariance,
    residual_degrees_of_freedom = residual_degrees_of_freedom,
    residual_variance_scale = residual_variance_scale,
    parameterization = list(
      coefficients = "multivariate_normal_covariance",
      residual_variance = "inverse_gamma_shape_scale")
  ), c("vb_lrm_prior_spec", "vb_prior_spec"))
}

.hierarchical_prior <- function(model, coefficient_mean,
                                coefficient_covariance,
                                coefficient_dispersion_df,
                                coefficient_dispersion_scale,
                                residual_degrees_of_freedom,
                                residual_variance_scale_shape,
                                residual_variance_scale_rate,
                                mixture_concentration = NULL) {
  coefficient_mean <- .validate_prior_vector(
    coefficient_mean, "coefficient_mean")
  p <- length(coefficient_mean)
  coefficient_covariance <- .spd_matrix(
    as.matrix(coefficient_covariance), p, "coefficient_covariance")
  coefficient_dispersion_scale <- .spd_matrix(
    as.matrix(coefficient_dispersion_scale), p,
    "coefficient_dispersion_scale")
  coefficient_dispersion_df <- .positive_scalar(
    coefficient_dispersion_df, "coefficient_dispersion_df")
  if (coefficient_dispersion_df <= p - 1L)
    stop("coefficient_dispersion_df must exceed coefficient_count - 1 for a proper inverse-Wishart prior")
  residual_degrees_of_freedom <- .positive_scalar(
    residual_degrees_of_freedom, "residual_degrees_of_freedom")
  residual_variance_scale_shape <- .positive_scalar(
    residual_variance_scale_shape, "residual_variance_scale_shape")
  residual_variance_scale_rate <- .positive_scalar(
    residual_variance_scale_rate, "residual_variance_scale_rate")
  out <- list(
    schema_version = "prior-spec-1",
    model = model,
    coefficient_mean = coefficient_mean,
    coefficient_covariance = coefficient_covariance,
    coefficient_dispersion_df = coefficient_dispersion_df,
    coefficient_dispersion_scale = coefficient_dispersion_scale,
    residual_degrees_of_freedom = residual_degrees_of_freedom,
    residual_variance_scale_shape = residual_variance_scale_shape,
    residual_variance_scale_rate = residual_variance_scale_rate,
    parameterization = list(
      population_coefficients = "multivariate_normal_covariance",
      coefficient_dispersion = "inverse_wishart_scale",
      local_residual_variance = "inverse_gamma_shape_scale",
      residual_variance_scale = "gamma_shape_rate")
  )
  if (model == "CHLRM") {
    mixture_concentration <- .validate_prior_vector(
      mixture_concentration, "mixture_concentration")
    if (any(mixture_concentration <= 0))
      stop("mixture_concentration entries must be positive")
    out$mixture_concentration <- mixture_concentration
    out$parameterization$mixture_weights <- "dirichlet_concentration"
  }
  .new_vb_spec(out, c(
    if (model == "HLRM") "vb_hlrm_prior_spec" else "vb_chlrm_prior_spec",
    "vb_prior_spec"))
}

hlrm_prior <- function(coefficient_mean, coefficient_covariance,
                       coefficient_dispersion_df,
                       coefficient_dispersion_scale,
                       residual_degrees_of_freedom,
                       residual_variance_scale_shape,
                       residual_variance_scale_rate) {
  .hierarchical_prior(
    "HLRM", coefficient_mean, coefficient_covariance,
    coefficient_dispersion_df, coefficient_dispersion_scale,
    residual_degrees_of_freedom, residual_variance_scale_shape,
    residual_variance_scale_rate)
}

chlrm_prior <- function(coefficient_mean, coefficient_covariance,
                        coefficient_dispersion_df,
                        coefficient_dispersion_scale,
                        residual_degrees_of_freedom,
                        residual_variance_scale_shape,
                        residual_variance_scale_rate,
                        mixture_concentration) {
  .hierarchical_prior(
    "CHLRM", coefficient_mean, coefficient_covariance,
    coefficient_dispersion_df, coefficient_dispersion_scale,
    residual_degrees_of_freedom, residual_variance_scale_shape,
    residual_variance_scale_rate, mixture_concentration)
}

.inverse_gamma_moments <- function(shape, scale) {
  mean_exists <- shape > 1
  variance_exists <- shape > 2
  list(
    mean = if (mean_exists) scale / (shape - 1) else NA_real_,
    variance = if (variance_exists)
      scale^2 / ((shape - 1)^2 * (shape - 2)) else NA_real_,
    mean_exists = mean_exists,
    variance_exists = variance_exists,
    parameterization = "inverse-gamma shape-scale")
}

.gamma_moments <- function(shape, rate) list(
  mean = shape / rate,
  variance = shape / rate^2,
  mean_exists = TRUE,
  variance_exists = TRUE,
  parameterization = "gamma shape-rate")

.dirichlet_moments <- function(concentration) {
  total <- sum(concentration)
  mean <- concentration / total
  covariance <- -tcrossprod(concentration) /
    (total^2 * (total + 1))
  diag(covariance) <- concentration * (total - concentration) /
    (total^2 * (total + 1))
  list(mean = mean, covariance = covariance,
       mean_exists = TRUE, covariance_exists = TRUE,
       parameterization = "Dirichlet concentration")
}

prior_moments <- function(prior) {
  if (!inherits(prior, "vb_prior_spec"))
    stop("prior must be created by lrm_prior(), hlrm_prior() or chlrm_prior()")
  coefficient_moments <- list(
    mean = prior$coefficient_mean,
    covariance = prior$coefficient_covariance,
    mean_exists = TRUE,
    covariance_exists = TRUE,
    parameterization = "multivariate normal covariance")
  if (prior$model == "LRM") {
    residual <- .inverse_gamma_moments(
      prior$residual_degrees_of_freedom / 2,
      prior$residual_degrees_of_freedom *
        prior$residual_variance_scale / 2)
    return(structure(list(
      schema_version = "prior-moments-1", model = "LRM",
      parameters = list(
        coefficients = coefficient_moments,
        residual_variance = residual)),
      class = c("vb_prior_moments", "list")))
  }

  p <- length(prior$coefficient_mean)
  covariance_mean_exists <- prior$coefficient_dispersion_df > p + 1
  covariance_mean <- if (covariance_mean_exists)
    prior$coefficient_dispersion_scale /
      (prior$coefficient_dispersion_df - p - 1) else
        matrix(NA_real_, p, p)
  population_covariance <- list(
    mean = covariance_mean,
    mean_exists = covariance_mean_exists,
    parameterization = "inverse-Wishart degrees-of-freedom and scale",
    note = "Elementwise variances are not summarized by this contract")
  variance_scale <- .gamma_moments(
    prior$residual_variance_scale_shape,
    prior$residual_variance_scale_rate)

  nu <- prior$residual_degrees_of_freedom
  local_variance_mean_exists <- nu > 2
  local_variance_variance_exists <- nu > 4
  expected_xi <- variance_scale$mean
  expected_xi_squared <- variance_scale$variance + expected_xi^2
  local_variance_mean <- if (local_variance_mean_exists)
    nu * expected_xi / (nu - 2) else NA_real_
  local_variance_second <- if (local_variance_variance_exists)
    nu^2 * expected_xi_squared / ((nu - 2) * (nu - 4)) else NA_real_
  local_variance <- list(
    mean = local_variance_mean,
    variance = if (local_variance_variance_exists)
      local_variance_second - local_variance_mean^2 else NA_real_,
    mean_exists = local_variance_mean_exists,
    variance_exists = local_variance_variance_exists,
    parameterization = paste0(
      "marginal under inverse-gamma residual variance conditional on ",
      "gamma variance_scale"))
  local_coefficients <- list(
    mean = prior$coefficient_mean,
    covariance = if (covariance_mean_exists)
      prior$coefficient_covariance + covariance_mean else
        matrix(NA_real_, p, p),
    mean_exists = TRUE,
    covariance_exists = covariance_mean_exists,
    parameterization = "marginal normal hierarchy")

  parameters <- list(
    population_coefficients = coefficient_moments,
    population_covariance = population_covariance,
    variance_scale = variance_scale)
  local_name <- if (prior$model == "HLRM")
    "group_coefficients" else "component_coefficients"
  variance_name <- if (prior$model == "HLRM")
    "group_variance" else "component_variance"
  parameters[[local_name]] <- local_coefficients
  parameters[[variance_name]] <- local_variance
  if (prior$model == "CHLRM") parameters$mixture_weights <-
    .dirichlet_moments(prior$mixture_concentration)
  structure(list(
    schema_version = "prior-moments-1", model = prior$model,
    parameters = parameters), class = c("vb_prior_moments", "list"))
}

.initialization_spec <- function(strategy, values = NULL, jitter = NULL) {
  fields <- list(schema_version = "initialization-spec-1",
                 strategy = strategy)
  if (!is.null(values)) fields$values <- values
  if (!is.null(jitter)) fields$jitter <- jitter
  .new_vb_spec(fields, "vb_initialization_spec")
}

initialization_from_prior <- function() .initialization_spec("from_prior")

initialization_ols <- function(jitter = 1e-6) {
  jitter <- .positive_scalar(jitter, "jitter")
  .initialization_spec("ols", jitter = jitter)
}

initialization_kmeans <- function(jitter = 1e-6) {
  jitter <- .positive_scalar(jitter, "jitter")
  .initialization_spec("kmeans", jitter = jitter)
}

initialization_values <- function(values) {
  if (!is.list(values) || !length(values) || is.null(names(values)) ||
      any(!nzchar(names(values))) || anyDuplicated(names(values)) ||
      any(vapply(values, is.null, logical(1L))))
    stop("values must be a nonempty named list without NULL entries")
  .initialization_spec("values", values = values)
}

mcmc_stopping <- function(min_post_warmup_draws = 1000L,
                          check_every = 500L, rhat_max = 1.01,
                          ess_bulk_min = 1000, ess_tail_min = 1000,
                          relative_mcse_mean_max = 0.05,
                          relative_mcse_sd_max = 0.10) {
  min_post_warmup_draws <- validate_integer_scalar(
    min_post_warmup_draws, "min_post_warmup_draws", 1L)
  check_every <- validate_integer_scalar(check_every, "check_every", 1L)
  numeric_thresholds <- c(
    rhat_max = rhat_max, ess_bulk_min = ess_bulk_min,
    ess_tail_min = ess_tail_min,
    relative_mcse_mean_max = relative_mcse_mean_max,
    relative_mcse_sd_max = relative_mcse_sd_max)
  if (any(!is.finite(numeric_thresholds)) || rhat_max < 1 ||
      ess_bulk_min <= 0 || ess_tail_min <= 0 ||
      relative_mcse_mean_max <= 0 || relative_mcse_sd_max <= 0)
    stop("MCMC thresholds require rhat_max >= 1 and positive ESS/relative MCSE limits")
  .new_vb_spec(list(
    schema_version = "mcmc-stopping-1",
    min_post_warmup_draws = min_post_warmup_draws,
    check_every = check_every,
    rhat_max = as.numeric(rhat_max),
    ess_bulk_min = as.numeric(ess_bulk_min),
    ess_tail_min = as.numeric(ess_tail_min),
    relative_mcse_mean_max = as.numeric(relative_mcse_mean_max),
    relative_mcse_sd_max = as.numeric(relative_mcse_sd_max)
  ), "vb_mcmc_stopping")
}

cavi_stopping <- function(objective_tolerance = 1e-8,
                          parameter_tolerance = NULL,
                          min_iterations = 1L, patience = 5L,
                          check_every = 1L) {
  if (!is.numeric(objective_tolerance) ||
      length(objective_tolerance) != 1L ||
      !is.finite(objective_tolerance) || objective_tolerance < 0)
    stop("objective_tolerance must be a nonnegative finite scalar")
  if (!is.null(parameter_tolerance) &&
      (!is.numeric(parameter_tolerance) || length(parameter_tolerance) != 1L ||
       !is.finite(parameter_tolerance) || parameter_tolerance < 0))
    stop("parameter_tolerance must be NULL or a nonnegative finite scalar")
  check_every <- validate_integer_scalar(check_every, "check_every", 1L)
  if (check_every != 1L)
    stop("CAVI currently evaluates its deterministic objective every iteration")
  .new_vb_spec(list(
    schema_version = "cavi-stopping-1",
    objective_tolerance = as.numeric(objective_tolerance),
    parameter_tolerance = parameter_tolerance,
    min_iterations = validate_integer_scalar(
      min_iterations, "min_iterations", 1L),
    patience = validate_integer_scalar(patience, "patience", 1L),
    check_every = check_every,
    smoothing_window = 1L
  ), "vb_cavi_stopping")
}

svi_stopping <- function(objective_tolerance = 1e-4,
                         parameter_tolerance = NULL,
                         min_iterations = 1000L, patience = 5L,
                         check_every = 50L, smoothing_window = 5L) {
  if (!is.numeric(objective_tolerance) ||
      length(objective_tolerance) != 1L ||
      !is.finite(objective_tolerance) || objective_tolerance < 0)
    stop("objective_tolerance must be a nonnegative finite scalar")
  if (!is.null(parameter_tolerance) &&
      (!is.numeric(parameter_tolerance) || length(parameter_tolerance) != 1L ||
       !is.finite(parameter_tolerance) || parameter_tolerance < 0))
    stop("parameter_tolerance must be NULL or a nonnegative finite scalar")
  .new_vb_spec(list(
    schema_version = "svi-stopping-1",
    objective_tolerance = as.numeric(objective_tolerance),
    parameter_tolerance = parameter_tolerance,
    min_iterations = validate_integer_scalar(
      min_iterations, "min_iterations", 1L),
    patience = validate_integer_scalar(patience, "patience", 1L),
    check_every = validate_integer_scalar(check_every, "check_every", 1L),
    smoothing_window = validate_integer_scalar(
      smoothing_window, "smoothing_window", 1L)
  ), "vb_svi_stopping")
}

power_schedule <- function(offset = 100, exponent = 0.7) {
  if (!is.numeric(offset) || length(offset) != 1L ||
      !is.finite(offset) || offset < 0)
    stop("offset must be a nonnegative finite scalar")
  if (!is.numeric(exponent) || length(exponent) != 1L ||
      !is.finite(exponent) || exponent <= 0.5 || exponent > 1)
    stop("exponent must satisfy 0.5 < exponent <= 1")
  .new_vb_spec(list(schema_version = "learning-rate-1",
                    family = "power",
                    offset = as.numeric(offset),
                    exponent = as.numeric(exponent)),
               "vb_learning_rate")
}

mcmc_control <- function(warmup = 1000L, post_warmup_draws = 1000L,
                         thinning = 1L, chains = 1L,
                         evaluation_draws = NULL,
                         stopping = NULL, diagnostic_scope = NULL,
                         max_memory_mb = NULL) {
  chains <- validate_integer_scalar(chains, "chains", 1L)
  if (chains == 1L && !is.null(stopping))
    stop("stopping must be NULL for one MCMC chain because R-hat is not defined")
  if (chains > 1L && !inherits(stopping, "vb_mcmc_stopping"))
    stop("multiple MCMC chains require stopping created by mcmc_stopping()")
  draws <- validate_integer_scalar(
    post_warmup_draws, "post_warmup_draws", 2L)
  evaluation_draws <- if (is.null(evaluation_draws)) min(500L, draws) else
    validate_integer_scalar(evaluation_draws, "evaluation_draws", 2L)
  if (evaluation_draws > draws)
    stop("evaluation_draws cannot exceed post_warmup_draws")
  thinning <- validate_integer_scalar(thinning, "thinning", 1L)
  if (chains > 1L && thinning != 1L)
    stop("multiple MCMC chains require thinning = 1 so ESS diagnoses autocorrelation")
  if (chains == 1L && !is.null(diagnostic_scope))
    stop("diagnostic_scope applies only to multiple MCMC chains")
  max_memory_mb <- .positive_finite_scalar(
    if (is.null(max_memory_mb)) 1024 else max_memory_mb, "max_memory_mb")
  if (chains > 1L) {
    if (stopping$min_post_warmup_draws > draws)
      stop("stopping$min_post_warmup_draws cannot exceed post_warmup_draws")
    diagnostic_scope <- if (is.null(diagnostic_scope)) "global" else
      match.arg(diagnostic_scope, c("global", "all"))
  }
  fields <- list(
    schema_version = "inference-control-1",
    method = "MCMC",
    warmup = validate_integer_scalar(warmup, "warmup", 0L),
    post_warmup_draws = draws,
    thinning = thinning,
    chains = chains,
    evaluation_draws = evaluation_draws,
    max_memory_mb = max_memory_mb)
  if (chains > 1L) {
    fields$stopping <- stopping
    fields$diagnostic_scope <- diagnostic_scope
  }
  .new_vb_spec(fields, c("vb_mcmc_control", "vb_inference_control"))
}

cavi_control <- function(max_iterations = 10000L,
                         evaluation_draws = 500L,
                         stopping = cavi_stopping(), starts = 1L,
                         require_converged = NULL, max_memory_mb = NULL) {
  if (!inherits(stopping, "vb_cavi_stopping"))
    stop("stopping must be created by cavi_stopping()")
  max_iterations <- validate_integer_scalar(
    max_iterations, "max_iterations", 1L)
  if (stopping$min_iterations > max_iterations)
    stop("stopping$min_iterations cannot exceed max_iterations")
  starts <- validate_integer_scalar(starts, "starts", 1L)
  if (starts == 1L && !is.null(require_converged))
    stop("require_converged applies only to multiple starts")
  max_memory_mb <- .positive_finite_scalar(
    if (is.null(max_memory_mb)) 1024 else max_memory_mb, "max_memory_mb")
  fields <- list(
    schema_version = "inference-control-1",
    method = "CAVI",
    max_iterations = max_iterations,
    evaluation_draws = validate_integer_scalar(
      evaluation_draws, "evaluation_draws", 2L),
    stopping = stopping,
    starts = starts,
    max_memory_mb = max_memory_mb)
  if (starts > 1L) {
    fields$require_converged <- .logical_scalar(
      if (is.null(require_converged)) FALSE else require_converged,
      "require_converged")
  }
  .new_vb_spec(fields, c("vb_cavi_control", "vb_inference_control"))
}

svi_control <- function(max_iterations = 50000L, minibatch_size,
                        learning_rate = power_schedule(),
                        evaluation_draws = 500L,
                        stopping = svi_stopping(), starts = 1L,
                        require_converged = NULL, max_memory_mb = NULL) {
  if (!inherits(learning_rate, "vb_learning_rate"))
    stop("learning_rate must be created by a learning-rate constructor")
  if (!inherits(stopping, "vb_svi_stopping"))
    stop("stopping must be created by svi_stopping()")
  max_iterations <- validate_integer_scalar(
    max_iterations, "max_iterations", 1L)
  if (stopping$min_iterations > max_iterations)
    stop("stopping$min_iterations cannot exceed max_iterations")
  starts <- validate_integer_scalar(starts, "starts", 1L)
  if (starts == 1L && !is.null(require_converged))
    stop("require_converged applies only to multiple starts")
  max_memory_mb <- .positive_finite_scalar(
    if (is.null(max_memory_mb)) 1024 else max_memory_mb, "max_memory_mb")
  fields <- list(
    schema_version = "inference-control-1",
    method = "SVI",
    max_iterations = max_iterations,
    minibatch_size = validate_integer_scalar(
      minibatch_size, "minibatch_size", 1L),
    learning_rate = learning_rate,
    evaluation_draws = validate_integer_scalar(
      evaluation_draws, "evaluation_draws", 2L),
    stopping = stopping,
    starts = starts,
    max_memory_mb = max_memory_mb)
  if (starts > 1L) {
    fields$require_converged <- .logical_scalar(
      if (is.null(require_converged)) FALSE else require_converged,
      "require_converged")
  }
  .new_vb_spec(fields, c("vb_svi_control", "vb_inference_control"))
}
