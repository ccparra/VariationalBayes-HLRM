#' @rdname inference-controls
#' @export
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

#' @rdname inference-controls
#' @export
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

#' @rdname inference-controls
#' @export
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

#' @rdname inference-controls
#' @export
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

#' @rdname inference-controls
#' @export
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

#' @rdname inference-controls
#' @export
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

#' @rdname inference-controls
#' @export
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
