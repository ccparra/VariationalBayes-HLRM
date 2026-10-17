.vi_convergence_result <- function(diagnostics, control) {
  converged <- isTRUE(diagnostics$converged)
  parameter_threshold <- control$stopping$parameter_tolerance
  metrics <- list(
    objective_relative_change = diagnostics$relative_change,
    stable_checks = diagnostics$stable_checks,
    smoothed_objective = diagnostics$smoothed_elbo)
  thresholds <- list(
    objective_tolerance = control$stopping$objective_tolerance,
    min_iterations = control$stopping$min_iterations,
    patience = control$stopping$patience,
    check_every = control$stopping$check_every,
    smoothing_window = control$stopping$smoothing_window)
  if (!is.null(parameter_threshold)) {
    metrics$parameter_relative_change <- diagnostics$parameter_change
    thresholds$parameter_tolerance <- parameter_threshold
  }
  new_convergence_result(
    control$method, if (converged) "converged" else "max_work",
    if (converged) {
      if (is.null(parameter_threshold)) "objective_stable" else
        "objective_and_parameters_stable"
    } else "max_iterations",
    diagnostics$iterations, nrow(diagnostics$history), metrics, thresholds,
    diagnostics$history, "iterations")
}

.single_chain_convergence_result <- function(control) {
  new_convergence_result(
    "MCMC", "max_work", "single_chain_fixed_budget",
    control$post_warmup_draws, 0L,
    metrics = list(), thresholds = list(), history = data.frame(),
    work_unit = "post_warmup_draws_per_chain")
}

.multichain_convergence_result <- function(orchestration, control) {
  variables <- orchestration$diagnostics$variables
  finite_max <- function(x) if (any(is.finite(x)))
    max(x[is.finite(x)]) else NA_real_
  finite_min <- function(x) if (any(is.finite(x)))
    min(x[is.finite(x)]) else NA_real_
  metrics <- list(
    max_rhat = finite_max(variables$rhat),
    min_ess_bulk = finite_min(variables$ess_bulk),
    min_ess_tail = finite_min(variables$ess_tail),
    max_relative_mcse_mean = finite_max(variables$relative_mcse_mean),
    max_relative_mcse_sd = finite_max(variables$relative_mcse_sd),
    variables = variables)
  stopping <- control$stopping
  new_convergence_result(
    "MCMC",
    if (isTRUE(orchestration$diagnostics$passed)) "converged" else "max_work",
    if (isTRUE(orchestration$diagnostics$passed)) "diagnostic_thresholds" else
      "max_post_warmup_draws",
    orchestration$run_settings$samples,
    nrow(orchestration$block_history), metrics,
    list(rhat_max = stopping$rhat_max,
         ess_bulk_min = stopping$ess_bulk_min,
         ess_tail_min = stopping$ess_tail_min,
         relative_mcse_mean_max = stopping$relative_mcse_mean_max,
         relative_mcse_sd_max = stopping$relative_mcse_sd_max,
         min_post_warmup_draws = stopping$min_post_warmup_draws,
         check_every = stopping$check_every),
    orchestration$block_history, "post_warmup_draws_per_chain")
}

.typed_vi_initial_states <- function(model, data, native_prior,
                                     native_control, initialization,
                                     expected_method, n_starts, seeds,
                                     lrm_ols, grouped_ols) {
  states <- vector("list", n_starts)
  for (start in seq_len(n_starts)) {
    settings <- native_control
    settings$seed <- seeds[start]
    backend <- .new_backend(
      model, tolower(expected_method), data, native_prior, settings)
    .apply_initialization(
      backend, model$model, expected_method, initialization,
      lrm_ols, grouped_ols, model$component_count)
    field <- .backend_state_field(model$model, expected_method)
    states[[start]] <- unserialize(serialize(backend[[field]]$params, NULL))
  }
  dimensions <- .vi_data_dimensions(
    model$model, data$y, data$X,
    if (model$model == "CHLRM") model$component_count else NULL)
  structure(list(
    schema_version = "vi-initial-states-1", model = model$model,
    dimensions = dimensions, states = states,
    metadata = data.frame(start = seq_len(n_starts),
                          initialization_seed = seeds),
    initialization = initialization$strategy
  ), class = "vi_initial_states")
}

.rebuild_vi_fit_result <- function(backend, attempt, model, method) {
  timing <- attempt$timing
  stochastic <- method == "SVI"
  state <- backend[[.backend_state_field(model, method)]]
  units <- if (model == "LRM") state$stats$n else state$stats$m
  new_fit_result(
    model, method,
    "variational parameters", backend$results, attempt$diagnostics,
    timing$total_seconds, timing$cumulative_seconds,
    evaluation_iterations = timing$evaluation_iterations,
    update_seconds = timing$update_seconds,
    objective_seconds = timing$objective_seconds,
    metadata = list(
      sampling_unit = if (!stochastic) "full_data" else
        if (model == "LRM") "observation" else "group",
      minibatch_size = if (stochastic)
        backend$run_settings$minibatch_m else units))
}

.finalize_typed_fit <- function(result, model, prior, initialization, control,
                                seed, backend, convergence,
                                orchestration = NULL, backends = NULL,
                                scaling = NULL, scaling_map = NULL,
                                resource_estimate = NULL) {
  result$schema_version <- "fit-result-2"
  result$method <- control$method
  result$diagnostics <- convergence
  result$specifications <- list(
    model = model, prior = prior, initialization = initialization,
    control = control)
  if (!is.null(scaling)) result$specifications$scaling <- scaling
  result$metadata$seed <- seed
  result$metadata$backend <- class(backend)[1L]
  result$metadata$parameter_scale <- "input"
  result$metadata$objective_scale <- if (is.null(scaling_map)) "input" else
    "standardized_backend"
  if (!is.null(resource_estimate))
    result$metadata$resource_estimate <- resource_estimate
  if (!is.null(orchestration)) {
    result$metadata$orchestration_schema <- orchestration$schema_version
    if (inherits(orchestration, "mcmc_multichain_fit")) {
      result$metadata$chains <- nrow(orchestration$chain_metadata)
      result$metadata$chain_metadata <- orchestration$chain_metadata
    } else {
      result$metadata$starts <- nrow(orchestration$summary)
      result$metadata$selected_start <- orchestration$selected_start
      result$metadata$start_summary <- orchestration$summary
    }
  }
  attr(result, "backend") <- backend
  if (!is.null(backends)) attr(result, "backends") <- backends
  if (!is.null(orchestration)) attr(result, "orchestration") <- orchestration
  if (!is.null(scaling_map)) {
    result$parameters <- restore_fit_parameters(
      result$parameters, model$model, control$method, scaling_map)
    attr(result, "scaling_map") <- scaling_map
  }
  class(result) <- c("vb_model_fit", "vb_fit_result", "list")
  result
}

