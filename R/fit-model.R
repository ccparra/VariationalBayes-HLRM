fit_model <- function(model, data, prior, method = c("mcmc", "cavi", "svi"),
                      initialization = initialization_from_prior(), control,
                      scaling = NULL, seed = NULL, verbose = TRUE) {
  method <- match.arg(tolower(method), c("mcmc", "cavi", "svi"))
  if (!inherits(model, "vb_model_spec"))
    stop("model must be created by lrm_model(), hlrm_model() or chlrm_model()")
  if (!is.null(scaling) && !inherits(scaling, "vb_data_scaling_spec"))
    stop("scaling must be NULL or created by data_scaling()")
  if (!inherits(initialization, "vb_initialization_spec"))
    stop("initialization must be created by an initialization constructor")
  if (!inherits(control, "vb_inference_control"))
    stop("control must be created by mcmc_control(), cavi_control() or svi_control()")
  expected_method <- .canonical_public_method(method)
  if (!identical(control$method, expected_method))
    stop("control class does not match method")
  input_validated <- .validate_fit_data(model, data)
  resource_estimate <- .estimate_fit_resources_validated(
    model, input_validated, control, data, scaling)
  if (!isTRUE(resource_estimate$within_budget))
    stop(sprintf(
      paste0("Estimated peak memory %.1f MB exceeds max_memory_mb %.1f. ",
             "Reduce draws/chains/starts or increase the explicit budget."),
      resource_estimate$estimated_peak_memory_mb,
      resource_estimate$max_memory_mb))
  .validate_prior_for_model(prior, model, input_validated$dimensions)
  scaling_map <- if (is.null(scaling)) NULL else fit_scaling(scaling, data)
  working_data <- if (is.null(scaling_map)) data else
    transform_data(data, scaling_map)
  working_prior <- if (is.null(scaling_map)) prior else
    transform_prior(prior, scaling_map)
  working_initialization <- if (is.null(scaling_map)) initialization else
    transform_initialization(
      initialization, scaling_map, model$model, expected_method)
  validated <- .validate_fit_data(model, working_data)
  .validate_prior_for_model(working_prior, model, validated$dimensions)
  if (method == "svi") {
    sampling_units <- if (model$model == "LRM")
      validated$dimensions$n else validated$dimensions$m
    sampling_name <- if (model$model == "LRM") "observations" else "groups"
    if (control$minibatch_size > sampling_units)
      stop("minibatch_size cannot exceed the number of ", sampling_name)
  }
  if (initialization$strategy == "kmeans" && model$model != "CHLRM")
    stop("initialization_kmeans() applies only to CHLRM")
  if (initialization$strategy == "ols" && model$model == "CHLRM")
    stop("CHLRM requires initialization_kmeans() rather than an unclustered OLS state")

  grouped_ols <- NULL
  lrm_ols <- NULL
  if (initialization$strategy %in% c("ols", "kmeans")) {
    if (model$model == "LRM")
      lrm_ols <- .lrm_ols_initialization(validated, initialization$jitter)
    else grouped_ols <- .group_ols_initialization(
      validated, initialization$jitter)
  }
  native_prior <- .prior_to_native(
    working_prior, model, validated$dimensions, grouped_ols)
  native_control <- .control_to_native(control, seed)
  if (model$model == "CHLRM" && method != "mcmc") {
    native_control$initialization <- if (initialization$strategy == "kmeans")
      "kmeans_ols" else "random_responsibilities"
  }

  prior_before <- serialize(prior, NULL)
  orchestration <- NULL
  backends <- NULL
  if (method == "mcmc" && control$chains > 1L) {
    if (is.null(seed)) stop("seed is required for reproducible multichain MCMC")
    initialize_chain <- function(backend, chain_index) {
      .apply_initialization(
        backend, model$model, expected_method, working_initialization,
        lrm_ols, grouped_ols, model$component_count)
    }
    stopping <- control$stopping
    orchestration <- run_mcmc_chains(
      model$model, validated$y, validated$X, native_prior, native_control,
      K = if (model$model == "CHLRM") model$component_count else NULL,
      n_chains = control$chains, seed = seed,
      diagnostic_scope = control$diagnostic_scope,
      rhat_max = stopping$rhat_max,
      ess_bulk_min = stopping$ess_bulk_min,
      ess_tail_min = stopping$ess_tail_min,
      relative_mcse_mean_max = stopping$relative_mcse_mean_max,
      relative_mcse_sd_max = stopping$relative_mcse_sd_max,
      max_memory_mb = control$max_memory_mb,
      block_size = min(stopping$check_every, control$post_warmup_draws),
      max_samples = control$post_warmup_draws,
      min_samples = stopping$min_post_warmup_draws,
      initialize = initialize_chain, verbose = verbose)
    backends <- attr(orchestration, "models", exact = TRUE)
    # The reporting backend owns combined draws. It must not alias the first
    # chain, whose terminal sampler/RNG state is needed for exact continuation.
    backend <- backends[[1L]]$clone(deep = TRUE)
    combined <- combine_mcmc_chains(orchestration, require_converged = FALSE)
    backend$chains <- combined
    backend$run_settings$samples <- if (model$model == "CHLRM")
      dim(combined$Betak)[1L] else nrow(combined$Beta)
    backend$fit_metadata <- list(
      total_seconds = sum(orchestration$chain_metadata$elapsed_seconds),
      warmup_iterations = control$warmup,
      stored_samples = backend$run_settings$samples,
      chains = control$chains)
    if (model$model == "CHLRM") {
      backend$relabeled_chains <- combined
      backend$label_reference <- orchestration$label_reference
    }
    backend_result <- new_fit_result(
      model$model, "MCMC", "posterior draws", combined,
      diagnostics = orchestration$diagnostics,
      total_seconds = backend$fit_metadata$total_seconds,
      cumulative_seconds = orchestration$block_history$elapsed_seconds,
      metadata = backend$fit_metadata)
    convergence <- .multichain_convergence_result(orchestration, control)
  } else if (method != "mcmc" && control$starts > 1L) {
    if (is.null(seed)) stop("seed is required for reproducible multistart inference")
    all_seeds <- derive_independent_seeds(2L * control$starts, seed)
    initial_states <- .typed_vi_initial_states(
      model, validated, native_prior, native_control, working_initialization,
      expected_method, control$starts, all_seeds[seq_len(control$starts)],
      lrm_ols, grouped_ols)
    orchestration <- run_vi_multistart(
      model$model, expected_method,
      validated$y, validated$X, native_prior, native_control,
      K = if (model$model == "CHLRM") model$component_count else NULL,
      n_starts = control$starts,
      fit_seeds = all_seeds[control$starts + seq_len(control$starts)],
      initial_states = initial_states,
      require_converged = control$require_converged,
      max_memory_mb = control$max_memory_mb, verbose = verbose)
    if (is.na(orchestration$selected_start))
      stop("No variational start satisfied the declared selection rule")
    backend <- selected_vi_model(orchestration)
    backends <- attr(orchestration, "models", exact = TRUE)
    selected_attempt <- orchestration$attempts[[orchestration$selected_start]]
    if (model$model == "CHLRM" &&
        !is.null(selected_attempt$aligned_results))
      backend$results <- selected_attempt$aligned_results
    backend_result <- .rebuild_vi_fit_result(
      backend, selected_attempt, model$model, expected_method)
    convergence <- .vi_convergence_result(backend$diagnostics, control)
  } else {
    backend <- .new_backend(
      model, method, validated, native_prior, native_control)
    .apply_initialization(
      backend, model$model, expected_method, working_initialization,
      lrm_ols, grouped_ols, model$component_count)
    backend_result <- backend$fit(verbose = verbose)
    convergence <- if (method == "mcmc")
      .single_chain_convergence_result(control) else
      .vi_convergence_result(backend$diagnostics, control)
  }
  if (!identical(prior_before, serialize(prior, NULL)))
    stop("Internal error: initialization modified the prior specification")
  .finalize_typed_fit(
    backend_result, model, prior, initialization, control, seed, backend,
    convergence, orchestration, backends, scaling, scaling_map,
    resource_estimate)
}

.fit_backend <- function(fit) {
  if (!inherits(fit, "vb_model_fit"))
    stop("fit must be produced by fit_model()")
  backend <- attr(fit, "backend", exact = TRUE)
  if (is.null(backend)) stop("fit does not contain its execution backend")
  backend
}
