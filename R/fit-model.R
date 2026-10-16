# Public typed construction facade. The existing R6 classes are execution
# backends; this is the only translation layer from semantic public names to
# their current equation-oriented notation.

.validate_fit_data <- function(model, data) {
  if (!inherits(model, "vb_model_spec"))
    stop("model must be created by lrm_model(), hlrm_model() or chlrm_model()")
  if (!is.list(data) || !all(c("X", "y") %in% names(data)))
    stop("data must be a list containing X and y")
  if (model$model == "LRM") {
    X <- as.matrix(data$X)
    y <- as.numeric(data$y)
    if (!nrow(X) || !ncol(X) || nrow(X) != length(y) ||
        any(!is.finite(X)) || any(!is.finite(y)))
      stop("LRM requires finite y and X with compatible nonzero dimensions")
    dimensions <- list(n = length(y), p = ncol(X))
  } else {
    X <- lapply(data$X, as.matrix)
    y <- lapply(data$y, as.numeric)
    dimensions <- validate_grouped_data(
      X, y, if (model$model == "CHLRM") model$component_count else NULL)
  }
  if (!is.null(model$coefficient_count) &&
      model$coefficient_count != dimensions$p)
    stop("model coefficient_count does not match the data")
  list(X = X, y = y, dimensions = dimensions)
}

.validate_prior_for_model <- function(prior, model, dimensions) {
  if (!inherits(prior, "vb_prior_spec"))
    stop("prior must be created by a model-specific prior constructor")
  if (!identical(prior$model, model$model))
    stop("prior and model must describe the same model family")
  if (length(prior$coefficient_mean) != dimensions$p)
    stop("prior coefficient dimension does not match the data")
  if (model$model == "CHLRM" &&
      length(prior$mixture_concentration) != model$component_count)
    stop("mixture_concentration must have one entry per component")
  invisible(prior)
}

.group_ols_initialization <- function(data, jitter) {
  m <- length(data$X)
  p <- ncol(data$X[[1L]])
  coefficients <- matrix(NA_real_, m, p)
  variances <- numeric(m)
  for (j in seq_len(m)) {
    precision <- crossprod(data$X[[j]]) + jitter * diag(p)
    coefficients[j, ] <- solve_spd(
      precision, crossprod(data$X[[j]], data$y[[j]]),
      paste0("regularized group XtX ", j))
    residual <- data$y[[j]] - drop(data$X[[j]] %*% coefficients[j, ])
    variances[j] <- max(mean(residual^2), jitter)
  }
  list(coefficients = coefficients, variances = variances)
}

.lrm_ols_initialization <- function(data, jitter) {
  p <- ncol(data$X)
  coefficient <- solve_spd(
    crossprod(data$X) + jitter * diag(p),
    crossprod(data$X, data$y), "regularized XtX")
  residual <- data$y - drop(data$X %*% coefficient)
  list(coefficients = as.numeric(coefficient),
       residual_variance = max(mean(residual^2), jitter))
}

.prior_to_native <- function(prior, model, dimensions,
                             grouped_ols = NULL) {
  if (model$model == "LRM") {
    return(validate_lrm_hyperparams(list(
      Beta0 = prior$coefficient_mean,
      Sigma0 = prior$coefficient_covariance,
      nu0 = prior$residual_degrees_of_freedom,
      sigma20 = prior$residual_variance_scale
    ), dimensions$p))
  }
  native <- list(
    Mu0 = prior$coefficient_mean,
    Lambda0 = prior$coefficient_covariance,
    n0 = prior$coefficient_dispersion_df,
    S0 = prior$coefficient_dispersion_scale,
    nu0 = prior$residual_degrees_of_freedom,
    a0 = prior$residual_variance_scale_shape,
    b0 = prior$residual_variance_scale_rate
  )
  if (!is.null(grouped_ols)) native$Betaj_ols <- grouped_ols$coefficients
  if (model$model == "CHLRM")
    native$alpha0 <- prior$mixture_concentration
  validate_hlrm_hyperparams(
    native, dimensions$p, dimensions$m, model$model,
    if (model$model == "CHLRM") model$component_count else NULL)
}

.control_to_native <- function(control, seed) {
  validate_rng_seed(seed)
  if (control$method == "MCMC") {
    return(list(samples = control$post_warmup_draws,
                burn_in = control$warmup,
                thin = control$thinning,
                samples_gof = control$evaluation_draws,
                seed = seed))
  }
  stopping <- control$stopping
  out <- list(
    iterations = control$max_iterations,
    epsilon = stopping$objective_tolerance,
    parameter_tolerance = stopping$parameter_tolerance,
    patience = stopping$patience,
    min_iter = stopping$min_iterations,
    elbo_every = stopping$check_every,
    elbo_window = stopping$smoothing_window,
    samples_gof = control$evaluation_draws,
    stochastic = identical(control$method, "SVI"),
    seed = seed)
  if (identical(control$method, "SVI")) {
    out$minibatch_m <- control$minibatch_size
    out$tau <- control$learning_rate$offset
    out$kappa <- control$learning_rate$exponent
  }
  out
}

.new_backend <- function(model, method, data, native_prior, native_control) {
  switch(model$model,
    LRM = if (method == "mcmc")
        mcmc_lrm$new(data$y, data$X, native_prior, native_control) else
        vi_lrm$new(data$y, data$X, native_prior, native_control),
    HLRM = if (method == "mcmc")
      mcmc_hlrm$new(data$y, data$X, native_prior, native_control) else
      vi_hlrm$new(data$y, data$X, native_prior, native_control),
    CHLRM = if (method == "mcmc")
      mcmc_chlrm$new(data$y, data$X, model$component_count,
                     native_prior, native_control) else
      vi_chlrm$new(data$y, data$X, model$component_count,
                   native_prior, native_control))
}

.backend_state_field <- function(model, method) {
  paste0(tolower(model), if (method == "MCMC") "_mcmc_object" else
    "_vi_object")
}

.initialization_name_map <- function(model, method) {
  variational <- method != "MCMC"
  switch(model,
    LRM = if (!variational) c(
      coefficients = "Beta",
      residual_variance = "sigma2") else c(
      coefficient_mean = "MuBeta",
      coefficient_covariance = "SigmaBeta",
      residual_variance_shape = "asigma2",
      residual_variance_rate = "bsigma2"),
    HLRM = if (!variational) c(
      group_coefficients = "Betaj",
      variance_scale = "xi2",
      group_residual_variances = "sigma2j",
      population_coefficients = "Beta",
      population_covariance = "Sigma") else c(
      group_coefficient_means = "MuBetaj",
      group_coefficient_covariances = "SigmaBetaj",
      population_coefficient_mean = "MuBeta",
      population_coefficient_covariance = "SigmaBeta",
      population_covariance_df = "nuS",
      population_covariance_scale = "SS",
      group_residual_variance_shapes = "aj",
      group_residual_variance_rates = "bj",
      variance_scale_shape = "axi2",
      variance_scale_rate = "bxi2"),
    CHLRM = if (!variational) c(
      mixture_weights = "omega",
      group_allocations = "z",
      population_coefficients = "Beta",
      population_covariance = "Sigma",
      variance_scale = "xi2",
      component_coefficients = "Betak",
      component_residual_variances = "sigma2k") else c(
      allocation_probabilities = "rho",
      mixture_concentration = "alpha",
      component_coefficient_means = "MuBetak",
      component_coefficient_covariances = "SigmaBetak",
      population_coefficient_mean = "MuBeta",
      population_coefficient_covariance = "SigmaBeta",
      population_covariance_df = "nuS",
      population_covariance_scale = "SS",
      component_residual_variance_shapes = "ak",
      component_residual_variance_rates = "bk",
      variance_scale_shape = "axi2",
      variance_scale_rate = "bxi2"),
    stop("Unknown model family"))
}

.validate_initialization_state <- function(values, model, method) {
  positive <- function(x, name) {
    if (any(x <= 0)) stop(name, " must contain positive values")
  }
  spd_list <- function(x, name) {
    for (index in seq_along(x))
      chol_spd(x[[index]], paste(name, index))
  }
  if (method == "MCMC") {
    if (model == "LRM") {
      positive(values$sigma2, "residual_variance")
    } else {
      chol_spd(values$Sigma, "population_covariance")
      positive(values$xi2, "variance_scale")
      if (model == "HLRM") {
        positive(values$sigma2j, "group_residual_variances")
      } else {
        positive(values$sigma2k, "component_residual_variances")
        positive(values$omega, "mixture_weights")
        if (abs(sum(values$omega) - 1) > 1e-10)
          stop("mixture_weights must sum to one")
        if (any(values$z != floor(values$z)) ||
            any(values$z < 1 | values$z > length(values$omega)))
          stop("group_allocations must contain valid component integers")
      }
    }
    return(invisible(values))
  }
  chol_spd(values$SigmaBeta, "population_coefficient_covariance")
  if (model == "LRM") {
    positive(values$asigma2, "residual_variance_shape")
    positive(values$bsigma2, "residual_variance_rate")
    return(invisible(values))
  }
  chol_spd(values$SS, "population_covariance_scale")
  if (values$nuS <= nrow(values$SS) - 1L)
    stop("population_covariance_df must exceed coefficient_count - 1")
  positive(values$axi2, "variance_scale_shape")
  positive(values$bxi2, "variance_scale_rate")
  if (model == "HLRM") {
    spd_list(values$SigmaBetaj, "group_coefficient_covariance")
    positive(values$aj, "group_residual_variance_shapes")
    positive(values$bj, "group_residual_variance_rates")
  } else {
    spd_list(values$SigmaBetak, "component_coefficient_covariance")
    positive(values$ak, "component_residual_variance_shapes")
    positive(values$bk, "component_residual_variance_rates")
    positive(values$alpha, "mixture_concentration")
    if (any(values$rho < 0) ||
        any(abs(rowSums(values$rho) - 1) > 1e-10))
      stop("allocation_probabilities must be nonnegative with rows summing to one")
  }
  invisible(values)
}

.apply_values_initialization <- function(backend, model, method, values) {
  field <- .backend_state_field(model, method)
  state <- backend[[field]]
  template <- state$params
  name_map <- .initialization_name_map(model, method)
  if (!setequal(names(values), names(name_map)))
    stop("initialization values must contain exactly these public fields: ",
         paste(names(name_map), collapse = ", "))
  values <- values[names(name_map)]
  names(values) <- unname(name_map)
  values <- values[names(template)]
  if (method != "MCMC") {
    .validate_vi_state(values, template)
  } else {
    candidate <- unlist(values, recursive = TRUE, use.names = FALSE)
    if (!is.numeric(candidate) || any(!is.finite(candidate)))
      stop("MCMC initialization values must contain only finite numeric parameters")
    if (!identical(.vi_state_signature(values), .vi_state_signature(template)))
      stop("MCMC initialization values have incompatible dimensions")
  }
  .validate_initialization_state(values, model, method)
  state$params <- unserialize(serialize(values, NULL))
  backend[[field]] <- state
  invisible(backend)
}

.apply_ols_initialization <- function(backend, model, method, ols) {
  field <- .backend_state_field(model, method)
  state <- backend[[field]]
  if (model == "LRM") {
    if (method == "MCMC") {
      state$params$Beta <- ols$coefficients
      state$params$sigma2 <- ols$residual_variance
    } else {
      state$params$MuBeta <- ols$coefficients
      shape <- max(state$params$asigma2, 2)
      state$params$asigma2 <- shape
      state$params$bsigma2 <- (shape - 1) * ols$residual_variance
    }
  } else if (model == "HLRM") {
    if (method == "MCMC") {
      state$params$Betaj <- ols$coefficients
      state$params$sigma2j <- ols$variances
      state$params$Beta <- colMeans(ols$coefficients)
      state$params$xi2 <- stats::median(ols$variances)
    } else {
      state$params$MuBetaj <- ols$coefficients
      state$params$MuBeta <- colMeans(ols$coefficients)
    }
  } else stop("initialization_ols() is not complete for CHLRM; use initialization_kmeans()")
  backend[[field]] <- state
  invisible(backend)
}

.apply_kmeans_initialization <- function(backend, model, method, ols,
                                         component_count) {
  if (model != "CHLRM")
    stop("initialization_kmeans() applies only to CHLRM")
  # The variational constructor consumes Betaj_ols directly. MCMC needs the equivalent
  # state written explicitly because its constructor otherwise draws from the
  # prior.
  if (method != "MCMC") return(invisible(backend))
  clustering <- tryCatch(
    stats::kmeans(ols$coefficients, centers = component_count),
    error = function(error) stop(
      "k-means initialization failed: ", conditionMessage(error)))
  state <- backend$chlrm_mcmc_object
  z <- as.integer(clustering$cluster)
  centers <- unname(clustering$centers)
  component_variances <- vapply(seq_len(component_count), function(k) {
    stats::median(ols$variances[z == k])
  }, numeric(1L))
  state$params$z <- z
  state$params$omega <- tabulate(z, nbins = component_count) / length(z)
  state$params$Betak <- centers
  state$params$sigma2k <- component_variances
  state$params$Beta <- colMeans(centers)
  if (component_count > 1L) {
    dispersion <- stats::cov(centers)
    if (is.null(dim(dispersion))) dispersion <- matrix(dispersion, 1L, 1L)
    state$params$Sigma <- symmetrize_numeric(
      dispersion + 1e-8 * diag(ncol(centers)),
      "k-means component dispersion")
  }
  state$params$xi2 <- stats::median(component_variances)
  backend$chlrm_mcmc_object <- state
  invisible(backend)
}

.canonical_public_method <- function(method) {
  switch(method, mcmc = "MCMC", cavi = "CAVI", svi = "SVI")
}

.apply_initialization <- function(backend, model, method, initialization,
                                  lrm_ols = NULL, grouped_ols = NULL,
                                  component_count = NULL) {
  if (initialization$strategy == "values")
    .apply_values_initialization(
      backend, model, method, initialization$values)
  if (initialization$strategy == "ols")
    .apply_ols_initialization(
      backend, model, method, if (model == "LRM") lrm_ols else grouped_ols)
  if (initialization$strategy == "kmeans")
    .apply_kmeans_initialization(
      backend, model, method, grouped_ols, component_count)
  invisible(backend)
}

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
