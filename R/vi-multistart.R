# Reproducible multistart orchestration for CAVI/SVI. Initial states are first-
# class serializable objects so both algorithms can start from the same q.

.clone_vi_state <- function(x) unserialize(serialize(x, NULL))

.vi_data_dimensions <- function(model, y, X, K = NULL) {
  if (model == "LRM") {
    X <- as.matrix(X); y <- as.numeric(y)
    if (!nrow(X) || !ncol(X) || length(y) != nrow(X) ||
        any(!is.finite(X)) || any(!is.finite(y)))
      stop("LRM requires finite y and X with compatible nonzero dimensions")
    return(list(n = length(y), p = ncol(X)))
  }
  dims <- validate_grouped_data(X, y, if (model == "CHLRM") K else NULL)
  if (model == "CHLRM" && is.null(K)) stop("K is required for CHLRM")
  list(m = dims$m, p = dims$p, nj = dims$nj,
       K = if (model == "CHLRM") as.integer(K) else NULL)
}

.new_vi_fit <- function(model, y, X, K, hyperparams, run_settings) {
  switch(model,
    LRM = vi_lrm$new(y, X, hyperparams, run_settings),
    HLRM = vi_hlrm$new(y, X, hyperparams, run_settings),
    CHLRM = vi_chlrm$new(y, X, K, hyperparams, run_settings))
}

.vi_object_name <- function(model) switch(model,
  LRM = "lrm_vi_object", HLRM = "hlrm_vi_object", CHLRM = "chlrm_vi_object")

.vi_state_signature <- function(x) {
  if (is.list(x) && is.null(dim(x))) {
    return(list(names = names(x), children = lapply(x, .vi_state_signature)))
  }
  list(type = typeof(x), class = class(x), dim = dim(x), length = length(x))
}

.validate_vi_state <- function(state, template) {
  if (!is.list(state) || !identical(.vi_state_signature(state),
                                    .vi_state_signature(template)))
    stop("An initial state is incompatible with the requested model dimensions")
  numeric_values <- unlist(state, recursive = TRUE, use.names = FALSE)
  if (!is.numeric(numeric_values) || any(!is.finite(numeric_values)))
    stop("Initial states must contain only finite numeric parameters")
  invisible(state)
}

.set_vi_state <- function(fit, model, state) {
  field <- .vi_object_name(model)
  object <- fit[[field]]
  .validate_vi_state(state, object$params)
  object$params <- .clone_vi_state(state)
  fit[[field]] <- object
  invisible(fit)
}

# Generate states independently from algorithm randomness.  Passing the same
# object to separate CAVI and SVI calls guarantees identical starting factors.
make_vi_initial_states <- function(model = c("LRM", "HLRM", "CHLRM"), y, X,
                                   hyperparams, run_settings, K = NULL,
                                   n_starts = 5L, seed = NULL, seeds = NULL) {
  model <- match.arg(model)
  n_starts <- validate_integer_scalar(n_starts, "n_starts", 1L)
  if (!is.list(hyperparams)) stop("hyperparams must be a list")
  dimensions <- .vi_data_dimensions(model, y, X, K)
  initialization_seeds <- derive_independent_seeds(n_starts, seed, seeds)
  states <- vector("list", n_starts)
  for (start in seq_len(n_starts)) {
    settings <- run_settings
    settings$seed <- initialization_seeds[start]
    scope <- inference_rng_scope(initialization_seeds[start])
    fit <- tryCatch(.new_vi_fit(model, y, X, K, hyperparams, settings),
                    finally = scope$restore())
    states[[start]] <- .clone_vi_state(fit[[.vi_object_name(model)]]$params)
  }
  out <- structure(list(
    schema_version = "vi-initial-states-1",
    model = model,
    dimensions = dimensions,
    states = states,
    metadata = data.frame(start = seq_len(n_starts),
                          initialization_seed = initialization_seeds),
    initialization = if (model == "CHLRM") run_settings$initialization else NULL
  ), class = "vi_initial_states")
  out
}

.validate_vi_initial_states <- function(initial_states, model, dimensions,
                                        n_starts = NULL) {
  if (!inherits(initial_states, "vi_initial_states") ||
      !identical(initial_states$model, model) ||
      !identical(initial_states$dimensions, dimensions) ||
      !is.list(initial_states$states) || !length(initial_states$states))
    stop("initial_states are incompatible with model/data dimensions")
  if (!is.null(n_starts) && length(initial_states$states) != n_starts)
    stop("n_starts does not match the supplied initial_states")
  invisible(initial_states)
}

.run_one_vi_start <- function(model, y, X, K, hyperparams, run_settings,
                              state, seed, verbose) {
  scope <- inference_rng_scope(seed)
  on.exit(scope$restore(), add = TRUE)
  settings <- run_settings
  settings$seed <- NULL
  fit <- .new_vi_fit(model, y, X, K, hyperparams, settings)
  .set_vi_state(fit, model, state)
  result <- fit$fit(verbose = verbose, seed = seed)
  list(model = fit, result = result)
}

.vi_reference_from_attempts <- function(attempts, selected_start) {
  successful <- which(vapply(attempts, function(x) identical(x$status, "complete"),
                             logical(1L)))
  results <- lapply(attempts[successful], `[[`, "results")
  K <- nrow(results[[1L]]$MuBetak); p <- ncol(results[[1L]]$MuBetak)
  draws <- list(
    Betak = array(NA_real_, c(length(results), K, p)),
    sigma2k = matrix(NA_real_, length(results), K),
    omega = matrix(NA_real_, length(results), K),
    logjoint = vapply(attempts[successful], `[[`, numeric(1L), "final_elbo")
  )
  for (s in seq_along(results)) {
    draws$Betak[s, , ] <- results[[s]]$MuBetak
    draws$sigma2k[s, ] <- exp(log(results[[s]]$bk) - digamma(results[[s]]$ak))
    draws$omega[s, ] <- results[[s]]$alpha / sum(results[[s]]$alpha)
  }
  reference_draw <- match(selected_start, successful)
  if (is.na(reference_draw)) reference_draw <- which.max(draws$logjoint)
  fit_label_reference_chlrm(draws, reference_draw = reference_draw)
}

# Fit every start, keep every outcome, and select the largest finite terminal
# ELBO among starts satisfying the declared convergence requirement.
run_vi_multistart <- function(model = c("LRM", "HLRM", "CHLRM"),
                              method = c("CAVI", "SVI"), y, X, hyperparams,
                              run_settings, K = NULL, n_starts = 5L,
                              seed = NULL, fit_seeds = NULL,
                              initial_states = NULL, label_reference = NULL,
                              require_converged = TRUE, max_memory_mb = 1024,
                              verbose = FALSE) {
  model <- match.arg(model); method <- match.arg(method)
  if (!is.logical(require_converged) || length(require_converged) != 1L ||
      is.na(require_converged)) stop("require_converged must be TRUE or FALSE")
  n_starts <- validate_integer_scalar(n_starts, "n_starts", 1L)
  dimensions <- .vi_data_dimensions(model, y, X, K)
  settings <- run_settings
  settings$stochastic <- identical(method, "SVI")
  vi_stopping_control(settings, if (model == "LRM") dimensions$n else dimensions$m)

  if (is.null(initial_states)) {
    all_seeds <- derive_independent_seeds(2L * n_starts, seed, NULL)
    initial_states <- make_vi_initial_states(model, y, X, hyperparams, settings,
      K, n_starts, seeds = all_seeds[seq_len(n_starts)])
    if (is.null(fit_seeds)) fit_seeds <- all_seeds[n_starts + seq_len(n_starts)]
  } else {
    .validate_vi_initial_states(initial_states, model, dimensions, n_starts)
  }
  if (is.null(fit_seeds)) fit_seeds <- derive_independent_seeds(n_starts, seed, NULL)
  fit_seeds <- derive_independent_seeds(n_starts, seeds = fit_seeds)

  estimated_mb <- 3 * sum(vapply(initial_states$states, object.size, numeric(1L))) / 1024^2
  if (!is.numeric(max_memory_mb) || length(max_memory_mb) != 1L ||
      !is.finite(max_memory_mb) || max_memory_mb <= 0)
    stop("max_memory_mb must be a positive finite scalar")
  if (estimated_mb > max_memory_mb)
    stop(sprintf("Estimated multistart storage %.1f MB exceeds max_memory_mb %.1f",
                 estimated_mb, max_memory_mb))

  attempts <- vector("list", n_starts)
  selected_models <- vector("list", n_starts)
  for (start in seq_len(n_starts)) {
    attempt <- tryCatch({
      fitted <- .run_one_vi_start(model, y, X, K, hyperparams, settings,
        initial_states$states[[start]], fit_seeds[start], verbose)
      final_elbo <- tail(fitted$result$objective$elbo[
        is.finite(fitted$result$objective$elbo)], 1L)
      if (!length(final_elbo)) stop("The fit produced no finite ELBO")
      selected_models[[start]] <- fitted$model
      list(status = "complete", error = NA_character_, results = fitted$model$results,
           diagnostics = fitted$result$diagnostics, timing = fitted$result$timing,
           final_elbo = as.numeric(final_elbo))
    }, error = function(e) list(
      status = "failed", error = conditionMessage(e), results = NULL,
      diagnostics = list(converged = FALSE, stop_reason = "error"),
      timing = NULL, final_elbo = NA_real_))
    attempts[[start]] <- attempt
  }

  complete <- vapply(attempts, function(x) identical(x$status, "complete"), logical(1L))
  converged <- vapply(attempts, function(x) isTRUE(x$diagnostics$converged), logical(1L))
  final_elbo <- vapply(attempts, `[[`, numeric(1L), "final_elbo")
  eligible <- complete & is.finite(final_elbo) & (!require_converged | converged)
  selected_start <- if (any(eligible)) which(eligible)[which.max(final_elbo[eligible])] else NA_integer_

  shared_reference <- label_reference
  successful <- which(complete & is.finite(final_elbo))
  if (model == "CHLRM" && length(successful)) {
    if (is.null(shared_reference))
      shared_reference <- .vi_reference_from_attempts(attempts, selected_start)
    .validate_chlrm_label_reference(shared_reference, dimensions$K, dimensions$p)
    for (start in successful) {
      attempts[[start]]$aligned_results <- align_to_reference_chlrm(
        attempts[[start]]$results, shared_reference)
      selected_models[[start]]$aligned_results <- attempts[[start]]$aligned_results
      selected_models[[start]]$label_reference <- shared_reference
    }
  }

  summary <- data.frame(
    start = seq_len(n_starts),
    initialization_seed = initial_states$metadata$initialization_seed,
    fit_seed = fit_seeds,
    status = vapply(attempts, `[[`, character(1L), "status"),
    converged = converged,
    stop_reason = vapply(attempts, function(x) x$diagnostics$stop_reason,
                         character(1L)),
    iterations = vapply(attempts, function(x)
      if (is.null(x$diagnostics$iterations)) NA_integer_ else
        as.integer(x$diagnostics$iterations), integer(1L)),
    final_elbo = final_elbo,
    elapsed_seconds = vapply(attempts, function(x)
      if (is.null(x$timing)) NA_real_ else x$timing$total_seconds, numeric(1L)),
    selected = seq_len(n_starts) == selected_start,
    error = vapply(attempts, `[[`, character(1L), "error")
  )

  out <- structure(list(
    schema_version = "vi-multistart-1",
    model = model,
    method = method,
    status = if (is.na(selected_start)) "no_valid_start" else "selected",
    require_converged = require_converged,
    attempts = attempts,
    summary = summary,
    selected_start = selected_start,
    selected_model = if (is.na(selected_start)) NULL else selected_models[[selected_start]],
    initial_states = initial_states,
    label_reference = shared_reference,
    run_settings = settings,
    estimated_storage_mb = estimated_mb
  ), class = "vi_multistart_fit")
  attr(out, "models") <- selected_models
  out
}

selected_vi_model <- function(fit) {
  if (!inherits(fit, "vi_multistart_fit")) stop("fit must be a vi_multistart_fit")
  if (is.null(fit$selected_model))
    stop("No start satisfied the declared convergence and finite-ELBO rule")
  fit$selected_model
}

print.vi_multistart_fit <- function(x, ...) {
  cat(sprintf("%s %s multistart: %d attempts; status %s\n",
              x$model, x$method, nrow(x$summary), x$status))
  if (!is.na(x$selected_start))
    cat(sprintf("Selected start %d with terminal ELBO %.6f\n",
                x$selected_start, x$summary$final_elbo[x$selected_start]))
  invisible(x)
}
