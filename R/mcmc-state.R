# Serializable state and transactional continuation for the three Gibbs
# samplers.  A continuation stores every new draw (thin = 1) and never repeats
# warm-up.

.mcmc_model_name <- function(fit) {
  if (inherits(fit, "mcmc_lrm")) return("LRM")
  if (inherits(fit, "mcmc_hlrm")) return("HLRM")
  if (inherits(fit, "mcmc_chlrm")) return("CHLRM")
  stop("fit must be an LRM, HLRM or CHLRM MCMC object")
}

.mcmc_object_name <- function(model) switch(model,
  LRM = "lrm_mcmc_object", HLRM = "hlrm_mcmc_object",
  CHLRM = "chlrm_mcmc_object")

.clone_mcmc_state <- function(x) unserialize(serialize(x, NULL))

.step_mcmc_sampler <- function(object, model) {
  switch(model,
    LRM = {
      object <- sample_coefficients_lrm(object)
      sample_residual_variance_lrm(object)
    },
    HLRM = {
      object <- sample_group_coefficients_hlrm(object)
      object <- sample_group_variances_hlrm(object)
      object <- sample_xi2_hlrm(object)
      object <- sample_global_coefficients_hlrm(object)
      sample_global_covariance_hlrm(object)
    },
    CHLRM = {
      object <- sample_allocations_chlrm(object)
      object <- sample_omega_chlrm(object)
      object <- sample_component_coefficients_chlrm(object)
      object <- sample_component_variances_chlrm(object)
      object <- sample_global_coefficients_chlrm(object)
      object <- sample_global_covariance_chlrm(object)
      sample_xi2_chlrm(object)
    })
}

.allocate_mcmc_draws <- function(model, stats, samples) {
  p <- stats$p
  switch(model,
    LRM = list(Beta = matrix(NA_real_, samples, p), sigma2 = numeric(samples),
               loglik = numeric(samples), logjoint = numeric(samples)),
    HLRM = list(Betaj = array(NA_real_, c(samples, stats$m, p)),
      sigma2j = matrix(NA_real_, samples, stats$m),
      Beta = matrix(NA_real_, samples, p),
      Sigma = array(NA_real_, c(samples, p, p)), xi2 = numeric(samples),
      loglik = numeric(samples), logjoint = numeric(samples)),
    CHLRM = list(Betaj = array(NA_real_, c(samples, stats$m, p)),
      sigma2j = matrix(NA_real_, samples, stats$m),
      z = matrix(NA_integer_, stats$m, samples),
      omega = matrix(NA_real_, samples, stats$K),
      Betak = array(NA_real_, c(samples, stats$K, p)),
      sigma2k = matrix(NA_real_, samples, stats$K),
      Beta = matrix(NA_real_, samples, p),
      Sigma = array(NA_real_, c(samples, p, p)), xi2 = numeric(samples),
      loglik = numeric(samples), logjoint = numeric(samples)))
}

.store_mcmc_draw <- function(draws, index, object, model) {
  params <- object$params
  if (model == "LRM") {
    draws$Beta[index, ] <- params$Beta
    draws$sigma2[index] <- params$sigma2
    draws$loglik[index] <- loglik_total_lrm(object$stats$y, object$stats$X,
                                             params$Beta, params$sigma2)
    draws$logjoint[index] <- logprior_lrm(object) + draws$loglik[index]
  } else if (model == "HLRM") {
    draws$Betaj[index, , ] <- params$Betaj
    draws$sigma2j[index, ] <- params$sigma2j
    draws$Beta[index, ] <- params$Beta
    draws$Sigma[index, , ] <- params$Sigma
    draws$xi2[index] <- params$xi2
    draws$loglik[index] <- loglik_total_hlrm(object$stats$y, object$stats$X,
                                              params$Betaj, params$sigma2j)
    draws$logjoint[index] <- logprior_hlrm(object) + draws$loglik[index]
  } else {
    group_coefficients <- params$Betak[params$z, , drop = FALSE]
    group_variances <- params$sigma2k[params$z]
    draws$Betaj[index, , ] <- group_coefficients
    draws$sigma2j[index, ] <- group_variances
    draws$z[, index] <- params$z
    draws$omega[index, ] <- params$omega
    draws$Betak[index, , ] <- params$Betak
    draws$sigma2k[index, ] <- params$sigma2k
    draws$Beta[index, ] <- params$Beta
    draws$Sigma[index, , ] <- params$Sigma
    draws$xi2[index] <- params$xi2
    draws$loglik[index] <- loglik_total_chlrm(object$stats$y, object$stats$X,
                                               group_coefficients,
                                               group_variances)
    draws$logjoint[index] <- logprior_chlrm(object) + draws$loglik[index]
  }
  draws
}

.sample_mcmc_block <- function(object, model, samples, verbose = FALSE) {
  draws <- .allocate_mcmc_draws(model, object$stats, samples)
  elapsed <- numeric(samples)
  start <- proc.time()[["elapsed"]]
  if (verbose) {
    progress <- txtProgressBar(min = 0, max = samples, width = 50,
                               char = "=", style = 3)
    on.exit(close(progress), add = TRUE)
  }
  for (iteration in seq_len(samples)) {
    object <- .step_mcmc_sampler(object, model)
    draws <- .store_mcmc_draw(draws, iteration, object, model)
    elapsed[iteration] <- proc.time()[["elapsed"]] - start
    if (verbose) setTxtProgressBar(progress, iteration)
  }
  list(object = object, draws = draws,
       timing = list(up_to_iter = elapsed, total = tail(elapsed, 1L)))
}

.bind_first_dimension <- function(objects) {
  first <- objects[[1L]]
  if (is.null(dim(first))) return(unlist(objects, use.names = FALSE))
  if (length(dim(first)) == 2L) return(do.call(rbind, objects))
  if (length(dim(first)) != 3L) stop("Unsupported chain field dimension")
  tail_dim <- dim(first)[-1L]
  if (!all(vapply(objects, function(x) identical(dim(x)[-1L], tail_dim), logical(1))))
    stop("Cannot combine incompatible chain fields")
  out <- array(NA_real_, c(sum(vapply(objects, function(x) dim(x)[1L], integer(1L))),
                           tail_dim))
  start <- 1L
  for (x in objects) {
    finish <- start + dim(x)[1L] - 1L
    out[start:finish, , ] <- x
    start <- finish + 1L
  }
  out
}

.combine_chain_draws <- function(chains, model) {
  fields <- names(chains[[1L]])
  out <- setNames(vector("list", length(fields)), fields)
  for (field in fields) {
    values <- lapply(chains, `[[`, field)
    out[[field]] <- if (model == "CHLRM" && field == "z")
      do.call(cbind, values) else .bind_first_dimension(values)
  }
  out
}

# Continue from the exact terminal sampler/RNG state.  The operation is
# transactional: existing chains and sampler parameters are replaced only
# after a complete successful block.
.resume_mcmc_backend <- function(fit, additional_samples, seed = NULL,
                                 verbose = FALSE) {
  model <- .mcmc_model_name(fit)
  additional_samples <- validate_integer_scalar(additional_samples,
                                                 "additional_samples", 1L)
  if (!identical(fit$fit_status, "fitted"))
    stop("Only a successfully fitted MCMC object can be resumed")
  thin <- if (is.null(fit$run_settings$thin)) 1L else fit$run_settings$thin
  if (!identical(as.integer(thin), 1L))
    stop("MCMC continuation requires thin = 1")
  if (is.null(seed) && is.null(fit$rng_state))
    stop("The fit has no saved RNG state; provide seed to start an explicit branch")
  validate_rng_seed(seed)
  field <- .mcmc_object_name(model)
  sampler <- .clone_mcmc_state(fit[[field]])
  scope <- inference_rng_scope(seed, if (is.null(seed)) fit$rng_state else NULL)
  on.exit(scope$restore(), add = TRUE)
  fit$fit_status <- "running"
  on.exit({ if (fit$fit_status == "running") fit$fit_status <- "failed" }, add = TRUE)

  block <- .sample_mcmc_block(sampler, model, additional_samples, verbose)
  new_rng_state <- rng_snapshot()
  combined <- .combine_chain_draws(list(fit$chains, block$draws), model)

  fit[[field]] <- block$object
  fit$chains <- combined
  fit$rng_state <- new_rng_state
  fit$run_settings$samples <- fit$run_settings$samples + additional_samples
  metadata <- fit$fit_metadata
  if (is.null(metadata)) metadata <- list(total_seconds = 0,
    warmup_iterations = fit$run_settings$burn_in,
    stored_samples = fit$run_settings$samples - additional_samples,
    extensions = list())
  extension <- list(samples = additional_samples, seconds = block$timing$total,
                    seed_override = seed)
  metadata$total_seconds <- metadata$total_seconds + block$timing$total
  metadata$stored_samples <- fit$run_settings$samples
  metadata$extensions[[length(metadata$extensions) + 1L]] <- extension
  fit$fit_metadata <- metadata
  fit$posterior_summary_result <- NULL
  if (model == "CHLRM") {
    fit$relabeled_chains <- list()
    fit$label_reference <- NULL
  }
  fit$fit_status <- "fitted"
  invisible(list(chains = block$draws, timing = block$timing,
                 total_stored_samples = fit$run_settings$samples))
}


.mcmc_extension_resource_estimate <- function(fit, control) {
  model <- fit$specifications$model
  backend <- .fit_backend(fit)
  data <- .mcmc_backend_data(backend, model$model)
  validated <- .validate_fit_data(model, data)
  scaling <- fit$specifications$scaling
  .estimate_fit_resources_validated(
    model, validated, control, data, scaling)
}

.new_extended_multichain <- function(backends, old_orchestration, model,
                                     control, extension_seeds) {
  chains <- lapply(backends, `[[`, "chains")
  elapsed <- vapply(backends, function(backend)
    backend$fit_metadata$total_seconds, numeric(1L))
  stopping <- control$stopping
  thresholds <- list(
    rhat = stopping$rhat_max,
    ess_bulk = stopping$ess_bulk_min,
    ess_tail = stopping$ess_tail_min,
    relative_mcse_mean = stopping$relative_mcse_mean_max,
    relative_mcse_sd = stopping$relative_mcse_sd_max)
  diagnosed <- .diagnose_mcmc_chain_set(
    chains, model, elapsed, control$diagnostic_scope, thresholds)
  samples <- .mcmc_draw_count(chains[[1L]], model)
  passed <- diagnosed$diagnostics$passed &&
    samples >= stopping$min_post_warmup_draws
  finite_max <- function(x) if (any(is.finite(x)))
    max(x[is.finite(x)]) else NA_real_
  finite_min <- function(x) if (any(is.finite(x)))
    min(x[is.finite(x)]) else NA_real_
  history <- old_orchestration$block_history
  history <- rbind(history, data.frame(
    block = nrow(history) + 1L,
    stored_samples = samples,
    status = if (samples < stopping$min_post_warmup_draws)
      "minimum_not_reached" else diagnosed$diagnostics$status,
    max_rhat = finite_max(diagnosed$diagnostics$variables$rhat),
    min_ess_bulk = finite_min(diagnosed$diagnostics$variables$ess_bulk),
    min_ess_tail = finite_min(diagnosed$diagnostics$variables$ess_tail),
    elapsed_seconds = sum(elapsed)))
  metadata <- old_orchestration$chain_metadata
  metadata$elapsed_seconds <- elapsed
  metadata$last_extension_seed <- if (all(vapply(
      extension_seeds, is.null, logical(1L)))) NA_integer_ else
    as.integer(unlist(extension_seeds))
  diagnostics <- diagnosed$diagnostics
  diagnostics$passed <- passed
  diagnostics$status <- if (passed) "passed" else
    if (samples < stopping$min_post_warmup_draws)
      "minimum_not_reached" else diagnostics$status
  dimensions <- .validate_multichain_data(
    model,
    backends[[1L]][[.mcmc_object_name(model)]]$stats$y,
    backends[[1L]][[.mcmc_object_name(model)]]$stats$X,
    if (model == "CHLRM")
      backends[[1L]][[.mcmc_object_name(model)]]$stats$K else NULL)
  component_count <- if (model == "CHLRM")
    backends[[1L]][[.mcmc_object_name(model)]]$stats$K else NULL
  out <- structure(list(
    schema_version = "mcmc-multichain-1",
    model = model,
    status = if (passed) "diagnosed" else "not_converged",
    chains = chains,
    aligned_chains = if (model == "CHLRM") diagnosed$chains else NULL,
    label_reference = diagnosed$reference,
    diagnostics = diagnostics,
    block_history = history,
    chain_metadata = metadata,
    run_settings = utils::modifyList(
      old_orchestration$run_settings, list(samples = samples, thin = 1L)),
    execution = list(
      incremental = TRUE, block_size = samples -
        old_orchestration$run_settings$samples,
      max_samples = samples, stop_reason = "explicit_extension"),
    estimated_storage_mb = .estimate_multichain_mb(
      model, dimensions, samples, length(backends),
      component_count,
      control$diagnostic_scope)),
    class = "mcmc_multichain_fit")
  attr(out, "models") <- backends
  out
}

# Continue a typed MCMC fit without mutating the input object. For multichain
# fits each chain resumes from its own sampler and RNG state, after which labels
# and convergence diagnostics are recomputed jointly.
#' @rdname orchestration
#' @export
resume_mcmc <- function(fit, additional_samples, seed = NULL,
                        verbose = FALSE) {
  .require_complete_mcmc_fit(fit)
  if (!inherits(fit, "vb_model_fit") || !identical(fit$method, "MCMC"))
    stop("fit must be an MCMC vb_model_fit produced by fit_model()")
  additional_samples <- validate_integer_scalar(
    additional_samples, "additional_samples", 1L)
  validate_rng_seed(seed)
  control <- unserialize(serialize(fit$specifications$control, NULL))
  old_samples <- control$post_warmup_draws
  control$post_warmup_draws <- old_samples + additional_samples
  estimate <- .mcmc_extension_resource_estimate(fit, control)
  if (!isTRUE(estimate$within_budget))
    stop(sprintf(
      "Estimated peak memory %.1f MB exceeds max_memory_mb %.1f",
      estimate$estimated_peak_memory_mb, estimate$max_memory_mb))

  stored_backends <- attr(fit, "backends", exact = TRUE)
  multichain <- !is.null(stored_backends)
  backends <- if (multichain) stored_backends else list(.fit_backend(fit))
  backends <- lapply(backends, function(backend) backend$clone(deep = TRUE))
  extension_seeds <- if (is.null(seed)) rep(list(NULL), length(backends)) else
    as.list(derive_independent_seeds(length(backends), seed))
  for (chain in seq_along(backends))
    .resume_mcmc_backend(
      backends[[chain]], additional_samples, extension_seeds[[chain]], verbose)

  model <- fit$model
  if (multichain) {
    old_orchestration <- attr(fit, "orchestration", exact = TRUE)
    if (is.null(old_orchestration))
      stop("Multichain fit does not retain orchestration state")
    orchestration <- .new_extended_multichain(
      backends, old_orchestration, model, control, extension_seeds)
    combined <- combine_mcmc_chains(orchestration, require_converged = FALSE)
    backend <- backends[[1L]]$clone(deep = TRUE)
    backend$chains <- combined
    backend$run_settings$samples <- .mcmc_draw_count(combined, model)
    backend$fit_metadata <- list(
      total_seconds = sum(orchestration$chain_metadata$elapsed_seconds),
      warmup_iterations = control$warmup,
      stored_samples = backend$run_settings$samples,
      chains = control$chains)
    if (model == "CHLRM") {
      backend$relabeled_chains <- combined
      backend$label_reference <- orchestration$label_reference
    }
    backend_result <- new_fit_result(
      model, "MCMC", "posterior draws", combined,
      diagnostics = orchestration$diagnostics,
      total_seconds = backend$fit_metadata$total_seconds,
      cumulative_seconds = orchestration$block_history$elapsed_seconds,
      metadata = backend$fit_metadata)
    convergence <- .multichain_convergence_result(orchestration, control)
  } else {
    orchestration <- NULL
    backend <- backends[[1L]]
    backend_result <- new_fit_result(
      model, "MCMC", "posterior draws", backend$chains,
      diagnostics = list(convergence = "not assessed for a single chain"),
      total_seconds = backend$fit_metadata$total_seconds,
      cumulative_seconds = numeric(), metadata = backend$fit_metadata)
    convergence <- .single_chain_convergence_result(control)
  }
  result <- .finalize_typed_fit(
    backend_result, fit$specifications$model, fit$specifications$prior,
    fit$specifications$initialization, control, fit$metadata$seed, backend,
    convergence, orchestration,
    if (multichain) backends else NULL,
    fit$specifications$scaling,
    attr(fit, "scaling_map", exact = TRUE), estimate)
  reference <- attr(fit, "label_reference", exact = TRUE)
  if (is.null(reference)) result else relabel_fit(result, reference)
}
n_samples_mcmc <- function(run_settings) {
  if (!is.null(run_settings$thin)) {
    run_settings$samples * run_settings$thin + run_settings$burn_in
  } else {
    run_settings$samples + run_settings$burn_in
  }
}

thin_chain_mcmc <- function(iteration, run_settings) {
  warmup <- if (is.null(run_settings$burn_in)) 0L else as.integer(run_settings$burn_in)
  thinning <- if (is.null(run_settings$thin)) 1L else as.integer(run_settings$thin)
  if (thinning < 1L) stop("run_settings$thin must be a positive integer")

  post_warmup <- as.integer(iteration) - warmup
  if (post_warmup < 1L || post_warmup %% thinning != 0L) return(NULL)
  post_warmup %/% thinning
}

validate_mcmc_settings <- function(run_settings) {
  required <- c("samples", "burn_in")
  missing <- setdiff(required, names(run_settings))
  if (length(missing)) {
    stop("Missing MCMC settings: ", paste(missing, collapse = ", "))
  }
  values <- unlist(
    run_settings[c("samples", "burn_in", "thin")],
    use.names = FALSE
  )
  values <- values[!is.na(values)]
  invalid <- any(!is.finite(values)) ||
    run_settings$samples < 1 ||
    run_settings$burn_in < 0 ||
    run_settings$samples != as.integer(run_settings$samples) ||
    run_settings$burn_in != as.integer(run_settings$burn_in) ||
    (!is.null(run_settings$thin) &&
      (run_settings$thin < 1 || run_settings$thin != as.integer(run_settings$thin)))
  if (invalid) {
    stop("samples/thin must be positive integers and burn_in a non-negative integer")
  }
  invisible(TRUE)
}
