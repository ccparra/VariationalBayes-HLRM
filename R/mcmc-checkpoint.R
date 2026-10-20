# Typed serialization and restoration of MCMC checkpoints.

.mcmc_backend_checkpoint <- function(fit) {
  model <- .mcmc_model_name(fit)
  if (!identical(fit$fit_status, "fitted"))
    stop("Only a successfully fitted MCMC object can be checkpointed")
  structure(list(
    schema_version = "mcmc-checkpoint-1",
    model = model,
    sampler = .clone_mcmc_state(fit[[.mcmc_object_name(model)]]),
    chains = .clone_mcmc_state(fit$chains),
    run_settings = .clone_mcmc_state(fit$run_settings),
    rng_state = .clone_mcmc_state(fit$rng_state),
    fit_metadata = .clone_mcmc_state(fit$fit_metadata)
  ), class = "mcmc_checkpoint")
}

.validate_mcmc_backend_checkpoint <- function(checkpoint) {
  required <- c("schema_version", "model", "sampler", "chains",
                "run_settings", "rng_state", "fit_metadata")
  if (!inherits(checkpoint, "mcmc_checkpoint") ||
      !identical(checkpoint$schema_version, "mcmc-checkpoint-1") ||
      any(!required %in% names(checkpoint)) ||
      !checkpoint$model %in% c("LRM", "HLRM", "CHLRM") ||
      !is.list(checkpoint$sampler) || !is.list(checkpoint$chains))
    stop("Invalid or unsupported MCMC checkpoint")
  invisible(checkpoint)
}

.restore_mcmc_backend_checkpoint <- function(checkpoint) {
  .validate_mcmc_backend_checkpoint(checkpoint)
  model <- checkpoint$model
  sampler <- checkpoint$sampler
  stats <- sampler$stats
  settings <- checkpoint$run_settings
  settings$seed <- NULL
  scope <- inference_rng_scope(1L)
  on.exit(scope$restore(), add = TRUE)
  fit <- switch(model,
    LRM = mcmc_lrm$new(stats$y, stats$X, sampler$hyperparams, settings),
    HLRM = mcmc_hlrm$new(stats$y, stats$X, sampler$hyperparams, settings),
    CHLRM = mcmc_chlrm$new(stats$y, stats$X, stats$K,
                           sampler$hyperparams, settings))
  fit[[.mcmc_object_name(model)]] <- .clone_mcmc_state(sampler)
  # Construction must not reseed, but the saved provenance settings must not
  # lose their original seed after the temporary constructor settings are used.
  fit$run_settings <- .clone_mcmc_state(checkpoint$run_settings)
  fit$chains <- .clone_mcmc_state(checkpoint$chains)
  fit$rng_state <- .clone_mcmc_state(checkpoint$rng_state)
  fit$fit_metadata <- .clone_mcmc_state(checkpoint$fit_metadata)
  fit$fit_status <- "fitted"
  fit
}

.mcmc_backend_data <- function(backend, model) {
  state <- backend[[.mcmc_object_name(model)]]
  list(X = state$stats$X, y = state$stats$y)
}

.mcmc_draw_count <- function(chains, model) {
  if (model == "CHLRM") dim(chains$Betak)[1L] else nrow(chains$Beta)
}

.typed_mcmc_checkpoint <- function(fit) {
  if (!inherits(fit, "vb_model_fit") || !identical(fit$method, "MCMC"))
    stop("fit must be an MCMC vb_model_fit produced by fit_model()")
  summary_backend <- .fit_backend(fit)
  chain_backends <- attr(fit, "backends", exact = TRUE)
  if (is.null(chain_backends)) chain_backends <- list(summary_backend)
  if (!all(vapply(chain_backends, inherits, logical(1L), "R6")))
    stop("fit does not retain valid MCMC chain backends")

  snapshot <- fit
  attr(snapshot, "backend") <- NULL
  attr(snapshot, "backends") <- NULL
  attr(snapshot, "orchestration") <- NULL
  orchestration <- attr(fit, "orchestration", exact = TRUE)
  if (!is.null(orchestration)) attr(orchestration, "models") <- NULL
  structure(list(
    schema_version = "typed-mcmc-checkpoint-1",
    model = fit$model,
    fit = snapshot,
    summary_backend = .mcmc_backend_checkpoint(summary_backend),
    chain_backends = lapply(chain_backends, .mcmc_backend_checkpoint),
    orchestration = orchestration),
    class = c("vb_mcmc_checkpoint", "list"))
}

.validate_typed_mcmc_checkpoint <- function(checkpoint) {
  required <- c("schema_version", "model", "fit", "summary_backend",
                "chain_backends", "orchestration")
  if (!inherits(checkpoint, "vb_mcmc_checkpoint") ||
      !identical(checkpoint$schema_version, "typed-mcmc-checkpoint-1") ||
      any(!required %in% names(checkpoint)) ||
      !inherits(checkpoint$fit, "vb_model_fit") ||
      !identical(checkpoint$fit$method, "MCMC") ||
      !checkpoint$model %in% c("LRM", "HLRM", "CHLRM") ||
      !length(checkpoint$chain_backends))
    stop("Invalid or unsupported typed MCMC checkpoint")
  .validate_mcmc_backend_checkpoint(checkpoint$summary_backend)
  lapply(checkpoint$chain_backends, .validate_mcmc_backend_checkpoint)
  models <- vapply(checkpoint$chain_backends, `[[`, character(1L), "model")
  if (!all(models == checkpoint$model) ||
      !identical(checkpoint$summary_backend$model, checkpoint$model))
    stop("Checkpoint contains inconsistent model families")
  invisible(checkpoint)
}

# Save only plain serializable state; mutable R6 environments are reconstructed
# on read and are never serialized as part of the public artifact.
#' @rdname orchestration
#' @export
save_mcmc_checkpoint <- function(fit, path, store = NULL, chunk_size = 5000L) {
  if (!is.character(path) || length(path) != 1L || !nzchar(path))
    stop("path must be one nonempty file path")
  if (file.exists(path))
    stop("Checkpoint path already exists; use a new versioned path")
  directory <- dirname(path)
  if (!dir.exists(directory)) dir.create(directory, recursive = TRUE)
  if (!is.null(store))
    return(.save_mcmc_checkpoint_store(fit, path, store, chunk_size))
  temporary <- tempfile(pattern = ".mcmc-checkpoint-", tmpdir = directory,
                        fileext = ".rds")
  on.exit(if (file.exists(temporary)) unlink(temporary), add = TRUE)
  saveRDS(.typed_mcmc_checkpoint(fit), temporary, version = 3)
  if (!file.rename(temporary, path))
    stop("Could not atomically publish checkpoint")
  invisible(normalizePath(path, mustWork = TRUE))
}

#' @rdname orchestration
#' @export
read_mcmc_checkpoint <- function(path) {
  if (!is.character(path) || length(path) != 1L || !nzchar(path))
    stop("path must be one nonempty file path")
  if (!file.exists(path)) stop("Checkpoint file does not exist")
  checkpoint <- readRDS(path)
  if (inherits(checkpoint, "vb_mcmc_store_checkpoint"))
    checkpoint <- .read_mcmc_checkpoint_store(checkpoint, path)
  .validate_typed_mcmc_checkpoint(checkpoint)
  fit <- checkpoint$fit
  summary_backend <- .restore_mcmc_backend_checkpoint(
    checkpoint$summary_backend)
  chain_backends <- lapply(
    checkpoint$chain_backends, .restore_mcmc_backend_checkpoint)
  attr(fit, "backend") <- summary_backend
  if (length(chain_backends) > 1L) {
    attr(fit, "backends") <- chain_backends
    orchestration <- checkpoint$orchestration
    if (is.null(orchestration) ||
        !inherits(orchestration, "mcmc_multichain_fit"))
      stop("Multichain checkpoint does not contain valid orchestration state")
    attr(orchestration, "models") <- chain_backends
    attr(fit, "orchestration") <- orchestration
  }
  fit
}

