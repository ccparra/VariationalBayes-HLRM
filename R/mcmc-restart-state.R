# Restart primitives contain current sampler/RNG state, never past draw arrays.
.restart_state_budget <- function(value) {
  if (!is.numeric(value) || length(value) != 1L || !is.finite(value) || value <= 0)
    stop("max_state_mb/max_memory_mb must be a positive finite scalar")
  value * 1024^2
}

.validate_restart_rng <- function(state) {
  if (!is.list(state) || !identical(names(state), c("kind", "seed")) ||
      !is.character(state$kind) || length(state$kind) != 3L ||
      state$kind[2L] != "Inversion" || state$kind[1L] == "user-supplied" ||
      !is.integer(state$seed) || length(state$seed) < 2L)
    stop("Restart requires a saved built-in RNG with Inversion normals")
  previous <- rng_snapshot(); on.exit(restore_rng_snapshot(previous), add = TRUE)
  do.call(RNGkind, as.list(state$kind)); set.seed(1L)
  expected <- rng_snapshot()$seed
  if (length(expected) != length(state$seed) || !identical(expected[1L], state$seed[1L]))
    stop("Invalid restart RNG encoding")
  invisible(NULL)
}

.validate_mcmc_restart_state <- function(state, max_state_mb = 256) {
  budget <- .restart_state_budget(max_state_mb)
  if (RNGkind()[2L] != "Inversion")
    stop("Restart primitives require the caller to use Inversion normals")
  if (!inherits(state, "vb_mcmc_restart_state") ||
      !identical(names(state), c("schema_version", "model", "context", "chains")) ||
      !identical(state$schema_version, "mcmc-restart-state-1") ||
      !is.character(state$model) || length(state$model) != 1L ||
      !state$model %in% c("LRM", "HLRM", "CHLRM") ||
      !is.list(state$context) || !identical(names(state$context), c("specifications", "scaling_map")) ||
      !is.list(state$chains) || !length(state$chains)) stop("Invalid MCMC restart state")
  specification <- state$context$specifications
  if (!is.list(specification) || !inherits(specification$model, "vb_model_spec") ||
      !identical(specification$model$model, state$model) ||
      !identical(specification$control$method, "MCMC") ||
      specification$control$chains != length(state$chains)) stop("Restart context disagrees with model/chains")
  .checkpoint_plain_state(state)
  if (as.numeric(object.size(state)) > budget) stop("Restart state exceeds max_state_mb")
  for (i in seq_along(state$chains)) {
    chain <- state$chains[[i]]
    if (!is.list(chain) || !identical(names(chain), c("chain", "sampler", "run_settings", "rng_state")) ||
        !identical(chain$chain, as.integer(i)) || !is.list(chain$sampler) ||
        !setequal(names(chain$sampler), c("stats", "params", "hyperparams")))
      stop("Invalid restart chain or sampler layout")
    validate_mcmc_settings(chain$run_settings)
    if (!identical(as.integer(chain$run_settings$thin), 1L)) stop("Restart requires thin = 1")
    .validate_restart_rng(chain$rng_state)
    stats <- chain$sampler$stats; params <- chain$sampler$params
    dims <- .validate_multichain_data(state$model, stats$y, stats$X,
      if (state$model == "CHLRM") stats$K else NULL)
    if (stats$p != dims$p || (state$model != "LRM" && stats$m != dims$m))
      stop("Restart sampler dimensions disagree with data")
    if ((!is.null(specification$model$coefficient_count) && specification$model$coefficient_count != dims$p) ||
        (state$model == "CHLRM" && specification$model$component_count != stats$K))
      stop("Restart context disagrees with sampler dimensions")
    p <- dims$p; m <- dims$m; K <- stats$K
    expected <- switch(state$model,
      LRM = list(Beta = c(p), sigma2 = 1L),
      HLRM = list(Betaj = c(m, p), sigma2j = c(m), Beta = c(p), Sigma = c(p, p), xi2 = 1L),
      CHLRM = list(omega = c(K), z = c(m), Beta = c(p), Sigma = c(p, p),
        xi2 = 1L, Betak = c(K, p), sigma2k = c(K)))
    if (!setequal(names(params), names(expected))) stop("Invalid restart parameter fields")
    for (name in names(expected)) {
      x <- params[[name]]; shape <- expected[[name]]
      if (!is.numeric(x) || any(!is.finite(x)) ||
          (length(shape) == 1L && (length(x) != shape || !is.null(dim(x)))) ||
          (length(shape) > 1L && !identical(dim(x), as.integer(shape))))
        stop("Invalid restart parameter dimensions: ", name)
    }
    .validate_initialization_state(params, state$model, "MCMC")
  }
  invisible(state)
}

.new_mcmc_restart_state <- function(model, context, chains, max_state_mb) {
  state <- structure(list(schema_version = "mcmc-restart-state-1", model = model,
    context = context, chains = chains), class = c("vb_mcmc_restart_state", "list"))
  .validate_mcmc_restart_state(state, max_state_mb)
  .clone_mcmc_state(state)
}

#' Extract or advance a history-free MCMC restart state
#'
#' These low-level operations preserve the existing Gibbs transition order and
#' saved per-chain RNG stream. They never repeat warm-up, append past histories,
#' relabel draws or certify convergence. The input is unchanged.
#'
#' @param fit A complete MCMC `vb_model_fit`, not a reporting view.
#' @param max_state_mb Guard on the plain restart state's in-memory size, in MiB.
#' @return `mcmc_restart_state()` returns a plain `vb_mcmc_restart_state` with
#'   current sampler/data/prior/RNG state, original settings and scaling context.
#'   It contains no past draws and is not a fitted model or posterior summary.
#'   Context records the extraction-time specifications; current retained counts
#'   are in each chain's `run_settings$samples`, which can differ during progress.
#' @details Restart primitives require Inversion normal generation for both the
#'   saved streams and caller. R's `.Random.seed` does not retain Box-Muller
#'   cached normal state. Built-in uniform generators, including L'Ecuyer-CMRG,
#'   are supported. Caller RNG state, including seed absence, is preserved.
#'
#'   Native block parameters remain on the sampler's working scale and retain
#'   raw component labels. The scaling context is carried unchanged. These
#'   primitives are not a complete disk-backed fitting/reporting interface;
#'   retain original histories separately and keep using the existing fit API
#'   for its posterior reports and scientific diagnostics.
#' @export
mcmc_restart_state <- function(fit, max_state_mb = 256) {
  .require_complete_mcmc_fit(fit)
  if (!inherits(fit, "vb_model_fit") || !identical(fit$method, "MCMC") ||
      !is.null(attr(fit, "reporting_view", exact = TRUE))) stop("Require a complete MCMC fit")
  backends <- attr(fit, "backends", exact = TRUE)
  if (is.null(backends)) backends <- list(.fit_backend(fit))
  chains <- lapply(seq_along(backends), function(i) {
    b <- backends[[i]]
    if (!identical(b$fit_status, "fitted") || .mcmc_model_name(b) != fit$model)
      stop("Require successfully fitted raw chain backends")
    list(chain = as.integer(i), sampler = b[[.mcmc_object_name(fit$model)]],
      run_settings = b$run_settings, rng_state = b$rng_state)
  })
  .new_mcmc_restart_state(fit$model, list(specifications = fit$specifications,
    scaling_map = attr(fit, "scaling_map", exact = TRUE)), chains, max_state_mb)
}

#' @rdname mcmc_restart_state
#' @param state A `vb_mcmc_restart_state`.
#' @param additional_samples Number of new unthinned retained draws for one chain.
#' @param chain One-based chain identity. Defaults to one for a single-chain state;
#'   must be supplied explicitly for a multichain state. Chains may progress
#'   separately; scientific diagnostics need a complete matched prefix.
#' @param max_memory_mb Analytical active-state/block memory guard, in MiB.
#'   This is not an operating-system RSS limit and does not use cumulative length.
#' @param verbose Show the existing sampler progress display.
#' @return `advance_mcmc_state()` returns a `vb_mcmc_draw_block` containing only
#'   the new native `draws`, their per-chain `draw_index`, next `state`, timings
#'   and a resource estimate. It is transactional on errors and does not write
#'   or discard histories. There is no seed override; continuation uses the
#'   original saved stream. Existing full-fit guards remain unchanged.
#' @export
advance_mcmc_state <- function(state, additional_samples, chain = NULL,
    max_memory_mb = 512, verbose = FALSE) {
  budget <- .restart_state_budget(max_memory_mb)
  .validate_mcmc_restart_state(state, max_state_mb = max_memory_mb)
  n <- validate_integer_scalar(additional_samples, "additional_samples", 1L)
  if (is.null(chain)) {
    if (length(state$chains) != 1L) stop("Specify chain for a multichain restart state")
    chain <- 1L
  }
  chain <- validate_integer_scalar(chain, "chain", 1L)
  if (chain > length(state$chains)) stop("Unknown restart chain")
  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) stop("verbose must be TRUE or FALSE")
  current <- state$chains[[chain]]; old <- current$run_settings$samples
  total <- as.double(old) + n
  if (total > .Machine$integer.max) stop("Cumulative retained iteration exceeds integer range")
  stats <- current$sampler$stats
  dimensions <- list(p = stats$p, m = if (state$model == "LRM") 1L else stats$m)
  payload <- 8 * as.double(n) * .mcmc_values_per_draw(state$model, dimensions, stats$K)
  reserve <- 1.35 * (3 * as.numeric(object.size(state)) + 2 * payload + 24 * n)
  if (!is.finite(reserve) || reserve > budget) stop("Active restart/block estimate exceeds max_memory_mb")
  next_state <- .clone_mcmc_state(state)
  scope <- inference_rng_scope(state = current$rng_state)
  on.exit(scope$restore(), add = TRUE)
  block <- .sample_mcmc_block(next_state$chains[[chain]]$sampler, state$model, n, verbose)
  next_state$chains[[chain]]$sampler <- block$object
  next_state$chains[[chain]]$rng_state <- rng_snapshot()
  next_state$chains[[chain]]$run_settings$samples <- as.integer(total)
  structure(list(schema_version = "mcmc-draw-block-1", model = state$model,
    chain = chain, draw_index = data.frame(chain = rep.int(chain, n),
      iteration = seq.int(as.integer(old) + 1L, as.integer(total))),
    draws = block$draws, state = next_state, timing = block$timing,
    resources = list(estimated_active_peak_mib = reserve / 1024^2,
      max_memory_mb = max_memory_mb, cumulative_history_allocated = FALSE)),
    class = c("vb_mcmc_draw_block", "list"))
}

.restart_file_magic <- "VBHLRM-MCMC-RESTART-1\n"

.publish_mcmc_restart_file <- function(temporary, path) {
  if (file.exists(path) || !file.rename(temporary, path)) stop("Could not publish restart file")
  invisible(NULL)
}

#' @rdname mcmc_restart_state
#' @param path A new restart file, or an existing restart/store-manifest file.
#' @return `save_mcmc_restart_state()` atomically writes a new checksum-protected
#'   state-only file and returns its path. It refuses overwrite. This is a
#'   dedicated versioned binary format; read it with `read_mcmc_restart_state()`.
#' @export
save_mcmc_restart_state <- function(state, path, max_state_mb = 256) {
  .validate_mcmc_restart_state(state, max_state_mb)
  if (!is.character(path) || length(path) != 1L || is.na(path) || !nzchar(path)) stop("Invalid restart path")
  if (file.exists(path)) stop("Restart path already exists")
  parent <- dirname(path); dir.create(parent, recursive = TRUE, showWarnings = FALSE)
  lock <- paste0(path, ".write-lock")
  if (!dir.create(lock, showWarnings = FALSE)) stop("Restart writer is locked")
  on.exit(unlink(lock, recursive = TRUE), add = TRUE)
  if (file.exists(path)) stop("Restart path already exists")
  payload <- tempfile(".restart-payload-", tmpdir = parent)
  temporary <- tempfile(".restart-file-", tmpdir = parent)
  input <- output <- NULL
  on.exit({try(close(input), silent = TRUE); try(close(output), silent = TRUE);
    unlink(c(payload, temporary))}, add = TRUE)
  saveRDS(state, payload, version = 3, compress = FALSE)
  input <- file(payload, "rb")
  output <- file(temporary, "wb")
  writeBin(charToRaw(paste0(.restart_file_magic, unname(tools::md5sum(payload)), "\n",
    format(file.info(payload)$size, scientific = FALSE, trim = TRUE), "\n")), output)
  repeat {
    bytes <- readBin(input, "raw", n = 1048576L)
    if (!length(bytes)) break
    writeBin(bytes, output)
  }
  close(output)
  .publish_mcmc_restart_file(temporary, path)
  invisible(normalizePath(path, mustWork = TRUE))
}

#' @rdname mcmc_restart_state
#' @param format `"restart"` reads the dedicated state-only format. Explicit
#'   `"checkpoint_store"` selectively reads sampler/settings/RNG/context nodes
#'   from a complete immutable object-store checkpoint without decoding histories.
#'   Legacy single-file full checkpoints are unsupported here; the default
#'   format rejects them before deserializing their payload.
#' @return `read_mcmc_restart_state()` returns the same plain restart state.
#'   Store mode validates selected nodes only, not historical draw chunks; it
#'   does not certify that the separate complete draw history is intact.
#' @export
read_mcmc_restart_state <- function(path, format = c("restart", "checkpoint_store"), max_state_mb = 256) {
  format <- match.arg(format); budget <- .restart_state_budget(max_state_mb)
  if (!is.character(path) || length(path) != 1L || is.na(path) || !file.exists(path)) stop("Restart file does not exist")
  if (format == "checkpoint_store") return(.read_restart_from_store(path, max_state_mb))
  input <- file(path, "rb"); on.exit(close(input), add = TRUE)
  if (!identical(readBin(input, "raw", n = nchar(.restart_file_magic)), charToRaw(.restart_file_magic)))
    stop("Unsupported restart format; full single-file checkpoints are not state-only files")
  checksum <- rawToChar(readBin(input, "raw", n = 33L))
  if (!grepl("^[0-9a-f]{32}\n$", checksum)) stop("Invalid restart checksum header")
  count <- raw()
  repeat {
    b <- readBin(input, "raw", n = 1L)
    if (!length(b) || length(count) > 20L) stop("Invalid restart size header")
    if (identical(b, as.raw(10))) break
    count <- c(count, b)
  }
  count_text <- rawToChar(count)
  if (!grepl("^[0-9]+$", count_text)) stop("Invalid restart size header")
  count <- as.numeric(count_text)
  if (!is.finite(count) || count > budget || count < 1 || file.info(path)$size - seek(input) != count)
    stop("Restart payload size differs or exceeds max_state_mb")
  temporary <- tempfile("restart-read-"); output <- NULL
  on.exit({try(close(output), silent = TRUE); unlink(temporary)}, add = TRUE)
  output <- file(temporary, "wb")
  remaining <- count
  while (remaining > 0) {
    b <- readBin(input, "raw", n = min(remaining, 1048576))
    if (!length(b)) stop("Truncated restart payload")
    writeBin(b, output); remaining <- remaining - length(b)
  }
  close(output)
  if (unname(tools::md5sum(temporary)) != substr(checksum, 1L, 32L)) stop("Restart payload checksum mismatch")
  state <- readRDS(temporary)
  .validate_mcmc_restart_state(state, max_state_mb)
  state
}
