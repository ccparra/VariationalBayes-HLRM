# Immutable history generations reference an existing complete checkpoint and
# new native draw blocks. Readers allocate selected outputs and one value chunk.
.mcmc_state_id <- function(state) {
  path <- tempfile('mcmc-state-id-'); on.exit(unlink(path), add = TRUE)
  saveRDS(state, path, version = 2, compress = FALSE)
  unname(tools::md5sum(path))
}
.history_child <- function(key, name, store, cache, decoded) {
  node <- .checkpoint_store_get(key, store, cache)
  if (node$kind != 'list' || is.null(node$attributes$names)) stop('Invalid named history metadata')
  labels <- .checkpoint_store_decode(node$attributes$names, store, cache, decoded)
  if (!is.character(labels) || length(labels) != length(node$children) || anyDuplicated(labels) || !name %in% labels)
    stop('Missing or duplicated history field: ', name)
  unname(node$children[match(name, labels)])
}
.history_relative_path <- function(target, directory) {
  a <- strsplit(normalizePath(target, winslash = '/'), '/', fixed = TRUE)[[1]]
  b <- strsplit(normalizePath(directory, winslash = '/'), '/', fixed = TRUE)[[1]]
  common <- 0L
  while (common < min(length(a), length(b)) && a[common + 1L] == b[common + 1L]) common <- common + 1L
  paste(c(rep('..', length(b) - common), a[seq.int(common + 1L, length(a))]), collapse = '/')
}
.history_checkpoint <- function(directory, baseline) {
  path <- file.path(directory, baseline$path)
  if (!file.exists(path) || unname(tools::md5sum(path)) != baseline$md5) stop('Baseline checkpoint manifest changed or missing')
  manifest <- readRDS(path); store <- .checkpoint_store_location(manifest, path)
  cache <- new.env(parent = emptyenv()); decoded <- new.env(parent = emptyenv())
  backends <- .checkpoint_store_get(.history_child(manifest$root, 'chain_backends', store, cache, decoded), store, cache)
  roots <- vapply(backends$children, .history_child, character(1L), name = 'chains',
    store = store, cache = cache, decoded = decoded)
  list(store = store, roots = unname(roots))
}
.history_snapshot_path <- function(directory, generation) file.path(directory, 'snapshots', sprintf('%08d.rds', generation))
.history_store <- function(directory) file.path(directory, 'history')
.history_config <- function(directory) {
  path <- file.path(.history_store(directory), 'store.rds')
  if (!file.exists(path)) stop('History store configuration missing')
  config <- readRDS(path)
  if (!identical(config, .checkpoint_store_config(config$chunk_size))) stop('Unsupported history store configuration')
  config
}
.history_manifest <- function(history) {
  if (!inherits(history, 'vb_mcmc_history') || !identical(names(history), c('directory', 'generation', 'manifest_md5')))
    stop('Require a history handle from create_mcmc_history()/read_mcmc_history()')
  path <- .history_snapshot_path(history$directory, history$generation)
  if (!file.exists(path) || file.info(path)$size > 16 * 1024^2 || unname(tools::md5sum(path)) != history$manifest_md5) stop('History snapshot checksum changed or missing')
  x <- readRDS(path)
  if (!identical(x$schema_version, 'mcmc-history-1') || !identical(x$generation, history$generation) ||
      !is.list(x$blocks) || !is.character(x$state_root) || !grepl('^[0-9a-f]{32}$', x$state_root) ||
      !is.character(x$state_id) || !grepl('^[0-9a-f]{32}$', x$state_id)) stop('Invalid history snapshot')
  .history_config(history$directory)
  if (!identical(x$config_md5, unname(tools::md5sum(file.path(.history_store(history$directory), 'store.rds')))))
    stop('History configuration checksum changed')
  counts <- x$baseline$counts; previous <- x$baseline$state_id
  if (!is.numeric(counts) || !length(counts) || anyNA(counts) || any(counts < 1 | counts != floor(counts)))
    stop('Invalid baseline chain counts')
  for (b in x$blocks) {
    if (!is.numeric(b$chain) || length(b$chain) != 1L || is.na(b$chain) ||
        b$chain < 1L || b$chain > length(counts) || b$chain != floor(b$chain) ||
        b$first != counts[b$chain] + 1L || b$last < b$first ||
        !identical(b$parent_state_id, previous) || !b$representation %in% c('native-1', 'chlrm-derived-groups-1') ||
        !is.character(b$root) || !grepl('^[0-9a-f]{32}$', b$root)) stop('Invalid history block index')
    counts[b$chain] <- b$last; previous <- b$next_state_id
  }
  if (!identical(as.integer(counts), x$counts) || !identical(previous, x$state_id)) stop('History frontier/state differs')
  x
}
.publish_mcmc_history_snapshot <- function(x, directory) {
  path <- .history_snapshot_path(directory, x$generation)
  temporary <- tempfile('.history-snapshot-', tmpdir = dirname(path))
  on.exit(unlink(temporary), add = TRUE)
  saveRDS(x, temporary, version = 3)
  if (file.exists(path) || !file.rename(temporary, path)) stop('Could not publish immutable history snapshot')
  invisible(NULL)
}

#' Immutable MCMC block history with bounded native draw access
#'
#' These operations retain an existing complete checkpoint as an immutable
#' baseline dependency and publish subsequent native blocks without rebuilding
#' past arrays. They neither fit a model nor assess convergence.
#'
#' @param checkpoint A complete immutable object-store checkpoint file.
#' @param directory A new history directory. It must not exist.
#' @param chunk_size Maximum iterations per newly encoded value chunk.
#' @param max_state_mb Guard on current restart state, in MiB.
#' @return `create_mcmc_history()` and `append_mcmc_history()` return an immutable
#'   `vb_mcmc_history` snapshot handle. Old handles remain readable. This is not
#'   a `vb_model_fit` or a posterior reporting/convergence certificate.
#' @details The baseline remains external and checksum-bound by relative path;
#'   retain it and its object store. To relocate, copy the entire dependency tree
#'   preserving relative paths. No existing checkpoint/history is deleted.
#'
#'   New CHLRM group-coefficient and group-variance draws are checked exactly
#'   against component draws/allocations, then encoded as a lossless derived
#'   representation. Native reads reconstruct those values on demand. Baseline
#'   native fields remain exactly as originally stored. Other native fields and
#'   models are stored without this transformation.
#'
#'   Writers use a single-writer lock and publish snapshots last. A failed
#'   publication leaves old snapshots valid and may leave unreferenced objects;
#'   no automatic garbage collection or sampling retry occurs. Reload the latest
#'   committed snapshot after an interrupted call. Stale append handles, duplicate
#'   blocks, gaps, wrong parent states and changes to other chains are rejected.
#' @export
create_mcmc_history <- function(checkpoint, directory, chunk_size = 1000L, max_state_mb = 256) {
  chunk_size <- validate_integer_scalar(chunk_size, 'chunk_size', 1L)
  if (!is.character(directory) || length(directory) != 1L || is.na(directory) || !nzchar(directory)) stop('Invalid history directory')
  if (file.exists(directory)) stop('History directory already exists')
  state <- read_mcmc_restart_state(checkpoint, 'checkpoint_store', max_state_mb)
  checkpoint <- normalizePath(checkpoint, mustWork = TRUE)
  if (!dir.create(directory, recursive = TRUE, showWarnings = FALSE)) stop('Could not create new history directory')
  directory <- normalizePath(directory)
  dir.create(file.path(directory, 'snapshots'))
  store <- .history_store(directory); dir.create(store); dir.create(file.path(store, 'objects'))
  config <- .checkpoint_store_config(chunk_size); saveRDS(config, file.path(store, 'store.rds'), version = 3)
  root <- .checkpoint_store_encode(state, store, config)
  x <- list(schema_version = 'mcmc-history-1', generation = 0L, model = state$model,
    baseline = list(path = .history_relative_path(checkpoint, directory), md5 = unname(tools::md5sum(checkpoint)),
      counts = vapply(state$chains, function(x) as.integer(x$run_settings$samples), integer(1L)),
      state_id = .mcmc_state_id(state)),
    counts = vapply(state$chains, function(x) as.integer(x$run_settings$samples), integer(1L)),
    config_md5 = unname(tools::md5sum(file.path(store, 'store.rds'))),
    state_root = root, state_id = .mcmc_state_id(state), blocks = list())
  .publish_mcmc_history_snapshot(x, directory)
  read_mcmc_history(directory)
}

#' @rdname create_mcmc_history
#' @param generation Optional nonnegative snapshot number; default is latest
#'   committed snapshot. Snapshot checksums bind a handle for subsequent reads.
#' @return `read_mcmc_history()` returns a snapshot handle; it does not restore a fit.
#' @export
read_mcmc_history <- function(directory, generation = NULL) {
  directory <- normalizePath(directory, mustWork = TRUE)
  files <- list.files(file.path(directory, 'snapshots'), pattern = '^[0-9]{8}[.]rds$')
  if (!length(files)) stop('No committed history snapshot')
  if (is.null(generation)) generation <- max(as.integer(sub('[.]rds$', '', files)))
  generation <- validate_integer_scalar(generation, 'generation', 0L)
  path <- .history_snapshot_path(directory, generation)
  if (!file.exists(path)) stop('Requested history generation does not exist')
  out <- structure(list(directory = directory, generation = generation,
    manifest_md5 = unname(tools::md5sum(path))), class = c('vb_mcmc_history', 'list'))
  .history_manifest(out)
  out
}

#' @rdname create_mcmc_history
#' @param history A `vb_mcmc_history` snapshot handle.
#' @return `mcmc_history_restart_state()` returns its exact current restart state,
#'   with no draw history restored.
#' @export
mcmc_history_restart_state <- function(history, max_state_mb = 256) {
  x <- .history_manifest(history); store <- .history_store(history$directory)
  cache <- new.env(parent = emptyenv()); decoded <- new.env(parent = emptyenv())
  # Existing selective restart output guard applies after this state-only decode;
  # this is not an RSS guarantee or an untrusted-deserialization sandbox.
  state <- .checkpoint_store_decode(x$state_root, store, cache, decoded)
  .validate_mcmc_restart_state(state, max_state_mb)
  if (.mcmc_state_id(state) != x$state_id) stop('History restart state identity changed')
  state
}

.history_validate_block <- function(block, state, chunk_size) {
  if (!inherits(block, 'vb_mcmc_draw_block') || !identical(block$schema_version, 'mcmc-draw-block-1') ||
      !identical(block$model, state$model) || !identical(block$parent_state_id, .mcmc_state_id(state)))
    stop('Block model/parent state differs from committed history')
  .validate_mcmc_restart_state(block$state)
  chain <- validate_integer_scalar(block$chain, 'chain', 1L)
  if (chain > length(state$chains) || !identical(block$state$context, state$context) ||
      !identical(block$state$model, state$model) || length(block$state$chains) != length(state$chains)) stop('Block context changed')
  for (i in setdiff(seq_along(state$chains), chain))
    if (!identical(block$state$chains[[i]], state$chains[[i]])) stop('Block changed another chain')
  old <- state$chains[[chain]]; next_chain <- block$state$chains[[chain]]
  n <- next_chain$run_settings$samples - old$run_settings$samples
  if (n < 1 || !identical(block$draw_index$chain, rep.int(chain, n)) ||
      !identical(block$draw_index$iteration, seq.int(as.integer(old$run_settings$samples) + 1L,
        as.integer(next_chain$run_settings$samples)))) stop('Block iterations have gaps, duplicates or wrong chain')
  settings <- old$run_settings; settings$samples <- as.integer(settings$samples + n)
  if (!identical(settings, next_chain$run_settings) ||
      !identical(old$sampler$stats, next_chain$sampler$stats) ||
      !identical(old$sampler$hyperparams, next_chain$sampler$hyperparams)) stop('Block changed sampler settings/data/prior')
  template <- .allocate_mcmc_draws(state$model, old$sampler$stats, 1L)
  if (!identical(names(block$draws), names(template))) stop('Native block fields changed')
  for (name in names(template)) {
    value <- block$draws[[name]]; expected <- template[[name]]
    axis <- if (state$model == 'CHLRM' && name == 'z') 2L else 1L
    shape <- dim(expected)
    if (is.null(shape)) valid <- is.null(dim(value)) && length(value) == n else {
      shape[axis] <- n; valid <- identical(dim(value), as.integer(shape))
    }
    if (!valid || !identical(typeof(value), typeof(expected))) stop('Native block shape/type changed: ', name)
    for (first in seq.int(1L, n, by = chunk_size)) {
      indices <- seq.int(first, min(n, first + chunk_size - 1L))
      if (any(!is.finite(.checkpoint_axis_subset(value, indices, axis)))) stop('Native block contains nonfinite draws')
    }
  }
  params <- next_chain$sampler$params
  for (name in names(params)) {
    value <- block$draws[[name]]
    actual <- if (name == 'z') value[, n] else if (is.null(dim(value))) value[n] else
      if (length(dim(value)) == 2L) value[n, ] else value[n, , ]
    if (!identical(as.vector(actual), as.vector(params[[name]])))
      stop('Terminal native draw disagrees with restart sampler: ', name)
  }
  if (state$model == 'CHLRM') {
    m <- old$sampler$stats$m; p <- old$sampler$stats$p
    for (j in seq_len(m)) {
      z <- block$draws$z[j, ]
      if (any(z < 1L | z > old$sampler$stats$K)) stop('Invalid native allocations')
      if (!identical(as.numeric(block$draws$sigma2j[, j]),
          as.numeric(block$draws$sigma2k[cbind(seq_len(n), z)]))) stop('Derived group variance differs')
      for (k in seq_len(p)) if (!identical(as.numeric(block$draws$Betaj[, j, k]),
          as.numeric(block$draws$Betak[cbind(seq_len(n), z, k)]))) stop('Derived group coefficient differs')
    }
  }
  invisible(n)
}

#' @rdname create_mcmc_history
#' @param block A new block returned by `advance_mcmc_state()` with its parent
#'   state identity. It must continue the latest committed snapshot exactly.
#' @return The append handle identifies the newly committed immutable generation.
#' @export
append_mcmc_history <- function(history, block) {
  directory <- history$directory; x <- .history_manifest(history)
  lock <- file.path(directory, '.write-lock')
  if (!dir.create(lock, showWarnings = FALSE)) stop('History writer is locked')
  on.exit(unlink(lock, recursive = TRUE), add = TRUE)
  if (read_mcmc_history(directory)$generation != history$generation) stop('Stale history snapshot cannot append')
  state <- mcmc_history_restart_state(history); config <- .history_config(directory)
  n <- .history_validate_block(block, state, config$chunk_size)
  draws <- block$draws
  representation <- 'native-1'
  if (x$model == 'CHLRM') {
    draws$Betaj <- NULL; draws$sigma2j <- NULL; representation <- 'chlrm-derived-groups-1'
  }
  store <- .history_store(directory)
  root <- .checkpoint_store_encode(draws, store, config)
  next_root <- .checkpoint_store_encode(block$state, store, config)
  record <- list(chain = block$chain, first = as.integer(state$chains[[block$chain]]$run_settings$samples) + 1L,
    last = as.integer(block$state$chains[[block$chain]]$run_settings$samples),
    root = root, representation = representation, parent_state_id = x$state_id,
    next_state_id = .mcmc_state_id(block$state))
  x$generation <- x$generation + 1L; x$state_root <- next_root; x$state_id <- record$next_state_id
  x$counts[block$chain] <- record$last
  x$blocks[[length(x$blocks) + 1L]] <- record
  .publish_mcmc_history_snapshot(x, directory)
  read_mcmc_history(directory, x$generation)
}

.history_atomic_subset <- function(key, indices, store, cache, decoded, max_block_mb,
                                  field_indices = NULL) {
  node <- .checkpoint_store_get(key, store, cache)
  if (!identical(node$kind, 'atomic') || !node$storage_type %in% c('integer', 'double') ||
      length(node$axis) != 1L || !node$axis %in% c(1L, 2L) ||
      length(node$chunks) != ceiling(node$count / node$chunk_size)) stop('Invalid native history layout')
  dimensions <- node$dimensions; axis <- node$axis
  if ((is.null(dimensions) && axis != 1L) || (!is.null(dimensions) &&
      (axis > length(dimensions) || prod(as.double(dimensions)) != node$length || dimensions[axis] != node$count)))
    stop('Native history dimensions disagree')
  width <- if (node$storage_type == 'integer') 4 else 8
  per_iteration <- if (is.null(dimensions)) 1 else prod(as.double(dimensions[-axis]))
  if (!is.null(field_indices)) {
    if (is.null(dimensions) || !is.list(field_indices) ||
        length(field_indices) != length(dimensions) - 1L)
      stop('field_indices must select non-iteration dimensions of an array field')
    non_iteration_axes <- setdiff(seq_along(dimensions), axis)
    for (i in seq_along(field_indices)) {
      selected <- field_indices[[i]]
      if (!is.null(selected) && (!is.numeric(selected) || !length(selected) ||
          anyNA(selected) || any(!is.finite(selected)) || any(selected != floor(selected)) ||
          any(selected < 1 | selected > dimensions[non_iteration_axes[i]]) ||
          anyDuplicated(selected)))
        stop('field_indices contain an invalid or duplicated coordinate')
    }
  } else if (!is.null(dimensions)) {
    field_indices <- rep(list(NULL), length(dimensions) - 1L)
  }
  blocks <- unique((indices - 1L) %/% node$chunk_size + 1L)
  output_dimensions <- dimensions
  if (!is.null(output_dimensions)) output_dimensions[axis] <- length(indices)
  if (!is.null(output_dimensions) && !is.null(field_indices)) {
    non_iteration_axes <- setdiff(seq_along(output_dimensions), axis)
    for (i in seq_along(field_indices)) if (!is.null(field_indices[[i]]))
      output_dimensions[non_iteration_axes[i]] <- length(field_indices[[i]])
  }
  output_per_iteration <- if (is.null(output_dimensions)) 1 else
    prod(as.double(output_dimensions[-axis]))
  out <- vector(node$storage_type, length(indices) * output_per_iteration)
  if (!is.null(output_dimensions)) dim(out) <- output_dimensions
  for (b in blocks) {
    first <- (b - 1L) * node$chunk_size + 1L; last <- min(node$count, b * node$chunk_size)
    expected <- (last - first + 1L) * per_iteration
    if (1.5 * expected * width / 1024^2 > max_block_mb) stop('Native value chunk exceeds max_block_mb')
    value <- .checkpoint_store_get(node$chunks[b], store, cache)
    if (!identical(value$kind, 'values') || typeof(value$values) != node$storage_type || length(value$values) != expected)
      stop('Native value chunk shape/type changed')
    block <- value$values
    shape <- dimensions; if (!is.null(shape)) {shape[axis] <- last - first + 1L; dim(block) <- shape}
    positions <- which(indices >= first & indices <= last)
    selected <- .checkpoint_axis_subset(block, indices[positions] - first + 1L, axis)
    if (!is.null(field_indices) && any(!vapply(field_indices, is.null, logical(1L)))) {
      args <- c(list(selected), rep(list(TRUE), length(dim(selected))), list(drop = FALSE))
      non_iteration_axes <- setdiff(seq_along(dim(selected)), axis)
      for (i in seq_along(field_indices)) if (!is.null(field_indices[[i]]))
        args[[non_iteration_axes[i] + 1L]] <- field_indices[[i]]
      selected <- do.call(`[`, args)
      rm(args)
    }
    if (is.null(dimensions)) out[positions] <- selected else {
      args <- c(list(out), rep(list(TRUE), length(dimensions)), list(value = selected))
      args[[axis + 1L]] <- positions; out <- do.call(`[<-`, args); rm(args)
    }
    rm(value, block, selected)
  }
  out
}
.history_read_fields <- function(root, fields, indices, store, cache, decoded,
    max_block_mb, representation, stats, field_indices = NULL) {
  derived <- representation == 'chlrm-derived-groups-1' & fields %in% c('Betaj', 'sigma2j')
  values <- lapply(fields[!derived], function(name) .history_atomic_subset(
    .history_child(root, name, store, cache, decoded), indices, store, cache, decoded,
    max_block_mb, field_indices[[name]]))
  names(values) <- fields[!derived]
  n <- length(indices)
  derived_values <- list()
  if (any(derived)) {
    requested_groups <- unique(unlist(lapply(fields[derived], function(name) {
      x <- field_indices[[name]]
      if (is.null(x) || is.null(x[[1L]])) seq_len(stats$m) else x[[1L]]
    }), use.names = FALSE))
    z_indices <- list(requested_groups)
    z <- .history_atomic_subset(.history_child(root, 'z', store, cache, decoded),
      indices, store, cache, decoded, max_block_mb, z_indices)
    if ('Betaj' %in% fields[derived]) {
      beta_indices <- field_indices[['Betaj']]
      coefficients <- if (is.null(beta_indices) || is.null(beta_indices[[2L]]))
        seq_len(stats$p) else beta_indices[[2L]]
      betak <- .history_atomic_subset(.history_child(root, 'Betak', store, cache, decoded),
        indices, store, cache, decoded, max_block_mb, list(NULL, coefficients))
      a <- array(NA_real_, c(n, length(requested_groups), length(coefficients)))
      for (j in seq_along(requested_groups)) for (k in seq_along(coefficients))
        a[, j, k] <- betak[cbind(seq_len(n), z[j, ], k)]
      derived_values$Betaj <- a
      rm(betak, a)
    }
    if ('sigma2j' %in% fields[derived]) {
      sigma2k <- .history_atomic_subset(.history_child(root, 'sigma2k', store, cache, decoded),
        indices, store, cache, decoded, max_block_mb)
      a <- matrix(NA_real_, n, length(requested_groups))
      for (j in seq_along(requested_groups))
        a[, j] <- sigma2k[cbind(seq_len(n), z[j, ])]
      derived_values$sigma2j <- a
      rm(sigma2k, a)
    }
    rm(z)
  }
  for (name in fields[derived]) {
    value <- derived_values[[name]]
    selected <- field_indices[[name]]
    if (!is.null(selected) && !is.null(selected[[1L]])) {
      positions <- match(selected[[1L]], requested_groups)
      args <- c(list(value), rep(list(TRUE), length(dim(value))), list(drop = FALSE))
      args[[3L]] <- positions
      value <- do.call(`[`, args)
      rm(args)
    }
    values[[name]] <- value
  }
  values[fields]
}

#' @rdname create_mcmc_history
#' @param chain One-based chain identity.
#' @param iterations Optional strictly increasing retained iteration identities.
#'   Default requests the full selected chain, subject to the output guard.
#' @param fields Optional unique native sampler fields. Defaults to every field.
#' @param field_indices Optional named list selecting non-iteration coordinates
#'   for requested fields. Each value is a list of coordinate vectors in the
#'   original field dimension order with the iteration axis omitted; `NULL`
#'   keeps that axis in full. This allows a single parameter/group trace to be
#'   read without allocating the complete native field.
#' @param max_output_mb Analytical requested-output guard, in MiB.
#' @param max_block_mb Guard for an individual decoded numeric value chunk, in
#'   MiB. Old baseline chunk sizes are preserved, so a large baseline field may
#'   require an explicit larger chunk budget. These are not RSS limits.
#' @return `mcmc_history_draws()` returns native values, shapes and types with
#'   explicit per-chain iteration identities. Labels and working scale remain
#'   native. Reads may span the baseline and multiple immutable blocks, and
#'   decode only chunks needed by selected fields/iterations (or exact derived
#'   dependencies). It does not align labels, sample draws, or certify mixing.
#' @export
mcmc_history_draws <- function(history, chain, iterations = NULL, fields = NULL,
    max_output_mb = 256, max_block_mb = 64, field_indices = NULL) {
  .restart_state_budget(max_output_mb); .restart_state_budget(max_block_mb)
  x <- .history_manifest(history); state <- mcmc_history_restart_state(history)
  chain <- validate_integer_scalar(chain, 'chain', 1L)
  if (chain > length(state$chains)) stop('Unknown history chain')
  available <- as.integer(state$chains[[chain]]$run_settings$samples)
  if (!is.null(iterations) && (!is.numeric(iterations) || !is.null(dim(iterations)) || !length(iterations) ||
      anyNA(iterations) || any(!is.finite(iterations)) || any(iterations < 1 | iterations > available) ||
      any(iterations != floor(iterations)) || any(diff(iterations) <= 0)))
    stop('iterations must be strictly increasing retained identities in range')
  n <- if (is.null(iterations)) available else length(iterations)
  stats <- state$chains[[chain]]$sampler$stats
  template <- .allocate_mcmc_draws(x$model, stats, 1L)
  if (is.null(fields)) fields <- names(template)
  if (!is.character(fields) || !length(fields) || anyNA(fields) || anyDuplicated(fields) || any(!fields %in% names(template)))
    stop('fields must contain unique native parameter names')
  if (is.null(field_indices)) field_indices <- list()
  if (!is.list(field_indices) || (length(field_indices) &&
      (is.null(names(field_indices)) || anyNA(names(field_indices)) ||
       anyDuplicated(names(field_indices)) || any(!names(field_indices) %in% fields))))
    stop('field_indices must be a named list for requested fields')
  for (name in names(field_indices)) {
    a <- template[[name]]
    axis <- if (x$model == 'CHLRM' && identical(name, 'z')) 2L else 1L
    d <- dim(a)
    if (is.null(d) || !is.list(field_indices[[name]]) ||
        length(field_indices[[name]]) != length(d) - 1L)
      stop('field_indices must contain one selector per non-iteration array dimension')
    non_iteration_axes <- setdiff(seq_along(d), axis)
    for (i in seq_along(field_indices[[name]])) {
      selected <- field_indices[[name]][[i]]
      if (!is.null(selected) && (!is.numeric(selected) || !length(selected) ||
          anyNA(selected) || any(!is.finite(selected)) || any(selected != floor(selected)) ||
          any(selected < 1 | selected > d[non_iteration_axes[i]]) || anyDuplicated(selected)))
        stop('field_indices contain an invalid or duplicated coordinate')
    }
  }
  selected_template <- lapply(fields, function(name) {
    a <- template[[name]]; selectors <- field_indices[[name]]
    if (!is.null(selectors)) {
      axis <- if (x$model == 'CHLRM' && identical(name, 'z')) 2L else 1L
      args <- c(list(a), rep(list(TRUE), length(dim(a))), list(drop = FALSE))
      non_iteration_axes <- setdiff(seq_along(dim(a)), axis)
      for (i in seq_along(selectors)) if (!is.null(selectors[[i]]))
        args[[non_iteration_axes[i] + 1L]] <- selectors[[i]]
      a <- do.call(`[`, args)
    }
    a
  }); names(selected_template) <- fields
  bytes <- sum(vapply(selected_template, function(a) length(a) * if (is.integer(a)) 4 else 8, numeric(1L)))
  estimated <- 1.5 * as.double(n) * (bytes + 8) / 1024^2
  if (!is.finite(estimated) || estimated > max_output_mb) stop('Native draws exceed max_output_mb')
  if (is.null(iterations)) iterations <- seq_len(available) else iterations <- as.integer(iterations)
  out <- Map(function(a, name) {
    shape <- dim(a); axis <- if (x$model == 'CHLRM' && identical(name, 'z')) 2L else 1L
    value <- vector(typeof(a), as.double(length(a)) * n)
    if (!is.null(shape)) {shape[axis] <- n; dim(value) <- shape}
    value
  }, selected_template, names(selected_template))
  names(out) <- fields
  cache <- new.env(parent = emptyenv()); decoded <- new.env(parent = emptyenv())
  segments <- list(list(first = 1L, last = x$baseline$counts[chain], representation = 'native-1'))
  baseline <- NULL
  for (b in x$blocks) if (b$chain == chain) segments[[length(segments) + 1L]] <- b
  covered <- integer(length(iterations))
  for (segment in segments) {
    positions <- which(iterations >= segment$first & iterations <= segment$last)
    if (!length(positions)) next
    if (is.null(segment$root)) {
      if (is.null(baseline)) baseline <- .history_checkpoint(history$directory, x$baseline)
      store <- baseline$store; root <- baseline$roots[chain]
      local_indices <- iterations[positions]
    } else {
      store <- .history_store(history$directory); root <- segment$root
      local_indices <- iterations[positions] - segment$first + 1L
    }
    cache <- new.env(parent = emptyenv()); decoded <- new.env(parent = emptyenv())
    part <- .history_read_fields(root, fields, local_indices, store, cache, decoded,
      max_block_mb, segment$representation, stats, field_indices)
    for (name in fields) {
      a <- out[[name]]; axis <- if (x$model == 'CHLRM' && name == 'z') 2L else 1L
      if (is.null(dim(a))) a[positions] <- part[[name]] else {
        args <- c(list(a), rep(list(TRUE), length(dim(a))), list(value = part[[name]]))
        args[[axis + 1L]] <- positions; a <- do.call(`[<-`, args); rm(args)
      }
      out[[name]] <- a
    }
    covered[positions] <- covered[positions] + 1L; rm(part)
  }
  if (any(covered != 1L)) stop('History index has missing or duplicate iterations')
  structure(list(model = x$model, draw_index = data.frame(chain = rep.int(chain, n), iteration = iterations),
    draws = out, metadata = list(estimated_output_mib = estimated,
      native_scale = TRUE, aligned = FALSE, generation = history$generation)),
    class = c('vb_mcmc_history_draws', 'list'))
}
