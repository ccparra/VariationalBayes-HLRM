# Optional immutable, content-addressed storage for complete typed checkpoints.
# This changes disk representation only: fitting and restoration remain in RAM.
.checkpoint_plain_state <- function(x) {
  if (!is.null(x) && !is.atomic(x) && !is.list(x))
    stop("Checkpoint stores require plain atomic/list state")
  if (isS4(x)) stop("Checkpoint stores do not support S4 state")
  if (is.list(x)) for (value in x) .checkpoint_plain_state(value)
  for (attribute in attributes(x)) .checkpoint_plain_state(attribute)
  invisible(NULL)
}

.checkpoint_store_config <- function(chunk_size) list(
  schema_version = "mcmc-object-store-1", chunk_size = chunk_size,
  inline_bytes = 4096L, checksum = "md5-file-bytes", serialization_version = 3L)

.checkpoint_store_put <- function(node, store) {
  temporary <- tempfile(".object-", tmpdir = file.path(store, "objects"))
  on.exit(if (file.exists(temporary)) unlink(temporary), add = TRUE)
  saveRDS(node, temporary, version = 3, compress = TRUE)
  key <- unname(tools::md5sum(temporary))
  target <- file.path(store, "objects", paste0(key, ".rds"))
  if (file.exists(target)) {
    if (!identical(unname(tools::md5sum(target)), key))
      stop("Existing checkpoint store object checksum mismatch")
  } else if (!file.rename(temporary, target)) stop("Could not publish checkpoint store object")
  key
}

.checkpoint_axis_subset <- function(x, indices, axis) {
  if (is.null(dim(x))) return(x[indices])
  args <- c(list(x), rep(list(TRUE), length(dim(x))), list(drop = FALSE))
  args[[axis + 1L]] <- indices
  do.call(`[`, args)
}

.checkpoint_store_attributes <- function(x, store, config) {
  lapply(attributes(x), .checkpoint_store_encode, store = store, config = config)
}

.checkpoint_store_encode <- function(x, store, config, name = "") {
  # Attributes are explicit nodes so repeated classes, names, data and model
  # specifications can be shared across backends and checkpoint generations.
  if (is.list(x)) {
    labels <- names(x)
    children <- vapply(seq_along(x), function(i)
      .checkpoint_store_encode(x[[i]], store, config,
        if (is.null(labels)) "" else labels[i]), character(1L))
    node <- list(kind = "list", children = children,
      attributes = .checkpoint_store_attributes(x, store, config))
  } else if (is.null(x) || (as.numeric(utils::object.size(x)) <= config$inline_bytes &&
              !name %in% c("Beta", "sigma2", "Betaj", "sigma2j", "z", "omega",
                "Betak", "sigma2k", "Sigma", "xi2", "loglik", "logjoint",
                "log_likelihood", "log_joint", "draws", "group_allocations", "permutations")) ||
             (is.object(x) && !is.matrix(x) && !is.array(x))) {
    node <- list(kind = "literal", value = x)
  } else {
    dimensions <- dim(x)
    axis <- if (!is.null(dimensions) && length(dimensions) >= 2L &&
      name %in% c("z", "group_allocations")) 2L else 1L
    count <- if (is.null(dimensions)) length(x) else dimensions[axis]
    starts <- if (count) seq.int(1L, count, by = config$chunk_size) else integer()
    chunks <- vapply(starts, function(first) {
      indices <- seq.int(first, min(as.double(first) + config$chunk_size - 1, count))
      values <- .checkpoint_axis_subset(x, indices, axis)
      attributes(values) <- NULL
      .checkpoint_store_put(list(kind = "values", values = values), store)
    }, character(1L))
    node <- list(kind = "atomic", storage_type = typeof(x), dimensions = dimensions,
      length = length(x), axis = axis, count = count, chunk_size = config$chunk_size,
      chunks = chunks, attributes = .checkpoint_store_attributes(x, store, config))
  }
  .checkpoint_store_put(node, store)
}

.checkpoint_store_get <- function(key, store, cache) {
  if (!is.character(key) || length(key) != 1L || is.na(key) ||
      !grepl("^[0-9a-f]{32}$", key)) stop("Invalid checkpoint store object reference")
  if (exists(key, cache, inherits = FALSE)) return(get(key, cache, inherits = FALSE))
  path <- file.path(store, "objects", paste0(key, ".rds"))
  if (!file.exists(path)) stop("Missing checkpoint store object: ", key)
  if (!identical(unname(tools::md5sum(path)), key))
    stop("Checkpoint store object checksum mismatch: ", key)
  node <- readRDS(path)
  if (!is.list(node) || !is.character(node$kind) || length(node$kind) != 1L)
    stop("Invalid checkpoint store node")
  assign(key, node, cache)
  node
}

.checkpoint_store_decode <- function(key, store, cache, decoded) {
  if (exists(key, decoded, inherits = FALSE)) return(get(key, decoded, inherits = FALSE))
  node <- .checkpoint_store_get(key, store, cache)
  if (identical(node$kind, "literal")) {
    value <- node$value
  } else if (identical(node$kind, "list")) {
    value <- lapply(node$children, .checkpoint_store_decode,
      store = store, cache = cache, decoded = decoded)
    attributes(value) <- if (!length(node$attributes)) NULL else
      lapply(node$attributes, .checkpoint_store_decode, store = store, cache = cache, decoded = decoded)
  } else if (identical(node$kind, "atomic")) {
    if (!node$storage_type %in% c("logical", "integer", "double", "complex", "character", "raw") ||
        !is.numeric(node$length) || length(node$length) != 1L ||
        !is.finite(node$length) || node$length < 0 ||
        node$length != floor(node$length) ||
        length(node$chunks) != ceiling(node$count / node$chunk_size))
      stop("Invalid checkpoint store atomic layout")
    value <- vector(node$storage_type, node$length)
    if (!is.null(node$dimensions)) dim(value) <- node$dimensions
    first <- 1L
    for (chunk in node$chunks) {
      block <- .checkpoint_store_get(chunk, store, cache)
      indices <- seq.int(first, min(as.double(first) + node$chunk_size - 1, node$count))
      dimensions <- node$dimensions
      if (!is.null(dimensions)) dimensions[node$axis] <- length(indices)
      expected <- if (is.null(dimensions)) length(indices) else prod(as.double(dimensions))
      if (!identical(block$kind, "values") || typeof(block$values) != node$storage_type ||
          length(block$values) != expected) stop("Invalid checkpoint store value block")
      if (is.null(dimensions)) {
        value[indices] <- block$values
      } else {
        args <- c(list(value), rep(list(TRUE), length(dimensions)), list(value = block$values))
        args[[node$axis + 1L]] <- indices
        value <- do.call(`[<-`, args)
      }
      first <- first + length(indices)
    }
    attributes(value) <- if (!length(node$attributes)) NULL else
      lapply(node$attributes, .checkpoint_store_decode, store = store, cache = cache, decoded = decoded)
  } else stop("Invalid checkpoint store node kind")
  assign(key, value, decoded)
  value
}

.save_mcmc_checkpoint_store <- function(fit, path, store, chunk_size) {
  if (!is.character(store) || length(store) != 1L || is.na(store) || !nzchar(store))
    stop("store must be one nonempty directory path")
  chunk_size <- validate_integer_scalar(chunk_size, "chunk_size", 1L)
  parent <- normalizePath(dirname(path), winslash = "/", mustWork = TRUE)
  # Resolve existing symlinks as well as relative paths. Keep the object store
  # within the manifest's directory tree so a complete bundle is relocatable.
  if (!dir.exists(store)) dir.create(store, recursive = TRUE)
  store <- normalizePath(store, winslash = "/", mustWork = TRUE)
  if (!startsWith(store, paste0(parent, "/")))
    stop("store must be a subdirectory of the checkpoint manifest directory")
  relative <- substring(store, nchar(parent) + 2L)
  lock <- file.path(store, ".write-lock")
  if (!dir.create(lock, showWarnings = FALSE))
    stop("Checkpoint store is locked; only one writer is supported")
  on.exit(unlink(lock, recursive = TRUE), add = TRUE)
  writeLines(c(paste0("pid=", Sys.getpid()), paste0("host=", Sys.info()[["nodename"]])),
    file.path(lock, "owner"))
  config <- .checkpoint_store_config(chunk_size)
  config_path <- file.path(store, "store.rds")
  if (file.exists(config_path)) {
    if (!identical(readRDS(config_path), config)) stop("Checkpoint store configuration differs")
  } else {
    if (length(setdiff(list.files(store, all.files = TRUE), c(".", "..", ".write-lock"))))
      stop("Refusing to initialize a nonempty checkpoint store")
    temporary <- tempfile(".config-", tmpdir = store)
    on.exit(if (file.exists(temporary)) unlink(temporary), add = TRUE)
    saveRDS(config, temporary, version = 3)
    if (!file.rename(temporary, config_path)) stop("Could not publish checkpoint store configuration")
  }
  dir.create(file.path(store, "objects"), showWarnings = FALSE)
  checkpoint <- .typed_mcmc_checkpoint(fit)
  .checkpoint_plain_state(checkpoint)
  root <- .checkpoint_store_encode(checkpoint, store, config)
  manifest <- structure(list(schema_version = "typed-mcmc-checkpoint-store-1",
    model = fit$model, store = relative,
    store_config_md5 = unname(tools::md5sum(config_path)), root = root),
    class = c("vb_mcmc_store_checkpoint", "list"))
  temporary_manifest <- tempfile(".checkpoint-manifest-", tmpdir = parent)
  on.exit(if (file.exists(temporary_manifest)) unlink(temporary_manifest), add = TRUE)
  saveRDS(manifest, temporary_manifest, version = 3)
  if (file.exists(path) || !file.rename(temporary_manifest, path))
    stop("Could not atomically publish checkpoint store manifest")
  invisible(normalizePath(path, mustWork = TRUE))
}

.read_mcmc_checkpoint_store <- function(manifest, path) {
  if (!identical(manifest$schema_version, "typed-mcmc-checkpoint-store-1") ||
      !is.character(manifest$store) || length(manifest$store) != 1L ||
      is.na(manifest$store) || !nzchar(manifest$store) ||
      grepl("^[/\\\\]|^[A-Za-z]:", manifest$store) ||
      any(strsplit(manifest$store, "[/\\\\]")[[1L]] == ".."))
    stop("Invalid checkpoint store manifest")
  parent <- normalizePath(dirname(path), winslash = "/", mustWork = TRUE)
  store <- normalizePath(file.path(parent, manifest$store), winslash = "/", mustWork = TRUE)
  if (!identical(store, parent) && !startsWith(store, paste0(parent, "/")))
    stop("Checkpoint store resolves outside the manifest directory tree")
  config_path <- file.path(store, "store.rds")
  if (!file.exists(config_path) ||
      !identical(unname(tools::md5sum(config_path)), manifest$store_config_md5))
    stop("Checkpoint store configuration checksum mismatch")
  config <- readRDS(config_path)
  if (!identical(config, .checkpoint_store_config(config$chunk_size)))
    stop("Unsupported checkpoint store configuration")
  checkpoint <- .checkpoint_store_decode(manifest$root, store,
    new.env(parent = emptyenv()), new.env(parent = emptyenv()))
  if (!identical(checkpoint$model, manifest$model)) stop("Checkpoint store model differs")
  checkpoint
}
