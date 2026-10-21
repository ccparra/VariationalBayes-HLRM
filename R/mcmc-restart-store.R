# Read only named current-state subtrees of an existing complete store manifest.
.read_restart_from_store <- function(path, max_state_mb) {
  if (file.info(path)$size > 16384) stop("Checkpoint store manifest exceeds metadata input guard")
  manifest <- readRDS(path)
  store <- .checkpoint_store_location(manifest, path)
  cache <- new.env(parent = emptyenv()); decoded <- new.env(parent = emptyenv())
  get <- function(key) .checkpoint_store_get(key, store, cache)
  decode <- function(key) .checkpoint_store_decode(key, store, cache, decoded)
  child <- function(key, name) {
    node <- get(key)
    if (!identical(node$kind, "list") || is.null(node$attributes$names)) stop("Invalid named restart metadata")
    labels <- decode(node$attributes$names)
    if (!is.character(labels) || length(labels) != length(node$children) || anyDuplicated(labels) || !name %in% labels)
      stop("Missing or duplicated restart metadata: ", name)
    node$children[match(name, labels)]
  }
  root <- get(manifest$root)
  if (!identical(root$kind, "list") ||
      !identical(decode(child(manifest$root, "schema_version")), "typed-mcmc-checkpoint-1"))
    stop("Require a complete typed checkpoint store")
  model <- decode(child(manifest$root, "model"))
  if (!identical(model, manifest$model)) stop("Checkpoint model differs")
  snapshot <- child(manifest$root, "fit"); fit_node <- get(snapshot)
  classes <- decode(fit_node$attributes$class)
  if (!"vb_model_fit" %in% classes || "vb_mcmc_reporting_view" %in% classes ||
      !is.null(fit_node$attributes$reporting_view) ||
      !identical(decode(child(snapshot, "method")), "MCMC")) stop("Require a complete original MCMC checkpoint")
  specifications <- child(snapshot, "specifications")
  scaling <- fit_node$attributes$scaling_map
  backends <- get(child(manifest$root, "chain_backends"))
  if (backends$kind != "list" || !length(backends$children)) stop("Missing restart chain backends")
  selected <- lapply(backends$children, function(key) {
    if (!identical(decode(child(key, "schema_version")), "mcmc-checkpoint-1") ||
        !identical(decode(child(key, "model")), model)) stop("Invalid restart backend metadata")
    c(sampler = child(key, "sampler"), run_settings = child(key, "run_settings"), rng_state = child(key, "rng_state"))
  })
  # Walk selected metadata, not value chunks or any history. This output-size
  # guard is distinct from whole-process RSS or adversarial R deserialization.
  visited <- new.env(parent = emptyenv()); total <- 0
  account <- function(key) {
    if (exists(key, visited, inherits = FALSE)) return(invisible(NULL))
    assign(key, TRUE, visited); node <- get(key)
    if (node$kind == "list") {
      total <<- total + 64 + 8 * length(node$children)
      for (k in node$children) account(k)
    } else if (node$kind == "literal") {
      total <<- total + as.numeric(object.size(node$value))
    } else if (node$kind == "atomic") {
      width <- switch(node$storage_type, logical = 4, integer = 4, double = 8,
        complex = 16, raw = 1, character = 8, stop("Unsupported restart atomic type"))
      if (!is.numeric(node$length) || length(node$length) != 1L || !is.finite(node$length) || node$length < 0)
        stop("Invalid restart atomic size")
      total <<- total + 64 + node$length * width
    } else stop("Invalid restart metadata kind")
    if (total > .restart_state_budget(max_state_mb)) stop("Selected restart nodes exceed max_state_mb")
    for (k in node$attributes) account(k)
    invisible(NULL)
  }
  keys <- c(specifications, scaling, unlist(selected, use.names = FALSE))
  for (key in keys) account(key)
  chains <- lapply(seq_along(selected), function(i) {
    keys <- selected[[i]]
    list(chain = as.integer(i), sampler = decode(keys[["sampler"]]),
      run_settings = decode(keys[["run_settings"]]), rng_state = decode(keys[["rng_state"]]))
  })
  .new_mcmc_restart_state(model, list(specifications = decode(specifications),
    scaling_map = if (is.null(scaling)) NULL else decode(scaling)), chains, max_state_mb)
}
