.checkpoint_store_fit <- local({
  cache <- new.env(parent = emptyenv())
  function(model = "CHLRM", chains = 2L) {
    key <- paste(model, chains)
    if (exists(key, cache, inherits = FALSE)) return(cache[[key]])
    X <- cbind(1, seq(-1, 1, length.out = 8))
    data <- list(X = rep(list(X), 6), y = lapply(seq_len(6), function(j)
      drop(X %*% c(j / 3, (-1)^j)) + sin(seq_len(8)) / 10))
    spec <- switch(model, LRM = lrm_model(2), HLRM = hlrm_model(2), CHLRM = chlrm_model(2, 2))
    prior <- switch(model,
      LRM = lrm_prior(c(0, 0), diag(2), 4, 1),
      HLRM = hlrm_prior(c(0, 0), diag(2), 4, diag(2), 4, 2, 2),
      CHLRM = chlrm_prior(c(0, 0), diag(2), 4, diag(2), 4, 2, 2, c(1, 1)))
    if (model == "LRM") data <- list(X = X, y = data$y[[1L]])
    cache[[key]] <- fit_model(spec, data, prior, "mcmc",
      if (model == "CHLRM") initialization_kmeans() else initialization_ols(),
      mcmc_control(warmup = 1L, post_warmup_draws = 8L, chains = chains,
        stopping = if (chains > 1L) mcmc_stopping(min_post_warmup_draws = 8L,
          check_every = 8L, rhat_max = 100, ess_bulk_min = .1, ess_tail_min = .1) else NULL,
        max_memory_mb = 50), seed = 9081, verbose = FALSE)
    cache[[key]]
  }
})
.checkpoint_store_files <- function(store) list.files(file.path(store, "objects"), full.names = TRUE)
.checkpoint_store_resume <- function(fit) {
  # These tiny fixtures establish exact state/trajectory identity, not usable
  # diagnostic precision. Muffle only posterior's known small-sample ESS cap.
  withCallingHandlers(resume_mcmc(fit, additional_samples = 4L), warning = function(w) {
    if (identical(trimws(conditionMessage(w)), "The ESS has been capped to avoid unstable estimates."))
      invokeRestart("muffleWarning")
  })
}

test_that("stored checkpoints reconstruct complete plain state for all model families", {
  directory <- tempfile("checkpoint-store-"); dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE))
  for (model in c("LRM", "HLRM", "CHLRM")) for (chains in c(1L, 2L)) {
    fit <- .checkpoint_store_fit(model, chains)
    parent <- file.path(directory, paste(model, chains)); dir.create(parent)
    path <- file.path(parent, "block-01.rds"); store <- file.path(parent, "history")
    scope <- inference_rng_scope(9082)
    before_rng <- .Random.seed
    before_fit <- .typed_mcmc_checkpoint(fit)
    save_mcmc_checkpoint(fit, path, store = store, chunk_size = 4L)
    actual <- read_mcmc_checkpoint(path)
    expect_identical(.typed_mcmc_checkpoint(actual), before_fit)
    expect_identical(.typed_mcmc_checkpoint(fit), before_fit)
    expect_identical(.Random.seed, before_rng)
    expect_s3_class(readRDS(path), "vb_mcmc_store_checkpoint")
    expect_false(dir.exists(file.path(store, ".write-lock")))
    scope$restore()
  }
})

test_that("duplicate snapshots reuse immutable objects and the single-file default stays unchanged", {
  fit <- .checkpoint_store_fit(); directory <- tempfile(); dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE))
  store <- file.path(directory, "history")
  first <- file.path(directory, "first.rds"); second <- file.path(directory, "second.rds")
  save_mcmc_checkpoint(fit, first, store, chunk_size = 4L)
  files <- .checkpoint_store_files(store); hashes <- unname(tools::md5sum(files))
  save_mcmc_checkpoint(fit, second, store, chunk_size = 4L)
  expect_identical(.checkpoint_store_files(store), files)
  expect_identical(unname(tools::md5sum(files)), hashes)
  expect_identical(readRDS(first), readRDS(second))
  legacy <- file.path(directory, "legacy.rds")
  save_mcmc_checkpoint(fit, legacy)
  expect_identical(readRDS(legacy), .typed_mcmc_checkpoint(fit))
  expect_identical(.typed_mcmc_checkpoint(read_mcmc_checkpoint(legacy)),
    .typed_mcmc_checkpoint(read_mcmc_checkpoint(second)))
  expect_error(save_mcmc_checkpoint(fit, first, store, 4L), "already exists")
})

test_that("stored terminal sampler and RNG state resume the same chains", {
  directory <- tempfile(); dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE))
  for (model in c("LRM", "HLRM", "CHLRM")) {
    fit <- .checkpoint_store_fit(model)
    parent <- file.path(directory, model); dir.create(parent)
    path <- file.path(parent, "first.rds"); store <- file.path(parent, "history")
    save_mcmc_checkpoint(fit, path, store, 4L)
    files <- .checkpoint_store_files(store); hashes <- unname(tools::md5sum(files))
    a <- .checkpoint_store_resume(fit)
    b <- .checkpoint_store_resume(read_mcmc_checkpoint(path))
    expect_identical(a$parameters, b$parameters)
    expect_identical(a$latent_variables, b$latent_variables)
    expect_identical(a$objective, b$objective)
    ca <- .typed_mcmc_checkpoint(a)$chain_backends
    cb <- .typed_mcmc_checkpoint(b)$chain_backends
    for (i in seq_along(ca)) for (name in c("sampler", "chains", "rng_state", "run_settings"))
      expect_identical(ca[[i]][[name]], cb[[i]][[name]])
    save_mcmc_checkpoint(a, file.path(parent, "second.rds"), store, 4L)
    expect_identical(unname(tools::md5sum(files)), hashes)
    expect_identical(.typed_mcmc_checkpoint(read_mcmc_checkpoint(path)), .typed_mcmc_checkpoint(fit))
    # A standalone encoding of generation two needs more objects than the
    # incremental additions to the shared store (old full blocks are reused).
    before <- length(files)
    added <- length(.checkpoint_store_files(store)) - before
    fresh <- file.path(parent, "fresh"); dir.create(fresh)
    save_mcmc_checkpoint(a, file.path(fresh, "second.rds"), file.path(fresh, "history"), 4L)
    expect_lt(added, length(.checkpoint_store_files(file.path(fresh, "history"))))
  }
})

test_that("store writing avoids whole-state clones and preserves existing object identities", {
  directory <- tempfile(); dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE))
  for (model in c("LRM", "HLRM", "CHLRM")) for (chains in c(1L, 2L)) {
    fit <- .checkpoint_store_fit(model, chains)
    expected <- .typed_mcmc_checkpoint(fit)
    parent <- file.path(directory, paste(model, chains)); dir.create(parent)
    path <- file.path(parent, "first.rds"); store <- file.path(parent, "history")
    # A deep-clone attempt is a regression even when tiny fixtures fit in RAM.
    local({
      testthat::local_mocked_bindings(.clone_mcmc_state = function(x)
        stop("unexpected whole-state clone"), .package = "VariationalBayesHLRM")
      save_mcmc_checkpoint(fit, path, store, 4L)
    })
    manifest <- readRDS(path)
    files <- .checkpoint_store_files(store); hashes <- unname(tools::md5sum(files))
    legacy_root <- .checkpoint_store_encode(expected, store, .checkpoint_store_config(4L))
    expect_identical(legacy_root, manifest$root)
    expect_identical(.checkpoint_store_files(store), files)
    expect_identical(unname(tools::md5sum(files)), hashes)
    expect_identical(.typed_mcmc_checkpoint(read_mcmc_checkpoint(path)), expected)
  }
})

test_that("a published store snapshot is isolated from subsequent source mutations", {
  directory <- tempfile(); dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE))
  original <- file.path(directory, "legacy.rds")
  save_mcmc_checkpoint(.checkpoint_store_fit(), original)
  fit <- read_mcmc_checkpoint(original) # Do not mutate cached fixture backends.
  expected <- .typed_mcmc_checkpoint(fit)
  path <- file.path(directory, "first.rds"); store <- file.path(directory, "history")
  save_mcmc_checkpoint(fit, path, store, 4L)
  fit$metadata$after_publication <- TRUE
  attr(fit, "backend")$rng_state[2L] <- 42L
  attr(fit, "backends")[[1L]]$chains$Betak[1L] <- 12345
  expect_false(identical(.typed_mcmc_checkpoint(fit), expected))
  expect_identical(.typed_mcmc_checkpoint(read_mcmc_checkpoint(path)), expected)
  # Mutable references in backend values or atomic attributes are rejected
  # before any object is encoded, even without the legacy cloning operation.
  attr(attr(fit, "backends")[[1L]]$chains$Betak, "unsupported") <- new.env()
  rejected <- file.path(directory, "rejected"); dir.create(rejected)
  expect_error(save_mcmc_checkpoint(fit, file.path(rejected, "bad.rds"),
    file.path(rejected, "history"), 4L), "plain atomic/list")
  expect_length(.checkpoint_store_files(file.path(rejected, "history")), 0L)
  expect_false(file.exists(file.path(rejected, "bad.rds")))
  expect_false(dir.exists(file.path(rejected, "history", ".write-lock")))
})

test_that("chunking preserves allocation axes, attributes and exact prefix blocks", {
  directory <- tempfile(); dir.create(directory); dir.create(file.path(directory, "objects"))
  on.exit(unlink(directory, recursive = TRUE))
  config <- .checkpoint_store_config(4L)
  state <- structure(list(z = matrix(rep(1:3, length.out = 50 * 16), 50, 16),
    draws = array(seq_len(16 * 10 * 3), c(16, 10, 3)),
    frame = data.frame(x = 1:16, label = rep(c("a", "b"), 8)),
    empty = list(NULL, structure(numeric(), names = character()))), class = "fixture")
  key <- .checkpoint_store_encode(state, directory, config)
  result <- .checkpoint_store_decode(key, directory, new.env(), new.env())
  expect_identical(result, state)
  objects <- lapply(.checkpoint_store_files(directory), readRDS)
  atomic <- Filter(function(x) identical(x$kind, "atomic"), objects)
  allocation <- Filter(function(x) identical(x$dimensions, c(50L, 16L)), atomic)
  expect_length(allocation, 1L)
  expect_identical(allocation[[1L]]$axis, 2L)
  chunk <- .checkpoint_store_get(allocation[[1L]]$chunks[1L], directory, new.env())
  expect_identical(chunk$values, as.integer(state$z[, 1:4]))
  extended <- state; extended$z <- cbind(state$z, matrix(2L, 50, 4))
  key2 <- .checkpoint_store_encode(extended, directory, config)
  objects <- lapply(.checkpoint_store_files(directory), readRDS)
  atomic2 <- Filter(function(x) identical(x$kind, "atomic") && identical(x$dimensions, c(50L, 20L)), objects)
  expect_identical(atomic2[[1L]]$chunks[1:4], allocation[[1L]]$chunks)
  expect_identical(.checkpoint_store_decode(key2, directory, new.env(), new.env()), extended)
})

test_that("stores reject competing writers, format drift and mutable state", {
  fit <- .checkpoint_store_fit(); directory <- tempfile(); dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE))
  store <- file.path(directory, "history")
  save_mcmc_checkpoint(fit, file.path(directory, "first.rds"), store, 4L)
  lock <- file.path(store, ".write-lock"); dir.create(lock)
  expect_error(save_mcmc_checkpoint(fit, file.path(directory, "locked.rds"), store, 4L), "locked")
  expect_true(dir.exists(lock)); unlink(lock, recursive = TRUE)
  expect_error(save_mcmc_checkpoint(fit, file.path(directory, "drift.rds"), store, 5L), "configuration differs")
  expect_false(file.exists(file.path(directory, "drift.rds")))
  expect_false(dir.exists(lock))
  expect_error(save_mcmc_checkpoint(fit, file.path(directory, "invalid.rds"), store, 0), "chunk_size")
  fit$metadata$unsupported <- new.env()
  expect_error(save_mcmc_checkpoint(fit, file.path(directory, "mutable.rds"), store, 4L), "plain atomic/list")
  expect_false(file.exists(file.path(directory, "mutable.rds")))
})

test_that("relocation works and missing or corrupted objects fail explicitly", {
  directory <- tempfile(); dir.create(directory); moved <- paste0(directory, "-moved")
  on.exit(unlink(c(directory, moved), recursive = TRUE))
  fit <- .checkpoint_store_fit(); path <- file.path(directory, "first.rds")
  save_mcmc_checkpoint(fit, path, file.path(directory, "history"), 4L)
  expect_true(file.rename(directory, moved))
  path <- file.path(moved, "first.rds")
  expect_identical(.typed_mcmc_checkpoint(read_mcmc_checkpoint(path)), .typed_mcmc_checkpoint(fit))
  manifest <- readRDS(path)
  object <- file.path(moved, "history", "objects", paste0(manifest$root, ".rds"))
  bytes <- readBin(object, "raw", n = file.info(object)$size)
  cat("changed", file = object, append = TRUE)
  expect_error(read_mcmc_checkpoint(path), "checksum mismatch")
  writeBin(bytes, object)
  unlink(object)
  expect_error(read_mcmc_checkpoint(path), "Missing checkpoint store object")
})

test_that("an interrupted store write preserves old manifests and never publishes a partial checkpoint", {
  directory <- tempfile(); dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE))
  fit <- .checkpoint_store_fit(); store <- file.path(directory, "history")
  path <- file.path(directory, "first.rds")
  save_mcmc_checkpoint(fit, path, store, 4L)
  files <- .checkpoint_store_files(store); before <- unname(tools::md5sum(files))
  changed <- fit; changed$metadata$branch_marker <- "new-generation-interruption"
  put <- .checkpoint_store_put
  testthat::local_mocked_bindings(.checkpoint_store_put = function(node, store) {
    key <- put(node, store)
    if (identical(node$value, "new-generation-interruption")) stop("injected interruption")
    key
  }, .package = "VariationalBayesHLRM")
  expect_error(save_mcmc_checkpoint(changed, file.path(directory, "partial.rds"), store, 4L),
    "injected interruption")
  expect_false(file.exists(file.path(directory, "partial.rds")))
  expect_false(dir.exists(file.path(store, ".write-lock")))
  expect_identical(unname(tools::md5sum(files)), before)
  expect_gt(length(.checkpoint_store_files(store)), length(files))
  expect_identical(.typed_mcmc_checkpoint(read_mcmc_checkpoint(path)), .typed_mcmc_checkpoint(fit))
})
