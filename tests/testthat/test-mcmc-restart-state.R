.restart_fixture <- local({
  cache <- new.env(parent = emptyenv())
  function(model = 'CHLRM', chains = 2L, scaling = FALSE) {
    key <- paste(model, chains, scaling, paste(RNGkind(), collapse = '-'))
    if (exists(key, cache, inherits = FALSE)) return(cache[[key]])
    X <- cbind(1, seq(-1, 1, length.out = 8))
    data <- list(X = rep(list(X), 6), y = lapply(seq_len(6), function(j)
      drop(X %*% c(j / 3, (-1)^j)) + sin(seq_len(8)) / 10))
    spec <- switch(model, LRM = lrm_model(2), HLRM = hlrm_model(2), CHLRM = chlrm_model(2, 2))
    prior <- switch(model, LRM = lrm_prior(c(0, 0), diag(2), 4, 1),
      HLRM = hlrm_prior(c(0, 0), diag(2), 4, diag(2), 4, 2, 2),
      CHLRM = chlrm_prior(c(0, 0), diag(2), 4, diag(2), 4, 2, 2, c(1, 1)))
    if (model == 'LRM') data <- list(X = X, y = data$y[[1L]])
    cache[[key]] <- fit_model(spec, data, prior, 'mcmc', initialization_from_prior(),
      mcmc_control(warmup = 2L, post_warmup_draws = 8L, chains = chains,
        stopping = if (chains > 1L) mcmc_stopping(min_post_warmup_draws = 8L, check_every = 8L,
          rhat_max = 100, ess_bulk_min = .1, ess_tail_min = .1,
          relative_mcse_mean_max = 1e10, relative_mcse_sd_max = 1e10) else NULL,
        max_memory_mb = 50), scaling = if (scaling) data_scaling() else NULL,
      seed = 70921, verbose = FALSE)
    cache[[key]]
  }
})
.restart_resume <- function(fit, n) withCallingHandlers(resume_mcmc(fit, n), warning = function(w) {
  if (identical(trimws(conditionMessage(w)), 'The ESS has been capped to avoid unstable estimates.')) invokeRestart('muffleWarning')
})
.restart_native_suffix <- function(draws, model, first) {
  lapply(names(draws), function(name) {
    x <- draws[[name]]; axis <- if (model == 'CHLRM' && name == 'z') 2L else 1L
    count <- if (is.null(dim(x))) length(x) else dim(x)[axis]
    .checkpoint_axis_subset(x, seq.int(first, count), axis)
  }) |> setNames(names(draws))
}

test_that('history-free blocks preserve every native field, sampler and RNG for all families', {
  for (model in c('LRM', 'HLRM', 'CHLRM')) for (chains in c(1L, 2L, 4L)) {
    fit <- .restart_fixture(model, chains)
    before <- .typed_mcmc_checkpoint(fit)
    state <- mcmc_restart_state(fit); original <- serialize(state, NULL)
    full <- .restart_resume(fit, 5L)
    backends <- attr(full, 'backends', exact = TRUE)
    if (is.null(backends)) backends <- list(.fit_backend(full))
    for (i in rev(seq_len(chains))) {
      first <- advance_mcmc_state(state, 3L, chain = i)
      second <- advance_mcmc_state(first$state, 2L, chain = i)
      expect_identical(.combine_chain_draws(list(first$draws, second$draws), model),
        .restart_native_suffix(backends[[i]]$chains, model, 9L))
      expect_identical(second$state$chains[[i]]$sampler, backends[[i]][[.mcmc_object_name(model)]])
      expect_identical(second$state$chains[[i]]$rng_state, backends[[i]]$rng_state)
      expect_identical(c(first$draw_index$iteration, second$draw_index$iteration), 9:13)
      expect_identical(first$draw_index$chain, rep.int(as.integer(i), 3L))
      expect_equal(second$state$chains[[i]]$run_settings$samples, 13L)
      expect_false(second$resources$cumulative_history_allocated)
      for (j in setdiff(seq_len(chains), i)) expect_identical(second$state$chains[[j]], state$chains[[j]])
    }
    expect_identical(serialize(state, NULL), original)
    expect_identical(.typed_mcmc_checkpoint(fit), before)
    expect_null(state$chains[[1]]$chains)
    expect_false(inherits(state, 'vb_model_fit'))
  }
})

test_that('working-scale blocks carry scaling context and preserve caller RNG including absent seed', {
  fit <- .restart_fixture('LRM', 1L, TRUE)
  state <- mcmc_restart_state(fit)
  expect_identical(state$context$scaling_map, attr(fit, 'scaling_map', exact = TRUE))
  ambient <- rng_snapshot(); on.exit(restore_rng_snapshot(ambient))
  for (kind in c('Mersenne-Twister', "L'Ecuyer-CMRG")) {
    RNGkind(kind, 'Inversion', 'Rejection'); set.seed(99)
    before <- rng_snapshot(); block <- advance_mcmc_state(state, 3L)
    expect_identical(rng_snapshot(), before)
    full <- .restart_resume(fit, 3L)
    expect_identical(block$draws, .restart_native_suffix(.fit_backend(full)$chains, 'LRM', 9L))
  }
  rm('.Random.seed', envir = .GlobalEnv)
  invisible(advance_mcmc_state(state, 1L))
  expect_false(exists('.Random.seed', .GlobalEnv, inherits = FALSE))
  RNGkind("L'Ecuyer-CMRG", 'Inversion', 'Rejection')
  # A genuine typed fit records the L'Ecuyer stream without a fabricated context.
  other <- .restart_fixture('LRM', 1L)
  genuine <- mcmc_restart_state(other)
  expect_identical(genuine$chains[[1]]$rng_state$kind[1], "L'Ecuyer-CMRG")
  before <- rng_snapshot(); result <- advance_mcmc_state(genuine, 2L)
  full <- .restart_resume(other, 2L)
  expect_identical(result$state$chains[[1]]$rng_state, .fit_backend(full)$rng_state)
  expect_identical(rng_snapshot(), before)
})

test_that('guards, malformed state and kernel errors preserve the source and RNG', {
  state <- mcmc_restart_state(.restart_fixture('CHLRM', 2L))
  before <- serialize(state, NULL); rng <- rng_snapshot()
  expect_error(advance_mcmc_state(state, 1L), 'Specify chain')
  expect_error(advance_mcmc_state(state, 1L, 3L), 'Unknown restart chain')
  expect_error(advance_mcmc_state(state, 1L, 1L, max_memory_mb = .00001), 'exceeds max_(state|memory)_mb')
  expect_error(advance_mcmc_state(state, 0L, 1L), 'additional_samples')
  bad <- state; bad$chains[[1]]$run_settings$thin <- 2L
  expect_error(advance_mcmc_state(bad, 1L, 1L), 'thin = 1')
  bad <- state; bad$chains[[1]]$rng_state$seed <- 1:2
  expect_error(advance_mcmc_state(bad, 1L, 1L), 'RNG encoding')
  bad <- state; bad$chains[[1]]$sampler$params$z[1] <- 999L
  expect_error(advance_mcmc_state(bad, 1L, 1L), 'valid component')
  bad <- state; bad$chains[[1]]$sampler$params$Sigma[1, 1] <- NA_real_
  expect_error(advance_mcmc_state(bad, 1L, 1L), 'parameter dimensions')
  bad <- state; bad$chains[[1]]$chains <- list(history = numeric(10))
  expect_error(advance_mcmc_state(bad, 1L, 1L), 'sampler layout')
  bad <- state; bad$chains[[1]]$run_settings$samples <- .Machine$integer.max
  expect_error(advance_mcmc_state(bad, 1L, 1L), 'integer range')
  local_mocked_bindings(.sample_mcmc_block = function(...) {runif(1); stop('deliberate kernel failure')})
  expect_error(advance_mcmc_state(state, 1L, 1L), 'deliberate kernel')
  expect_identical(serialize(state, NULL), before)
  expect_identical(rng_snapshot(), rng)
})

test_that('dedicated restart files are exact, atomic, checked and relocatable', {
  state <- mcmc_restart_state(.restart_fixture())
  directory <- tempfile(); dir.create(directory); on.exit(unlink(directory, recursive = TRUE))
  path <- file.path(directory, 'one.state')
  save_mcmc_restart_state(state, path)
  expect_identical(read_mcmc_restart_state(path), state)
  expect_false(dir.exists(paste0(path, '.write-lock')))
  expect_error(save_mcmc_restart_state(state, path), 'already exists')
  copied <- file.path(directory, 'two.state'); file.copy(path, copied)
  expect_identical(read_mcmc_restart_state(copied), state)
  bytes <- readBin(path, 'raw', n = file.info(path)$size)
  bytes[length(bytes)] <- as.raw(bitwXor(as.integer(bytes[length(bytes)]), 1L))
  writeBin(bytes, copied)
  expect_error(read_mcmc_restart_state(copied), 'checksum mismatch')
  expect_error(read_mcmc_restart_state(path, max_state_mb = .00001), 'exceeds max_state_mb')
  legacy <- file.path(directory, 'legacy.rds'); save_mcmc_checkpoint(.restart_fixture(), legacy)
  expect_error(read_mcmc_restart_state(legacy), 'Unsupported restart format')
  lock <- paste0(file.path(directory, 'locked.state'), '.write-lock'); dir.create(lock)
  expect_error(save_mcmc_restart_state(state, file.path(directory, 'locked.state')), 'locked')
  expect_true(dir.exists(lock))
})

test_that('store restart reading never decodes native histories and validates selected nodes', {
  directory <- tempfile(); dir.create(directory); on.exit(unlink(directory, recursive = TRUE))
  for (model in c('LRM', 'HLRM', 'CHLRM')) for (n in c(1L, 2L)) {
    fit <- .restart_fixture(model, n, model == 'LRM')
    parent <- file.path(directory, paste(model, n)); dir.create(parent)
    path <- file.path(parent, 'original.rds'); store <- file.path(parent, 'history')
    save_mcmc_checkpoint(fit, path, store, 4L)
    files <- list.files(store, full.names = TRUE, recursive = TRUE)
    before <- unname(tools::md5sum(files)); rng <- rng_snapshot()
    expect_identical(read_mcmc_restart_state(path, 'checkpoint_store'), mcmc_restart_state(fit))
    expect_identical(unname(tools::md5sum(files)), before)
    expect_identical(rng_snapshot(), rng)
    expect_error(read_mcmc_restart_state(path, 'checkpoint_store', max_state_mb = .00001), 'exceed')
  }
  # Remove only a native historical value chunk. State recovery must still
  # succeed, whereas complete-fit restoration must report corruption/missing data.
  path <- file.path(directory, 'CHLRM 2/original.rds'); manifest <- readRDS(path)
  store <- file.path(dirname(path), manifest$store)
  objects <- list.files(file.path(store, 'objects'), full.names = TRUE)
  candidates <- objects[vapply(objects, function(f) {
    x <- readRDS(f)
    identical(x$kind, 'values') && is.double(x$values) && length(x$values) == 4L * 6L * 2L
  }, logical(1L))]
  expect_gt(length(candidates), 0L)
  original <- read_mcmc_restart_state(path, 'checkpoint_store')
  unlink(candidates[1L])
  expect_identical(read_mcmc_restart_state(path, 'checkpoint_store'), original)
  expect_error(read_mcmc_checkpoint(path), 'Missing checkpoint store object')
  # Current-state corruption is fatal, with no source mutation or sampling.
  objects <- list.files(file.path(store, 'objects'), full.names = TRUE)
  current <- objects[vapply(objects, function(f) {
    x <- readRDS(f); identical(x$kind, 'literal') && is.integer(x$value) && length(x$value) > 500L
  }, logical(1L))]
  if (!length(current)) {
    # RNG seed can instead be encoded as an atomic integer chunk.
    current <- objects[vapply(objects, function(f) {
      x <- readRDS(f); identical(x$kind, 'values') && is.integer(x$values) && length(x$values) > 100L
    }, logical(1L))]
  }
  expect_gt(length(current), 0L)
  writeBin(as.raw(1:4), current[1L])
  expect_error(read_mcmc_restart_state(path, 'checkpoint_store'), 'checksum mismatch')
})

test_that('restart primitives reject reporting views and unsupported cached normal RNGs', {
  fit <- .restart_fixture()
  view <- mcmc_reporting_view(fit)
  expect_error(mcmc_restart_state(view), 'reporting view')
  class(view) <- class(fit)
  expect_error(mcmc_restart_state(view), 'complete MCMC fit')
  state <- mcmc_restart_state(fit)
  original <- rng_snapshot(); on.exit(restore_rng_snapshot(original))
  RNGkind(normal.kind = 'Box-Muller')
  expect_error(advance_mcmc_state(state, 1L, 1L), 'caller.*Inversion')
})


test_that('failed restart publication preserves previous versions and releases owned files', {
  state <- mcmc_restart_state(.restart_fixture('LRM', 1L))
  directory <- tempfile(); dir.create(directory); on.exit(unlink(directory, recursive = TRUE))
  first <- file.path(directory, 'first.state'); next_path <- file.path(directory, 'next.state')
  save_mcmc_restart_state(state, first)
  before <- unname(tools::md5sum(first)); files <- list.files(directory, all.files = TRUE)
  local_mocked_bindings(.publish_mcmc_restart_file = function(...) stop('deliberate publication failure'))
  expect_error(save_mcmc_restart_state(state, next_path), 'deliberate publication')
  expect_identical(unname(tools::md5sum(first)), before)
  expect_identical(read_mcmc_restart_state(first), state)
  expect_identical(list.files(directory, all.files = TRUE), files)
  expect_false(file.exists(next_path))
})
