test_that("LRM checkpoint continuation equals one uninterrupted chain", {
  d <- .test_lrm_simulation(20, c(1, -0.5), 1, matrix(c(-1, 1), 1, 2), seed = 901)
  h <- hyperparams_lrm(d$train, "OLS")
  short <- mcmc_lrm$new(d$train$y, d$train$X, h,
    list(samples = 4, burn_in = 2, thin = 1, seed = 902))
  short$fit(FALSE)
  checkpoint <- .mcmc_backend_checkpoint(short)
  .resume_mcmc_backend(short, 3)

  full <- mcmc_lrm$new(d$train$y, d$train$X, h,
    list(samples = 7, burn_in = 2, thin = 1, seed = 902))
  full$fit(FALSE)
  expect_identical(short$chains, full$chains)
  expect_equal(short$fit_metadata$stored_samples, 7)
  expect_length(short$fit_metadata$extensions, 1)

  restored <- .restore_mcmc_backend_checkpoint(checkpoint)
  .resume_mcmc_backend(restored, 3)
  expect_identical(restored$chains, full$chains)
  expect_identical(restored[["lrm_mcmc_object"]]$params,
                   full[["lrm_mcmc_object"]]$params)
})

test_that("HLRM and CHLRM continuation preserve their full sampler state", {
  dh <- .test_hlrm_hierarchy(rep(4, 3), c(0, 1), diag(2), 1,
                                seed = 910)
  hh <- hyperparams_hlrm(dh$train, "OLS")
  a <- mcmc_hlrm$new(dh$train$y, dh$train$X, hh,
    list(samples = 3, burn_in = 2, thin = 1, seed = 911))
  a$fit(FALSE); .resume_mcmc_backend(a, 2)
  b <- mcmc_hlrm$new(dh$train$y, dh$train$X, hh,
    list(samples = 5, burn_in = 2, thin = 1, seed = 911))
  b$fit(FALSE)
  expect_identical(a$chains, b$chains)
  expect_identical(a$hlrm_mcmc_object$params, b$hlrm_mcmc_object$params)

  dc <- .test_chlrm_simulation(4, 5, 2, 5, seed = 912)
  hc <- .test_chlrm_native_prior(5, 2)
  x <- mcmc_chlrm$new(dc$train$y, dc$train$X, 2, hc,
    list(samples = 3, burn_in = 2, thin = 1, seed = 913))
  x$fit(FALSE); x$summarize_posterior()
  .resume_mcmc_backend(x, 2)
  expect_length(x$relabeled_chains, 0)
  expect_null(x$label_reference)
  z <- mcmc_chlrm$new(dc$train$y, dc$train$X, 2, hc,
    list(samples = 5, burn_in = 2, thin = 1, seed = 913))
  z$fit(FALSE)
  expect_identical(x$chains, z$chains)
  expect_identical(x$chlrm_mcmc_object$params, z$chlrm_mcmc_object$params)
})

test_that("typed checkpoint files are atomic plain-state round trips", {
  d <- .test_lrm_simulation(12, c(0, 1), 1, matrix(c(-1, 1), 1, 2), seed = 920)
  prior <- lrm_prior(c(0, 0), 10 * diag(2), 4, 1)
  fit <- fit_model(
    lrm_model(2), d$train, prior, "mcmc",
    initialization_from_prior(),
    mcmc_control(warmup = 1, post_warmup_draws = 3,
                 evaluation_draws = 2),
    seed = 921, verbose = FALSE)
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  saved <- save_mcmc_checkpoint(fit, path)
  expect_true(file.exists(saved))
  restored <- read_mcmc_checkpoint(path)
  expect_s3_class(restored, "vb_model_fit")
  expect_identical(restored$parameters, fit$parameters)
  expect_identical(restored$metadata, fit$metadata)

  continued <- resume_mcmc(restored, 3)
  expect_equal(continued$specifications$control$post_warmup_draws, 6)
  expect_equal(dim(posterior_draws(continued)$parameters$coefficients), c(6, 2))
  expect_equal(dim(posterior_draws(fit)$parameters$coefficients), c(3, 2))
})

test_that("typed multichain continuation preserves every individual chain", {
  d <- .test_lrm_simulation(
    18, c(.5, -1), 1, matrix(c(-1, 1), 1, 2), seed = 925)
  stopping <- mcmc_stopping(
    rhat_max = 100, ess_bulk_min = .1, ess_tail_min = .1,
    relative_mcse_mean_max = 100, relative_mcse_sd_max = 100,
    min_post_warmup_draws = 2, check_every = 2)
  fit <- suppressWarnings(fit_model(
    lrm_model(2), d$train,
    lrm_prior(c(0, 0), 10 * diag(2), 4, 1), "mcmc",
    initialization_from_prior(),
    mcmc_control(
      warmup = 1, post_warmup_draws = 4, chains = 2,
      evaluation_draws = 2, stopping = stopping,
      max_memory_mb = 50),
    seed = 926, verbose = FALSE))

  raw <- attr(fit, "backends", exact = TRUE)
  expect_equal(vapply(raw, function(x) nrow(x$chains$Beta), integer(1L)),
               c(4L, 4L))
  expect_equal(nrow(fit$parameters$coefficients$draws), 8)

  continued <- suppressWarnings(resume_mcmc(fit, 2))
  continued_raw <- attr(continued, "backends", exact = TRUE)
  expect_equal(vapply(continued_raw, function(x) nrow(x$chains$Beta),
                      integer(1L)), c(6L, 6L))
  expect_equal(nrow(continued$parameters$coefficients$draws), 12)
  expect_equal(nrow(fit$parameters$coefficients$draws), 8)
  expect_equal(tail(fit_orchestration(continued)$block_history$stored_samples,
                    1L), 6)

  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  save_mcmc_checkpoint(continued, path)
  restored <- read_mcmc_checkpoint(path)
  expect_equal(nrow(restored$parameters$coefficients$draws), 12)
  expect_equal(vapply(attr(restored, "backends", exact = TRUE),
                      function(x) nrow(x$chains$Beta), integer(1L)),
               c(6L, 6L))
})

test_that("typed continuation covers HLRM and CHLRM quantities", {
  hierarchical <- .test_hlrm_hierarchy(
    rep(4, 3), c(0, 1), diag(2), 1, seed = 927)
  hlrm_fit <- fit_model(
    hlrm_model(2), hierarchical$train,
    hlrm_prior(c(0, 0), 10 * diag(2), 4, diag(2), 4, 2, 2),
    "mcmc", initialization_from_prior(),
    mcmc_control(warmup = 1, post_warmup_draws = 3,
                 evaluation_draws = 2),
    seed = 928, verbose = FALSE)
  hlrm_extended <- resume_mcmc(hlrm_fit, 2)
  expect_equal(dim(hlrm_extended$parameters$group_coefficients$draws),
               c(5, 3, 2))
  expect_equal(dim(hlrm_fit$parameters$group_coefficients$draws),
               c(3, 3, 2))

  clustered <- .test_chlrm_simulation(4, 5, 2, 2, seed = 929)
  chlrm_fit <- fit_model(
    chlrm_model(2, 2), clustered$train,
    chlrm_prior(c(0, 0), 10 * diag(2), 4, diag(2), 4, 2, 2,
                c(1, 1)),
    "mcmc", initialization_from_prior(),
    mcmc_control(warmup = 1, post_warmup_draws = 3,
                 evaluation_draws = 2),
    seed = 930, verbose = FALSE)
  chlrm_extended <- resume_mcmc(chlrm_fit, 2)
  expect_equal(dim(chlrm_extended$parameters$component_coefficients$draws),
               c(5, 2, 2))
  expect_equal(dim(chlrm_extended$latent_variables$group_allocations),
               c(4, 5))
})

test_that("typed continuation preserves the original scale transformation", {
  d <- .test_lrm_simulation(
    20, c(10, -2), 4, matrix(c(100, 200), 1, 2), seed = 932)
  prior <- lrm_prior(c(0, 0), 100 * diag(2), 4, 4)
  short <- fit_model(
    lrm_model(2), d$train, prior, "mcmc", initialization_from_prior(),
    mcmc_control(warmup = 2, post_warmup_draws = 3,
                 evaluation_draws = 2),
    scaling = data_scaling(), seed = 933, verbose = FALSE)
  continued <- resume_mcmc(short, 2)
  full <- fit_model(
    lrm_model(2), d$train, prior, "mcmc", initialization_from_prior(),
    mcmc_control(warmup = 2, post_warmup_draws = 5,
                 evaluation_draws = 2),
    scaling = data_scaling(), seed = 933, verbose = FALSE)

  expect_identical(continued$parameters, full$parameters)
  expect_s3_class(attr(continued, "scaling_map", exact = TRUE),
                  "vb_data_scaling_map")
})

test_that("resume validates state and never resumes thinned chains", {
  d <- .test_lrm_simulation(12, c(0, 1), 1, matrix(c(-1, 1), 1, 2), seed = 930)
  h <- hyperparams_lrm(d$train, "OLS")
  fit <- mcmc_lrm$new(d$train$y, d$train$X, h,
    list(samples = 3, burn_in = 1, thin = 2, seed = 931))
  fit$fit(FALSE)
  expect_error(.resume_mcmc_backend(fit, 2), "thin = 1")
  expect_error(.resume_mcmc_backend(fit, 0), "additional_samples")
  expect_error(.restore_mcmc_backend_checkpoint(list()), "Invalid")
})

test_that("incremental multichain records every diagnostic block and stops safely", {
  d <- .test_lrm_simulation(25, c(1, 0.5), 1, matrix(c(-1, 1), 1, 2), seed = 940)
  h <- hyperparams_lrm(d$train, "OLS")
  fit <- suppressWarnings(run_mcmc_chains(
    "LRM", d$train$y, d$train$X, h,
    list(samples = 99, burn_in = 3, thin = 1),
    n_chains = 2, seed = 941, block_size = 8, max_samples = 16,
    rhat_max = 1.01, ess_bulk_min = 1e6, ess_tail_min = 1e6))
  expect_true(fit$execution$incremental)
  expect_identical(fit$execution$stop_reason, "max_samples")
  expect_equal(fit$block_history$stored_samples, c(8, 16))
  expect_equal(fit$run_settings$samples, 16)
  expect_error(combine_mcmc_chains(fit), "did not satisfy")
  expect_equal(nrow(combine_mcmc_chains(fit, FALSE)$Beta), 32)

  early <- suppressWarnings(run_mcmc_chains(
    "LRM", d$train$y, d$train$X, h,
    list(samples = 99, burn_in = 3, thin = 1),
    n_chains = 2, seed = 942, block_size = 8, max_samples = 24,
    rhat_max = 100, ess_bulk_min = 0.1, ess_tail_min = 0.1))
  expect_identical(early$execution$stop_reason, "diagnostic_thresholds")
  expect_equal(nrow(early$block_history), 1)
})
