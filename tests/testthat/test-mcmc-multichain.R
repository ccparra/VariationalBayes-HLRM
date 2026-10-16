test_that("multichain LRM is reproducible, diagnosed and guarded", {
  d <- .test_lrm_simulation(25, c(1, -0.5), 1, matrix(c(-1, 1), 1, 2), seed = 701)
  h <- hyperparams_lrm(d$train, "OLS")
  settings <- list(samples = 20, burn_in = 5, thin = 1)
  args <- list(model = "LRM", y = d$train$y, X = d$train$X, hyperparams = h,
               run_settings = settings, n_chains = 2, seed = 702,
               rhat_max = 10, ess_bulk_min = 0.1, ess_tail_min = 0.1)
  a <- suppressWarnings(do.call(run_mcmc_chains, args))
  b <- suppressWarnings(do.call(run_mcmc_chains, args))

  expect_s3_class(a, "mcmc_multichain_fit")
  expect_identical(a$chains, b$chains)
  expect_identical(a$chain_metadata$seed, b$chain_metadata$seed)
  expect_true(all(c("rhat", "ess_bulk", "ess_tail", "mcse_mean", "mcse_sd",
                    "relative_mcse_mean", "relative_mcse_sd",
                    "ess_per_second", "passed") %in%
                    names(a$diagnostics$variables)))
  expect_match(a$diagnostics$implementation, "posterior")
  expect_true(a$diagnostics$passed)

  combined <- combine_mcmc_chains(a)
  expect_s3_class(combined, "mcmc_combined_draws")
  expect_equal(dim(combined$Beta), c(40, 2))
  expect_equal(nrow(attr(combined, "draw_origin")), 40)
  expect_equal(unique(attr(combined, "draw_origin")$chain), 1:2)

  strict <- a
  strict$diagnostics$passed <- FALSE
  expect_error(combine_mcmc_chains(strict), "did not satisfy")
  expect_equal(nrow(combine_mcmc_chains(strict, FALSE)$Beta), 40)
})

test_that("multichain rejects thinning, duplicate seeds and unsafe storage", {
  d <- .test_lrm_simulation(12, c(0, 1), 1, matrix(c(-1, 1), 1, 2), seed = 710)
  h <- hyperparams_lrm(d$train, "OLS")
  base <- list(model = "LRM", y = d$train$y, X = d$train$X,
               hyperparams = h,
               run_settings = list(samples = 8, burn_in = 2, thin = 1),
               n_chains = 2, seed = 711)
  thin <- base; thin$run_settings$thin <- 2
  expect_error(do.call(run_mcmc_chains, thin), "thin = 1")
  duplicate <- base; duplicate$seed <- NULL; duplicate$seeds <- c(1, 1)
  expect_error(do.call(run_mcmc_chains, duplicate), "distinct")
  memory <- base; memory$max_memory_mb <- 1e-12
  expect_error(do.call(run_mcmc_chains, memory), "exceeds")
})

test_that("HLRM and CHLRM share the same multichain contract", {
  dh <- .test_hlrm_hierarchy(rep(5, 3), c(0, 0.5), diag(2), 1,
                                seed = 720)
  hh <- hyperparams_hlrm(dh$train, "OLS")
  settings <- list(samples = 12, burn_in = 3, thin = 1)
  fh <- suppressWarnings(run_mcmc_chains(
    "HLRM", dh$train$y, dh$train$X, hh, settings,
    n_chains = 2, seed = 721,
    rhat_max = 100, ess_bulk_min = 0.1, ess_tail_min = 0.1))
  expect_true(fh$diagnostics$passed)
  expect_null(fh$label_reference)
  expect_equal(dim(combine_mcmc_chains(fh)$Betaj), c(24, 3, 2))

  dc <- .test_chlrm_simulation(4, 5, 2, 5, seed = 722)
  hc <- .test_chlrm_native_prior(5, 2)
  fc <- suppressWarnings(run_mcmc_chains(
    "CHLRM", dc$train$y, dc$train$X, hc, settings, K = 2,
    n_chains = 2, seed = 723, rhat_max = 100,
    ess_bulk_min = 0.1, ess_tail_min = 0.1))
  expect_true(fc$diagnostics$passed)
  expect_s3_class(fc$label_reference, "chlrm_label_reference")
  expect_length(fc$aligned_chains, 2)
  expect_true(all(vapply(fc$aligned_chains, function(x)
    identical(attr(x, "label_reference"), fc$label_reference), logical(1))))
  combined <- combine_mcmc_chains(fc)
  expect_equal(dim(combined$Betak), c(24, 2, 5))
  expect_identical(attr(combined, "label_reference"), fc$label_reference)
})
