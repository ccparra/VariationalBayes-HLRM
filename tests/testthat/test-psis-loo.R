test_that("PSIS-LOO uses converged MCMC draws and preserves its estimand", {
  skip_if_not_installed("loo")
  scope <- inference_rng_scope(1500); on.exit(scope$restore())
  X <- cbind(intercept = 1, x = stats::rnorm(20))
  data <- list(X = X, y = drop(X %*% c(2, -1) + stats::rnorm(20)))
  prior <- lrm_prior(c(0, 0), diag(10, 2), 4, 2)
  stopping <- mcmc_stopping(
    min_post_warmup_draws = 50, check_every = 50,
    rhat_max = 100, ess_bulk_min = .1, ess_tail_min = .1,
    relative_mcse_mean_max = 100, relative_mcse_sd_max = 100)
  fit <- fit_model(
    lrm_model(2), data, prior, "mcmc", initialization_ols(),
    mcmc_control(
      warmup = 10, post_warmup_draws = 100, chains = 2,
      evaluation_draws = 80, stopping = stopping, max_memory_mb = 50),
    scaling = data_scaling(), seed = 1501, verbose = FALSE)
  expect_identical(fit$diagnostics$status, "converged")
  result <- suppressWarnings(psis_loo(fit, n_draws = 80, seed = 1502))
  expect_s3_class(result, "vb_psis_loo")
  expect_true(loo::is.loo(result$loo_object))
  expect_named(result$estimates,
               c("quantity", "Estimate", "SE"))
  expect_identical(result$unit, "observation")
  expect_identical(result$conditioning, "parameter_conditional")
  expect_identical(result$metadata$density_scale, "input")
  expect_equal(length(result$diagnostics$pareto_k), 20)
  expect_true(result$metadata$draws_per_chain_for_relative_efficiency >= 4)
  expect_error(psis_loo(fit, n_draws = 6, seed = 1503),
               "at least four selected draws from every chain")

  variational <- fit_model(
    lrm_model(2), data, prior, "cavi", initialization_ols(),
    cavi_control(3, evaluation_draws = 4), seed = 1504, verbose = FALSE)
  expect_error(psis_loo(variational), "requires MCMC")
})
