test_that("posterior draws convert to standard draws_df without losing values", {
  X <- cbind(intercept = 1, predictor = seq(-1, 1, length.out = 16))
  y <- drop(X %*% c(1, 2)) + sin(seq_len(16)) / 10
  fit <- fit_model(
    lrm_model(2), list(X = X, y = y),
    lrm_prior(c(0, 0), 10 * diag(2), 4, 1), "cavi",
    initialization_ols(), cavi_control(4, evaluation_draws = 5),
    seed = 730, verbose = FALSE)
  semantic <- posterior_draws(fit, n_draws = 5, seed = 731)
  converted <- posterior::as_draws(semantic)
  expect_s3_class(converted, "draws_df")
  expect_true(all(c("coefficients[1]", "coefficients[2]",
                    "residual_variance") %in%
                  posterior::variables(converted)))
  expect_equal(converted$`coefficients[1]`,
               semantic$parameters$coefficients[, 1])
  direct <- posterior::as_draws(fit, n_draws = 5, seed = 731)
  expect_equal(as.data.frame(direct), as.data.frame(converted))
})

test_that("draws conversion preserves MCMC chain and iteration provenance", {
  X <- cbind(1, seq(-1, 1, length.out = 12))
  y <- drop(X %*% c(.5, -1)) + cos(seq_len(12)) / 10
  fit <- fit_model(
    lrm_model(2), list(X = X, y = y),
    lrm_prior(c(0, 0), diag(2), 4, 1), "mcmc",
    initialization_ols(),
    mcmc_control(
      warmup = 1, post_warmup_draws = 6, chains = 2,
      evaluation_draws = 4,
      stopping = mcmc_stopping(
        min_post_warmup_draws = 4, check_every = 3,
        rhat_max = 100, ess_bulk_min = .1, ess_tail_min = .1),
      max_memory_mb = 50),
    seed = 732, verbose = FALSE)
  semantic <- posterior_draws(fit, n_draws = 6, seed = 733)
  converted <- posterior::as_draws(semantic)
  frame <- as.data.frame(converted)
  expect_equal(frame$.chain, semantic$draw_index$chain)
  expect_equal(frame$source_iteration, semantic$draw_index$iteration)
  expect_equal(frame$source_stored_draw, semantic$draw_index$stored_draw)
})
