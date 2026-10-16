test_that("standard grouped PPC implements T1 through T5 exactly", {
  y <- list(c(1, 3, 5), c(2, 4, 8))
  X <- list(cbind(1, c(-1, 0, 1)), cbind(1, c(-1, 0, 1)))
  replicated <- list(c(2, 3, 7), c(1, 5, 9))
  check <- pp_check_stream(
    y, 2, function(s) replicated, "test", seed = 700,
    X_list = X, discrepancies = ppc_discrepancies())

  flat <- unlist(y, use.names = FALSE)
  fitted <- unlist(Map(function(response, design)
    qr.fitted(qr(design, tol = 1e-10), response), y, X), use.names = FALSE)
  expected <- c(
    mean = mean(flat),
    sd = sd(flat),
    between_group_mean_variance = var(vapply(y, mean, numeric(1))),
    mean_within_group_variance = mean(vapply(y, var, numeric(1))),
    ols_mse = mean((flat - fitted)^2))
  expect_equal(check$t0, expected)
  expect_identical(names(check$repdata), names(expected))
  expect_equal(check$ppp,
    colMeans(sweep(as.matrix(check$repdata), 2, expected, `>=`)))
})

test_that("PPC discrepancy applicability is explicit", {
  expect_identical(ppc_discrepancies(c("T1", "T5"))$statistics,
                   c("mean", "ols_mse"))
  expect_error(ppc_discrepancies("made_up"), "Unknown")
  expect_error(pp_check_stream(
    list(1:3), 2, function(s) list(1:3), "test",
    discrepancies = ppc_discrepancies("T3")), "at least two groups")
  expect_error(pp_check_stream(
    list(1, 2:3), 2, function(s) list(1, 2:3), "test",
    discrepancies = ppc_discrepancies("T4")), "at least two observations")
})

test_that("public PPC reports dependence-aware or explicitly unavailable MCSE", {
  X <- cbind(1, seq(-1, 1, length.out = 14))
  y <- drop(X %*% c(1, 2)) + sin(seq_len(14)) / 10
  fit <- fit_model(
    lrm_model(2), list(X = X, y = y),
    lrm_prior(c(0, 0), diag(2), 4, 1), "mcmc",
    initialization_ols(),
    mcmc_control(post_warmup_draws = 8, warmup = 2,
                 evaluation_draws = 8),
    seed = 701, verbose = FALSE)
  check <- posterior_predictive_check(fit, n_draws = 8, seed = 702)
  expect_identical(check$discrepancies$mcse_method,
                   "indicator_sequence_ess_mean")
  expect_true(all(is.finite(check$discrepancies$effective_sample_size)))
  expect_true(all(is.finite(check$discrepancies$mcse)))

  short <- posterior_predictive_check(fit, n_draws = 2, seed = 703)
  expect_match(short$discrepancies$mcse_method, "unavailable")
  expect_true(all(is.na(short$discrepancies$mcse)))
})
