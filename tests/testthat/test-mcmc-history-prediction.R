history_prediction_fixture <- function(directory) {
  X <- cbind(1, seq(-1, 1, length.out = 8))
  data <- list(X = rep(list(X), 4L),
    y = lapply(1:4, function(j) seq(-1, 1, length.out = 8) + j / 3))
  fit <- fit_model(chlrm_model(2L, 2L), data,
    chlrm_prior(c(0, 0), diag(2), 4, diag(2), 2, 2, 2, c(0.5, 0.5)),
    method = "mcmc", initialization = initialization_from_prior(),
    control = mcmc_control(warmup = 2L, post_warmup_draws = 8L,
      chains = 1L, evaluation_draws = 2L, stopping = NULL, max_memory_mb = 50),
    seed = 72031L, verbose = FALSE)
  checkpoint <- file.path(directory, "baseline.rds")
  save_mcmc_checkpoint(fit, checkpoint, file.path(directory, "baseline-store"),
    chunk_size = 3L)
  history <- create_mcmc_history(checkpoint, file.path(directory, "archive"),
    chunk_size = 2L)
  list(fit = fit, history = history, X = X)
}

test_that("typed predictions from persisted CHLRM histories equal fit predictions", {
  directory <- tempfile(); dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  fixture <- history_prediction_fixture(directory)
  fit <- fixture$fit; history <- fixture$history
  history_draws <- posterior_draws_from_mcmc_history(
    fit, list(history), chains = 1L, chain_ids = 1L, iterations = 1:8,
    max_output_mb = 16)
  fit_draws <- posterior_draws(fit, draw_indices = 1:8, include_latent = TRUE,
    max_output_mb = 16)

  expect_equal(history_draws$parameters$component_coefficients,
    fit_draws$parameters$component_coefficients, tolerance = 0)
  expect_equal(history_draws$parameters$component_variances,
    fit_draws$parameters$component_variances, tolerance = 0)
  expect_equal(history_draws$parameters$mixture_weights,
    fit_draws$parameters$mixture_weights, tolerance = 0)
  expect_equal(history_draws$latent_variables$group_allocations,
    fit_draws$latent_variables$group_allocations, tolerance = 0)
  expect_identical(history_draws$draw_index$iteration, 1:8)
  selected_history <- posterior_draws_from_mcmc_history(fit, list(history),
    n_draws = 4, seed = 7791, max_output_mb = 16)
  selected_fit <- posterior_draws(fit, n_draws = 4, seed = 7791,
    include_latent = TRUE, max_output_mb = 16)
  expect_identical(selected_history$draw_index$stored_draw,
    selected_fit$draw_index$stored_draw)
  expect_identical(selected_history$parameters$component_coefficients,
    selected_fit$parameters$component_coefficients)

  newdata <- list(fixture$X[1:3, , drop = FALSE])
  fit_mean <- predict(fit, newdata, group_ids = 1L, n_draws = 8,
    return_draws = TRUE, max_output_mb = 16)
  history_mean <- predict(history_draws, newdata, group_ids = 1L,
    return_draws = TRUE, max_output_mb = 16)
  expect_identical(history_mean$summary, fit_mean$summary)
  expect_identical(history_mean$draws, fit_mean$draws)

  fit_response <- posterior_predict(fit, newdata, group_ids = 1L,
    n_draws = 8, seed = 4381, return_draws = TRUE, max_output_mb = 16)
  history_response <- posterior_predict(history_draws, newdata,
    group_ids = 1L, seed = 4381, return_draws = TRUE, max_output_mb = 16)
  expect_identical(history_response$summary, fit_response$summary)
  expect_identical(history_response$draws, fit_response$draws)

  first_half <- posterior_draws_from_mcmc_history(fit, list(history),
    iterations = 1:4, max_output_mb = 16)
  second_half <- posterior_draws_from_mcmc_history(fit, list(history),
    iterations = 5:8, max_output_mb = 16)
  set.seed(4381)
  streamed_response <- lapply(list(first_half, second_half), function(block)
    posterior_predict(block, newdata, group_ids = 1L, return_draws = TRUE,
      max_output_mb = 16)$draws[[1L]])
  expect_identical(do.call(cbind, streamed_response), fit_response$draws[[1L]])

  held_out <- list(X = list(fixture$X[1:3, , drop = FALSE]),
    y = list(c(-0.4, 0.1, 0.7)))
  fit_performance <- predictive_performance(fit, held_out,
    predictive_target = "observed_groups", predictive_unit = "observation",
    group_ids = 1L, n_draws = 8, seed = 4382, max_output_mb = 16)
  history_performance <- predictive_performance_from_draws(fit, held_out,
    history_draws, predictive_target = "observed_groups",
    predictive_unit = "observation", group_ids = 1L, seed = 4382,
    max_output_mb = 16)
  expect_identical(history_performance$pointwise, fit_performance$pointwise)
  expect_identical(history_performance$log_score_pointwise,
    fit_performance$log_score_pointwise)
  expect_identical(history_performance$statistics, fit_performance$statistics)
  direct_history_performance <- predictive_performance_from_mcmc_history(
    fit, held_out, list(history), chains = 1L, chain_ids = 1L,
    n_draws = 8, seed = 4382, group_ids = 1L, max_output_mb = 16)
  expect_identical(direct_history_performance$pointwise,
    fit_performance$pointwise)
  expect_identical(direct_history_performance$log_score_pointwise,
    fit_performance$log_score_pointwise)
  fit_sampled_performance <- predictive_performance(fit, held_out,
    predictive_target = "observed_groups", predictive_unit = "observation",
    group_ids = 1L, n_draws = 4, seed = 4383, max_output_mb = 16)
  sampled_history_performance <- predictive_performance_from_mcmc_history(
    fit, held_out, list(history), chains = 1L, chain_ids = 1L,
    n_draws = 4, seed = 4383, group_ids = 1L, max_output_mb = 16)
  expect_identical(sampled_history_performance$pointwise,
    fit_sampled_performance$pointwise)
  expect_identical(sampled_history_performance$log_score_pointwise,
    fit_sampled_performance$log_score_pointwise)
  expect_identical(sampled_history_performance$statistics,
    fit_sampled_performance$statistics)
})

test_that("history posterior draws enforce explicit memory and source contracts", {
  directory <- tempfile(); dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  fixture <- history_prediction_fixture(directory)
  fit <- fixture$fit; history <- fixture$history
  expect_error(posterior_draws_from_mcmc_history(fit, list(history),
    iterations = 1:8, max_output_mb = 1e-8), "max_output_mb")
  expect_error(posterior_draws_from_mcmc_history(fit, list(history),
    iterations = c(2L, 1L)), "increasing retained identities")
  expect_error(posterior_draws_from_mcmc_history(fit, list(history, history),
    chains = c(1L, 1L), chain_ids = c(1L, 1L)), "unique positive integer")
  selected <- posterior_draws_from_mcmc_history(fit, list(history, history),
    chains = c(1L, 1L), chain_ids = c(1L, 2L),
    iterations = list(c(1L, 3L), c(2L, 4L, 6L)), max_output_mb = 16)
  expect_identical(selected$draw_index$chain, c(1L, 1L, 2L, 2L, 2L))
  expect_identical(selected$draw_index$iteration, c(1L, 3L, 2L, 4L, 6L))
})
