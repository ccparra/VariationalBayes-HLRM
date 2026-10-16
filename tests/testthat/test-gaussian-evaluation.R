test_that("Gaussian fit, PPC and predictive criteria are separate contracts", {
  X <- list(matrix(1, 2, 1), matrix(1, 3, 1))
  y <- list(c(0, 1), c(1, 2, 3))
  draws <- list(Betaj = array(c(.1, .3, 1.5, 1.8), c(2, 2, 1)),
                sigma2j = matrix(c(1, 1.2, .5, .8), 2))
  point_beta <- matrix(c(.2, 1.65), 2, 1)
  ll <- sapply(1:2, function(s) unlist(lapply(1:2, function(j)
    dnorm(y[[j]], draws$Betaj[s, j, 1], sqrt(draws$sigma2j[s, j]), log = TRUE))))

  fit <- evaluate_in_sample_fit_gaussian(y, X, point_beta, "HLRM", "MCMC")
  ppc <- posterior_predictive_check_gaussian(y, X, draws, "HLRM", "MCMC", seed = 1)
  criteria <- predictive_criteria_gaussian(
    y, X, draws, "HLRM", "MCMC")

  expect_s3_class(fit, "vb_in_sample_fit")
  expect_s3_class(ppc, "vb_posterior_predictive_check")
  expect_s3_class(criteria, "vb_predictive_criteria")
  expect_equal(fit$statistics$mse,
    mean((unlist(y) - c(.2, .2, 1.65, 1.65, 1.65))^2))
  expect_equal(criteria$criteria$waic_observation_conditional$WAIC,
               waic_from_loglik(ll)$WAIC)
  expect_equal(criteria$criteria$waic_group_conditional$WAIC,
               waic_from_loglik(rowsum(ll, c(1, 1, 2, 2, 2)))$WAIC)
  expect_identical(criteria$criteria$waic_observation_conditional$estimand$unit,
                   "observation")
  expect_identical(criteria$criteria$waic_group_conditional$estimand$unit,
                   "group")
  expect_identical(names(criteria$criteria),
                   c("waic_observation_conditional", "waic_group_conditional"))
  expect_error(predictive_criteria_gaussian(
    y, X, draws, "HLRM", "CAVI", max_output_mb = 1e-8), "exceed")
  draws$sigma2j[1, 1] <- 0
  expect_error(predictive_criteria_gaussian(y, X, draws, "HLRM", "CAVI"),
               "positive")
})

test_that("Gaussian PPC streams the same canonical discrepancies", {
  y <- list(1, c(2, 3)); X <- list(matrix(1), matrix(1, 2, 1))
  draws <- list(Betaj = array(c(1, 1, 2, 2), c(2, 2, 1)),
                sigma2j = matrix(1, 2, 2))
  set.seed(66); before <- rng_snapshot()
  a <- pp_check_gaussian(y, X, draws, seed = 67)
  b <- pp_check_gaussian(y, X, draws, seed = 67, return_yrep = TRUE)
  expect_identical(a, b[setdiff(names(b), "y_rep")])
  expect_identical(rng_snapshot(), before)
  expect_equal(b$repdata$mean, colMeans(b$y_rep))
  expect_identical(names(a$repdata),
    c("mean", "sd", "between_group_mean_variance", "ols_mse"))
  expect_false("mean_within_group_variance" %in% names(a$repdata))
  expect_error(pp_check_gaussian(
    y, X, draws, return_yrep = TRUE, max_output_mb = 1e-8), "exceed")
})

test_that("LRM MCMC, CAVI and SVI expose the evaluation responsibilities", {
  X <- cbind(1, seq(-1, 1, length.out = 12))
  y <- 1 + 2 * X[, 2] + sin(seq_len(12)) / 3
  hyper <- hyperparams_lrm(list(X = X, y = y), "OLS")
  models <- list(
    mcmc_lrm$new(y, X, hyper,
      list(samples = 6, burn_in = 2, thin = 1, samples_gof = 4)),
    vi_lrm$new(y, X, hyper,
      list(iterations = 8, epsilon = NULL, samples_gof = 4)),
    vi_lrm$new(y, X, hyper,
      list(iterations = 8, epsilon = NULL, stochastic = TRUE,
           minibatch_m = 4, elbo_every = 2, min_iter = 1,
           samples_gof = 4, seed = 100)))
  for (model in models) {
    expect_error(model$evaluate_in_sample_fit(), "Call fit")
    expect_error(model$posterior_predictive_check(), "Call fit")
    expect_error(model$predictive_criteria(), "Call fit")
    model$fit(FALSE)
    before <- rng_snapshot()
    fit <- model$evaluate_in_sample_fit()
    ppc <- model$posterior_predictive_check(4, seed = 102)
    criteria <- model$predictive_criteria(4, seed = 102)
    retained <- model$predictive_criteria(
      4, seed = 102, return_log_likelihood = TRUE)
    expect_identical(criteria, model$predictive_criteria(4, seed = 102))
    expect_identical(rng_snapshot(), before)
    expect_identical(fit$schema_version, "in-sample-fit-1")
    expect_identical(ppc$schema_version, "posterior-predictive-check-2")
    expect_identical(criteria$schema_version, "predictive-criteria-2")
    expect_true(is.finite(criteria$criteria$waic_observation$WAIC))
    expect_s3_class(retained$log_likelihood$observation,
                    "vb_log_likelihood")
    grouped <- model$predictive_criteria(
      4, seed = 102, predictive_unit = "group",
      groups = rep(1:3, each = 4))
    expect_true(is.finite(grouped$criteria$waic_group$WAIC))
    expect_identical(grouped$criteria$waic_group$estimand$conditioning,
                     "parameter_conditional")
    expect_error(model$predictive_criteria(1), "integer >= 2")
  }
  expect_error(models[[1]]$predictive_criteria(7), "exceeds")
  expect_error(models[[1]]$predictive_criteria(
    4, predictive_target = "new_groups"), "observed_data")
})

test_that("HLRM selects conditional or marginal likelihood internally", {
  set.seed(110)
  X <- rep(list(cbind(1, seq(-1, 1, length.out = 8))), 3)
  y <- lapply(seq_along(X), function(j)
    drop(X[[j]] %*% c(j, .5)) + rnorm(8, sd = .3))
  hyper <- hyperparams_hlrm(list(X = X, y = y), "OLS")
  models <- list(
    mcmc_hlrm$new(y, X, hyper,
      list(samples = 5, burn_in = 2, thin = 1, samples_gof = 3)),
    vi_hlrm$new(y, X, hyper,
      list(iterations = 8, epsilon = NULL, samples_gof = 3)),
    vi_hlrm$new(y, X, hyper,
      list(iterations = 8, epsilon = NULL, stochastic = TRUE,
           minibatch_m = 2, samples_gof = 3)))
  for (i in seq_along(models)) {
    model <- models[[i]]; model$fit(FALSE)
    fit <- model$evaluate_in_sample_fit()
    ppc <- model$posterior_predictive_check(3, 112)
    conditional <- model$predictive_criteria(3, 112, "observed_groups")
    retained <- model$predictive_criteria(
      3, 112, "observed_groups", return_log_likelihood = TRUE)
    marginal <- model$predictive_criteria(
      3, 112, "new_groups", marginal_control = list(tolerance = 1e-5))
    expect_equal(c(fit$method, ppc$method, conditional$method, marginal$method),
                 rep(c("MCMC", "CAVI", "SVI")[i], 4))
    expect_identical(
      conditional$criteria$waic_group_conditional$estimand$conditioning,
      "parameter_conditional")
    expect_identical(
      marginal$criteria$waic_group_marginal$estimand$conditioning,
      "group_marginal")
    expect_true(is.finite(conditional$criteria$waic_group_conditional$WAIC))
    expect_s3_class(retained$log_likelihood$group,
                    "vb_log_likelihood")
    expect_true(is.finite(marginal$criteria$waic_group_marginal$WAIC))
    expect_identical(marginal, model$predictive_criteria(
      3, 112, "new_groups", marginal_control = list(tolerance = 1e-5)))
    expect_error(model$predictive_criteria(
      3, 112, "new_groups",
      marginal_control = list(max_evaluations = 1)), "max_evaluations")
  }
})

test_that("singleton dimensions and LRM/HLRM likelihood equivalence are preserved", {
  y <- list(c(1, 2, 3)); X <- list(matrix(1, 3, 1))
  draws <- list(Betaj = array(c(1, 2), c(2, 1, 1)),
                sigma2j = matrix(c(1, 2), 2, 1))
  lrm <- predictive_criteria_gaussian(y, X, draws, "LRM", "MCMC")
  hlrm <- predictive_criteria_gaussian(y, X, draws, "HLRM", "MCMC")
  expect_equal(lrm$criteria$waic_observation$WAIC,
               hlrm$criteria$waic_observation_conditional$WAIC)
  expect_identical(lrm$criteria$waic_observation$estimand$conditioning,
                   hlrm$criteria$waic_observation_conditional$estimand$conditioning)
  expect_true(is.na(hlrm$criteria$waic_group_conditional$se))
  grouped <- predictive_criteria_gaussian(
    y, X, draws, "LRM", "MCMC", predictive_unit = "group",
    groups = c("a", "b", "b"))
  ll <- sapply(1:2, function(s)
    dnorm(y[[1]], draws$Betaj[s, 1, 1],
          sqrt(draws$sigma2j[s, 1]), log = TRUE))
  expect_equal(grouped$criteria$waic_group$WAIC,
               waic_from_loglik(rowsum(ll, c(1, 2, 2)))$WAIC)
})
