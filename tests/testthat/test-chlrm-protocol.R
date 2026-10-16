# Load independently of the order in which testthat executes test files.

test_that("simulation is reproducible and preserves a caller's seeded RNG", {
  set.seed(431)
  before <- .Random.seed
  a <- .test_chlrm_simulation(4, c(8, 9, 10, 11), 2, 5,
                               test_group_sizes = 3, seed = 52)
  expect_identical(.Random.seed, before)
  b <- .test_chlrm_simulation(4, c(8, 9, 10, 11), 2, 5,
                               test_group_sizes = 3, seed = 52)
  expect_identical(a, b)
  expect_equal(lengths(a$train$y), 8:11)
  expect_equal(lengths(a$test$y), rep(3, 4))
  expect_true(all(vapply(a$train$X, function(x) all(x[, 1] == 1), logical(1))))
  expect_false(identical(a$train$X[[1]][1:3, ], a$test$X[[1]]))
  expect_equal(a$truth$mixture_weights, c(.5, .5))
  expect_s3_class(a, "vb_simulation_result")
  expect_identical(a$schema_version, "simulation-result-2")
})

test_that("simulation obeys covariate covariance and conditional noise variance", {
  V <- matrix(c(1, .6, .6, 2), 2)
  a <- .test_chlrm_simulation(50, 100, 1, 3,
    covariate_covariance = V, component_variances = 2,
    component_coefficients = matrix(c(2, 1, -1), 1), seed = 81)
  X <- do.call(rbind, a$train$X)
  residual <- unlist(a$train$y) -
    drop(X %*% t(a$truth$component_coefficients))
  expect_equal(unname(cov(X[, -1])), V, tolerance = .08)
  expect_lt(abs(mean(residual)), .08)
  expect_lt(abs(var(residual) - 2), .15)
  expect_true(all(a$truth$allocations == 1))
  single <- .test_chlrm_simulation(1, 2, 1, 1,
                                    component_coefficients = matrix(3), seed = 1)
  expect_equal(single$train$X[[1]], matrix(1, 2, 1))
  expect_error(.test_chlrm_simulation(2, 0, 2, 5), "observations_per_group")
  expect_error(.test_chlrm_simulation(2, 3, 2, 5,
                                      mixture_weights = c(1, 1)),
               "mixture_weights")
  expect_error(.test_chlrm_simulation(2, 3, 2, 5,
                                      covariate_covariance = diag(3)),
               "covariance dimension")
  expect_error(.test_chlrm_simulation(2, 3, 2, 5,
                                      covariate_covariance = -diag(4)))
})

test_that("marginal likelihood integrates one allocation per group", {
  X <- list(matrix(1, 2, 1), matrix(1, 1, 1))
  y <- list(c(-1, 1), 2)
  draws <- list(Betak = array(c(-1, -2, 1, 2), c(2, 2, 1)),
                sigma2k = rbind(c(1, 1), c(.7, 1.4)),
                omega = rbind(c(.3, .7), c(.6, .4)))
  ll <- loglik_group_marginal_chlrm(y, X, draws)
  manual <- log(.3 * prod(dnorm(y[[1]], -1)) + .7 * prod(dnorm(y[[1]], 1)))
  observation_mixture <- sum(log(.3 * dnorm(y[[1]], -1) + .7 * dnorm(y[[1]], 1)))
  expect_equal(ll[1, 1], manual)
  expect_gt(abs(ll[1, 1] - observation_mixture), .1)
  permuted <- draws
  # Change labels in just one draw, jointly for all component parameters.
  permuted$Betak[2, , ] <- draws$Betak[2, 2:1, ]
  permuted$sigma2k[2, ] <- draws$sigma2k[2, 2:1]
  permuted$omega[2, ] <- draws$omega[2, 2:1]
  expect_equal(loglik_group_marginal_chlrm(y, X, permuted), ll)
  w <- waic_from_loglik(ll)
  expected_lppd <- sum(log(rowMeans(exp(ll))))
  expected_penalty <- sum(apply(ll, 1, var))
  expect_equal(w$WAIC, -2 * (expected_lppd - expected_penalty))
  expect_equal(waic_from_loglik(loglik_group_marginal_chlrm(y, X, permuted)), w)
  expect_error(loglik_group_marginal_chlrm(y, X, draws, max_output_mb = 1e-8),
               "exceeds")
  draws$omega[1, ] <- c(0, 1)
  expect_equal(loglik_group_marginal_chlrm(y, X, draws)[1, 1],
               sum(dnorm(y[[1]], 1, log = TRUE)))
  draws$sigma2k[1, 1] <- 0
  expect_error(loglik_group_marginal_chlrm(y, X, draws), "strictly positive")
})

test_that("marginal likelihood handles singleton dimensions and underflow", {
  draws <- list(Betak = array(0, c(1, 1, 1)), sigma2k = matrix(1), omega = matrix(1))
  ll <- loglik_group_marginal_chlrm(list(rep(50, 10)), list(matrix(1, 10, 1)), draws)
  expect_equal(dim(ll), c(1, 1))
  expect_equal(ll[1, 1], 10 * dnorm(50, log = TRUE))
  expect_true(is.finite(ll[1, 1]))
})

test_that("protocol prior and marginal likelihood integrate with MCMC CAVI and SVI", {
  d <- .test_chlrm_simulation(6, 10, 2, 5, seed = 80)
  prior <- .test_chlrm_native_prior(5, 2)
  for (stochastic in c(FALSE, TRUE)) {
    set.seed(91)
    fit <- vi_chlrm$new(d$train$y, d$train$X, 2, prior,
                       list(iterations = 12, epsilon = 0, stochastic = stochastic,
                            minibatch_m = 3, tau = 10, kappa = .7))
    fit$fit(verbose = FALSE)
    expect_true(all(is.finite(fit$results$elbo)))
    draws <- vi_sample_draws_chlrm(fit$chlrm_vi_object, fit$results, 10, seed = 92,
                                   include_allocations = FALSE)
    expect_setequal(names(draws), c("Betak", "sigma2k", "omega"))
    expect_equal(rowSums(draws$omega), rep(1, 10))
    ll <- loglik_group_marginal_chlrm(d$train$y, d$train$X, draws)
    expect_equal(dim(ll), c(6, 10))
    expect_true(all(is.finite(ll)))
    expect_true(is.finite(waic_from_loglik(ll)$WAIC))
    full <- vi_sample_draws_chlrm(fit$chlrm_vi_object, fit$results, 10, seed = 92)
    expect_identical(full$Betak, draws$Betak)
    expect_identical(full$sigma2k, draws$sigma2k)
    expect_lt(as.numeric(object.size(draws)), as.numeric(object.size(full)))
    many <- vi_sample_draws_chlrm(fit$chlrm_vi_object, fit$results, 5000, seed = 17,
                                  include_allocations = FALSE)
    alpha <- fit$results$alpha
    expect_equal(colMeans(many$omega), alpha / sum(alpha), tolerance = .02)
  }
  set.seed(93)
  fit <- mcmc_chlrm$new(d$train$y, d$train$X, 2, prior,
                       list(samples = 5, burn_in = 2, thin = 1))
  fit$fit(verbose = FALSE)
  ll <- loglik_group_marginal_chlrm(d$train$y, d$train$X, fit$chains)
  expect_equal(dim(ll), c(6, 5))
  expect_true(all(is.finite(ll)))
})
