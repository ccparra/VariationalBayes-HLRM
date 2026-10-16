test_that("thinning is indexed relative to burn-in", {
  settings <- list(samples = 4, burn_in = 3, thin = 2)
  expect_equal(n_samples_mcmc(settings), 11)
  expect_equal(vapply(seq_len(11), function(b) {
    i <- thin_chain_mcmc(b, settings)
    if (is.null(i)) NA_integer_ else i
  }, integer(1)), c(NA, NA, NA, NA, 1, NA, 2, NA, 3, NA, 4))
  expect_equal(thin_chain_mcmc(4, list(burn_in = 3, thin = NULL)), 1)
  expect_error(thin_chain_mcmc(2, list(burn_in = 0, thin = 0)))
})

test_that("HLRM expected residual uses covariance only once", {
  object <- list(
    hyperparams = list(nu0 = 2),
    params = list(
      axi2 = 2, bxi2 = 2,
      MuBetaj = matrix(c(1, 2), 1),
      SigmaBetaj = list(diag(c(0.2, 0.3))),
      aj = 1, bj = 1
    ),
    stats = list(m = 1, nj = 2, X = list(diag(2)), y = list(c(2, 4)), XtX = list(diag(2)))
  )
  updated <- update_q_group_variances_hlrm(object)
  expected_ssr <- sum((c(2, 4) - c(1, 2))^2) + 0.5
  expect_equal(updated$params$bj, 0.5 * (2 + expected_ssr))
})

test_that("HLRM xi conditional uses inverse variances", {
  object <- list(
    hyperparams = list(a0 = 2, b0 = 3, nu0 = 4),
    params = list(sigma2j = c(2, 5)),
    stats = list(m = 2)
  )
  set.seed(11)
  draws <- replicate(20000, sample_xi2_hlrm(object)$params$xi2)
  target <- (2 + 4) / (3 + 2 * sum(1 / c(2, 5)))
  expect_equal(mean(draws), target, tolerance = 0.02)
})

test_that("HLRM log-prior matches its hierarchical factorization", {
  object <- list(
    hyperparams = list(Mu0 = c(0, 0), Lambda0 = diag(2), n0 = 4,
                       S0 = diag(2), nu0 = 3, a0 = 2, b0 = 1),
    params = list(Betaj = rbind(c(1, 0), c(0, 1)), Beta = c(.2, -.1),
                  Sigma = matrix(c(1, .1, .1, 1.2), 2), sigma2j = c(.8, 1.3), xi2 = .9)
  )
  manual <- sum(apply(object$params$Betaj, 1, mvtnorm::dmvnorm,
                      mean = object$params$Beta, sigma = object$params$Sigma, log = TRUE)) +
    mvtnorm::dmvnorm(object$params$Beta, object$hyperparams$Mu0,
                     object$hyperparams$Lambda0, log = TRUE) +
    LaplacesDemon::dinvwishart(object$params$Sigma, 4, diag(2), log = TRUE) +
    sum(dgamma(1 / object$params$sigma2j, 1.5, rate = 1.5 * .9, log = TRUE) -
          2 * log(object$params$sigma2j)) + dgamma(.9, 2, rate = 1, log = TRUE)
  expect_equal(logprior_hlrm(object), as.numeric(manual))
})

test_that("CHLRM relabeling jointly permutes all component state", {
  chains <- list(
    Betak = array(c(-2, 1, 3, 4, 3, 4, -2, 1), dim = c(2, 2, 2)),
    sigma2k = rbind(c(1, 4), c(4, 1)),
    omega = rbind(c(.3, .7), c(.7, .3)),
    z = cbind(c(1L, 2L, 1L), c(2L, 1L, 2L))
  )
  # Explicit construction avoids ambiguity in R's column-major array fill.
  chains$Betak[1, , ] <- rbind(c(-2, 1), c(3, 4))
  chains$Betak[2, , ] <- rbind(c(3, 4), c(-2, 1))
  reference <- fit_label_reference_chlrm(chains, reference_draw = 1)
  rel <- align_to_reference_chlrm(chains, reference)
  expect_equal(rel$Betak[1, , ], rel$Betak[2, , ])
  expect_equal(rel$sigma2k[1, ], rel$sigma2k[2, ])
  expect_equal(rel$omega[1, ], rel$omega[2, ])
  expect_equal(rel$z[, 1], rel$z[, 2])
  expect_equal(incidence_matrix_mcmc_chlrm(t(chains$z)),
               incidence_matrix_mcmc_chlrm(t(rel$z)))
})

test_that("short HLRM and CHLRM chains are complete and finite", {
  dh <- .test_hlrm_simulation(rep(10, 3),
    rbind(c(1, 2), c(-1, .5), c(2, -1)), rep(1, 3), seed = 8)
  hh <- hyperparams_hlrm(dh$train, "OLS")
  mh <- mcmc_hlrm$new(dh$train$y, dh$train$X, hh,
    list(samples = 8, burn_in = 3, thin = 2, samples_gof = 3))
  rh <- mh$fit(FALSE, seed = 8)
  expect_length(rh$objective$log_likelihood, 8)
  expect_true(all(is.finite(rh$objective$log_joint)))

  dc <- .test_chlrm_simulation(6, 8, 2, 2,
    component_coefficients = rbind(c(-2, 1), c(3, -1)),
    component_variances = c(1, 1.5),
    allocations = c(1, 1, 1, 2, 2, 2), seed = 9)
  hc <- hyperparams_chlrm(dc$train, 2, "OLS")
  mc <- mcmc_chlrm$new(dc$train$y, dc$train$X, 2, hc,
                       list(samples = 8, burn_in = 3, thin = 2, samples_gof = 3))
  rc <- mc$fit(FALSE, seed = 9)
  expect_length(rc$objective$log_likelihood, 8)
  expect_true(all(is.finite(rc$objective$log_joint)))
  expect_equal(mc$chains$Betaj[8, , ],
               mc$chains$Betak[8, mc$chains$z[, 8], ])
  expect_false(any(c("Betaj", "sigma2j", "nk") %in%
                   names(mc$chlrm_mcmc_object$params)))
})

test_that("CHLRM SVI remains finite and recovers separated components", {
  truth <- rbind(c(-3, .8), c(2, -.5))
  z <- rep(1:2, each = 5)
  d <- .test_chlrm_simulation(10, 15, 2, 2,
    component_coefficients = truth, component_variances = c(.8, 1.1),
    allocations = z, seed = 15)
  h <- hyperparams_chlrm(d$train, 2, "OLS")
  model <- vi_chlrm$new(d$train$y, d$train$X, 2, h,
                        list(iterations = 500, epsilon = NULL, stochastic = TRUE,
                             minibatch_m = 5, tau = 30, kappa = .7, samples_gof = 5))
  fit <- model$fit(FALSE)
  hat <- model$summarize_posterior()
  expect_true(all(is.finite(fit$objective$elbo)))
  expect_true(all(is.finite(hat$parameters$component_coefficients$mean)))
  expect_true(all(is.finite(hat$parameters$component_variances$mean)))
  expect_equal(length(unique(hat$allocations$labels)), 2)
})

test_that("CHLRM sufficient-statistic SSE equals the raw-data formula", {
  set.seed(19)
  X <- lapply(c(7, 9, 6), function(n) cbind(1, matrix(rnorm(2 * n), n, 2)))
  y <- lapply(X, function(Xj) rnorm(nrow(Xj)))
  XtX <- lapply(X, crossprod); Xty <- Map(crossprod, X, y)
  stats <- list(X = X, y = y, XtX = XtX, Xty = Xty,
                Xty_mat = t(vapply(Xty, as.numeric, numeric(3))),
                XtX_mat = vapply(XtX, as.numeric, numeric(9)),
                yty = vapply(y, function(v) drop(crossprod(v)), numeric(1)))
  mu <- c(.4, -.2, .7); Sigma <- diag(c(.2, .3, .1))
  optimized <- .expected_sse_batch_chlrm(stats, 1:3, mu, Sigma)
  raw <- vapply(1:3, function(j) {
    residual <- y[[j]] - drop(X[[j]] %*% mu)
    drop(crossprod(residual)) + sum(XtX[[j]] * Sigma)
  }, numeric(1))
  expect_equal(optimized, raw, tolerance = 1e-12)
})

test_that("SVI evaluates the full ELBO at configured spacing and final iteration", {
  d <- .test_hlrm_simulation(rep(8, 6),
    matrix(rep(c(1, .2), 6), 6, 2, byrow = TRUE), rep(1, 6), seed = 20)
  h <- hyperparams_hlrm(d$train, "OLS")
  model <- vi_hlrm$new(d$train$y, d$train$X, h,
                       list(iterations = 10, epsilon = NULL, stochastic = TRUE,
                            minibatch_m = 2, tau = 20, kappa = .7,
                            elbo_every = 3, samples_gof = 2))
  fit <- model$fit(FALSE)
  expect_equal(fit$objective$evaluation_iterations, c(3L, 6L, 9L, 10L))
  expect_length(fit$objective$elbo, 4)
  expect_true(fit$timing$objective_seconds >= 0)
  expect_true(fit$timing$update_seconds >= 0)
})

test_that("stable utilities respect limiting cases and WAIC definition", {
  expect_equal(logmeanexp(c(-Inf, -Inf)), -Inf)
  expect_equal(logsumexp(numeric()), -Inf)
  expect_true(is.na(r2(rep(1, 3), rep(1, 3))))
  ll <- matrix(c(-1, -2, -1.5, -2.5, -0.5, -3), nrow = 2)
  w <- waic_from_loglik(ll)
  expect_equal(w$p_waic, sum(apply(ll, 1, var)))
  expect_equal(w$WAIC, -2 * (w$lppd - w$p_waic))
  expect_error(logdet_spd(matrix(c(1, 2, 2, 1), 2)), "positive definite")
})

test_that("hierarchical simulators follow dimensions and positive variance support", {
  dh <- .test_hlrm_hierarchy(4:8, c(1, -.5),
    matrix(c(1, .2, .2, .5), 2), variance_scale = 1.2, seed = 71)
  expect_equal(dim(dh$truth$group_coefficients), c(5, 2))
  expect_true(all(dh$truth$group_variances > 0))
  expect_equal(lengths(dh$train$y), 4:8)

  dc <- .test_chlrm_hierarchy(rep(5, 9), 3, c(0, 1), diag(2),
    variance_scale = 1, mixture_weights = c(.2, .3, .5), seed = 72)
  expect_equal(dim(dc$truth$component_coefficients), c(3, 2))
  expect_true(all(dc$truth$allocations %in% 1:3))
  expect_true(all(dc$truth$component_variances > 0))
})

test_that("MCMC stores the complete hierarchy", {
  d <- .test_hlrm_hierarchy(rep(8, 3), c(0, 1), diag(2),
                               variance_scale = 1, seed = 80)
  h <- hyperparams_hlrm(d$train, "OLS")
  model <- mcmc_hlrm$new(d$train$y, d$train$X, h,
                         list(samples = 5, burn_in = 2, thin = 1, samples_gof = 3))
  fit <- model$fit(FALSE, seed = 81)
  expect_equal(dim(fit$parameters$population_coefficients$draws), c(5, 2))
  expect_equal(dim(fit$parameters$population_covariance$draws), c(5, 2, 2))
  expect_length(fit$parameters$variance_scale$draws, 5)
  summary <- model$summarize_posterior()
  expect_true(all(is.finite(summary$parameters$population_coefficients$mean)))
  expect_true(all(is.finite(summary$parameters$population_covariance$mean)))
  expect_true(all(is.finite(summary$parameters$variance_scale$mean)))
})
