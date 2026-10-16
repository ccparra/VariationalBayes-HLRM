test_that("variational settings fail early and distinguish CAVI and SVI defaults", {
  vi <- vi_stopping_control(list(iterations = 10), 4)
  svi <- vi_stopping_control(list(iterations = 10, stochastic = TRUE), 4)
  expect_equal(c(vi$every, vi$window, vi$patience, vi$epsilon), c(1, 1, 5, 1e-8))
  expect_equal(c(svi$every, svi$window, svi$patience, svi$min_iter, svi$epsilon),
               c(50, 5, 5, 1000, 1e-4))
  expect_null(vi_stopping_control(list(iterations = 10, epsilon = NULL), 4)$epsilon)
  for (bad in list(0, NA_real_, 2.5, c(2, 3)))
    expect_error(vi_stopping_control(list(iterations = bad), 4), "iterations")
  expect_error(vi_stopping_control(list(iterations = 10, stochastic = 1), 4), "stochastic")
  expect_error(vi_stopping_control(list(iterations = 10, epsilon = -1), 4), "epsilon")
  expect_error(vi_stopping_control(list(iterations = 10, elbo_every = 2), 4), "CAVI")
  expect_error(vi_stopping_control(list(iterations = 10, minibatch_m = 5), 4), "minibatch")
  expect_error(vi_stopping_control(list(iterations = 10, kappa = .5), 4), "kappa")
  expect_error(vi_stopping_control(list(iterations = 10, tau = -1), 4), "tau")
})

test_that("CAVI patience counts consecutive changes and resets after a large change", {
  control <- vi_stopping_control(list(iterations = 20, epsilon = .01), 4)
  state <- new_vi_monitor()
  values <- c(-100, -100, -100, -50, rep(-50, 5))
  for (i in seq_along(values)) {
    state <- update_vi_monitor(state, values[i], i, control)
    if (i < length(values)) expect_false(state$converged)
  }
  expect_true(state$converged)
  expect_equal(state$stable, 5)
  state <- update_vi_monitor(new_vi_monitor(), -2, 1, control)
  state <- update_vi_monitor(state, -1.5, 2, control)
  expect_equal(state$relative_change, .5 / 3)
  expect_error(update_vi_monitor(state, -Inf, 3, control), "Non-finite ELBO")
})

test_that("SVI waits for the window, minimum iterations and five eligible checkpoints", {
  control <- vi_stopping_control(list(iterations = 2000, stochastic = TRUE), 4)
  state <- new_vi_monitor()
  for (i in seq(50, 1200, by = 50)) {
    state <- update_vi_monitor(state, -100, i, control)
    if (i < 1200) expect_false(state$converged)
    if (i < 1000) expect_equal(state$stable, 0)
  }
  expect_true(state$converged)
  expect_equal(state$stable, 5)
  expect_equal(state$smoothed_elbo, -100)
  state <- new_vi_monitor()
  control$min_iter <- 1
  for (i in seq_len(5)) state <- update_vi_monitor(state, i, i * 50, control)
  expect_equal(state$smoothed_elbo, 3)
  expect_true(is.na(state$relative_change))
  state <- update_vi_monitor(state, 6, 300, control)
  expect_equal(state$smoothed_elbo, 4)
  expect_equal(state$relative_change, 1 / 4)
  for (epsilon in list(NULL, 0)) {
    control$epsilon <- epsilon
    state <- new_vi_monitor()
    for (i in 1:20) state <- update_vi_monitor(state, -100, i * 50, control)
    expect_false(state$converged)
  }
})

test_that("isolated RNG restores caller state on success and errors, including absent seed", {
  original <- rng_snapshot()
  on.exit(restore_rng_snapshot(original))
  seeded_operation <- function(fail = FALSE) {
    scope <- inference_rng_scope(123)
    on.exit(scope$restore())
    x <- rnorm(3)
    if (fail) stop("deliberate failure")
    x
  }
  set.seed(900)
  before <- rng_snapshot()
  expect_identical(seeded_operation(), seeded_operation())
  expect_identical(rng_snapshot(), before)
  expect_error(seeded_operation(TRUE), "deliberate")
  expect_identical(rng_snapshot(), before)
  rm(".Random.seed", envir = .GlobalEnv)
  seeded_operation()
  expect_false(exists(".Random.seed", .GlobalEnv, inherits = FALSE))
  expect_error(inference_rng_scope(-1), "seed")
})

test_that("seed controls constructor and fit, survives serialization and isolates callers", {
  d <- .test_chlrm_simulation(6, 8, 2, 5, seed = 110)
  prior <- .test_chlrm_native_prior(5, 2)
  for (stochastic in c(FALSE, TRUE)) {
    settings <- list(iterations = 12, stochastic = stochastic, seed = 111,
                     epsilon = NULL, minibatch_m = 3, elbo_every = if (stochastic) 5 else 1)
    set.seed(190)
    before <- rng_snapshot()
    a <- vi_chlrm$new(d$train$y, d$train$X, 2, prior, settings)
    expect_identical(rng_snapshot(), before)
    b <- unserialize(serialize(a, NULL))
    expect_error(a$summarize_posterior(), "Call fit")
    expect_error(a$predictive_criteria(), "Call fit")
    fit <- a$fit(FALSE)
    expect_identical(rng_snapshot(), before)
    runif(20)
    before <- rng_snapshot()
    b$fit(FALSE)
    expect_identical(a$results, b$results)
    expect_identical(a$rng_state, b$rng_state)
    expect_identical(rng_snapshot(), before)
    expect_identical(a$fit_status, "fitted")
    expect_identical(a$predictive_criteria(3, seed = 901)$schema_version,
                     "predictive-criteria-2")
    expect_false(fit$diagnostics$converged)
    expect_identical(fit$diagnostics$stop_reason, "max_iterations")
    expect_equal(fit$diagnostics$iterations, 12)
    if (stochastic) expect_equal(fit$objective$evaluation_iterations,
                                 c(5, 10, 12))
    expect_error(a$fit(FALSE), "already been fitted")
  }
  settings <- list(samples = 4, burn_in = 2, thin = 1, seed = 112)
  before <- rng_snapshot()
  a <- mcmc_chlrm$new(d$train$y, d$train$X, 2, prior, settings)
  b <- mcmc_chlrm$new(d$train$y, d$train$X, 2, prior, settings)
  expect_identical(rng_snapshot(), before)
  expect_error(a$summarize_posterior(), "Call fit")
  expect_error(a$predictive_criteria(), "Call fit")
  a$fit(FALSE)
  runif(10)
  before <- rng_snapshot()
  b$fit(FALSE)
  expect_identical(a$chains, b$chains)
  expect_identical(rng_snapshot(), before)
  expect_identical(a$fit_status, "fitted")
  expect_identical(a$predictive_criteria(3, seed = 901)$schema_version,
                   "predictive-criteria-2")
  expect_error(a$fit(FALSE), "already been fitted")
})

test_that("same seed yields matching CAVI/SVI initial states and full-batch updates", {
  d <- .test_chlrm_simulation(4, 8, 2, 5, seed = 121)
  make_model <- function(stochastic) vi_chlrm$new(d$train$y, d$train$X, 2,
    .test_chlrm_native_prior(5, 2), list(iterations = 10, stochastic = stochastic, seed = 122))
  a <- make_model(FALSE)
  b <- make_model(TRUE)
  expect_identical(a$chlrm_vi_object$params, b$chlrm_vi_object$params)
  x <- update_q_allocations_chlrm(a$chlrm_vi_object)
  y <- update_q_allocations_chlrm(b$chlrm_vi_object, batch_idx = 4:1)
  for (f in list(update_q_component_coefficients_chlrm,
                 update_q_component_variances_chlrm,
                 update_q_mixture_weights_chlrm)) {
    x <- f(x)
    y <- f(y, batch_idx = 4:1, learning_rate = 1)
  }
  for (f in list(update_q_global_coefficients_chlrm,
                 update_q_global_covariance_chlrm,
                 update_q_xi2_chlrm)) {
    x <- f(x)
    y <- f(y, learning_rate = 1)
  }
  expect_equal(x$params, y$params, tolerance = 1e-12)
})

test_that("fit reports tolerance stop, refuses partial restart and restores RNG on failure", {
  d <- .test_chlrm_simulation(4, 8, 2, 5, seed = 131)
  model <- function() vi_chlrm$new(d$train$y, d$train$X, 2,
    .test_chlrm_native_prior(5, 2), list(iterations = 10, seed = 132, epsilon = 1e10))
  a <- model()
  fit <- a$fit(FALSE)
  expect_true(fit$diagnostics$converged)
  expect_equal(fit$diagnostics$iterations, 6)
  expect_identical(fit$diagnostics$stop_reason, "objective_tolerance")
  svi <- vi_chlrm$new(d$train$y, d$train$X, 2, .test_chlrm_native_prior(5, 2),
    list(iterations = 12, seed = 132, stochastic = TRUE, epsilon = 1e10,
         elbo_every = 2, elbo_window = 2, min_iter = 6, patience = 2))
  svi_fit <- svi$fit(FALSE)
  expect_true(svi_fit$diagnostics$converged)
  expect_equal(svi_fit$diagnostics$iterations, 8)
  expect_equal(svi_fit$objective$evaluation_iterations, c(2, 4, 6, 8))
  bad <- model()
  bad$chlrm_vi_object$params$MuBetak[1, 1] <- NA_real_
  before <- rng_snapshot()
  expect_error(bad$fit(FALSE))
  expect_identical(bad$fit_status, "failed")
  expect_identical(bad$diagnostics$stop_reason, "error")
  expect_identical(rng_snapshot(), before)
  expect_error(bad$fit(FALSE), "already been fitted or failed")
})

test_that("parameter stability is enforced when explicitly requested", {
  control <- vi_stopping_control(list(
    iterations = 10, epsilon = 1, parameter_tolerance = 0.01,
    min_iter = 1, patience = 2), 1)
  state <- new_vi_monitor()
  state <- update_vi_monitor(state, -10, 1, control, list(mean = 0))
  state <- update_vi_monitor(state, -10, 2, control, list(mean = 1))
  state <- update_vi_monitor(state, -10, 3, control, list(mean = 2))
  expect_false(state$converged)
  expect_gt(state$parameter_change, control$parameter_tolerance)
  state <- update_vi_monitor(state, -10, 4, control, list(mean = 2))
  state <- update_vi_monitor(state, -10, 5, control, list(mean = 2))
  expect_true(state$converged)
  expect_identical(vi_convergence_reason(control, TRUE),
                   "objective_and_parameter_tolerance")
})

test_that("variational stopping coordinates use finite natural parameters", {
  lrm <- vi_convergence_parameters(list(
    MuBeta = c(1, -1), SigmaBeta = diag(c(2, 3)),
    asigma2 = 4, bsigma2 = 5), "LRM")
  expect_equal(lrm[1:2], c(.5, -1 / 3))
  expect_equal(tail(lrm, 2), c(-5, -5))
  expect_true(all(is.finite(lrm)))
})
