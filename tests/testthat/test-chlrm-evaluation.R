evaluation_fixture <- function() list(
  X = list(matrix(1, 2, 1), matrix(1, 3, 1)),
  y = list(c(-2, -1), c(2, 3, 4)),
  draws = list(Betak = array(c(-2, -1, 2, 3), c(2, 2, 1)),
               sigma2k = rbind(c(.7, 1.2), c(.9, 1.5)),
               omega = rbind(c(.4, .6), c(.3, .7))),
  z = matrix(c(1, 2, 1, 2), 2))

test_that("compact alignment supports permutations and singleton dimensions", {
  d <- evaluation_fixture()
  switched <- d$draws
  switched$Betak[2, , ] <- d$draws$Betak[2, 2:1, ]
  switched$sigma2k[2, ] <- d$draws$sigma2k[2, 2:1]
  switched$omega[2, ] <- d$draws$omega[2, 2:1]
  reference <- fit_label_reference_chlrm(switched, reference_draw = 1)
  aligned <- align_to_reference_chlrm(switched, reference)
  expect_equal(aligned$Betak, d$draws$Betak)
  expect_equal(aligned$sigma2k, d$draws$sigma2k)
  expect_equal(aligned$omega, d$draws$omega)
  expect_null(aligned$z)
  one <- list(Betak = array(c(1, 2), c(2, 1, 1)),
              sigma2k = matrix(c(.5, 1), 2), omega = matrix(1, 2))
  one_reference <- fit_label_reference_chlrm(one, reference_draw = 1)
  expect_equal(align_to_reference_chlrm(one, one_reference)$Betak, one$Betak)
})

test_that("streaming CHLRM PPC matches stored replicates", {
  d <- evaluation_fixture()
  set.seed(98); before <- rng_snapshot()
  compact <- pp_check_chlrm(d$y, d$X, d$draws, z = d$z, seed = 99)
  full <- pp_check_chlrm(
    d$y, d$X, d$draws, z = d$z, seed = 99, return_yrep = TRUE)
  expect_identical(before, rng_snapshot())
  expect_identical(compact, full[setdiff(names(full), "y_rep")])
  expect_equal(dim(full$y_rep), c(5, 2))
  expect_equal(full$repdata$mean, colMeans(full$y_rep))
  expect_error(pp_check_chlrm(
    d$y, d$X, d$draws, z = d$z, return_yrep = TRUE,
    max_output_mb = 1e-8), "exceed")
  expect_error(pp_check_chlrm(d$y, d$X, d$draws), "exactly one")
})

test_that("CHLRM evaluation responsibilities have independent schemas", {
  d <- evaluation_fixture()
  fit <- evaluate_in_sample_fit_chlrm(
    d$y, d$X, d$draws, "MCMC", z = d$z)
  ppc <- posterior_predictive_check_chlrm(
    d$y, d$X, d$draws, "MCMC", z = d$z, seed = 10)
  criteria <- predictive_criteria_chlrm(
    d$y, d$X, d$draws, "MCMC")
  expect_identical(fit$schema_version, "in-sample-fit-1")
  expect_identical(ppc$schema_version, "posterior-predictive-check-2")
  expect_identical(criteria$schema_version, "predictive-criteria-2")
  expect_equal(fit$statistics$mse,
    mse(unlist(d$y), c(-1.5, -1.5, 2.5, 2.5, 2.5)))
  manual <- waic_from_loglik(log_likelihood(
    "CHLRM", d$y, d$X, d$draws,
    unit = "group", conditioning = "group_marginal"))
  expect_equal(criteria$criteria$waic_group_marginal$WAIC, manual$WAIC)
  expect_identical(
    criteria$criteria$waic_group_marginal$estimand$conditioning,
    "group_marginal")
  expect_identical(names(criteria$criteria), "waic_group_marginal")
  expect_error(predictive_criteria_chlrm(
    d$y, d$X, d$draws, "MCMC", max_output_mb = 1e-8), "exceed")
})

test_that("CHLRM CAVI and SVI evaluation is reproducible and non-mutating", {
  d <- .test_chlrm_simulation(4, 8, 2, 5, seed = 150)
  for (stochastic in c(FALSE, TRUE)) {
    model <- vi_chlrm$new(
      d$train$y, d$train$X, 2, .test_chlrm_native_prior(5, 2),
      list(iterations = 10, stochastic = stochastic, seed = 151,
           samples_gof = 8))
    model$fit(FALSE)
    state <- serialize(model$chlrm_vi_object, NULL)
    before <- rng_snapshot()
    fit <- model$evaluate_in_sample_fit()
    ppc <- model$posterior_predictive_check(seed = 152)
    criteria <- model$predictive_criteria(seed = 152)
    retained <- model$predictive_criteria(
      seed = 152, return_log_likelihood = TRUE)
    expect_identical(ppc, model$posterior_predictive_check(seed = 152))
    expect_identical(criteria, model$predictive_criteria(seed = 152))
    expect_identical(before, rng_snapshot())
    expect_identical(state, serialize(model$chlrm_vi_object, NULL))
    expect_equal(fit$method, if (stochastic) "SVI" else "CAVI")
    expect_true(is.finite(criteria$criteria$waic_group_marginal$WAIC))
    expect_s3_class(retained$log_likelihood$group,
                    "vb_log_likelihood")
    expect_error(model$predictive_criteria(1), "integer >= 2")
    expect_error(model$posterior_predictive_check(
      8, max_output_mb = 1e-8), "exceed")
  }
})

test_that("CHLRM MCMC criteria do not mutate chains", {
  d <- .test_chlrm_simulation(4, 8, 2, 5, seed = 150)
  model <- mcmc_chlrm$new(
    d$train$y, d$train$X, 2, .test_chlrm_native_prior(5, 2),
    list(samples = 4, burn_in = 2, thin = 1, seed = 160,
         samples_gof = 3))
  model$fit(FALSE)
  state <- serialize(model$chains, NULL)
  criteria <- model$predictive_criteria(seed = 161)
  retained <- model$predictive_criteria(
    seed = 161, return_log_likelihood = TRUE)
  expect_identical(criteria, model$predictive_criteria(seed = 161))
  expect_identical(state, serialize(model$chains, NULL))
  expect_true(is.finite(criteria$criteria$waic_group_marginal$WAIC))
  expect_s3_class(retained$log_likelihood$group,
                  "vb_log_likelihood")
  expect_error(model$predictive_criteria(5), "exceeds")
  expect_error(model$predictive_criteria(
    3, predictive_target = "observed_groups"), "new_groups")
})

test_that("intercept-only single-component fits support every inference engine", {
  d <- .test_chlrm_simulation(
    3, 4, 1, 1, component_coefficients = matrix(2), seed = 170)
  prior <- .test_chlrm_native_prior(1, 1)
  for (stochastic in c(FALSE, TRUE)) {
    model <- vi_chlrm$new(d$train$y, d$train$X, 1, prior,
      list(iterations = 4, epsilon = NULL, stochastic = stochastic,
           minibatch_m = 1, seed = 171))
    model$fit(FALSE)
    expect_true(is.finite(
      model$predictive_criteria(5, 172)$criteria$waic_group_marginal$WAIC))
    expect_equal(dim(model$chlrm_vi_object$stats$Xty_mat), c(3, 1))
  }
  model <- mcmc_chlrm$new(d$train$y, d$train$X, 1, prior,
    list(samples = 4, burn_in = 1, thin = 1, seed = 171))
  model$fit(FALSE)
  expect_true(is.finite(
    model$predictive_criteria(4, 172)$criteria$waic_group_marginal$WAIC))
})
