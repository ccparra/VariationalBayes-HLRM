.fit_relabeling_fixture <- function(seed = 1701) {
  .test_chlrm_simulation(
    group_count = 6, group_sizes = 5, component_count = 2,
    coefficient_count = 2,
    component_coefficients = rbind(c(-1, 1), c(2, -1)),
    component_variances = c(.5, .8), seed = seed)
}

.fit_relabeling_prior <- function() chlrm_prior(
  c(0, 0), 10 * diag(2), 4, diag(2), 4, 2, 2, c(1, 1))

.fit_relabeling_cavi <- function(data, seed, scaling = NULL) fit_model(
  chlrm_model(2), data, .fit_relabeling_prior(), "cavi",
  initialization_kmeans(),
  cavi_control(
    max_iterations = 8, evaluation_draws = 4, starts = 2,
    require_converged = FALSE),
  scaling = scaling, seed = seed, verbose = FALSE)

test_that("relabel_fit applies one reference consistently and purely", {
  simulation <- .fit_relabeling_fixture()
  fit <- .fit_relabeling_cavi(simulation$train, 1702)
  reference <- fit_orchestration(fit)$label_reference
  baseline <- relabel_fit(fit, reference)
  permuted_reference <- reference
  permuted_reference$eta <- reference$eta[2:1, , drop = FALSE]
  before <- serialize(fit$parameters, NULL)

  relabeled <- relabel_fit(fit, permuted_reference)
  expect_identical(serialize(fit$parameters, NULL), before)
  expect_s3_class(relabeled, "vb_model_fit")
  expect_equal(
    relabeled$parameters$component_coefficients$mean,
    baseline$parameters$component_coefficients$mean[2:1, , drop = FALSE])
  expect_equal(
    relabeled$parameters$component_variances$shape,
    baseline$parameters$component_variances$shape[2:1])
  expect_equal(
    relabeled$parameters$mixture_weights$concentration,
    baseline$parameters$mixture_weights$concentration[2:1])
  expect_equal(
    relabeled$latent_variables$allocation_probabilities,
    baseline$latent_variables$allocation_probabilities[, 2:1, drop = FALSE])
  expect_equal(responsibilities(relabeled),
               responsibilities(baseline)[, 2:1, drop = FALSE],
               ignore_attr = TRUE)
  expect_equal(co_clustering(relabeled), co_clustering(baseline))
  expect_equal(
    predict(relabeled, simulation$train$X, n_draws = 4, seed = 1703),
    predict(baseline, simulation$train$X, n_draws = 4, seed = 1703))
  expect_identical(
    fit_orchestration(relabeled)$label_reference, permuted_reference)

  repeated <- relabel_fit(relabeled, permuted_reference)
  expect_identical(repeated$parameters, relabeled$parameters)
  expect_identical(repeated$latent_variables, relabeled$latent_variables)
  expect_equal(responsibilities(repeated), responsibilities(relabeled))
})

test_that("relabel_fit validates model dimensions and scaling coordinates", {
  simulation <- .fit_relabeling_fixture(1710)
  unscaled <- .fit_relabeling_cavi(simulation$train, 1711)
  unscaled_reference <- fit_orchestration(unscaled)$label_reference
  scaled <- .fit_relabeling_cavi(
    simulation$train, 1712, scaling = data_scaling())
  scaled_reference <- fit_orchestration(scaled)$label_reference

  scaled_relabeled <- relabel_fit(scaled, scaled_reference)
  expect_identical(scaled_relabeled$metadata$parameter_scale, "input")
  expect_equal(
    summarize_posterior(scaled_relabeled)$parameters$component_coefficients$mean,
    scaled_relabeled$parameters$component_coefficients$mean)
  expect_error(relabel_fit(scaled, unscaled_reference),
               "different coordinate scales")

  wrong_dimension <- unscaled_reference
  wrong_dimension$eta <- cbind(wrong_dimension$eta, 0)
  wrong_dimension$scale <- c(wrong_dimension$scale, 1)
  expect_error(relabel_fit(unscaled, wrong_dimension),
               "different coefficient dimensions")

  lrm_data <- .test_lrm_simulation(
    10, c(0, 1), 1, matrix(c(-1, 1), 1, 2), seed = 1713)
  lrm_fit <- fit_model(
    lrm_model(2), lrm_data$train,
    lrm_prior(c(0, 0), diag(2), 4, 1), "cavi",
    initialization_ols(), cavi_control(3, evaluation_draws = 2),
    seed = 1714, verbose = FALSE)
  expect_error(relabel_fit(lrm_fit, unscaled_reference), "only to CHLRM")
})

test_that("MCMC continuation preserves an externally applied label reference", {
  simulation <- .fit_relabeling_fixture(1720)
  reference_fit <- .fit_relabeling_cavi(simulation$train, 1721)
  reference <- fit_orchestration(reference_fit)$label_reference
  fit <- fit_model(
    chlrm_model(2), simulation$train, .fit_relabeling_prior(), "mcmc",
    initialization_from_prior(),
    mcmc_control(warmup = 1, post_warmup_draws = 3,
                 evaluation_draws = 2),
    seed = 1722, verbose = FALSE)
  relabeled <- relabel_fit(fit, reference)
  continued <- resume_mcmc(relabeled, additional_samples = 2)

  expect_identical(attr(continued, "label_reference", exact = TRUE), reference)
  expect_equal(
    dim(continued$parameters$component_coefficients$draws), c(5, 2, 2))
  expect_equal(dim(continued$latent_variables$group_allocations), c(6, 5))
  expect_identical(
    summarize_posterior(continued)$relabeling$algorithm,
    reference$algorithm)
})
