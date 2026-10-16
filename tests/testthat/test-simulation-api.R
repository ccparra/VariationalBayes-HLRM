test_that("typed simulation specifications reject ambiguous requests", {
  expect_error(conditional_dgp(lrm_model(2), list(coefficients = c(0, 1))),
               "Missing conditional")
  expect_error(hierarchical_dgp(lrm_model(2), list(a = 1)),
               "only to HLRM")
  expect_error(prior_predictive_dgp(
    lrm_model(2), hlrm_prior(c(0, 0), diag(2), 4, diag(2), 4, 2, 2)),
    "model-matched")
  expect_error(simulation_sample(sample_size = 10, group_count = 2,
                                 observations_per_group = 5),
               "either sample_size")
  expect_error(simulation_design("mixed", binary_probabilities = 1),
               "strictly between")
})

test_that("conditional DGP is uniform across LRM, HLRM and CHLRM", {
  design <- simulation_design(
    "mixed", mean = 2, sd = .5, binary_probabilities = .35,
    binary_position = "first",
    column_names = c("intercept", "binary", "numeric"))
  grouped_sample <- simulation_sample(
    group_count = 4, observations_per_group = 20,
    test_observations_per_group = 5)
  lrm <- simulate_data(
    conditional_dgp(lrm_model(3), list(
      coefficients = c(1, .5, -1), residual_variance = .7)),
    design, simulation_sample(sample_size = 80, test_sample_size = 20),
    seed = 1301)
  hlrm <- simulate_data(
    conditional_dgp(hlrm_model(3), list(
      group_coefficients = matrix(rep(c(1, .5, -1), 4), 4, byrow = TRUE),
      group_variances = rep(.7, 4))),
    design, grouped_sample, seed = 1302)
  chlrm <- simulate_data(
    conditional_dgp(chlrm_model(2, 3), list(
      component_coefficients = rbind(c(1, .5, -1), c(-1, -.5, 1)),
      component_variances = c(.7, 1.2), mixture_weights = c(.4, .6),
      allocations = c(1, 2, 1, 2))),
    design, grouped_sample, seed = 1303)
  simulations <- list(lrm, hlrm, chlrm)
  expect_identical(vapply(simulations, `[[`, character(1), "schema_version"),
                   rep("simulation-result-2", 3))
  expect_identical(vapply(simulations, `[[`, character(1), "model"),
                   c("LRM", "HLRM", "CHLRM"))
  expect_true(all(vapply(simulations, function(value)
    all(c("train", "test", "truth", "conditioning", "parameter_source",
          "design_source") %in% names(value)), logical(1L))))
  expect_equal(unname(colMeans(do.call(rbind, hlrm$train$X))[-1]), c(.35, 2),
               tolerance = .15)
  expect_identical(colnames(lrm$train$X),
                   c("intercept", "binary", "numeric"))
  expect_identical(chlrm$truth$allocations, c(1L, 2L, 1L, 2L))
  expect_equal(chlrm$truth$realized_mixture_proportions, c(.5, .5))
  expect_false("population_covariance" %in% names(chlrm$truth))
})

test_that("hierarchical DGP reproduces its population moments", {
  group_count <- 4000L
  population_covariance <- matrix(c(1, .3, .3, .5), 2)
  parameters <- list(
    population_coefficients = c(1, -1),
    population_covariance = population_covariance,
    variance_scale = 2,
    residual_degrees_of_freedom = 6)
  simulation <- simulate_data(
    hierarchical_dgp(hlrm_model(2), parameters),
    simulation_design("independent_normal"),
    simulation_sample(group_count = group_count, observations_per_group = 1),
    seed = 1310)
  expect_equal(colMeans(simulation$truth$group_coefficients), c(1, -1),
               tolerance = .05)
  expect_equal(stats::cov(simulation$truth$group_coefficients),
               population_covariance, tolerance = .06)
  expect_equal(mean(simulation$truth$group_variances), 3, tolerance = .15)
  expect_identical(simulation$parameter_source$local, "hierarchy")
})

test_that("prior-predictive DGP uses the same prior parameterization", {
  prior <- chlrm_prior(
    c(0, 0), diag(2), 5, diag(2), 6, 2, 2, c(2, 3))
  sample <- simulation_sample(
    group_count = 5, observations_per_group = 3,
    test_observations_per_group = 2, test_group_mode = "new",
    test_group_count = 3)
  set.seed(1320)
  before <- rng_snapshot()
  simulation <- simulate_data(
    prior_predictive_dgp(chlrm_model(2, 2), prior),
    simulation_design("multivariate_normal", covariance = matrix(1)),
    sample, seed = 1321)
  expect_identical(rng_snapshot(), before)
  expect_identical(simulation, simulate_data(
    prior_predictive_dgp(chlrm_model(2, 2), prior),
    simulation_design("multivariate_normal", covariance = matrix(1)),
    sample, seed = 1321))
  expect_equal(sum(simulation$truth$mixture_weights), 1)
  expect_equal(dim(simulation$truth$component_coefficients), c(2, 2))
  expect_equal(length(simulation$truth$test_allocations), 3)
  expect_identical(simulation$parameter_source$global, "prior")
  expect_identical(simulation$conditioning$test_units, "new_groups")
  expect_equal(anyDuplicated(names(simulation$truth)), 0L)

  lrm_prior_spec <- lrm_prior(c(0, 0), diag(2), 6, 1.5)
  lrm <- simulate_data(
    prior_predictive_dgp(lrm_model(2), lrm_prior_spec),
    simulation_design("independent_normal"),
    simulation_sample(sample_size = 8), seed = 1322)
  expect_true(all(c("coefficients", "residual_variance") %in%
                  names(lrm$truth)))
  expect_identical(lrm$parameter_source$global, "prior")

  hlrm_prior_spec <- hlrm_prior(
    c(0, 0), diag(2), 5, diag(2), 6, 2, 2)
  hlrm <- simulate_data(
    prior_predictive_dgp(hlrm_model(2), hlrm_prior_spec),
    simulation_design("independent_normal"),
    simulation_sample(group_count = 3, observations_per_group = 4),
    seed = 1323)
  expect_true(all(c("population_coefficients", "population_covariance",
                    "variance_scale", "group_coefficients",
                    "group_variances") %in% names(hlrm$truth)))
  expect_identical(hlrm$parameter_source$global, "prior")
  expect_identical(hlrm$parameter_source$local, "hierarchy")
})

test_that("generated responses agree with the canonical conditional likelihood", {
  sample <- simulation_sample(group_count = 3, observations_per_group = 4)
  parameters <- list(
    component_coefficients = rbind(c(0, 1), c(2, -1)),
    component_variances = c(.5, 1.5), mixture_weights = c(.5, .5),
    allocations = c(1, 2, 1))
  simulation <- simulate_data(
    conditional_dgp(chlrm_model(2, 2), parameters),
    simulation_design("independent_normal"), sample, seed = 1330,
    retain_observation_truth = TRUE)
  group_coefficients <- parameters$component_coefficients[
    parameters$allocations, , drop = FALSE]
  group_variances <- parameters$component_variances[parameters$allocations]
  expected <- sum(unlist(Map(function(y, X, beta, variance)
    stats::dnorm(y, drop(X %*% beta), sqrt(variance), log = TRUE),
    simulation$train$y, simulation$train$X,
    split(group_coefficients, row(group_coefficients)), group_variances)))
  observed <- loglik_total_chlrm(
    simulation$train$y, simulation$train$X,
    group_coefficients, group_variances)
  expect_equal(observed, expected)
  expect_equal(unlist(simulation$truth$train_conditional_mean),
               unlist(Map(function(X, beta) drop(X %*% beta),
                          simulation$train$X,
                          split(group_coefficients, row(group_coefficients)))))
})

test_that("fixed designs and conditional new-group restrictions are explicit", {
  X <- list(cbind(1, c(-1, 0, 1)), cbind(1, c(2, 3)))
  design <- simulation_design("fixed", fixed_train = X)
  sample <- simulation_sample(group_count = 2,
                              observations_per_group = c(3, 2))
  parameters <- list(
    group_coefficients = rbind(c(0, 1), c(1, -1)),
    group_variances = c(1, 1))
  simulation <- simulate_data(
    conditional_dgp(hlrm_model(2), parameters), design, sample, seed = 1340)
  expect_identical(simulation$train$X, X)
  new_group_sample <- simulation_sample(
    group_count = 2, observations_per_group = 3,
    test_observations_per_group = 2, test_group_mode = "new")
  expect_error(simulate_data(
    conditional_dgp(hlrm_model(2), parameters),
    simulation_design("independent_normal"), new_group_sample, seed = 1341),
    "no population distribution")
})

test_that("tidy_truth is uniform and excludes observation noise by default", {
  simulation <- simulate_data(
    conditional_dgp(chlrm_model(2, 2), list(
      component_coefficients = rbind(c(0, 1), c(2, -1)),
      component_variances = c(.5, 1.5), mixture_weights = c(.4, .6),
      allocations = c(1, 2, 1))),
    simulation_design("independent_normal"),
    simulation_sample(group_count = 3, observations_per_group = 4),
    seed = 1350, retain_observation_truth = TRUE)
  compact <- tidy_truth(simulation)
  complete <- tidy_truth(simulation, include_observation_level = TRUE)
  expect_s3_class(compact, "vb_tidy_truth")
  expect_named(compact,
               c("quantity", "index_1", "index_2", "index_3", "value", "source"))
  expect_false(any(grepl("conditional_mean|residuals", compact$quantity)))
  expect_true(any(grepl("conditional_mean|residuals", complete$quantity)))
  expect_equal(nrow(compact[compact$quantity == "component_coefficients", ]), 4)
  expect_identical(unique(compact$source[compact$quantity == "allocations"]),
                   "conditioned")
})

test_that("mixed designs preserve covariance, ordering and binary support", {
  covariance <- matrix(c(1, .65, .65, 2), 2)
  simulation <- simulate_data(
    conditional_dgp(lrm_model(4), list(
      coefficients = c(0, 1, -1, .5), residual_variance = 1)),
    simulation_design(
      "mixed", mean = c(1, -2), covariance = covariance,
      binary_probabilities = .3, binary_position = "last"),
    simulation_sample(sample_size = 5000), seed = 1360)
  X <- simulation$train$X
  expect_equal(colMeans(X[, 2:3, drop = FALSE]), c(1, -2), tolerance = .05)
  expect_equal(stats::cov(X[, 2:3, drop = FALSE]), covariance, tolerance = .08)
  expect_true(all(X[, 4] %in% 0:1))
  expect_equal(mean(X[, 4]), .3, tolerance = .03)
  residual <- simulation$train$y - drop(X %*% c(0, 1, -1, .5))
  expect_equal(mean(residual), 0, tolerance = .04)
  expect_equal(stats::var(residual), 1, tolerance = .05)
})

test_that("degenerate and fixed designs preserve their declared schema", {
  intercept_only <- simulate_data(
    conditional_dgp(lrm_model(1), list(
      coefficients = 2, residual_variance = .5)),
    simulation_design(
      "grid", bounds = matrix(numeric(), 0, 2),
      column_names = "intercept"),
    simulation_sample(sample_size = 5), seed = 1365)
  expect_equal(intercept_only$train$X, matrix(1, 5, 1,
    dimnames = list(NULL, "intercept")))

  fixed <- simulate_data(
    conditional_dgp(lrm_model(2), list(
      coefficients = c(0, 1), residual_variance = 1)),
    simulation_design(
      "fixed", fixed_train = cbind(1, 1:3),
      column_names = c("intercept", "x")),
    simulation_sample(sample_size = 3), seed = 1366)
  expect_identical(colnames(fixed$train$X), c("intercept", "x"))
})

test_that("simulation output budgets and optional observation truth are explicit", {
  generator <- conditional_dgp(lrm_model(2), list(
    coefficients = c(0, 1), residual_variance = 1))
  design <- simulation_design("independent_normal")
  sample <- simulation_sample(sample_size = 100)
  compact <- simulate_data(generator, design, sample, seed = 1370)
  retained <- simulate_data(
    generator, design, sample, seed = 1370,
    retain_observation_truth = TRUE)
  expect_false("train_residuals" %in% names(compact$truth))
  expect_true(all(c("train_conditional_mean", "train_residuals") %in%
                  names(retained$truth)))
  expect_identical(compact$train, retained$train)
  expect_error(simulate_data(
    generator, design, simulation_sample(sample_size = 1e6), seed = 1371,
    max_output_mb = .1), "exceeding max_output_mb")
})
