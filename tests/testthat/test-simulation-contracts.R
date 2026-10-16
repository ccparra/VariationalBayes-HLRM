test_that("all simulators share one versioned public container", {
  lrm <- .test_lrm_simulation(10, c(1, -0.5), 1, c(-2, 2), seed = 1201)
  hlrm <- .test_hlrm_hierarchy(rep(5, 3), c(0, 1), diag(2), 1,
                                  seed = 1202)
  chlrm <- .test_chlrm_simulation(4, 6, 2, 5, seed = 1203)

  simulations <- list(lrm, hlrm, chlrm)
  expect_true(all(vapply(simulations, inherits, logical(1),
                         "vb_simulation_result")))
  expect_identical(vapply(simulations, `[[`, character(1), "schema_version"),
                   rep("simulation-result-2", 3))
  expect_identical(vapply(simulations, `[[`, character(1), "model"),
                   c("LRM", "HLRM", "CHLRM"))
  expect_true(all(vapply(simulations, function(x)
    all(c("train", "truth", "configuration") %in% names(x)), logical(1))))
  expect_false(any(vapply(simulations, function(x)
    any(c("real_parameters", "cluster_info") %in% names(x)), logical(1))))
  expect_false("test" %in% names(lrm))
  expect_false("test" %in% names(hlrm))
  expect_true("test" %in% names(chlrm))
})

test_that("simulation truth uses model semantics instead of implementation symbols", {
  lrm <- .test_lrm_simulation(8, c(1, 2), .5, c(-1, 1), seed = 1210)
  expect_true(all(c("coefficients", "residual_variance") %in%
                  names(lrm$truth)))
  expect_false(any(c("train_conditional_mean", "train_residuals") %in%
                   names(lrm$truth)))

  hlrm <- .test_hlrm_hierarchy(rep(4, 2), c(0, 1), diag(2), 1,
                                  seed = 1211)
  expect_true(all(c("population_coefficients", "population_covariance",
                    "variance_scale", "residual_degrees_of_freedom",
                    "group_coefficients", "group_variances") %in%
                  names(hlrm$truth)))

  chlrm <- .test_chlrm_hierarchy(rep(4, 5), 2, c(0, 1), diag(2), 1,
                                    seed = 1212)
  expect_true(all(c("population_coefficients", "population_covariance",
                    "variance_scale", "residual_degrees_of_freedom",
                    "component_coefficients", "component_variances",
                    "mixture_weights", "allocations") %in%
                  names(chlrm$truth)))
})

test_that("seeded simulation is reproducible without changing caller RNG", {
  set.seed(1220)
  before <- rng_snapshot()
  a <- .test_lrm_simulation(9, c(1, .2), 1, c(-1, 1), seed = 1221)
  expect_identical(rng_snapshot(), before)
  expect_identical(a,
    .test_lrm_simulation(9, c(1, .2), 1, c(-1, 1), seed = 1221))

  b <- .test_hlrm_simulation(rep(4, 2), rbind(c(0, 1), c(1, -1)), c(1, 2),
                     seed = 1222)
  expect_identical(rng_snapshot(), before)
  expect_identical(b,
    .test_hlrm_simulation(rep(4, 2), rbind(c(0, 1), c(1, -1)), c(1, 2),
                   seed = 1222))
})

test_that("conditional HLRM covariates are explicit rather than group-indexed", {
  coefficients <- rbind(c(0, 1, -1), c(1, .5, .2))
  simulation <- .test_hlrm_simulation(
    group_sizes = c(1000, 1000), group_coefficients = coefficients,
    group_variances = c(1, 1), binary_probabilities = .35,
    covariate_mean = 2, covariate_sd = .5, seed = 1225)
  X <- do.call(rbind, simulation$train$X)
  expect_true(all(X[, 2] %in% 0:1))
  expect_equal(mean(X[, 2]), .35, tolerance = .04)
  expect_equal(mean(X[, 3]), 2, tolerance = .04)
  expect_equal(simulation$configuration$design$binary_probabilities, .35)
  expect_error(.test_hlrm_simulation(c(4, 4), coefficients, c(1, 1),
                              binary_probabilities = 1),
               "strictly between")
})

test_that("historical simulator APIs are not left callable", {
  expect_false(exists("simulate_lrm", mode = "function"))
  expect_false(exists("simulate_hlrm", mode = "function"))
  expect_false(exists("simulate_hlrm_hierarchy", mode = "function"))
  expect_false(exists("simulate_chlrm_protocol", mode = "function"))
  expect_false(exists("simulate_chlrm_hierarchy", mode = "function"))
  expect_false(exists("sim_lrm", mode = "function"))
  expect_false(exists("sim_hlrm", mode = "function"))
  expect_false(exists("sim_hlrm_hierarchical", mode = "function"))
  expect_false(exists("sim_chlrm", mode = "function"))
  expect_false(exists("sim_chlrm_hierarchical", mode = "function"))
})
