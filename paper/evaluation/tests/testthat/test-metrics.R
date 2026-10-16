test_that("ARI is label invariant and agrees with a manual contingency example", {
  truth <- c(1, 1, 2, 2)
  expect_equal(adjusted_rand_index(truth, c("b", "b", "a", "a")), 1)
  # Contingency [[1,1],[1,1]] gives index -1/2.
  expect_equal(adjusted_rand_index(truth, c(1, 2, 1, 2)), -.5)
  expect_equal(adjusted_rand_index(1, "x"), 1)
})

test_that("co-clustering error and normalized entropy match limiting cases", {
  truth <- c(1, 1, 2)
  actual <- outer(truth, truth, `==`) * 1
  expect_equal(co_clustering_error(truth, actual), 0)
  wrong <- matrix(1, 3, 3)
  expect_equal(co_clustering_error(truth, wrong), 2 / 3)
  expect_equal(normalized_classification_entropy(diag(3)), 0)
  expect_equal(normalized_classification_entropy(
    matrix(1 / 3, 2, 3)), 1)
  expect_equal(normalized_classification_entropy(matrix(1, 4, 1)), 0)
})

test_that("recovery and block coverage use the review definitions", {
  scalar <- scalar_recovery_metrics(
    c(1, 2, 3), truth = 2, lower = c(0, 1, 2.5), upper = c(2, 3, 3.5))
  expect_equal(scalar$bias, 0)
  expect_equal(scalar$rmse, sqrt(2 / 3))
  expect_equal(scalar$coverage, 2 / 3)
  expect_equal(scalar$mean_width, 5 / 3)

  block <- block_interval_coverage(
    c(0, 1), rbind(c(-1, 0), c(.1, 0)), rbind(c(1, 2), c(1, .5)))
  expect_equal(block$per_replicate$coverage, c(1, 0))
  expect_equal(block$coverage, .5)
  expect_equal(block$coverage_mcse, .5)
})

test_that("posterior fidelity is zero for the reference itself", {
  reference <- cbind(a = 1:6, b = c(2, 1, 4, 3, 6, 5))
  result <- posterior_fidelity(reference, reference)
  expect_equal(result$mean_distance, 0)
  expect_equal(result$mean_relative_sd_difference, 0)
  expect_equal(result$mean_sd_ratio, 1)
  expect_equal(result$correlation_frobenius, 0)
  expect_equal(result$maximum_off_diagonal_correlation_error, 0)
})

test_that("selection summaries and predictive regret use paired replicates", {
  expect_equal(predictive_regret(c(-3, -1, -2), 3), 1)
  result <- selection_recovery(
    c(3, 3, 2, 4), truth = 3, candidate_values = 1:4,
    regret = c(0, .2, .1, .3))
  expect_equal(result$exact_recovery, .5)
  expect_equal(result$mean_absolute_error, .5)
  expect_equal(result$distribution$count, c(0, 1, 2, 1))
  expect_equal(result$mean_predictive_regret, .15)
})

test_that("metric validators reject incoherent inputs", {
  expect_error(adjusted_rand_index(1:2, 1:3), "equal length")
  expect_error(co_clustering_error(1:3, diag(2)), "m by m")
  expect_error(normalized_classification_entropy(matrix(.2, 2, 2)),
               "row-stochastic")
  expect_error(scalar_recovery_metrics(1:2, 1, lower = 0:1),
               "supplied together")
  expect_error(posterior_fidelity(matrix(1, 3, 1), matrix(1, 3, 1)),
               "positive finite")
  expect_error(predictive_regret(1:3, 4), "identify")
})

test_that("E1 aggregation reports means and Monte Carlo standard errors", {
  metrics <- data.frame(
    scenario = rep("S1", 4), method = rep(c("MCMC", "CAVI"), each = 2),
    ari = c(1, .8, .7, .5), cpe = c(0, .1, .2, .3),
    normalized_entropy = c(.1, .2, .3, .4), mlpd = c(-1, -3, -2, -4),
    coefficient_rmse_per_dataset = c(1, 2, 3, 4),
    coefficient_coverage_per_dataset = c(.9, 1, .8, .9),
    variance_rmse_per_dataset = c(1, 2, 3, 4),
    variance_coverage_per_dataset = c(.9, 1, .8, .9),
    weight_rmse_per_dataset = c(.1, .2, .3, .4),
    test_rmse = 1:4, test_mae = 1:4,
    predictive_coverage = c(1, .5, .5, 0),
    mean_interval_width = c(2, 2, 3, 3), fit_seconds = c(4, 6, 1, 3))
  result <- aggregate_e1_metrics(metrics)
  expect_equal(result$method, c("MCMC", "CAVI"))
  expect_equal(result$ari_mean, c(.9, .6))
  expect_equal(result$mlpd_mean, c(-2, -3))
  expect_equal(result$fit_seconds_mcse, c(1, 1))
})

test_that("aligned parameter aggregation preserves coordinates and missing means", {
  recovery <- data.frame(
    scenario = "S1", method = "CAVI", replicate = rep(1:2, each = 2),
    parameter = rep(c("component_coefficients", "component_variances"), 2),
    index_1 = 1L, index_2 = c(1L, NA, 1L, NA), truth = c(2, 1, 2, 1),
    estimate = c(1, NA, 3, NA), sd = c(.5, NA, .5, NA),
    error = c(-1, NA, 1, NA), covered = c(TRUE, TRUE, FALSE, TRUE),
    interval_width = c(2, 3, 2, 3))
  result <- aggregate_e1_parameter_recovery(recovery)
  coefficient <- result[result$parameter == "component_coefficients", ]
  variance <- result[result$parameter == "component_variances", ]
  expect_equal(coefficient$bias, 0)
  expect_equal(coefficient$rmse, 1)
  expect_equal(coefficient$coverage, .5)
  expect_equal(coefficient$replicates, 2)
  expect_equal(variance$finite_posterior_means, 0)
  expect_true(is.na(variance$bias))
  expect_equal(variance$coverage, 1)
})

test_that("truth alignment recovers a global component permutation", {
  truth_coefficients <- rbind(c(-2, 1), c(0, -1), c(3, .5))
  truth_variances <- c(.5, 1, 2)
  permutation <- c(3, 1, 2)
  aligned <- match_components_to_truth(
    truth_coefficients[permutation, ], truth_variances[permutation],
    truth_coefficients, truth_variances)
  expect_identical(aligned$estimated_index_for_truth, c(2L, 3L, 1L))
  expect_equal(aligned$total_cost, 0)
})

test_that("clustering stability averages all pairs of completed starts", {
  classifications <- rbind(
    data.frame(start = 1, group = 1:4, component = c(1, 1, 2, 2),
               assignment_probability = 1),
    data.frame(start = 2, group = 1:4, component = c(2, 2, 1, 1),
               assignment_probability = 1),
    data.frame(start = 3, group = 1:4, component = c(1, 2, 1, 2),
               assignment_probability = 1))
  orchestration <- structure(list(
    model = "CHLRM", type = "multistart",
    start_classifications = classifications),
    class = c("vb_fit_orchestration", "list"))
  result <- clustering_stability(orchestration)
  expect_equal(result$pairs, 3)
  expect_equal(result$stability, mean(c(1, -.5, -.5)))
})
