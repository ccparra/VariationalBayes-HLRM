test_that("SBC ranks randomize ties and preserve caller RNG", {
  set.seed(2010)
  before <- .test_rng_state()
  truth <- c(theta = 1)
  draws <- matrix(c(.5, 1, 1, 2), ncol = 1,
                  dimnames = list(NULL, "theta"))
  rank <- sbc_ranks(truth, draws, seed = 2011)
  expect_identical(.test_rng_state(), before)
  expect_true(rank$rank %in% 1:3)
  expect_equal(rank$ties, 2)
  expect_true(rank$normalized_rank > 0 && rank$normalized_rank < 1)
})

test_that("CHLRM SBC inputs are invariant to component permutations", {
  simulation <- simulate_data(
    conditional_dgp(chlrm_model(3, 2), list(
      component_coefficients = rbind(c(-2, 1), c(0, -1), c(3, .5)),
      component_variances = c(.5, 1, 2),
      mixture_weights = c(.2, .3, .5))),
    simulation_design("independent_normal"),
    simulation_sample(group_count = 8, observations_per_group = 3),
    seed = 2020)
  coefficients <- array(rep(simulation$truth$component_coefficients, 4),
                        dim = c(4, 3, 2))
  # Rebuild explicitly because array recycling does not define draw ordering.
  for (draw in 1:4) coefficients[draw, , ] <-
    simulation$truth$component_coefficients
  variances <- matrix(rep(simulation$truth$component_variances, each = 4), 4, 3)
  weights <- matrix(rep(simulation$truth$mixture_weights, each = 4), 4, 3)
  posterior <- structure(list(
    schema_version = "posterior-draws-1", model = "CHLRM", method = "CAVI",
    parameters = list(
      population_coefficients = matrix(0, 4, 2),
      population_covariance = array(diag(2), c(4, 2, 2)),
      variance_scale = rep(1, 4),
      component_coefficients = coefficients,
      component_variances = variances,
      mixture_weights = weights)),
    class = c("vb_posterior_draws", "list"))
  baseline <- chlrm_invariant_sbc_inputs(simulation, posterior)
  permutation <- c(3, 1, 2)
  posterior$parameters$component_coefficients <-
    coefficients[, permutation, , drop = FALSE]
  posterior$parameters$component_variances <- variances[, permutation, drop = FALSE]
  posterior$parameters$mixture_weights <- weights[, permutation, drop = FALSE]
  permuted <- chlrm_invariant_sbc_inputs(simulation, posterior)
  expect_equal(permuted, baseline)
  expect_identical(baseline$metadata$estimand, "label_invariant")
})

test_that("prior-predictive simulation, fitting and SBC ranks connect end to end", {
  model <- lrm_model(2)
  prior <- lrm_prior(c(0, 0), diag(2), 6, 1)
  simulation <- simulate_data(
    prior_predictive_dgp(model, prior),
    simulation_design("independent_normal"),
    simulation_sample(sample_size = 20), seed = 2030)
  fit <- fit_model(
    model, simulation$train, prior, "cavi", initialization_from_prior(),
    cavi_control(
      max_iterations = 5, evaluation_draws = 20,
      stopping = cavi_stopping(
        objective_tolerance = 1e20, parameter_tolerance = 1e20,
        min_iterations = 2, patience = 1)),
    seed = 2031, verbose = FALSE)
  posterior <- posterior_draws(fit, n_draws = 20, seed = 2032)
  truth <- simulation$truth$coefficients
  names(truth) <- paste0("coefficients[", seq_along(truth), "]")
  draws <- posterior$parameters$coefficients
  colnames(draws) <- names(truth)
  ranks <- sbc_ranks(truth, draws, seed = 2033)
  expect_equal(nrow(ranks), 2)
  expect_true(all(ranks$rank >= 0 & ranks$rank <= 20))
  expect_identical(fit$specifications$prior, prior)
})

test_that("review coefficient geometry has the requested minimum separation", {
  for (K in c(2L, 3L, 6L, 10L)) {
    coefficients <- make_review_component_coefficients(K, 10L, 1.5)
    expect_equal(min(as.numeric(stats::dist(coefficients))), 1.5,
                 tolerance = 1e-12)
    expect_equal(coefficients[, 6:10, drop = FALSE], matrix(0, K, 5))
  }
})

test_that("E1 catalog reproduces all nine review scenarios", {
  catalog <- e1_scenario_catalog()
  expect_identical(catalog$scenario,
                   c("S1", "S2", "S3", "S4", "S5", "S6",
                     "S7a", "S7b", "S7c"))
  expect_equal(catalog$component_count, c(3, 3, 3, 3, 3, 3, 3, 6, 10))
  expect_equal(catalog$group_count, c(50, 50, 50, 50, 80, 100, 90, 180, 300))
  expect_equal(catalog$coefficient_count, c(5, 5, 5, 5, 5, 20, 10, 10, 10))

  set.seed(2200); before <- .test_rng_state()
  scenarios <- lapply(catalog$scenario, scenario_e1, replicate = 2L)
  expect_identical(.test_rng_state(), before)
  expect_true(all(vapply(scenarios, inherits, logical(1L),
                         "paper_e1_scenario")))
  expect_true(all(vapply(scenarios, function(x)
    identical(x$sample$test_group_mode, "same"), logical(1L))))
  expect_true(all(vapply(scenarios, function(x)
    all(x$sample$test_observations_per_group == 20L), logical(1L))))
  expect_true(all(scenarios[[4L]]$sample$observations_per_group %in% 5:50))
  expect_true(all(scenarios[[5L]]$sample$observations_per_group %in% 5:50))
  expect_equal(scenarios[[3L]]$design$covariance,
               .8^abs(outer(1:4, 1:4, `-`)))
  expect_equal(scenarios[[6L]]$design$covariance,
               .5^abs(outer(1:19, 1:19, `-`)))
  expect_equal(scenarios[[5L]]$generator$parameters$mixture_weights,
               c(.75, .20, .05))
  expect_equal(scenarios[[9L]]$generator$parameters$mixture_weights,
               rep(.1, 10))
})

test_that("an E1 scenario generates independent test responses for the same groups", {
  scenario <- scenario_e1("S1", replicate = 3L)
  simulation <- simulate_data(
    scenario$generator, scenario$design, scenario$sample,
    seed = scenario$seed)
  expect_equal(length(simulation$train$y), 50)
  expect_equal(lengths(simulation$train$y), rep(20, 50))
  expect_equal(lengths(simulation$test$y), rep(20, 50))
  expect_identical(simulation$conditioning$test_units, "same_groups")
  expect_false("test_allocations" %in% names(simulation$truth))
  expect_false(identical(simulation$train$X, simulation$test$X))
  expect_identical(scenario$prior$mixture_concentration, rep(1 / 3, 3))
})

test_that("E1 task expansion is deterministic and profile sizes are explicit", {
  smoke <- expand_e1_tasks("smoke")
  pilot <- expand_e1_tasks("pilot")
  full <- expand_e1_tasks("full")
  expect_equal(nrow(smoke), 6)
  expect_equal(nrow(pilot), 9 * 5 * 3)
  expect_equal(nrow(full), 9 * 200 * 3)
  expect_identical(anyDuplicated(full$task_id), 0L)
  expect_identical(anyDuplicated(unique(full[c(
    "scenario", "replicate", "dataset_seed")])), 0L)
  expect_identical(smoke, expand_e1_tasks("smoke"))
  expect_true(all(!smoke$confirmatory))
  expect_true(all(full$confirmatory))
})

test_that("one smoke task crosses simulation, inference and evaluation contracts", {
  task <- expand_e1_tasks("smoke")[3L, , drop = FALSE]
  result <- fit_e1_task(task)
  expect_s3_class(result, "paper_task_result")
  expect_identical(result$schema_version, "paper-task-result-4")
  expect_identical(result$fit$method, task$method)
  expect_s3_class(result$parameters, "vb_tidy_parameters")
  expect_s3_class(result$evaluation$in_sample, "vb_in_sample_fit")
  expect_s3_class(result$evaluation$predictive, "vb_predictive_criteria")
  expect_s3_class(result$evaluation$ppc, "vb_posterior_predictive_check")
  expect_s3_class(result$evaluation$held_out, "vb_predictive_performance")
  expect_true(result$allocations$metrics$adjusted_rand_index >= -1)
  expect_true(result$allocations$metrics$adjusted_rand_index <= 1)
  expect_true(result$allocations$metrics$co_clustering_error >= 0)
  expect_true(result$allocations$metrics$co_clustering_error <= 1)
  expect_true(result$allocations$metrics$normalized_entropy >= 0)
  expect_true(result$allocations$metrics$normalized_entropy <= 1)
  expect_s3_class(result$recovery, "paper_chlrm_recovery")
  expect_setequal(result$recovery$summary$parameter,
                  c("component_coefficients", "component_variances",
                    "mixture_weights"))
  expect_s3_class(result$fit$metadata$resource_estimate,
                  "vb_fit_resource_estimate")
  expect_true(result$fit$metadata$resource_estimate$within_budget)
  expect_equal(unname(rowSums(result$allocations$co_clustering > -Inf)),
               rep(50, 50))
  metric <- e1_task_metric_row(result)
  expect_identical(metric$task_id, task$task_id)
  expect_equal(metric$mlpd,
               result$evaluation$held_out$statistics$mean_log_predictive_density)
})

test_that("resource summaries keep runtime, memory and disk distinct", {
  manifest <- data.frame(
    method = c("MCMC", "MCMC", "CAVI"),
    fit_seconds = c(2, 4, 1),
    estimated_peak_memory_mb = c(10, 12, 5),
    estimated_evaluation_storage_mb = c(3, 4, 2),
    result_file_mb = c(.2, .3, .1))
  summary <- summarize_e1_resources(manifest)
  expect_equal(summary$method, c("MCMC", "CAVI"))
  expect_equal(summary$mean_fit_seconds, c(3, 1))
  expect_equal(summary$max_estimated_peak_memory_mb, c(12, 5))
  expect_equal(summary$total_result_file_mb, c(.5, .1))
})

test_that("E2 through E5 protocols reproduce the review targets", {
  e2 <- e2_protocol("full")
  expect_identical(e2$data_generating$minimum_distance, c(1.5, .5))
  expect_equal(length(e2$execution$replicates), 100)
  expect_identical(e2$inference$candidate_components, 1:8)
  expect_identical(e2$inference$cavi_starts, 10L)

  e3 <- e3_protocol("full")
  expect_equal(nrow(e3$full_grid), 36)
  expect_equal(length(e3$target_seeds), 20)
  expect_equal(nrow(e3$execution$grid) * length(e3$execution$seeds), 720)

  e4 <- e4_protocol("full")
  expect_identical(e4$nested_group_counts,
                   c(50L, 200L, 1000L, 5000L))
  expect_equal(length(e4$execution$replicates), 20)

  e5 <- e5_protocol("full")
  expect_identical(e5$selection$initial_components, 1:12)
  expect_identical(e5$selection$cavi_starts, 20L)
  expect_identical(e5$overfitted_mcmc_components, 23L)
  expect_identical(e5$waic_conditioning, "group_marginal")
})

test_that("smoke protocols remain non-confirmatory and computationally bounded", {
  expect_equal(e2_protocol("smoke")$execution$candidate_components, c(2, 3))
  expect_equal(nrow(e3_protocol("smoke")$execution$grid), 1)
  expect_equal(e4_protocol("smoke")$execution$nested_group_counts, c(10, 20))
  expect_equal(e5_protocol("smoke")$execution$evaluation_draws, 4)
})

test_that("E2 smoke tasks cover fixed K fits and overfitted MCMC", {
  tasks <- expand_e2_tasks("smoke")
  expect_equal(nrow(tasks), 10)
  expect_equal(sum(tasks$task_type == "fixed_k"), 8)
  expect_equal(sum(tasks$task_type == "overfitted_mcmc"), 2)
  expect_identical(anyDuplicated(tasks$task_id), 0L)
  scenario <- scenario_e2(1.5, profile = "smoke")
  expect_equal(scenario$sample$group_count, 12)
  expect_equal(scenario$sample$observations_per_group, rep(4, 12))
})

test_that("E3 expansion separates one CAVI reference from the SVI factorial", {
  smoke <- expand_e3_tasks("smoke")
  pilot <- expand_e3_tasks("pilot")
  full <- expand_e3_tasks("full")
  expect_equal(nrow(smoke), 2)
  expect_equal(sum(smoke$method == "CAVI"), 1)
  expect_equal(nrow(pilot), 109)
  expect_equal(nrow(full), 721)
  expect_equal(sum(full$method == "SVI"), 720)
  scenario <- scenario_e3("smoke")
  expect_equal(scenario$sample$group_count, 30)
  expect_true(all(scenario$sample$observations_per_group %in% 10:30))
})

test_that("E3 controls materialize the review stopping rules uniformly", {
  tasks <- expand_e3_tasks("pilot")
  cavi_task <- tasks[tasks$method == "CAVI", , drop = FALSE]
  svi_task <- tasks[tasks$method == "SVI", , drop = FALSE][1L, , drop = FALSE]
  cavi <- e3_inference_control(cavi_task, 1000L)
  svi <- e3_inference_control(svi_task, 1000L)

  expect_identical(cavi$method, "CAVI")
  expect_identical(cavi$max_iterations, 10000L)
  expect_identical(cavi$starts, 5L)
  expect_equal(cavi$stopping$objective_tolerance, 1e-8)
  expect_identical(cavi$stopping$patience, 5L)
  expect_identical(svi$method, "SVI")
  expect_identical(svi$max_iterations, 50000L)
  expect_identical(svi$starts, 1L)
  expect_identical(svi$stopping$min_iterations, 1000L)
  expect_identical(svi$stopping$check_every, 50L)
  expect_identical(svi$stopping$smoothing_window, 5L)
  expect_equal(svi$learning_rate$offset, svi_task$offset)
  expect_equal(svi$learning_rate$exponent, svi_task$exponent)
})

test_that("E3 pilot screen is bounded and preserves factorial corner cases", {
  screen <- expand_e3_pilot_stage("screen")
  factorial <- expand_e3_pilot_stage("factorial")
  expect_equal(nrow(screen), 10)
  expect_equal(sum(screen$method == "CAVI"), 1)
  expect_equal(sum(screen$method == "SVI"), 9)
  expect_equal(nrow(factorial), 109)
  svi <- screen[screen$method == "SVI", , drop = FALSE]
  expect_true(all(c(25L, 200L) %in% svi$minibatch_size))
  expect_true(all(c(.6, .8) %in% svi$exponent))
  expect_true(all(c(1, 100) %in% svi$offset))
  expect_true(any(svi$minibatch_size == 100L & svi$exponent == .7 &
                  svi$offset == 10))
})

test_that("E3 pilot and full fits remain guarded outside frozen runners", {
  pilot_task <- expand_e3_pilot_stage("screen")[1L, , drop = FALSE]
  full_task <- expand_e3_tasks("full")[1L, , drop = FALSE]
  expect_error(fit_e3_task(pilot_task), "validated frozen configuration")
  expect_error(fit_e3_task(full_task, pilot_authorized = TRUE),
               "full tasks remain blocked")
  expect_error(validate_e3_frozen_configuration(list()),
               "freeze_e3_configuration")
})

test_that("E4 smoke uses paired nested datasets", {
  tasks <- expand_e4_tasks("smoke")
  expect_equal(nrow(tasks), 5)
  expect_equal(sum(tasks$method == "MCMC"), 1)
  expect_identical(sort(unique(tasks$group_count)), c(10L, 20L))
  expect_equal(length(unique(tasks$dataset_seed)), 1)
  scenario <- scenario_e4("smoke")
  expect_equal(scenario$sample$group_count, 20)
  expect_true(all(scenario$sample$observations_per_group %in% 10:30))
})

test_that("E5 smoke covers fixed K and overfitted MCMC tasks", {
  tasks <- expand_e5_tasks("smoke")
  expect_equal(nrow(tasks), 7)
  expect_equal(sum(tasks$task_type == "fixed_k"), 6)
  expect_equal(sum(tasks$task_type == "overfitted_mcmc"), 1)
  expect_equal(e5_protocol("smoke")$execution$overfitted_mcmc_components, 4)
})
