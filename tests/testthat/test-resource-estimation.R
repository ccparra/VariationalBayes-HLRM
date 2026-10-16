.resource_test_data <- function(seed = 1600) {
  scope <- inference_rng_scope(seed)
  on.exit(scope$restore())
  X <- cbind(intercept = 1, predictor = stats::rnorm(30))
  lrm <- list(X = X, y = drop(X %*% c(1, 2) + stats::rnorm(30, sd = .5)))
  groups <- lapply(seq_len(6), function(j) {
    Xj <- cbind(intercept = 1, predictor = stats::rnorm(8))
    coefficient <- if (j <= 3) c(0, 1) else c(2, -1)
    list(X = Xj, y = drop(Xj %*% coefficient + stats::rnorm(8, sd = .5)))
  })
  list(lrm = lrm, grouped = list(
    X = lapply(groups, `[[`, "X"), y = lapply(groups, `[[`, "y")))
}

.resource_test_priors <- function() {
  list(
    lrm = lrm_prior(c(0, 0), diag(10, 2), 4, 1),
    hlrm = hlrm_prior(c(0, 0), diag(10, 2), 4, diag(2), 4, 2, 2),
    chlrm = chlrm_prior(c(0, 0), diag(10, 2), 4, diag(2), 4, 2, 2,
                        c(1, 1)))
}

test_that("resource estimates expose one uniform preflight contract", {
  data <- .resource_test_data(seed = 1600)
  lrm <- estimate_fit_resources(
    lrm_model(2), data$lrm,
    mcmc_control(warmup = 1, post_warmup_draws = 10,
                 evaluation_draws = 4, max_memory_mb = 100))

  expect_s3_class(lrm, "vb_fit_resource_estimate")
  expect_identical(lrm$schema_version, "fit-resource-estimate-1")
  expect_identical(lrm$model, "LRM")
  expect_identical(lrm$method, "MCMC")
  expect_true(lrm$within_budget)
  expect_gt(lrm$estimated_peak_memory_mb, 0)
  expect_gt(lrm$estimated_evaluation_storage_mb, 0)
  expect_named(lrm$components, c("component", "estimated_mb"))
  expect_output(print(lrm), "within budget")

  scaled <- estimate_fit_resources(
    lrm_model(2), data$lrm,
    cavi_control(max_iterations = 2, evaluation_draws = 2,
      stopping = cavi_stopping(min_iterations = 1), max_memory_mb = 100),
    scaling = data_scaling())
  unscaled <- estimate_fit_resources(
    lrm_model(2), data$lrm,
    cavi_control(max_iterations = 2, evaluation_draws = 2,
      stopping = cavi_stopping(min_iterations = 1), max_memory_mb = 100))
  expect_gt(scaled$estimated_peak_memory_mb,
            unscaled$estimated_peak_memory_mb)
})

test_that("resource estimates increase with retained work and model structure", {
  data <- .resource_test_data(seed = 1601)
  one <- estimate_fit_resources(
    hlrm_model(2), data$grouped,
    cavi_control(max_iterations = 5, evaluation_draws = 2,
      stopping = cavi_stopping(min_iterations = 1), starts = 1))
  three <- estimate_fit_resources(
    hlrm_model(2), data$grouped,
    cavi_control(max_iterations = 5, evaluation_draws = 2,
      stopping = cavi_stopping(min_iterations = 1), starts = 3,
      require_converged = FALSE))
  clustered <- estimate_fit_resources(
    chlrm_model(5, 2), data$grouped,
    cavi_control(max_iterations = 5, evaluation_draws = 2,
      stopping = cavi_stopping(min_iterations = 1), starts = 3,
      require_converged = FALSE))

  expect_gt(three$estimated_peak_memory_mb, one$estimated_peak_memory_mb)
  expect_gt(clustered$estimated_peak_memory_mb,
            three$estimated_peak_memory_mb)

  short <- estimate_fit_resources(
    lrm_model(2), data$lrm,
    mcmc_control(warmup = 0, post_warmup_draws = 10,
                 evaluation_draws = 2))
  long <- estimate_fit_resources(
    lrm_model(2), data$lrm,
    mcmc_control(warmup = 0, post_warmup_draws = 100,
                 evaluation_draws = 2))
  expect_gt(long$estimated_peak_memory_mb, short$estimated_peak_memory_mb)
})

test_that("fit_model enforces the preflight budget for a single run", {
  data <- .resource_test_data(seed = 1602)$lrm
  prior <- .resource_test_priors()$lrm
  control <- cavi_control(
    max_iterations = 2, evaluation_draws = 2,
    stopping = cavi_stopping(min_iterations = 1), max_memory_mb = 1e-12)
  estimate <- estimate_fit_resources(lrm_model(2), data, control)

  expect_false(estimate$within_budget)
  expect_error(
    fit_model(lrm_model(2), data, prior, "cavi", initialization_ols(),
              control, seed = 1602, verbose = FALSE),
    "Estimated peak memory.*exceeds")
})

test_that("successful fits retain their resource estimate", {
  data <- .resource_test_data(seed = 1603)$lrm
  fit <- fit_model(
    lrm_model(2), data, .resource_test_priors()$lrm, "cavi",
    initialization_ols(),
    cavi_control(max_iterations = 2, evaluation_draws = 2,
      stopping = cavi_stopping(min_iterations = 1)),
    seed = 1603, verbose = FALSE)
  estimate <- fit$metadata$resource_estimate
  expect_s3_class(estimate, "vb_fit_resource_estimate")
  expect_true(estimate$within_budget)
  expect_identical(estimate$method, "CAVI")
})
