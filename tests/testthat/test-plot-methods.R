.with_plot_device <- function(code) {
  path <- tempfile(fileext = ".pdf")
  grDevices::pdf(path)
  on.exit({
    grDevices::dev.off()
    unlink(path)
  })
  force(code)
}

test_that("base plot methods cover fitting, prediction, PPC and likelihood", {
  simulation <- .test_lrm_simulation(
    18, c(1, 2), .5, matrix(c(-1, 1), 1, 2), seed = 1801)
  fit <- fit_model(
    lrm_model(2), simulation$train,
    lrm_prior(c(0, 0), 10 * diag(2), 4, 1), "cavi",
    initialization_ols(), cavi_control(8, evaluation_draws = 6),
    seed = 1802, verbose = FALSE)
  prediction <- posterior_predict(
    fit, simulation$train$X, n_draws = 6, seed = 1803)
  check <- posterior_predictive_check(fit, n_draws = 6, seed = 1804)
  draws <- posterior_draws(fit, n_draws = 20, seed = 1806)
  likelihood <- predictive_criteria(
    fit, n_draws = 6, seed = 1805,
    return_log_likelihood = TRUE)$log_likelihood$observation

  plotted <- .with_plot_device(list(
    simulation = plot(simulation),
    convergence = plot(fit, type = "convergence"),
    coefficients = plot(fit, type = "coefficients"),
    densities = plot(draws, parameter = "coefficients"),
    prediction = plot(prediction, observed = simulation$train$y),
    ppc = plot(check),
    log_likelihood = plot(likelihood)))

  expect_s3_class(plotted$simulation, "data.frame")
  expect_named(plotted$convergence, c("step", "value", "metric"))
  expect_true(all(c("estimate", "lower", "upper", "label") %in%
                  names(plotted$coefficients)))
  expect_named(plotted$densities, c("draw", "parameter", "value"))
  expect_equal(nrow(plotted$densities), 40)
  expect_identical(plotted$prediction, prediction$summary)
  expect_named(plotted$ppc, c(
    "replicated", "observed", "posterior_predictive_p_value"))
  expect_equal(nrow(plotted$log_likelihood), nrow(likelihood$values))
})

test_that("CHLRM and orchestration plots expose mixture and timing data", {
  simulation <- .test_chlrm_simulation(
    5, 5, 2, 2,
    component_coefficients = rbind(c(-1, 1), c(2, -1)), seed = 1810)
  fit <- fit_model(
    chlrm_model(2, 2), simulation$train,
    chlrm_prior(c(0, 0), 10 * diag(2), 4, diag(2), 4, 2, 2, c(1, 1)),
    "cavi",
    initialization_kmeans(),
    cavi_control(6, evaluation_draws = 4, starts = 2,
                 require_converged = FALSE),
    seed = 1811, verbose = FALSE)
  allocation <- .with_plot_device(
    plot(fit, type = "responsibilities"))
  timing <- .with_plot_device(plot(fit_orchestration(fit)))

  expect_equal(allocation, responsibilities(fit))
  expect_true(all(c("start", "elapsed_seconds", "selected") %in%
                  names(timing)))
  expect_error(plot(
    fit_model(
      lrm_model(2),
      list(X = cbind(1, 1:4), y = 1:4),
      lrm_prior(c(0, 0), diag(2), 4, 1), "mcmc",
      initialization_from_prior(),
      mcmc_control(warmup = 1, post_warmup_draws = 2,
                   evaluation_draws = 2),
      seed = 1812, verbose = FALSE),
    type = "convergence"), "no convergence history")
})
