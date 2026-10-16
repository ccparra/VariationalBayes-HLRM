.scoring_data <- function(seed = 1300) {
  scope <- inference_rng_scope(seed); on.exit(scope$restore())
  X <- cbind(intercept = 1, x = stats::rnorm(24))
  lrm <- list(X = X, y = drop(X %*% c(3, 1.5) + stats::rnorm(24)))
  groups <- lapply(seq_len(5), function(group) {
    Xj <- cbind(intercept = 1, x = stats::rnorm(7))
    list(X = Xj, y = drop(Xj %*% c(3, 1.5) + stats::rnorm(7)))
  })
  list(
    lrm = lrm,
    grouped = list(X = lapply(groups, `[[`, "X"),
                   y = lapply(groups, `[[`, "y")))
}

test_that("LRM held-out scoring equals its posterior draw calculation", {
  data <- .scoring_data()
  fit <- fit_model(
    lrm_model(2), data$lrm,
    lrm_prior(c(0, 0), diag(20, 2), 4, 2), "cavi",
    initialization_ols(), cavi_control(4, evaluation_draws = 4),
    scaling = data_scaling(), seed = 1301, verbose = FALSE)
  held_out <- list(X = data$lrm$X[1:4, , drop = FALSE],
                   y = data$lrm$y[1:4])
  sampled <- posterior_draws(fit, n_draws = 4, seed = 1302)
  expected <- vapply(seq_len(4), function(draw)
    stats::dnorm(
      held_out$y,
      drop(held_out$X %*% sampled$parameters$coefficients[draw, ]),
      sqrt(sampled$parameters$residual_variance[draw]), log = TRUE),
    numeric(4))
  likelihood <- predictive_log_likelihood(
    fit, held_out, n_draws = 4, seed = 1302)
  expect_s3_class(likelihood, "vb_log_likelihood")
  expect_equal(likelihood$values, expected)
  expect_identical(likelihood$metadata$data_role, "held_out")
  expect_identical(likelihood$metadata$density_scale, "input")

  score <- predictive_score(
    fit, held_out, n_draws = 4, seed = 1302,
    return_log_likelihood = TRUE)
  pointwise <- apply(expected, 1L, logmeanexp)
  expect_s3_class(score, "vb_predictive_score")
  expect_equal(score$statistics$elpd, sum(pointwise))
  expect_equal(score$statistics$mean_log_predictive_density,
               mean(pointwise))
  expect_equal(score$pointwise$log_predictive_density, pointwise)
  expect_s3_class(score$log_likelihood, "vb_log_likelihood")
})

test_that("HLRM scoring distinguishes observed from new groups", {
  data <- .scoring_data()
  fit <- fit_model(
    hlrm_model(2), data$grouped,
    hlrm_prior(c(0, 0), diag(20, 2), 4, diag(2), 4, 2, 2),
    "cavi", initialization_ols(), cavi_control(4, evaluation_draws = 3),
    seed = 1310, verbose = FALSE)
  held_out <- list(
    X = lapply(data$grouped$X[1:2], function(X) X[1:2, , drop = FALSE]),
    y = lapply(data$grouped$y[1:2], function(y) y[1:2]))
  observed <- predictive_log_likelihood(
    fit, held_out, predictive_target = "observed_groups",
    group_ids = 1:2, n_draws = 3, seed = 1311)
  expect_identical(observed$unit, "observation")
  expect_identical(observed$conditioning, "parameter_conditional")
  expect_equal(dim(observed$values), c(4, 3))

  new_group <- predictive_log_likelihood(
    fit, held_out, predictive_target = "new_groups",
    n_draws = 3, seed = 1311)
  expect_identical(new_group$unit, "group")
  expect_identical(new_group$conditioning, "group_marginal")
  expect_equal(dim(new_group$values), c(2, 3))
  expect_error(predictive_log_likelihood(
    fit, held_out, predictive_target = "new_groups",
    predictive_unit = "observation", n_draws = 3),
    "requires predictive_unit")
  expect_error(predictive_log_likelihood(
    fit, held_out, predictive_target = "observed_groups",
    n_draws = 3), "requires one valid")
})

test_that("CHLRM held-out scoring uses observed allocations or new-group mixture", {
  data <- .scoring_data()
  fit <- fit_model(
    chlrm_model(2, 2), data$grouped,
    chlrm_prior(c(0, 0), diag(20, 2), 4, diag(2), 4, 2, 2,
                c(1, 1)),
    "cavi", initialization_kmeans(),
    cavi_control(4, evaluation_draws = 3), seed = 1320, verbose = FALSE)
  held_out <- list(
    X = lapply(data$grouped$X[1:2], function(X) X[1:3, , drop = FALSE]),
    y = lapply(data$grouped$y[1:2], function(y) y[1:3]))
  observed <- predictive_log_likelihood(
    fit, held_out, predictive_target = "observed_groups",
    predictive_unit = "group", group_ids = 1:2,
    n_draws = 3, seed = 1321)
  expect_identical(observed$conditioning, "parameter_conditional")
  expect_equal(dim(observed$values), c(2, 3))
  expect_true(all(is.finite(observed$values)))

  new_group <- predictive_log_likelihood(
    fit, held_out, predictive_target = "new_groups",
    n_draws = 3, seed = 1321)
  expect_identical(new_group$conditioning, "group_marginal")
  expect_equal(dim(new_group$values), c(2, 3))
  expect_true(all(is.finite(new_group$values)))
})

test_that("LRM grouped held-out units are explicit", {
  data <- .scoring_data()
  fit <- fit_model(
    lrm_model(2), data$lrm,
    lrm_prior(c(0, 0), diag(20, 2), 4, 2), "cavi",
    initialization_ols(), cavi_control(3, evaluation_draws = 3),
    seed = 1330, verbose = FALSE)
  labels <- rep(c("a", "b"), each = 3)
  held_out <- list(X = data$lrm$X[1:6, , drop = FALSE],
                   y = data$lrm$y[1:6])
  grouped <- predictive_log_likelihood(
    fit, held_out, predictive_unit = "group", groups = labels,
    n_draws = 3, seed = 1331)
  expect_equal(dim(grouped$values), c(2, 3))
  observation <- predictive_log_likelihood(
    fit, held_out, n_draws = 3, seed = 1331)
  expect_equal(unname(grouped$values),
               unname(rowsum(
                 observation$values, match(labels, unique(labels)),
                 reorder = FALSE)))
})

test_that("predictive performance combines test error intervals and log score", {
  data <- .scoring_data()
  fit <- fit_model(
    lrm_model(2), data$lrm,
    lrm_prior(c(0, 0), diag(20, 2), 4, 2), "cavi",
    initialization_ols(), cavi_control(4, evaluation_draws = 4),
    scaling = data_scaling(), seed = 1340, verbose = FALSE)
  held_out <- list(X = data$lrm$X[1:6, , drop = FALSE],
                   y = data$lrm$y[1:6])
  performance <- predictive_performance(
    fit, held_out, n_draws = 4, level = .8, seed = 1341)
  expect_s3_class(performance, "vb_predictive_performance")
  expect_identical(performance$schema_version, "predictive-performance-1")
  expect_equal(performance$statistics$rmse,
               sqrt(mean(performance$pointwise$squared_error)))
  expect_equal(performance$statistics$mae,
               mean(performance$pointwise$absolute_error))
  expect_equal(performance$statistics$empirical_coverage,
               mean(performance$pointwise$covered))
  expect_equal(performance$statistics$mean_interval_width,
               mean(performance$pointwise$interval_width))
  expect_equal(nrow(performance$pointwise), 6)
  expect_equal(nrow(performance$log_score_pointwise), 6)
  expect_identical(performance$metadata$prediction_scale, "input")
  expect_error(predictive_performance(fit, held_out), "seed is required")
})

test_that("grouped predictive performance defaults transparently to new groups", {
  data <- .scoring_data()
  fit <- fit_model(
    chlrm_model(2, 2), data$grouped,
    chlrm_prior(c(0, 0), diag(20, 2), 4, diag(2), 4, 2, 2,
                c(1, 1)),
    "cavi", initialization_from_prior(),
    cavi_control(3, evaluation_draws = 3), seed = 1350, verbose = FALSE)
  held_out <- list(
    X = lapply(data$grouped$X[1:2], function(X) X[1:2, , drop = FALSE]),
    y = lapply(data$grouped$y[1:2], function(y) y[1:2]))
  performance <- predictive_performance(
    fit, held_out, n_draws = 3, seed = 1351)
  expect_identical(performance$predictive_target, "new_groups")
  expect_true(all(performance$pointwise$group_status == "new"))
  expect_equal(performance$statistics$observations, 4)
  expect_equal(performance$statistics$log_score_units, 2)
  expect_equal(nrow(performance$log_score_pointwise), 2)
})
