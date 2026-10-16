.cv_data <- function(seed = 1400) {
  scope <- inference_rng_scope(seed); on.exit(scope$restore())
  X <- cbind(intercept = 1, x = stats::rnorm(24))
  lrm <- list(X = X, y = drop(X %*% c(4, -1) + stats::rnorm(24)))
  groups <- lapply(seq_len(6), function(group) {
    Xj <- cbind(intercept = 1, x = stats::rnorm(6))
    list(X = Xj, y = drop(Xj %*% c(4, -1) + stats::rnorm(6)))
  })
  list(lrm = lrm,
       grouped = list(X = lapply(groups, `[[`, "X"),
                      y = lapply(groups, `[[`, "y")))
}

test_that("cross-validation plans are reproducible and preserve every unit", {
  data <- .cv_data()
  set.seed(77); before <- .Random.seed
  first <- cross_validation_plan(
    data$lrm, "observation_kfold", folds = 4, seed = 1401)
  second <- cross_validation_plan(
    data$lrm, "observation_kfold", folds = 4, seed = 1401)
  expect_identical(.Random.seed, before)
  expect_identical(first$assignments, second$assignments)
  expect_equal(tabulate(first$assignments, 4), rep(6, 4))
  held_out <- unlist(lapply(seq_len(4), function(fold)
    cross_validation_split(first, data$lrm, fold)$scoring$original_observation))
  expect_equal(sort(held_out), seq_len(24))

  within <- cross_validation_plan(
    data$grouped, "within_group_kfold", folds = 3, seed = 1402)
  expect_identical(within$predictive_target, "observed_groups")
  expect_true(all(vapply(within$assignments, function(x)
    identical(sort(unique(x)), 1:3), logical(1L))))
  split <- cross_validation_split(within, data$grouped, 1)
  expect_equal(lengths(split$training$y), rep(4L, 6))
  expect_equal(lengths(split$assessment$y), rep(2L, 6))

  logo <- cross_validation_plan(
    data$grouped, "leave_one_group_out", seed = 1403)
  expect_equal(logo$fold_count, 6L)
  expect_identical(logo$predictive_target, "new_groups")
  expect_false("observation_groups" %in% names(logo))
})

test_that("ungrouped group folds retain labels and aggregate scoring units", {
  data <- .cv_data()
  labels <- rep(letters[1:6], each = 4)
  plan <- cross_validation_plan(
    data$lrm, "group_kfold", folds = 3, groups = labels, seed = 1410)
  split <- cross_validation_split(plan, data$lrm, 1)
  expect_true(all(split$scoring$groups %in%
    plan$assignments$group[plan$assignments$fold == 1]))
  expect_equal(nrow(split$training$X) + nrow(split$assessment$X), 24)
})

test_that("bounded LRM CV learns scaling independently inside every fold", {
  data <- .cv_data()
  plan <- cross_validation_plan(
    data$lrm, "observation_kfold", folds = 3, seed = 1420)
  result <- cross_validate_model(
    lrm_model(2), data$lrm,
    lrm_prior(c(0, 0), diag(20, 2), 4, 2), "cavi",
    initialization_ols(), cavi_control(3, evaluation_draws = 3), plan,
    scaling = data_scaling(), n_draws = 3, seed = 1421,
    max_fits = 3, retain_fits = TRUE)
  expect_s3_class(result, "vb_cross_validation_result")
  expect_equal(nrow(result$folds), 3)
  expect_equal(nrow(result$pointwise), 24)
  expect_equal(result$statistics$elpd, sum(result$pointwise$log_predictive_density))
  expect_equal(result$statistics$rmse,
               sqrt(mean(result$prediction_pointwise$squared_error)))
  expect_equal(nrow(result$prediction_pointwise), 24)
  for (fold in seq_len(3)) {
    split <- cross_validation_split(plan, data$lrm, fold)
    expect_equal(scaling_map(result$fits[[fold]])$response_center,
                 mean(split$training$y))
  }
  expect_error(cross_validate_model(
    lrm_model(2), data$lrm,
    lrm_prior(c(0, 0), diag(20, 2), 4, 2), "cavi",
    initialization_ols(), cavi_control(2, evaluation_draws = 2), plan,
    seed = 1422, max_fits = 2), "exceeds max_fits")

  score_only <- cross_validate_model(
    lrm_model(2), data$lrm,
    lrm_prior(c(0, 0), diag(20, 2), 4, 2), "cavi",
    initialization_ols(), cavi_control(3, evaluation_draws = 3), plan,
    scaling = data_scaling(), n_draws = 3, seed = 1421,
    max_fits = 3, compute_performance = FALSE)
  expect_equal(score_only$pointwise, result$pointwise)
  expect_false("prediction_pointwise" %in% names(score_only))
})

test_that("grouped CV supports within-group HLRM and new-group CHLRM targets", {
  data <- .cv_data()
  control <- cavi_control(3, evaluation_draws = 2)
  within <- cross_validation_plan(
    data$grouped, "within_group_kfold", folds = 2, seed = 1430)
  hlrm <- cross_validate_model(
    hlrm_model(2), data$grouped,
    hlrm_prior(c(0, 0), diag(20, 2), 4, diag(2), 4, 2, 2),
    "cavi", initialization_ols(), control, within,
    n_draws = 2, seed = 1431, max_fits = 2)
  expect_identical(hlrm$predictive_target, "observed_groups")
  expect_identical(hlrm$unit, "observation")
  expect_equal(nrow(hlrm$pointwise), 36)
  expect_equal(nrow(hlrm$prediction_pointwise), 36)

  group_plan <- cross_validation_plan(
    data$grouped, "group_kfold", folds = 2, seed = 1432)
  chlrm <- cross_validate_model(
    chlrm_model(2, 2), data$grouped,
    chlrm_prior(c(0, 0), diag(20, 2), 4, diag(2), 4, 2, 2,
                c(1, 1)),
    "cavi", initialization_from_prior(), control, group_plan,
    n_draws = 2, seed = 1433, max_fits = 2)
  expect_identical(chlrm$predictive_target, "new_groups")
  expect_identical(chlrm$unit, "group")
  expect_equal(nrow(chlrm$pointwise), 6)
  expect_equal(sort(chlrm$pointwise$original_group), 1:6)
  expect_equal(nrow(chlrm$prediction_pointwise), 36)
  expect_equal(sort(unique(chlrm$prediction_pointwise$original_group)), 1:6)
})

test_that("invalid CV schemes fail before fitting", {
  data <- .cv_data()
  expect_error(cross_validation_plan(
    data$grouped, "observation_kfold", seed = 1440), "Grouped data require")
  expect_error(cross_validation_plan(
    data$lrm, "within_group_kfold", seed = 1440), "requires grouped")
  expect_error(cross_validation_plan(
    data$lrm, "group_kfold", seed = 1440), "requires one nonmissing")
  expect_error(cross_validation_plan(
    data$grouped, "within_group_kfold", folds = 7, seed = 1440),
    "cannot exceed")
})
