# Reproducible cross-validation plans and a bounded sequential runner. Scenario
# grids, checkpointing and parallel scheduling remain responsibilities of the
# paper evaluation suite.

.balanced_fold_assignment <- function(count, folds) {
  sample(rep(seq_len(folds), length.out = count), count, replace = FALSE)
}

#' @rdname predictive-evaluation
#' @export
cross_validation_plan <- function(
    data, scheme = c("observation_kfold", "within_group_kfold",
                     "group_kfold", "leave_one_group_out"),
    folds = 5L, groups = NULL, seed) {
  scheme <- match.arg(scheme)
  if (missing(seed) || is.null(seed))
    stop("seed is required for a reproducible cross-validation plan")
  validate_rng_seed(seed)
  grouped <- is.list(data$X) && !is.data.frame(data$X)
  parts <- .scaling_data_parts(data)
  if (grouped && scheme == "observation_kfold")
    stop("Grouped data require within_group_kfold, group_kfold or leave_one_group_out")
  if (!grouped && scheme == "within_group_kfold")
    stop("within_group_kfold requires grouped data")
  if (!grouped && scheme %in% c("group_kfold", "leave_one_group_out")) {
    if (is.null(groups) || !is.atomic(groups) ||
        length(groups) != parts$group_sizes || anyNA(groups))
      stop("Ungrouped group-level CV requires one nonmissing group label per row")
    labels <- unique(groups)
    group_count <- length(labels)
  } else if (grouped) {
    labels <- seq_along(parts$X)
    group_count <- length(labels)
  }

  if (scheme == "leave_one_group_out") {
    folds <- group_count
  } else {
    folds <- validate_integer_scalar(folds, "folds", 2L)
  }
  units <- switch(scheme,
    observation_kfold = parts$group_sizes,
    within_group_kfold = min(parts$group_sizes),
    group_kfold = group_count,
    leave_one_group_out = group_count)
  if (folds > units)
    stop("folds cannot exceed the available units in the selected scheme")

  scope <- inference_rng_scope(seed); on.exit(scope$restore(), add = TRUE)
  assignments <- switch(scheme,
    observation_kfold = .balanced_fold_assignment(parts$group_sizes, folds),
    within_group_kfold = lapply(
      parts$group_sizes, .balanced_fold_assignment, folds = folds),
    group_kfold = data.frame(
      group = labels,
      fold = .balanced_fold_assignment(group_count, folds),
      stringsAsFactors = FALSE),
    leave_one_group_out = data.frame(
      group = labels, fold = seq_len(group_count),
      stringsAsFactors = FALSE))
  predictive_target <- if (!grouped) "observed_data" else if (
    scheme == "within_group_kfold") "observed_groups" else "new_groups"
  predictive_unit <- if (scheme %in% c("group_kfold", "leave_one_group_out"))
    "group" else "observation"
  out <- list(
    schema_version = "cross-validation-plan-1",
    data_structure = if (grouped) "grouped" else "ungrouped",
    scheme = scheme,
    fold_count = as.integer(folds),
    predictive_target = predictive_target,
    predictive_unit = predictive_unit,
    assignments = assignments,
    data_signature = list(
      observations = sum(parts$group_sizes),
      coefficient_count = parts$p,
      group_sizes = as.integer(parts$group_sizes)),
    metadata = list(seed = seed))
  if (!grouped && scheme %in% c("group_kfold", "leave_one_group_out"))
    out$observation_groups <- groups
  structure(out, class = c("vb_cross_validation_plan", "list"))
}

.subset_ungrouped_data <- function(data, indices) {
  out <- data
  out$X <- as.matrix(data$X)[indices, , drop = FALSE]
  out$y <- as.numeric(data$y)[indices]
  out
}

.subset_grouped_rows <- function(data, indices) {
  out <- data
  out$X <- Map(function(X, rows) as.matrix(X)[rows, , drop = FALSE],
               data$X, indices)
  out$y <- Map(function(y, rows) as.numeric(y)[rows], data$y, indices)
  out
}

#' @rdname predictive-evaluation
#' @export
cross_validation_split <- function(plan, data, fold) {
  if (!inherits(plan, "vb_cross_validation_plan"))
    stop("plan must be created by cross_validation_plan()")
  fold <- validate_integer_scalar(fold, "fold", 1L)
  if (fold > plan$fold_count) stop("fold exceeds the plan's fold_count")
  parts <- .scaling_data_parts(data)
  signature <- list(
    observations = sum(parts$group_sizes), coefficient_count = parts$p,
    group_sizes = as.integer(parts$group_sizes))
  if (!identical(signature, plan$data_signature))
    stop("data dimensions differ from the cross-validation plan")

  if (plan$scheme == "observation_kfold") {
    assessment <- which(plan$assignments == fold)
    training <- which(plan$assignments != fold)
    return(list(
      training = .subset_ungrouped_data(data, training),
      assessment = .subset_ungrouped_data(data, assessment),
      scoring = list(original_observation = assessment)))
  }
  if (plan$scheme == "within_group_kfold") {
    assessment <- lapply(plan$assignments, function(x) which(x == fold))
    training <- lapply(plan$assignments, function(x) which(x != fold))
    return(list(
      training = .subset_grouped_rows(data, training),
      assessment = .subset_grouped_rows(data, assessment),
      scoring = list(
        group_ids = seq_along(data$X),
        original_group = rep.int(seq_along(assessment), lengths(assessment)),
        original_observation = unlist(assessment, use.names = FALSE))))
  }

  held_labels <- plan$assignments$group[plan$assignments$fold == fold]
  if (plan$data_structure == "ungrouped") {
    assessment <- which(plan$observation_groups %in% held_labels)
    training <- which(!plan$observation_groups %in% held_labels)
    return(list(
      training = .subset_ungrouped_data(data, training),
      assessment = .subset_ungrouped_data(data, assessment),
      scoring = list(
        groups = plan$observation_groups[assessment],
        original_observation = assessment)))
  }
  assessment_groups <- as.integer(held_labels)
  training_groups <- setdiff(seq_along(data$X), assessment_groups)
  training_data <- data
  training_data$X <- data$X[training_groups]
  training_data$y <- data$y[training_groups]
  assessment_data <- data
  assessment_data$X <- data$X[assessment_groups]
  assessment_data$y <- data$y[assessment_groups]
  list(
    training = training_data,
    assessment = assessment_data,
    scoring = list(
      group_ids = rep.int(NA_integer_, length(assessment_groups)),
      original_group = assessment_groups))
}

#' @rdname predictive-evaluation
#' @export
cross_validate_model <- function(
    model, data, prior, method = c("mcmc", "cavi", "svi"),
    initialization = initialization_from_prior(), control, plan,
    scaling = NULL, n_draws = NULL, seed, max_fits = 20L,
    max_output_mb = 256, retain_fits = FALSE, compute_performance = TRUE,
    prediction_level = 0.95, verbose = FALSE) {
  method <- match.arg(tolower(method), c("mcmc", "cavi", "svi"))
  if (!inherits(plan, "vb_cross_validation_plan"))
    stop("plan must be created by cross_validation_plan()")
  if (missing(seed) || is.null(seed))
    stop("seed is required for reproducible cross-validation")
  validate_rng_seed(seed)
  max_fits <- validate_integer_scalar(max_fits, "max_fits", 1L)
  if (plan$fold_count > max_fits)
    stop("Cross-validation exceeds max_fits; raise the explicit budget to proceed")
  if (!is.logical(retain_fits) || length(retain_fits) != 1L ||
      is.na(retain_fits)) stop("retain_fits must be TRUE or FALSE")
  if (!is.logical(compute_performance) || length(compute_performance) != 1L ||
      is.na(compute_performance))
    stop("compute_performance must be TRUE or FALSE")
  if (compute_performance) validate_prediction_level(prediction_level)
  seeds <- matrix(
    derive_independent_seeds(3L * plan$fold_count, seed),
    nrow = plan$fold_count, ncol = 3L,
    dimnames = list(NULL, c("fit", "score", "performance")))
  fold_statistics <- vector("list", plan$fold_count)
  pointwise <- vector("list", plan$fold_count)
  prediction_pointwise <- if (compute_performance)
    vector("list", plan$fold_count) else NULL
  fits <- if (retain_fits) vector("list", plan$fold_count) else NULL

  for (fold in seq_len(plan$fold_count)) {
    split <- cross_validation_split(plan, data, fold)
    fitted <- fit_model(
      model, split$training, prior, method, initialization, control,
      scaling = scaling, seed = seeds[fold, "fit"], verbose = verbose)
    score_arguments <- list(
      fit = fitted, newdata = split$assessment,
      predictive_target = plan$predictive_target,
      predictive_unit = plan$predictive_unit,
      n_draws = n_draws, seed = seeds[fold, "score"],
      max_output_mb = max_output_mb)
    if (plan$scheme == "within_group_kfold")
      score_arguments$group_ids <- split$scoring$group_ids
    if (!is.null(split$scoring$groups))
      score_arguments$groups <- split$scoring$groups
    score <- do.call(predictive_score, score_arguments)
    score_statistics <- score$statistics
    score_pointwise <- score$pointwise
    if (compute_performance) {
      performance_arguments <- score_arguments
      performance_arguments$level <- prediction_level
      performance_arguments$seed <- seeds[fold, "performance"]
      performance <- do.call(predictive_performance, performance_arguments)
    }
    fold_statistics[[fold]] <- data.frame(
      fold = fold, elpd = score_statistics$elpd,
      mean_log_predictive_density =
        score_statistics$mean_log_predictive_density,
      predictive_units = score_statistics$predictive_units,
      rmse = if (compute_performance) performance$statistics$rmse else NA_real_,
      mae = if (compute_performance) performance$statistics$mae else NA_real_,
      empirical_coverage = if (compute_performance)
        performance$statistics$empirical_coverage else NA_real_,
      mean_interval_width = if (compute_performance)
        performance$statistics$mean_interval_width else NA_real_)
    fold_pointwise <- score_pointwise
    fold_pointwise$fold <- fold
    if (!is.null(split$scoring$original_group))
      fold_pointwise$original_group <- split$scoring$original_group
    if (!is.null(split$scoring$original_observation))
      fold_pointwise$original_observation <-
        split$scoring$original_observation
    pointwise[[fold]] <- fold_pointwise
    if (compute_performance) {
      prediction <- performance$pointwise
      prediction$fold <- fold
      if (!is.null(split$scoring$original_group)) {
        if (length(split$scoring$original_group) == nrow(prediction)) {
          prediction$original_group <- split$scoring$original_group
        } else {
          prediction$original_group <-
            split$scoring$original_group[prediction$group]
        }
      }
      if (!is.null(split$scoring$original_observation))
        prediction$original_observation <-
          split$scoring$original_observation
      prediction_pointwise[[fold]] <- prediction
    }
    if (retain_fits) fits[[fold]] <- fitted
  }
  new_cross_validation_result(
    model$model, .canonical_public_method(method), plan,
    do.call(rbind, fold_statistics), do.call(rbind, pointwise),
    as.data.frame(seeds), fits,
    if (compute_performance) do.call(rbind, prediction_pointwise) else NULL)
}
