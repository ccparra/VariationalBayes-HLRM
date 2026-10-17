# Exact held-out log scoring from posterior or variational parameter draws.
# This module scores new responses; unlike predictive_criteria(), it does not
# reuse the observations that produced the fit.

.predictive_scoring_data <- function(fit, newdata) {
  if (!is.list(newdata) || !all(c("X", "y") %in% names(newdata)))
    stop("newdata must be a list containing X and y")
  if (fit$model == "LRM") {
    X <- as.matrix(newdata$X); y <- as.numeric(newdata$y)
    if (!is.numeric(X) || !nrow(X) || nrow(X) != length(y) ||
        any(!is.finite(X)) || any(!is.finite(y)))
      stop("LRM newdata require finite X and matching y")
    return(list(X = X, y = y, m = 1L, p = ncol(X), nj = length(y)))
  }
  X <- lapply(newdata$X, as.matrix)
  y <- lapply(newdata$y, as.numeric)
  info <- validate_grouped_data(X, y)
  list(X = X, y = y, m = info$m, p = info$p, nj = info$nj)
}

.validate_score_group_ids <- function(group_ids, group_count, fitted_groups,
                                      predictive_target) {
  if (predictive_target == "new_groups") {
    if (!is.null(group_ids) &&
        !(length(group_ids) == group_count && all(is.na(group_ids))))
      stop("new_groups requires group_ids = NULL or all NA")
    return(rep.int(NA_integer_, group_count))
  }
  if (!is.numeric(group_ids) || length(group_ids) != group_count ||
      any(!is.finite(group_ids)) || any(group_ids != floor(group_ids)) ||
      any(group_ids < 1L | group_ids > fitted_groups))
    stop("observed_groups requires one valid fitted group_id per data group")
  as.integer(group_ids)
}

.chlrm_observed_log_likelihood <- function(y, X, draws, allocations,
                                           group_ids, unit,
                                           max_output_mb) {
  info <- validate_grouped_data(X, y)
  dimensions <- validate_component_draws_chlrm(draws, info$p)
  if (!is.matrix(allocations) || nrow(allocations) != dimensions$S ||
      any(!is.finite(allocations)) || any(allocations != floor(allocations)) ||
      any(allocations < 1L | allocations > dimensions$K) ||
      any(group_ids > ncol(allocations)))
    stop("Invalid CHLRM allocation draws for observed groups")
  rows <- if (unit == "observation") sum(info$nj) else info$m
  validate_evaluation_size(
    dimensions$S, max_output_mb, as.double(rows) * dimensions$S)
  values <- matrix(NA_real_, rows, dimensions$S)
  starts <- cumsum(c(1L, head(info$nj, -1L)))
  for (group in seq_len(info$m)) for (draw in seq_len(dimensions$S)) {
    component <- allocations[draw, group_ids[group]]
    contribution <- stats::dnorm(
      y[[group]], drop(X[[group]] %*% draws$Betak[draw, component, ]),
      sqrt(draws$sigma2k[draw, component]), log = TRUE)
    if (unit == "observation") {
      indices <- starts[group]:(starts[group] + info$nj[group] - 1L)
      values[indices, draw] <- contribution
    } else {
      values[group, draw] <- sum(contribution)
    }
  }
  new_log_likelihood_result(
    values, "CHLRM", unit, "parameter_conditional",
    if (unit == "observation") .observation_ids(info$nj) else
      .group_ids(seq_len(info$m)),
    metadata = list(
      conditioned_on = c("component_parameters", "observed_group_allocation")))
}

#' Held-out pointwise log likelihood
#'
#' Draws are taken from a fitted object, already restored to the input scale.
#' Grouped models require an explicit predictive target: observed groups use
#' their learned local state, while new groups integrate that state out.
#' @rdname predictive-evaluation
#' @export
predictive_log_likelihood <- function(
    fit, newdata, predictive_target = NULL, predictive_unit = NULL,
    group_ids = NULL, groups = NULL, n_draws = NULL, seed = NULL,
    max_output_mb = 256) {
  if (!inherits(fit, "vb_model_fit"))
    stop("fit must be produced by fit_model()")
  data <- .predictive_scoring_data(fit, newdata)
  expected_p <- fit$specifications$model$coefficient_count
  if (!is.null(expected_p) && data$p != expected_p)
    stop("newdata coefficient count differs from the fitted model")
  if (is.null(predictive_target)) predictive_target <- if (fit$model == "LRM")
    "observed_data" else "new_groups"
  allowed_target <- if (fit$model == "LRM") "observed_data" else
    c("observed_groups", "new_groups")
  predictive_target <- match.arg(predictive_target, allowed_target)
  if (is.null(predictive_unit)) predictive_unit <-
    if (predictive_target == "new_groups") "group" else "observation"
  predictive_unit <- match.arg(predictive_unit, c("observation", "group"))
  if (predictive_target == "new_groups" && predictive_unit != "group")
    stop("new_groups requires predictive_unit = 'group'")

  include_latent <- fit$model == "CHLRM" &&
    predictive_target == "observed_groups"
  sampled <- posterior_draws(
    fit, n_draws = n_draws, seed = seed, include_latent = include_latent,
    max_output_mb = max_output_mb)
  if (sampled$metadata$draw_count < 2L)
    stop("Predictive scoring requires at least two posterior draws")
  parameters <- sampled$parameters

  if (fit$model == "LRM") {
    if (predictive_unit == "group" && is.null(groups))
      stop("groups are required for grouped LRM scoring")
    if (predictive_unit == "observation" && !is.null(groups))
      stop("groups must be NULL for observation-level LRM scoring")
    likelihood <- log_likelihood(
      "LRM", data$y, data$X,
      list(Beta = parameters$coefficients,
           sigma2 = parameters$residual_variance),
      unit = predictive_unit, conditioning = "parameter_conditional",
      groups = groups, max_output_mb = max_output_mb)
  } else if (fit$model == "HLRM") {
    fitted_groups <- dim(parameters$group_coefficients)[2L]
    ids <- .validate_score_group_ids(
      group_ids, data$m, fitted_groups, predictive_target)
    if (predictive_target == "observed_groups") {
      likelihood <- log_likelihood(
        "HLRM", data$y, data$X,
        list(
          Betaj = parameters$group_coefficients[, ids, , drop = FALSE],
          sigma2j = parameters$group_variances[, ids, drop = FALSE]),
        unit = predictive_unit, conditioning = "parameter_conditional",
        max_output_mb = max_output_mb)
    } else {
      likelihood <- log_likelihood(
        "HLRM", data$y, data$X,
        list(Beta = parameters$population_coefficients,
             Sigma = parameters$population_covariance,
             xi2 = parameters$variance_scale),
        unit = "group", conditioning = "group_marginal",
        nu0 = fit$specifications$prior$residual_degrees_of_freedom,
        max_output_mb = max_output_mb)
    }
  } else {
    fitted_groups <- if (fit$method == "MCMC")
      nrow(fit$latent_variables$group_allocations) else
      nrow(fit$latent_variables$allocation_probabilities)
    ids <- .validate_score_group_ids(
      group_ids, data$m, fitted_groups, predictive_target)
    component_draws <- list(
      Betak = parameters$component_coefficients,
      sigma2k = parameters$component_variances,
      omega = parameters$mixture_weights)
    likelihood <- if (predictive_target == "observed_groups")
      .chlrm_observed_log_likelihood(
        data$y, data$X, component_draws,
        sampled$latent_variables$group_allocations, ids,
        predictive_unit, max_output_mb) else
      log_likelihood(
        "CHLRM", data$y, data$X, component_draws,
        unit = "group", conditioning = "group_marginal",
        max_output_mb = max_output_mb)
  }
  likelihood$metadata$predictive_target <- predictive_target
  likelihood$metadata$data_role <- "held_out"
  likelihood$metadata$density_scale <- "input"
  likelihood
}

#' @rdname predictive-evaluation
#' @export
predictive_score <- function(
    fit, newdata, predictive_target = NULL, predictive_unit = NULL,
    group_ids = NULL, groups = NULL, n_draws = NULL, seed = NULL,
    return_log_likelihood = FALSE, max_output_mb = 256) {
  likelihood <- predictive_log_likelihood(
    fit, newdata, predictive_target, predictive_unit, group_ids, groups,
    n_draws, seed, max_output_mb)
  new_predictive_score_result(fit, likelihood, return_log_likelihood)
}

#' @rdname predictive-evaluation
#' @export
predictive_performance <- function(
    fit, newdata, predictive_target = NULL, predictive_unit = NULL,
    group_ids = NULL, groups = NULL, n_draws = NULL, level = 0.95,
    seed, max_output_mb = 256) {
  if (missing(seed) || is.null(seed))
    stop("seed is required for reproducible predictive performance")
  validate_rng_seed(seed)
  data <- .predictive_scoring_data(fit, newdata)
  if (is.null(predictive_target)) predictive_target <- if (fit$model == "LRM")
    "observed_data" else "new_groups"
  predictive_target <- match.arg(
    predictive_target, if (fit$model == "LRM") "observed_data" else
      c("observed_groups", "new_groups"))
  if (fit$model != "LRM" && predictive_target == "new_groups") {
    if (is.null(group_ids)) group_ids <- rep.int(NA_integer_, data$m)
    if (length(group_ids) != data$m || !all(is.na(group_ids)))
      stop("new_groups requires group_ids = NULL or all NA")
  }
  if (fit$model != "LRM" && predictive_target == "observed_groups" &&
      is.null(group_ids))
    stop("observed_groups requires fitted group_ids")
  streams <- derive_independent_seeds(3L, seed)
  prediction_arguments <- list(
    object = fit, newdata = newdata$X, n_draws = n_draws,
    level = level, seed = streams[1L], return_draws = FALSE,
    max_output_mb = max_output_mb)
  if (fit$model != "LRM") prediction_arguments$group_ids <- group_ids
  mean_prediction <- do.call(predict, prediction_arguments)
  names(prediction_arguments)[names(prediction_arguments) == "object"] <- "fit"
  prediction_arguments$seed <- streams[2L]
  response_prediction <- do.call(posterior_predict, prediction_arguments)
  score <- predictive_score(
    fit, newdata, predictive_target, predictive_unit, group_ids, groups,
    n_draws, streams[3L], return_log_likelihood = FALSE,
    max_output_mb = max_output_mb)
  observed <- if (fit$model == "LRM") data$y else
    unlist(data$y, use.names = FALSE)
  new_predictive_performance_result(
    fit, predictive_target, level, observed,
    mean_prediction, response_prediction, score)
}
# Vehtari et al.'s variance form of WAIC. Rows are predictive units and
# columns are posterior draws; grouping observations changes the estimand.
#' @rdname predictive-evaluation
#' @export
waic_from_loglik <- function(log_lik) {
  contract <- if (inherits(log_lik, "vb_log_likelihood")) log_lik else NULL
  values <- as.matrix(log_lik)
  if (!is.numeric(values) || !nrow(values) || any(!is.finite(values))) {
    stop("WAIC requires a nonempty finite numeric log-likelihood matrix")
  }
  if (ncol(values) < 2L) {
    stop("WAIC requires at least two posterior draws")
  }
  pointwise_lppd <- apply(values, 1L, logmeanexp)
  pointwise_penalty <- apply(values, 1L, stats::var)
  result <- summarize_waic_pointwise(data.frame(
    lppd = pointwise_lppd,
    p_waic = pointwise_penalty
  ))
  if (!is.null(contract)) {
    result$estimand <- log_likelihood_estimand(contract)
    result$predictive_units <- contract$predictive_units
  }
  result
}
