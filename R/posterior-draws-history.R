# Convert exact native CHLRM history ranges into the typed posterior-draw contract.
# This keeps history persistence separate from prediction while allowing callers
# to process large histories in bounded, ordered blocks.
#' Posterior draws from a persisted CHLRM history
#'
#' Read selected retained iterations from an immutable MCMC history and expose
#' them through the same typed draw contract as [posterior_draws()]. The fitting
#' object supplies model and scaling metadata; native history values are not
#' modified. Supply one or more history handles, with one chain identity per
#' handle. Use successive iteration ranges to process large histories in blocks.
#'
#' @param fit A fitted CHLRM MCMC object, used for model and scaling metadata.
#' @param histories A non-empty list of `vb_mcmc_history` handles.
#' @param chains Integer chain identity within each corresponding history.
#' @param chain_ids Integer identities assigned to the resulting draw index.
#' @param iterations Strictly increasing retained iteration identities. Supply
#'   one shared vector or a list with one vector per history. Defaults to all
#'   retained iterations in every history. Cannot be combined with `n_draws`.
#' @param n_draws Optional number of draws sampled without replacement from the
#'   chain-major concatenation of all supplied histories.
#' @param seed Optional seed used only when `n_draws` requests random selection.
#' @param max_output_mb Positive limit for the estimated typed draw allocation.
#' @return A `vb_posterior_draws` object with component parameters and group
#'   allocations in the requested chain-major and iteration order.
#' @export
posterior_draws_from_mcmc_history <- function(fit, histories, chains = 1L,
    chain_ids = seq_along(histories), iterations = NULL, n_draws = NULL,
    seed = NULL, max_output_mb = 256) {
  if (!inherits(fit, "vb_model_fit") || fit$model != "CHLRM" || fit$method != "MCMC")
    stop("fit must be a fitted CHLRM MCMC object")
  if (!is.list(histories) || !length(histories) ||
      !all(vapply(histories, inherits, logical(1L), "vb_mcmc_history")))
    stop("histories must be a non-empty list of MCMC history handles")
  if (length(chains) != length(histories) || anyNA(chains) ||
      any(!is.finite(chains)) || any(chains < 1L) || any(chains != floor(chains)))
    stop("chains must identify one valid chain per history")
  if (length(chain_ids) != length(histories) || anyNA(chain_ids) ||
      any(!is.finite(chain_ids)) || any(chain_ids < 1L) ||
      any(chain_ids != floor(chain_ids)) || anyDuplicated(chain_ids))
    stop("chain_ids must be unique positive integer identities")
  if (!is.numeric(max_output_mb) || length(max_output_mb) != 1L ||
      !is.finite(max_output_mb) || max_output_mb <= 0)
    stop("max_output_mb must be a positive finite scalar")

  states <- lapply(histories, mcmc_history_restart_state)
  available <- vapply(seq_along(states), function(i) {
    if (states[[i]]$model != "CHLRM" || chains[i] > length(states[[i]]$chains))
      stop("History model or chain identity does not match fit")
    as.integer(states[[i]]$chains[[chains[i]]]$run_settings$samples)
  }, integer(1L))
  if (length(unique(available)) != 1L) stop("Histories must have equal retained lengths")
  n <- available[1L]
  total_available <- length(histories) * n
  validate_rng_seed(seed)
  if (!is.null(n_draws)) {
    if (!is.null(iterations)) stop("Specify either iterations or n_draws, not both")
    request <- .validate_posterior_draw_request(n_draws, max_output_mb)
    if (request$n_draws > total_available)
      stop("n_draws cannot exceed the combined retained history length")
    scope <- inference_rng_scope(seed); on.exit(scope$restore(), add = TRUE)
    selected <- prediction_draw_indices(total_available, request$n_draws)
    iterations <- lapply(seq_along(histories), function(i) {
      first <- (i - 1L) * n + 1L; last <- i * n
      as.integer(selected[selected >= first & selected <= last] - first + 1L)
    })
  } else if (is.null(iterations)) iterations <- rep(list(seq_len(n)), length(histories))
  else if (!is.list(iterations)) iterations <- rep(list(iterations), length(histories))
  if (length(iterations) != length(histories) ||
      !all(vapply(iterations, function(index) is.numeric(index) &&
        (!length(index) || (!anyNA(index) && all(is.finite(index)) &&
          all(index >= 1L & index <= n) && all(index == floor(index)) &&
          all(diff(index) > 0))), logical(1L))))
    stop("iterations must contain increasing retained identities for every history")
  iterations <- lapply(iterations, as.integer)

  template <- mcmc_history_draws(histories[[1L]], chains[1L], iterations = 1L,
    fields = c("Betak", "sigma2k", "omega", "z"), max_output_mb = 32)
  dims <- lapply(template$draws, dim)
  k <- dims$Betak[2L]; p <- dims$Betak[3L]; m <- dims$z[1L]
  if (fit$specifications$model$component_count != k ||
      fit$specifications$model$coefficient_count != p ||
      dim(fit$latent_variables$group_allocations)[1L] != m)
    stop("History dimensions do not match fit")
  counts <- lengths(iterations)
  total <- sum(counts)
  if (total < 1L) stop("At least one retained draw must be selected")
  bytes_per_draw <- 8 * (k * p + 2 * k) + 4 * m
  estimated <- 1.5 * as.double(total) * (bytes_per_draw + 8) / 1024^2
  if (!is.finite(estimated) || estimated > max_output_mb)
    stop("Posterior history draws exceed max_output_mb; request a smaller iteration block")

  pieces <- lapply(seq_along(histories), function(i) {
    if (!counts[i]) return(NULL)
    mcmc_history_draws(histories[[i]], chains[i], iterations = iterations[[i]],
      fields = c("Betak", "sigma2k", "omega", "z"),
      max_output_mb = max_output_mb, max_block_mb = max_output_mb)
  })
  component <- array(NA_real_, c(total, k, p))
  variance <- matrix(NA_real_, total, k)
  weights <- matrix(NA_real_, total, k)
  allocation <- matrix(NA_integer_, total, m)
  for (i in seq_along(pieces)) {
    count <- counts[i]
    if (!count) next
    if (!identical(dim(pieces[[i]]$draws$Betak), c(count, k, p)) ||
        !identical(dim(pieces[[i]]$draws$sigma2k), c(count, k)) ||
        !identical(dim(pieces[[i]]$draws$omega), c(count, k)) ||
        !identical(dim(pieces[[i]]$draws$z), c(m, count)))
      stop("History parameter dimensions differ across sources")
    rows <- (sum(counts[seq_len(i - 1L)]) + 1L):(sum(counts[seq_len(i - 1L)]) + count)
    component[rows, , ] <- pieces[[i]]$draws$Betak
    variance[rows, ] <- pieces[[i]]$draws$sigma2k
    weights[rows, ] <- pieces[[i]]$draws$omega
    allocation[rows, ] <- t(pieces[[i]]$draws$z)
  }
  map <- .fit_scaling_map(fit)
  if (!is.null(map)) {
    component <- .transform_coefficient_array(component, map, "to_input")
    variance <- restore_residual_variance(variance, map)
  }
  parameters <- list(component_coefficients = component,
    component_variances = variance, mixture_weights = weights)
  source_index <- unlist(lapply(seq_along(histories), function(i)
    (i - 1L) * n + iterations[[i]]), use.names = FALSE)
  draw_index <- data.frame(draw = seq_len(total), stored_draw = as.integer(source_index),
    chain = rep(as.integer(chain_ids), times = counts),
    iteration = as.integer(unlist(iterations, use.names = FALSE)))
  sampled <- list(parameters = parameters,
    latent_variables = list(group_allocations = allocation))
  out <- new_posterior_draws_result(fit, sampled, draw_index, seed = seed,
    estimated_storage_mb = estimated,
    parameter_names = names(parameters))
  out$metadata$draw_source <- "persisted_mcmc_history"
  out
}

.posterior_prediction_from_draws <- function(draws, newdata, n_draws, seed,
    group_ids, level, return_draws, max_output_mb, include_noise) {
  if (!inherits(draws, "vb_posterior_draws") || draws$model != "CHLRM" ||
      draws$method != "MCMC" || is.null(draws$latent_variables$group_allocations))
    stop("Require CHLRM MCMC posterior draws with group allocations")
  required <- c("component_coefficients", "component_variances", "mixture_weights")
  if (!all(required %in% names(draws$parameters)))
    stop("Posterior draws are missing component parameters")
  available <- nrow(draws$parameters$component_coefficients)
  count <- prediction_draw_count(n_draws, available, available)
  validate_rng_seed(seed)
  scope <- inference_rng_scope(seed); on.exit(scope$restore(), add = TRUE)
  indices <- prediction_draw_indices(available, count)
  parameters <- draws$parameters
  selected <- list(Betak = parameters$component_coefficients[indices, , , drop = FALSE],
    sigma2k = parameters$component_variances[indices, , drop = FALSE],
    omega = parameters$mixture_weights[indices, , drop = FALSE])
  z <- t(draws$latent_variables$group_allocations[indices, , drop = FALSE])
  ans <- posterior_prediction_chlrm(newdata, selected, group_ids, "MCMC", z = z,
    include_noise = include_noise, level = level, return_draws = return_draws,
    max_output_mb = max_output_mb)
  if ("stored_draw" %in% names(draws$draw_index))
    ans$metadata$draw_indices <- draws$draw_index$stored_draw[indices]
  ans$metadata$chain_ids <- draws$draw_index$chain[indices]
  ans$metadata$iterations <- draws$draw_index$iteration[indices]
  ans
}

#' @rdname predictive-evaluation
#' @export
predict.vb_posterior_draws <- function(object, newdata, n_draws = NULL,
    seed = NULL, group_ids = NULL, level = 0.95, return_draws = FALSE,
    max_output_mb = 256, ...) {
  .posterior_prediction_from_draws(object, newdata, n_draws, seed, group_ids,
    level, return_draws, max_output_mb, include_noise = FALSE)
}

#' @rdname predictive-evaluation
#' @export
posterior_predict.vb_posterior_draws <- function(fit, newdata, n_draws = NULL,
    seed = NULL, group_ids = NULL, level = 0.95, return_draws = TRUE,
    max_output_mb = 256, ...) {
  .posterior_prediction_from_draws(fit, newdata, n_draws, seed, group_ids,
    level, return_draws, max_output_mb, include_noise = TRUE)
}


#' Held-out predictive performance directly from persisted histories
#'
#' Sample a reproducible number of draws from one or more persisted CHLRM
#' histories, then compute held-out predictive performance using package
#' prediction and scoring kernels. No model fitting or history mutation occurs.
#'
#' @param fit A fitted CHLRM MCMC object used as model and scaling metadata.
#' @param newdata Grouped held-out data containing `X` and `y`.
#' @param histories A non-empty list of `vb_mcmc_history` handles.
#' @param chains Integer chain identity within each history.
#' @param chain_ids Integer chain labels for returned provenance.
#' @param n_draws Number of posterior draws sampled without replacement.
#' @param seed Seed for draw selection and posterior response simulation.
#' @param predictive_target The supported target, `observed_groups`.
#' @param predictive_unit Predictive unit: `observation` or `group`.
#' @param group_ids Fitted group identity for every held-out group.
#' @param level Pointwise predictive interval probability.
#' @param max_output_mb Maximum estimated allocation for draws and evaluation.
#' @return A `vb_predictive_performance` result.
#' @export
predictive_performance_from_mcmc_history <- function(fit, newdata, histories,
    chains = 1L, chain_ids = seq_along(histories), n_draws, seed,
    predictive_target = "observed_groups", predictive_unit = "observation",
    group_ids, level = 0.95, max_output_mb = 256) {
  if (missing(n_draws) || missing(seed) || missing(group_ids))
    stop("n_draws, seed, and group_ids are required")
  streams <- derive_independent_seeds(3L, seed)
  sources <- lapply(seq_len(3L), function(i) posterior_draws_from_mcmc_history(
    fit, histories, chains, chain_ids, n_draws = n_draws, seed = streams[i],
    max_output_mb = max_output_mb))
  states <- lapply(histories, mcmc_history_restart_state)
  available <- vapply(seq_along(histories), function(i)
    states[[i]]$chains[[chains[i]]]$run_settings$samples, integer(1L))
  if (length(unique(available)) != 1L) stop("Histories must have equal retained lengths")
  response_scope <- inference_rng_scope(streams[2L])
  on.exit(response_scope$restore(), add = TRUE)
  invisible(prediction_draw_indices(length(histories) * available[1L], n_draws))
  .predictive_performance_from_draw_sets(fit, newdata, sources[[1L]],
    sources[[2L]], sources[[3L]], predictive_target, predictive_unit,
    group_ids, NULL, level, max_output_mb)
}

#' Held-out predictive performance from typed posterior draws
#'
#' Evaluate held-out CHLRM responses from draws already extracted from a fit or
#' persisted history. The supplied draw object is used for all three summaries;
#' use [predictive_performance_from_mcmc_history()] to reproduce the independent
#' seeded draw subsets used by [predictive_performance()].
#'
#' @param fit A fitted CHLRM MCMC object used as model metadata.
#' @param newdata A grouped list with design matrices `X` and responses `y`.
#' @param posterior_draws A `vb_posterior_draws` object containing CHLRM
#'   component parameters and group allocations.
#' @param predictive_target The supported target, `observed_groups`.
#' @param predictive_unit Predictive unit: `observation` or `group`.
#' @param group_ids Fitted group identity for each held-out group.
#' @param seed Seed for posterior response simulation.
#' @param level Pointwise predictive interval probability.
#' @param max_output_mb Maximum estimated prediction/scoring output size.
#' @return A `vb_predictive_performance` result using the supplied draws.
#' @export
predictive_performance_from_draws <- function(fit, newdata, posterior_draws,
    predictive_target = "observed_groups", predictive_unit = "observation",
    group_ids, seed, level = 0.95, max_output_mb = 256) {
  if (missing(seed) || is.null(seed))
    stop("seed is required for reproducible predictive performance")
  streams <- derive_independent_seeds(3L, seed)
  .predictive_performance_from_draw_sets(fit, newdata, posterior_draws,
    posterior_draws, posterior_draws, predictive_target, predictive_unit,
    group_ids, streams[2L], level, max_output_mb)
}

.predictive_performance_from_draw_sets <- function(fit, newdata, mean_draws,
    response_draws, score_draws, predictive_target, predictive_unit,
    group_ids, response_seed, level, max_output_mb) {
  if (!inherits(fit, "vb_model_fit") || fit$model != "CHLRM" || fit$method != "MCMC")
    stop("fit must be a fitted CHLRM MCMC object")
  sources <- list(mean_draws, response_draws, score_draws)
  if (!all(vapply(sources, function(x) inherits(x, "vb_posterior_draws") &&
      x$model == "CHLRM" && x$method == "MCMC" &&
      !is.null(x$latent_variables$group_allocations), logical(1L))))
    stop("All draw sources must contain CHLRM MCMC component and allocation draws")
  predictive_target <- match.arg(predictive_target, "observed_groups")
  predictive_unit <- match.arg(predictive_unit, c("observation", "group"))
  if (is.null(group_ids)) stop("group_ids are required for observed-group prediction")
  validate_rng_seed(response_seed)
  data <- .predictive_scoring_data(fit, newdata)
  expected_p <- fit$specifications$model$coefficient_count
  if (data$p != expected_p) stop("newdata coefficient count differs from the fitted model")
  ids <- .validate_score_group_ids(group_ids, data$m,
    ncol(mean_draws$latent_variables$group_allocations), predictive_target)
  X <- lapply(newdata$X, as.matrix)
  mean_prediction <- predict(mean_draws, X, group_ids = ids,
    level = level, return_draws = FALSE, max_output_mb = max_output_mb)
  response_prediction <- posterior_predict(response_draws, X, group_ids = ids,
    seed = response_seed, level = level, return_draws = FALSE,
    max_output_mb = max_output_mb)
  parameters <- score_draws$parameters
  allocations <- score_draws$latent_variables$group_allocations
  likelihood <- .chlrm_observed_log_likelihood(data$y, data$X,
    list(Betak = parameters$component_coefficients,
      sigma2k = parameters$component_variances,
      omega = parameters$mixture_weights), allocations, ids,
    predictive_unit, max_output_mb)
  likelihood$metadata$predictive_target <- predictive_target
  likelihood$metadata$data_role <- "held_out"
  likelihood$metadata$density_scale <- "input"
  score <- new_predictive_score_result(fit, likelihood)
  score$metadata$draw_source <- score_draws$metadata$draw_source %||% "typed_posterior_draws"
  result <- new_predictive_performance_result(fit, predictive_target, level,
    unlist(data$y, use.names = FALSE), mean_prediction, response_prediction, score)
  result$metadata$draw_source <- score$metadata$draw_source
  result$metadata$draw_counts <- c(mean = nrow(mean_draws$draw_index),
    response = nrow(response_draws$draw_index), scoring = nrow(score_draws$draw_index))
  result
}
