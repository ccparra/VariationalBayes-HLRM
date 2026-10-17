# Public operations for vb_model_fit objects. Mutable R6 engines remain an
# implementation detail reached through one dispatcher. When fitting used an
# affine scaling map, engines stay in numerical model units while every public
# operation is restored to the input-data units here.

.invoke_fit_backend <- function(fit, operation, arguments = list()) {
  if (!inherits(fit, "vb_model_fit"))
    stop("fit must be produced by fit_model()")
  backend <- .fit_backend(fit)
  method <- backend[[operation]]
  if (!is.function(method))
    stop("The selected backend does not implement ", operation, "()")
  result <- do.call(method, arguments)
  if (is.list(result) && "method" %in% names(result))
    result$method <- fit$method
  result
}

.fit_scaling_map <- function(fit) attr(fit, "scaling_map", exact = TRUE)

.summarize_canonical_parameter <- function(parameter, method, level) {
  if (method == "MCMC") {
    dimensions <- dim(parameter$draws)
    margin <- if (is.null(dimensions) || length(dimensions) == 1L) NULL else
      seq.int(2L, length(dimensions))
    return(summarize_draws_parameter(parameter$draws, margin, level))
  }
  switch(parameter$distribution,
    "normal variational factor" = summarize_normal_parameter(
      parameter$mean, parameter$covariance, level),
    "independent normal variational factors" = summarize_independent_normals(
      parameter$mean, parameter$covariance, level),
    "inverse-gamma variational factor" = summarize_inverse_gamma_parameter(
      parameter$shape, parameter$scale, level),
    "gamma variational factor" = summarize_gamma_parameter(
      parameter$shape, parameter$rate, level),
    "Dirichlet variational factor" = summarize_dirichlet_parameter(
      parameter$concentration, level),
    "inverse-Wishart variational factor" = summarize_inverse_wishart_mean(
      parameter$df, parameter$scale, level),
    stop("Unknown canonical parameter distribution: ", parameter$distribution))
}

.transform_prediction_X <- function(newdata, fit, map) {
  if (fit$model == "LRM") return(.validate_scaling_X(newdata, map))
  if (!is.list(newdata) || !length(newdata))
    stop("Grouped prediction requires a nonempty list of design matrices")
  lapply(newdata, .validate_scaling_X, map = map)
}

.restore_prediction_result <- function(result, map) {
  location <- intersect(c("mean", "median", "lower", "upper"),
                        names(result$summary))
  result$summary[location] <- lapply(result$summary[location], function(value)
    map$response_center + map$response_scale * value)
  if ("sd" %in% names(result$summary))
    result$summary$sd <- map$response_scale * result$summary$sd
  if (!is.null(result$draws))
    result$draws <- lapply(result$draws, function(value)
      map$response_center + map$response_scale * value)
  result$metadata$output_scale <- "input"
  result
}

.restore_in_sample_result <- function(result, map) {
  result$statistics$mse <- result$statistics$mse * map$response_scale^2
  if (!is.null(result$fitted_values))
    result$fitted_values <- lapply(result$fitted_values, function(value)
      map$response_center + map$response_scale * value)
  if (!is.null(result$residuals))
    result$residuals <- lapply(result$residuals, function(value)
      map$response_scale * value)
  result$metadata$output_scale <- "input"
  result
}

.restore_ppc_result <- function(result, map) {
  location <- "mean"
  spread <- "sd"
  squared <- c("between_group_mean_variance",
               "mean_within_group_variance", "ols_mse")
  present <- intersect(location, names(result$discrepancies$repdata))
  result$discrepancies$repdata[present] <- lapply(
    result$discrepancies$repdata[present], function(value)
      map$response_center + map$response_scale * value)
  present <- intersect(spread, names(result$discrepancies$repdata))
  result$discrepancies$repdata[present] <- lapply(
    result$discrepancies$repdata[present], function(value)
      map$response_scale * value)
  present <- intersect(squared, names(result$discrepancies$repdata))
  result$discrepancies$repdata[present] <- lapply(
    result$discrepancies$repdata[present], function(value)
      map$response_scale^2 * value)
  t0 <- result$discrepancies$t0
  t0[intersect(names(t0), location)] <- map$response_center +
    map$response_scale * t0[intersect(names(t0), location)]
  t0[intersect(names(t0), spread)] <-
    map$response_scale * t0[intersect(names(t0), spread)]
  t0[intersect(names(t0), squared)] <-
    map$response_scale^2 * t0[intersect(names(t0), squared)]
  result$discrepancies$t0 <- t0
  if (!is.null(result$y_rep))
    result$y_rep <- map$response_center + map$response_scale * result$y_rep
  result$metadata$output_scale <- "input"
  result
}

.ppc_effective_mcse <- function(fit, result) {
  discrepancies <- result$discrepancies
  indicators <- sweep(as.matrix(discrepancies$repdata), 2L,
                      discrepancies$t0, `>=`)
  S <- nrow(indicators)
  if (fit$method != "MCMC") {
    discrepancies$effective_sample_size <-
      stats::setNames(rep.int(S, ncol(indicators)), colnames(indicators))
    discrepancies$mcse <- discrepancies$mcse_iid
    discrepancies$mcse_method <- "independent_variational_draws"
    result$discrepancies <- discrepancies
    return(result)
  }
  if (!requireNamespace("posterior", quietly = TRUE))
    stop("Package 'posterior' is required for MCMC PPC MCSE")
  indices <- result$metadata$draw_indices
  if (is.null(indices) || length(indices) != S)
    stop("MCMC PPC metadata must retain one draw index per replicate")
  origin <- attr(.fit_backend(fit)$chains, "draw_origin", exact = TRUE)
  if (is.null(origin)) {
    selected <- data.frame(chain = rep.int(1L, S), iteration = indices)
  } else {
    selected <- origin[indices, c("chain", "iteration"), drop = FALSE]
  }
  rows_by_chain <- split(seq_len(S), selected$chain)
  rows_by_chain <- lapply(rows_by_chain, function(rows)
    rows[order(selected$iteration[rows])])
  balanced <- min(lengths(rows_by_chain))
  effective <- stats::setNames(rep.int(NA_real_, ncol(indicators)),
                               colnames(indicators))
  if (balanced >= 4L) {
    array <- array(NA_real_, c(balanced, length(rows_by_chain), ncol(indicators)),
                   dimnames = list(iteration = NULL, chain = NULL,
                                   variable = colnames(indicators)))
    for (chain in seq_along(rows_by_chain))
      array[, chain, ] <- indicators[rows_by_chain[[chain]][seq_len(balanced)], ,
                                     drop = FALSE]
    constant <- vapply(seq_len(ncol(indicators)), function(j)
      length(unique(indicators[, j])) == 1L, logical(1L))
    if (any(!constant)) {
      estimates <- suppressWarnings(as.data.frame(posterior::summarise_draws(
        posterior::as_draws_array(array[, , !constant, drop = FALSE]),
        "ess_mean")))
      effective[estimates$variable] <- estimates$ess_mean
    }
    effective[constant] <- balanced * length(rows_by_chain)
  }
  valid <- is.finite(effective) & effective > 0
  mcse <- stats::setNames(rep.int(NA_real_, length(effective)), names(effective))
  mcse[valid] <- sqrt(discrepancies$ppp[valid] *
                        (1 - discrepancies$ppp[valid]) / effective[valid])
  discrepancies$effective_sample_size <- effective
  discrepancies$mcse <- mcse
  discrepancies$mcse_method <- if (balanced >= 4L)
    "indicator_sequence_ess_mean" else
    "unavailable_fewer_than_four_selected_draws_per_chain"
  result$discrepancies <- discrepancies
  result$metadata$mcse_balanced_draws_per_chain <- balanced
  result$metadata$mcse_chains <- length(rows_by_chain)
  result
}

.log_likelihood_unit_counts <- function(log_likelihood, fit, arguments, map) {
  if (log_likelihood$unit == "observation")
    return(rep.int(1, nrow(log_likelihood$values)))
  if (fit$model == "LRM") {
    groups <- arguments$groups
    if (is.null(groups))
      stop("Internal error: grouped LRM likelihood has no group labels")
    labels <- unique(groups)
    return(as.numeric(tabulate(match(groups, labels), nbins = length(labels))))
  }
  counts <- map$training$group_sizes
  if (length(counts) != nrow(log_likelihood$values))
    stop("Internal error: scaling group sizes do not match log likelihood")
  as.numeric(counts)
}

.restore_log_likelihood_result <- function(log_likelihood, fit, arguments,
                                           map) {
  counts <- .log_likelihood_unit_counts(log_likelihood, fit, arguments, map)
  log_likelihood$values <- sweep(
    log_likelihood$values, 1L, counts * log(map$response_scale), `-`)
  log_likelihood$metadata$density_scale <- "input"
  log_likelihood$metadata$jacobian_log_adjustment_per_unit <-
    -counts * log(map$response_scale)
  log_likelihood
}

.restore_predictive_criteria_result <- function(result, fit, arguments, map,
                                                return_log_likelihood) {
  restored <- lapply(result$log_likelihood, .restore_log_likelihood_result,
                     fit = fit, arguments = arguments, map = map)
  for (name in names(result$criteria)) {
    likelihood_name <- if (grepl("observation", name, fixed = TRUE))
      "observation" else "group"
    result$criteria[[name]] <- waic_from_loglik(restored[[likelihood_name]])
  }
  result$metadata$density_scale <- "input"
  if (return_log_likelihood) result$log_likelihood <- restored else
    result$log_likelihood <- NULL
  result
}

#' @rdname fitting-and-draws
#' @export
summarize_posterior <- function(fit, ...) {
  UseMethod("summarize_posterior")
}

#' @rdname fitting-and-draws
#' @export
summarize_posterior.vb_model_fit <- function(fit, ...) {
  result <- .invoke_fit_backend(fit, "summarize_posterior", list(...))
  if (!is.null(.fit_scaling_map(fit))) {
    result$parameters <- lapply(
      fit$parameters, .summarize_canonical_parameter,
      method = fit$method, level = result$interval$level)
  }
  result
}

#' @rdname predictive-evaluation
#' @export
evaluate_in_sample_fit <- function(fit, ...) {
  UseMethod("evaluate_in_sample_fit")
}

#' @rdname predictive-evaluation
#' @export
evaluate_in_sample_fit.vb_model_fit <- function(fit, ...) {
  result <- .invoke_fit_backend(fit, "evaluate_in_sample_fit", list(...))
  map <- .fit_scaling_map(fit)
  if (is.null(map)) result else .restore_in_sample_result(result, map)
}

#' @rdname predictive-evaluation
#' @export
posterior_predictive_check <- function(fit, ...) {
  UseMethod("posterior_predictive_check")
}

#' @rdname predictive-evaluation
#' @export
posterior_predictive_check.vb_model_fit <- function(fit, ...) {
  result <- .invoke_fit_backend(fit, "posterior_predictive_check", list(...))
  result <- .ppc_effective_mcse(fit, result)
  map <- .fit_scaling_map(fit)
  if (is.null(map)) result else .restore_ppc_result(result, map)
}

#' @rdname predictive-evaluation
#' @export
predictive_criteria <- function(fit, ...) {
  UseMethod("predictive_criteria")
}

#' @rdname predictive-evaluation
#' @export
predictive_criteria.vb_model_fit <- function(fit, ...) {
  arguments <- list(...)
  map <- .fit_scaling_map(fit)
  if (is.null(map))
    return(.invoke_fit_backend(fit, "predictive_criteria", arguments))
  requested <- if (is.null(arguments$return_log_likelihood)) FALSE else
    arguments$return_log_likelihood
  if (!is.logical(requested) || length(requested) != 1L || is.na(requested))
    return(.invoke_fit_backend(fit, "predictive_criteria", arguments))
  backend_arguments <- arguments
  backend_arguments$return_log_likelihood <- TRUE
  result <- .invoke_fit_backend(fit, "predictive_criteria", backend_arguments)
  .restore_predictive_criteria_result(
    result, fit, arguments, map, return_log_likelihood = requested)
}

#' @rdname predictive-evaluation
#' @export
predict.vb_model_fit <- function(object, newdata, ...) {
  map <- .fit_scaling_map(object)
  backend_data <- if (is.null(map)) newdata else
    .transform_prediction_X(newdata, object, map)
  result <- .invoke_fit_backend(
    object, "predict", c(list(backend_data), list(...)))
  if (is.null(map)) result else .restore_prediction_result(result, map)
}

#' @rdname predictive-evaluation
#' @export
posterior_predict <- function(fit, newdata, ...) {
  UseMethod("posterior_predict")
}

#' @rdname predictive-evaluation
#' @export
posterior_predict.vb_model_fit <- function(fit, newdata, ...) {
  map <- .fit_scaling_map(fit)
  backend_data <- if (is.null(map)) newdata else
    .transform_prediction_X(newdata, fit, map)
  result <- .invoke_fit_backend(
    fit, "posterior_predict", c(list(backend_data), list(...)))
  if (is.null(map)) result else .restore_prediction_result(result, map)
}

#' @rdname fitting-and-draws
#' @export
fit_orchestration <- function(fit) {
  if (!inherits(fit, "vb_model_fit"))
    stop("fit must be produced by fit_model()")
  orchestration <- attr(fit, "orchestration", exact = TRUE)
  if (is.null(orchestration))
    stop("This fit used one chain or one variational start")
  if (inherits(orchestration, "mcmc_multichain_fit")) {
    out <- list(
      schema_version = "fit-orchestration-1",
      model = fit$model, method = fit$method, type = "multichain",
      status = orchestration$status,
      chain_metadata = orchestration$chain_metadata,
      diagnostics = orchestration$diagnostics,
      block_history = orchestration$block_history,
      estimated_storage_mb = orchestration$estimated_storage_mb)
    if (!is.null(orchestration$label_reference))
      out$label_reference <- orchestration$label_reference
  } else if (inherits(orchestration, "vi_multistart_fit")) {
    out <- list(
      schema_version = "fit-orchestration-1",
      model = fit$model, method = fit$method, type = "multistart",
      status = orchestration$status,
      summary = orchestration$summary,
      selected_start = orchestration$selected_start,
      require_converged = orchestration$require_converged,
      estimated_storage_mb = orchestration$estimated_storage_mb)
    if (!is.null(orchestration$label_reference))
      out$label_reference <- orchestration$label_reference
    if (fit$model == "CHLRM") {
      rows <- lapply(seq_along(orchestration$attempts), function(start) {
        attempt <- orchestration$attempts[[start]]
        if (!identical(attempt$status, "complete")) return(NULL)
        result <- if (!is.null(attempt$aligned_results))
          attempt$aligned_results else attempt$results
        probabilities <- result$rho
        labels <- max.col(probabilities, ties.method = "first")
        data.frame(
          start = start, group = seq_len(nrow(probabilities)),
          component = labels,
          assignment_probability = probabilities[
            cbind(seq_len(nrow(probabilities)), labels)],
          stringsAsFactors = FALSE)
      })
      out$start_classifications <- do.call(rbind, rows)
      rownames(out$start_classifications) <- NULL
    }
  } else {
    stop("Unknown orchestration object stored in fit")
  }
  structure(out, class = c("vb_fit_orchestration", "list"))
}
