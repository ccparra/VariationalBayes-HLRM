# Public posterior-draw extraction. The first dimension of every parameter is
# always the draw index. MCMC draws retain their chain/iteration provenance;
# CAVI/SVI draws are explicitly samples from the variational approximation.

.validate_posterior_draw_request <- function(n_draws, max_output_mb) {
  if (!is.null(n_draws))
    n_draws <- validate_integer_scalar(n_draws, "n_draws", 1L)
  if (!is.numeric(max_output_mb) || length(max_output_mb) != 1L ||
      !is.finite(max_output_mb) || max_output_mb <= 0)
    stop("max_output_mb must be a positive finite scalar")
  list(n_draws = n_draws, max_output_mb = as.numeric(max_output_mb))
}

.draw_first_subset <- function(x, indices) {
  dimensions <- dim(x)
  if (is.null(dimensions)) return(x[indices])
  arguments <- c(list(x, indices), rep(list(TRUE), length(dimensions) - 1L),
                 list(drop = FALSE))
  do.call(`[`, arguments)
}

.posterior_draw_storage <- function(fit, n_draws, include_latent,
                                    parameters = NULL) {
  if (!is.null(parameters)) {
    # Only selected MCMC blocks are copied. Include the four provenance columns
    # even for latent-only requests; this remains an output estimate, not RSS.
    values_per_draw <- 4 + sum(vapply(fit$parameters[parameters], function(x) {
      dimensions <- dim(x$draws)
      if (is.null(dimensions)) 1 else prod(as.double(dimensions[-1L]))
    }, numeric(1L)))
    if (include_latent)
      values_per_draw <- values_per_draw +
        nrow(fit$latent_variables$group_allocations)
    return(1.5 * 8 * values_per_draw * n_draws / 1024^2)
  }
  model <- fit$model
  specification <- fit$specifications$model
  p <- specification$coefficient_count
  if (is.null(p)) {
    coefficient_draws <- fit$parameters[[if (model == "LRM")
      "coefficients" else "population_coefficients"]]
    p <- if (fit$method == "MCMC") ncol(coefficient_draws$draws) else
      length(coefficient_draws$mean)
  }
  values_per_draw <- switch(model,
    LRM = p + 1,
    HLRM = {
      local <- fit$parameters$group_coefficients
      m <- if (fit$method == "MCMC") dim(local$draws)[2L] else
        nrow(local$mean)
      p + p * p + 1 + m * p + m
    },
    CHLRM = {
      components <- fit$parameters$component_coefficients
      K <- if (fit$method == "MCMC") dim(components$draws)[2L] else
        nrow(components$mean)
      m <- if (include_latent) {
        latent <- fit$latent_variables
        if (fit$method == "MCMC") nrow(latent$group_allocations) else
          nrow(latent$allocation_probabilities)
      } else 0
      p + p * p + 1 + K * p + K + K + m
    })
  # A modest overhead allowance covers list/data-frame bookkeeping and
  # temporary matrices created by the random generators.
  1.5 * 8 * as.double(values_per_draw) * n_draws / 1024^2
}

.sample_inverse_wishart_draws <- function(n_draws, df, scale) {
  p <- nrow(scale)
  out <- array(NA_real_, dim = c(n_draws, p, p))
  for (draw in seq_len(n_draws))
    out[draw, , ] <- LaplacesDemon::rinvwishart(nu = df, S = scale)
  out
}

.sample_independent_normals <- function(n_draws, means, covariances) {
  count <- nrow(means)
  p <- ncol(means)
  out <- array(NA_real_, dim = c(n_draws, count, p))
  for (index in seq_len(count))
    out[, index, ] <- mvtnorm::rmvnorm(
      n_draws, means[index, ], covariances[[index]], method = "chol")
  out
}

.sample_variational_parameters <- function(fit, n_draws, include_latent) {
  parameters <- fit$parameters
  coefficient_name <- if (fit$model == "LRM")
    "coefficients" else "population_coefficients"
  population <- parameters[[coefficient_name]]
  draws <- stats::setNames(list(mvtnorm::rmvnorm(
    n_draws, population$mean, population$covariance, method = "chol")),
    coefficient_name)

  if (fit$model == "LRM") {
    residual <- parameters$residual_variance
    draws$residual_variance <- 1 / stats::rgamma(
      n_draws, shape = residual$shape, rate = residual$scale)
    return(list(parameters = draws))
  }

  covariance <- parameters$population_covariance
  draws$population_covariance <- .sample_inverse_wishart_draws(
    n_draws, covariance$df, covariance$scale)
  variance_scale <- parameters$variance_scale
  draws$variance_scale <- stats::rgamma(
    n_draws, shape = variance_scale$shape, rate = variance_scale$rate)

  if (fit$model == "HLRM") {
    coefficients <- parameters$group_coefficients
    residuals <- parameters$group_variances
    draws$group_coefficients <- .sample_independent_normals(
      n_draws, coefficients$mean, coefficients$covariance)
    draws$group_variances <- vapply(seq_along(residuals$shape), function(j)
      1 / stats::rgamma(n_draws, residuals$shape[j], rate = residuals$scale[j]),
      numeric(n_draws))
    return(list(parameters = draws))
  }

  coefficients <- parameters$component_coefficients
  residuals <- parameters$component_variances
  weights <- parameters$mixture_weights
  draws$component_coefficients <- .sample_independent_normals(
    n_draws, coefficients$mean, coefficients$covariance)
  draws$component_variances <- vapply(seq_along(residuals$shape), function(k)
    1 / stats::rgamma(n_draws, residuals$shape[k], rate = residuals$scale[k]),
    numeric(n_draws))
  draws$mixture_weights <- gtools::rdirichlet(
    n_draws, weights$concentration)

  out <- list(parameters = draws)
  if (include_latent) {
    probabilities <- fit$latent_variables$allocation_probabilities
    allocations <- vapply(seq_len(nrow(probabilities)), function(j)
      sample.int(ncol(probabilities), n_draws, replace = TRUE,
                 prob = probabilities[j, ]), integer(n_draws))
    out$latent_variables <- list(group_allocations = allocations)
  }
  out
}

.extract_mcmc_parameters <- function(fit, indices, include_latent,
                                     parameters = names(fit$parameters)) {
  draws <- lapply(fit$parameters[parameters], function(parameter)
    .draw_first_subset(parameter$draws, indices))
  out <- list(parameters = draws)
  if (fit$model == "CHLRM" && include_latent) {
    # Stored internally as group x draw; public draws always put draw first.
    allocations <- fit$latent_variables$group_allocations[, indices,
                                                           drop = FALSE]
    out$latent_variables <- list(group_allocations = t(allocations))
  }
  out
}

#' @rdname fitting-and-draws
#' @export
posterior_draws <- function(fit, n_draws = NULL, seed = NULL,
                            include_latent = FALSE,
                            max_output_mb = 256, draw_indices = NULL,
                            parameters = NULL) {
  UseMethod("posterior_draws")
}

#' @rdname fitting-and-draws
#' @export
posterior_draws.vb_model_fit <- function(fit, n_draws = NULL, seed = NULL,
                                         include_latent = FALSE,
                                         max_output_mb = 256,
                                         draw_indices = NULL,
                                         parameters = NULL) {
  request <- .validate_posterior_draw_request(n_draws, max_output_mb)
  validate_rng_seed(seed)
  if (!is.logical(include_latent) || length(include_latent) != 1L ||
      is.na(include_latent)) stop("include_latent must be TRUE or FALSE")
  if (include_latent && fit$model != "CHLRM")
    stop("Latent allocation draws apply only to CHLRM")
  if (fit$method != "MCMC" &&
      (!is.null(draw_indices) || !is.null(parameters)))
    stop("draw_indices and parameters selection apply only to stored MCMC draws")
  parameter_names <- parameter_names_for_model(fit$model)
  if (!is.null(parameters)) {
    if (!is.character(parameters) || anyNA(parameters) ||
        anyDuplicated(parameters) || any(!parameters %in% parameter_names))
      stop("parameters must contain unique applicable parameter names")
    parameter_names <- parameter_names[parameter_names %in% parameters]
    if (!length(parameter_names) && !include_latent)
      stop("Select at least one parameter or include latent allocations")
  }

  if (fit$method == "MCMC") {
    coefficient_name <- if (fit$model == "LRM")
      "coefficients" else "population_coefficients"
    available <- nrow(fit$parameters[[coefficient_name]]$draws)
    if (!is.null(draw_indices)) {
      if (!is.null(request$n_draws))
        stop("draw_indices and n_draws are mutually exclusive")
      if (!is.numeric(draw_indices) || !is.null(dim(draw_indices)) ||
          !length(draw_indices) || anyNA(draw_indices) ||
          any(!is.finite(draw_indices)) ||
          any(draw_indices != floor(draw_indices)) ||
          any(draw_indices < 1 | draw_indices > available) ||
          any(diff(draw_indices) <= 0))
        stop("draw_indices must be strictly increasing stored MCMC indices in range")
      draw_indices <- as.integer(draw_indices)
    }
    count <- if (!is.null(draw_indices)) length(draw_indices) else
      if (is.null(request$n_draws)) available else request$n_draws
    if (count > available)
      stop("n_draws cannot exceed the stored MCMC draws")
  } else {
    count <- if (is.null(request$n_draws))
      fit$specifications$control$evaluation_draws else request$n_draws
  }
  estimated <- .posterior_draw_storage(
    fit, count, include_latent,
    parameters = if (is.null(parameters) && is.null(draw_indices)) NULL else
      parameter_names)
  if (estimated > request$max_output_mb)
    stop("Posterior draws exceed max_output_mb; reduce n_draws or increase the explicit budget")

  scope <- inference_rng_scope(seed)
  on.exit(scope$restore(), add = TRUE)
  if (fit$method == "MCMC") {
    indices <- if (is.null(draw_indices))
      prediction_draw_indices(available, count) else draw_indices
    sampled <- .extract_mcmc_parameters(
      fit, indices, include_latent, parameter_names)
    origin <- attr(.fit_backend(fit)$chains, "draw_origin", exact = TRUE)
    if (is.null(origin))
      origin <- data.frame(chain = rep.int(1L, available),
                           iteration = seq_len(available))
    draw_index <- data.frame(
      draw = seq_len(count), stored_draw = indices,
      chain = origin$chain[indices], iteration = origin$iteration[indices])
  } else {
    sampled <- .sample_variational_parameters(fit, count, include_latent)
    draw_index <- data.frame(draw = seq_len(count))
  }
  new_posterior_draws_result(
    fit, sampled, draw_index, seed, estimated, parameter_names)
}

.flatten_draw_parameter <- function(value, name, draw_count) {
  dimensions <- dim(value)
  if (is.null(dimensions)) {
    if (length(value) != draw_count)
      stop("Every posterior parameter must put the draw axis first")
    out <- matrix(as.numeric(value), ncol = 1L)
    colnames(out) <- name
    return(out)
  }
  if (dimensions[1L] != draw_count)
    stop("Every posterior parameter must put the draw axis first")
  out <- matrix(value, nrow = draw_count)
  if (length(dimensions) == 2L) {
    suffix <- as.character(seq_len(dimensions[2L]))
  } else {
    indices <- expand.grid(
      lapply(dimensions[-1L], seq_len), KEEP.OUT.ATTRS = FALSE)
    suffix <- apply(indices, 1L, paste, collapse = ",")
  }
  colnames(out) <- paste0(name, "[", suffix, "]")
  out
}

# Interoperability with posterior/bayesplot. This conversion is lossless with
# respect to parameter values and MCMC chain/iteration provenance; semantic
# tensor dimensions are encoded in standard bracketed variable names.
#' @rdname fitting-and-draws
#' @export
as_draws.vb_posterior_draws <- function(x, include_latent = FALSE, ...) {
  if (!inherits(x, "vb_posterior_draws"))
    stop("x must be produced by posterior_draws()")
  if (!is.logical(include_latent) || length(include_latent) != 1L ||
      is.na(include_latent)) stop("include_latent must be TRUE or FALSE")
  draw_count <- x$metadata$draw_count
  blocks <- Map(
    .flatten_draw_parameter, x$parameters, names(x$parameters),
    MoreArgs = list(draw_count = draw_count))
  if (include_latent) {
    if (is.null(x$latent_variables))
      stop("This posterior draw object has no latent variables")
    latent <- Map(
      .flatten_draw_parameter, x$latent_variables,
      names(x$latent_variables), MoreArgs = list(draw_count = draw_count))
    blocks <- c(blocks, latent)
  }
  if (!length(blocks))
    stop("No parameter draws selected; use include_latent = TRUE for latent-only draws")
  values <- as.data.frame(do.call(cbind, blocks), check.names = FALSE)
  index <- x$draw_index
  is_mcmc <- "chain" %in% names(index)
  if (is_mcmc) {
    # draws_df requires a rectangular sequential .iteration index. Keep the
    # original (possibly subsampled and gapped) provenance as ordinary columns.
    values$source_iteration <- index$iteration
    values$source_stored_draw <- index$stored_draw
  }
  values$.chain <- if (is_mcmc) index$chain else 1L
  values$.iteration <- if ("iteration" %in% names(index))
    index$iteration else index$draw
  values$.draw <- seq_len(draw_count)
  posterior::as_draws_df(values)
}

#' @rdname fitting-and-draws
#' @export
as_draws.vb_model_fit <- function(x, n_draws = NULL, seed = NULL,
                                  include_latent = FALSE,
                                  max_output_mb = 256,
                                  draw_indices = NULL, parameters = NULL, ...) {
  posterior::as_draws(posterior_draws(
    x, n_draws = n_draws, seed = seed, include_latent = include_latent,
    max_output_mb = max_output_mb, draw_indices = draw_indices,
    parameters = parameters), include_latent = include_latent)
}
