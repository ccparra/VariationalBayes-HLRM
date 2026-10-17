.dgp_dimension <- function(generator) {
  type <- generator$type; model <- generator$model; p <- generator$parameters
  coefficient_count <- if (type == "prior_predictive")
    length(generator$prior$coefficient_mean) else if (type == "conditional")
      switch(model$model,
        LRM = length(p$coefficients),
        HLRM = ncol(as.matrix(p$group_coefficients)),
        CHLRM = ncol(as.matrix(p$component_coefficients))) else
          length(p$population_coefficients)
  if (!is.null(model$coefficient_count) &&
      model$coefficient_count != coefficient_count)
    stop("DGP coefficient dimension does not match the model specification")
  coefficient_count
}

.validate_conditional_parameters <- function(generator, sample, p) {
  parameters <- generator$parameters; model <- generator$model$model
  if (model == "LRM") {
    parameters$coefficients <- as.numeric(parameters$coefficients)
    if (length(parameters$coefficients) != p ||
        any(!is.finite(parameters$coefficients)))
      stop("coefficients must be a finite coefficient vector")
    parameters$residual_variance <- .positive_scalar(
      parameters$residual_variance, "residual_variance")
  } else if (model == "HLRM") {
    parameters$group_coefficients <- as.matrix(parameters$group_coefficients)
    if (!identical(dim(parameters$group_coefficients),
                   c(sample$group_count, p)) ||
        any(!is.finite(parameters$group_coefficients)))
      stop("group_coefficients must be group_count by coefficient_count")
    variances <- parameters$group_variances
    if (length(variances) == 1L) variances <- rep(variances, sample$group_count)
    if (!is.numeric(variances) || length(variances) != sample$group_count ||
        any(!is.finite(variances)) || any(variances <= 0))
      stop("group_variances must be positive with one value per group")
    parameters$group_variances <- as.numeric(variances)
    if (!is.null(sample$test_group_mode) && sample$test_group_mode == "new")
      stop("A conditional HLRM has no population distribution for new groups")
  } else {
    K <- generator$model$component_count
    parameters$component_coefficients <- as.matrix(parameters$component_coefficients)
    if (!identical(dim(parameters$component_coefficients), c(K, p)) ||
        any(!is.finite(parameters$component_coefficients)))
      stop("component_coefficients must be component_count by coefficient_count")
    variances <- parameters$component_variances
    if (length(variances) == 1L) variances <- rep(variances, K)
    if (!is.numeric(variances) || length(variances) != K ||
        any(!is.finite(variances)) || any(variances <= 0))
      stop("component_variances must be positive with one value per component")
    parameters$component_variances <- as.numeric(variances)
    parameters$mixture_weights <- .validate_probability_vector(
      parameters$mixture_weights, K, "mixture_weights")
    if (!is.null(parameters$allocations)) {
      parameters$allocations <- validate_integer_vector(
        parameters$allocations, "allocations", sample$group_count)
      if (any(parameters$allocations > K))
        stop("allocations cannot exceed component_count")
    }
  }
  parameters
}

.validate_population_parameters <- function(parameters, model, p) {
  parameters$population_coefficients <- as.numeric(
    parameters$population_coefficients)
  if (length(parameters$population_coefficients) != p ||
      any(!is.finite(parameters$population_coefficients)))
    stop("population_coefficients must be a finite coefficient vector")
  parameters$population_covariance <- .spd_matrix(
    as.matrix(parameters$population_covariance), p, "population_covariance")
  parameters$variance_scale <- .positive_scalar(
    parameters$variance_scale, "variance_scale")
  parameters$residual_degrees_of_freedom <- .positive_scalar(
    parameters$residual_degrees_of_freedom, "residual_degrees_of_freedom")
  if (model$model == "CHLRM")
    parameters$mixture_weights <- .validate_probability_vector(
      parameters$mixture_weights, model$component_count, "mixture_weights")
  parameters
}

.draw_prior_parameters <- function(model, prior) {
  p <- length(prior$coefficient_mean)
  if (model$model == "LRM") return(list(
    coefficients = drop(mvtnorm::rmvnorm(
      1L, prior$coefficient_mean, prior$coefficient_covariance,
      method = "chol")),
    residual_variance = 1 / stats::rgamma(
      1L, prior$residual_degrees_of_freedom / 2,
      rate = prior$residual_degrees_of_freedom *
        prior$residual_variance_scale / 2)))
  out <- list(
    population_coefficients = drop(mvtnorm::rmvnorm(
      1L, prior$coefficient_mean, prior$coefficient_covariance,
      method = "chol")),
    population_covariance = LaplacesDemon::rinvwishart(
      nu = prior$coefficient_dispersion_df,
      S = prior$coefficient_dispersion_scale),
    variance_scale = stats::rgamma(
      1L, prior$residual_variance_scale_shape,
      rate = prior$residual_variance_scale_rate),
    residual_degrees_of_freedom = prior$residual_degrees_of_freedom)
  if (model$model == "CHLRM") out$mixture_weights <- drop(
    gtools::rdirichlet(1L, prior$mixture_concentration))
  out
}

.draw_group_truth <- function(parameters, count) {
  list(
    group_coefficients = unname(mvtnorm::rmvnorm(
      count, parameters$population_coefficients,
      parameters$population_covariance, method = "chol")),
    group_variances = 1 / stats::rgamma(
      count, parameters$residual_degrees_of_freedom / 2,
      rate = parameters$residual_degrees_of_freedom *
        parameters$variance_scale / 2))
}

.draw_component_truth <- function(parameters, model) {
  K <- model$component_count
  list(
    component_coefficients = unname(mvtnorm::rmvnorm(
      K, parameters$population_coefficients,
      parameters$population_covariance, method = "chol")),
    component_variances = 1 / stats::rgamma(
      K, parameters$residual_degrees_of_freedom / 2,
      rate = parameters$residual_degrees_of_freedom *
        parameters$variance_scale / 2))
}

.simulate_lrm_split <- function(design, n, p, parameters, split) {
  X <- .design_matrix(design, n, p, split)
  conditional_mean <- gaussian_mean_draw(X, parameters$coefficients)
  response <- gaussian_response_draw(
    X, parameters$coefficients, parameters$residual_variance)
  residuals <- response - conditional_mean
  list(data = list(X = X, y = response),
       conditional_mean = conditional_mean, residuals = residuals)
}

.simulate_grouped_split <- function(design, sizes, p, coefficients, variances,
                                    split) {
  count <- length(sizes)
  X <- lapply(seq_len(count), function(group)
    .design_matrix(design, sizes[group], p, split, group))
  conditional_mean <- lapply(seq_len(count), function(group)
    gaussian_mean_draw(X[[group]], coefficients[group, ]))
  y <- lapply(seq_len(count), function(group) gaussian_response_draw(
    X[[group]], coefficients[group, ], variances[group]))
  residuals <- Map(`-`, y, conditional_mean)
  list(data = list(X = X, y = y, group = lapply(
    seq_len(count), function(group) rep.int(group, sizes[group]))),
    conditional_mean = conditional_mean, residuals = residuals)
}

