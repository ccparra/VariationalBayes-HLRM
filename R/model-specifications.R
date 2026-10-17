# Typed, model-level specifications. These objects contain no data-dependent
# estimates and no backend notation; translation to the current kernels lives
# in fit-model.R.

.new_vb_spec <- function(fields, classes) {
  structure(fields, class = c(classes, "vb_spec", "list"))
}

.optional_count <- function(x, name) {
  if (is.null(x)) return(NULL)
  validate_integer_scalar(x, name, 1L)
}

.logical_scalar <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x))
    stop(name, " must be TRUE or FALSE")
  x
}

.positive_finite_scalar <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x <= 0)
    stop(name, " must be a positive finite scalar")
  as.numeric(x)
}

#' @rdname model-specifications
#' @export
lrm_model <- function(coefficient_count = NULL) {
  .new_vb_spec(list(
    schema_version = "model-spec-1",
    model = "LRM",
    coefficient_count = .optional_count(coefficient_count, "coefficient_count"),
    likelihood = "gaussian_linear_regression"
  ), c("vb_lrm_model_spec", "vb_model_spec"))
}

#' @rdname model-specifications
#' @export
hlrm_model <- function(coefficient_count = NULL) {
  .new_vb_spec(list(
    schema_version = "model-spec-1",
    model = "HLRM",
    coefficient_count = .optional_count(coefficient_count, "coefficient_count"),
    likelihood = "gaussian_hierarchical_linear_regression"
  ), c("vb_hlrm_model_spec", "vb_model_spec"))
}

#' @rdname model-specifications
#' @export
chlrm_model <- function(component_count, coefficient_count = NULL) {
  .new_vb_spec(list(
    schema_version = "model-spec-1",
    model = "CHLRM",
    component_count = validate_integer_scalar(
      component_count, "component_count", 1L),
    coefficient_count = .optional_count(coefficient_count, "coefficient_count"),
    likelihood = "clustered_gaussian_hierarchical_linear_regression"
  ), c("vb_chlrm_model_spec", "vb_model_spec"))
}

.validate_prior_vector <- function(x, name) {
  x <- as.numeric(x)
  if (!length(x) || any(!is.finite(x)))
    stop(name, " must be a nonempty finite numeric vector")
  x
}

#' @rdname model-specifications
#' @export
lrm_prior <- function(coefficient_mean, coefficient_covariance,
                      residual_degrees_of_freedom,
                      residual_variance_scale) {
  coefficient_mean <- .validate_prior_vector(
    coefficient_mean, "coefficient_mean")
  p <- length(coefficient_mean)
  coefficient_covariance <- .spd_matrix(
    as.matrix(coefficient_covariance), p, "coefficient_covariance")
  residual_degrees_of_freedom <- .positive_scalar(
    residual_degrees_of_freedom, "residual_degrees_of_freedom")
  residual_variance_scale <- .positive_scalar(
    residual_variance_scale, "residual_variance_scale")
  .new_vb_spec(list(
    schema_version = "prior-spec-1",
    model = "LRM",
    coefficient_mean = coefficient_mean,
    coefficient_covariance = coefficient_covariance,
    residual_degrees_of_freedom = residual_degrees_of_freedom,
    residual_variance_scale = residual_variance_scale,
    parameterization = list(
      coefficients = "multivariate_normal_covariance",
      residual_variance = "inverse_gamma_shape_scale")
  ), c("vb_lrm_prior_spec", "vb_prior_spec"))
}

.hierarchical_prior <- function(model, coefficient_mean,
                                coefficient_covariance,
                                coefficient_dispersion_df,
                                coefficient_dispersion_scale,
                                residual_degrees_of_freedom,
                                residual_variance_scale_shape,
                                residual_variance_scale_rate,
                                mixture_concentration = NULL) {
  coefficient_mean <- .validate_prior_vector(
    coefficient_mean, "coefficient_mean")
  p <- length(coefficient_mean)
  coefficient_covariance <- .spd_matrix(
    as.matrix(coefficient_covariance), p, "coefficient_covariance")
  coefficient_dispersion_scale <- .spd_matrix(
    as.matrix(coefficient_dispersion_scale), p,
    "coefficient_dispersion_scale")
  coefficient_dispersion_df <- .positive_scalar(
    coefficient_dispersion_df, "coefficient_dispersion_df")
  if (coefficient_dispersion_df <= p - 1L)
    stop("coefficient_dispersion_df must exceed coefficient_count - 1 for a proper inverse-Wishart prior")
  residual_degrees_of_freedom <- .positive_scalar(
    residual_degrees_of_freedom, "residual_degrees_of_freedom")
  residual_variance_scale_shape <- .positive_scalar(
    residual_variance_scale_shape, "residual_variance_scale_shape")
  residual_variance_scale_rate <- .positive_scalar(
    residual_variance_scale_rate, "residual_variance_scale_rate")
  out <- list(
    schema_version = "prior-spec-1",
    model = model,
    coefficient_mean = coefficient_mean,
    coefficient_covariance = coefficient_covariance,
    coefficient_dispersion_df = coefficient_dispersion_df,
    coefficient_dispersion_scale = coefficient_dispersion_scale,
    residual_degrees_of_freedom = residual_degrees_of_freedom,
    residual_variance_scale_shape = residual_variance_scale_shape,
    residual_variance_scale_rate = residual_variance_scale_rate,
    parameterization = list(
      population_coefficients = "multivariate_normal_covariance",
      coefficient_dispersion = "inverse_wishart_scale",
      local_residual_variance = "inverse_gamma_shape_scale",
      residual_variance_scale = "gamma_shape_rate")
  )
  if (model == "CHLRM") {
    mixture_concentration <- .validate_prior_vector(
      mixture_concentration, "mixture_concentration")
    if (any(mixture_concentration <= 0))
      stop("mixture_concentration entries must be positive")
    out$mixture_concentration <- mixture_concentration
    out$parameterization$mixture_weights <- "dirichlet_concentration"
  }
  .new_vb_spec(out, c(
    if (model == "HLRM") "vb_hlrm_prior_spec" else "vb_chlrm_prior_spec",
    "vb_prior_spec"))
}

#' @rdname model-specifications
#' @export
hlrm_prior <- function(coefficient_mean, coefficient_covariance,
                       coefficient_dispersion_df,
                       coefficient_dispersion_scale,
                       residual_degrees_of_freedom,
                       residual_variance_scale_shape,
                       residual_variance_scale_rate) {
  .hierarchical_prior(
    "HLRM", coefficient_mean, coefficient_covariance,
    coefficient_dispersion_df, coefficient_dispersion_scale,
    residual_degrees_of_freedom, residual_variance_scale_shape,
    residual_variance_scale_rate)
}

#' @rdname model-specifications
#' @export
chlrm_prior <- function(coefficient_mean, coefficient_covariance,
                        coefficient_dispersion_df,
                        coefficient_dispersion_scale,
                        residual_degrees_of_freedom,
                        residual_variance_scale_shape,
                        residual_variance_scale_rate,
                        mixture_concentration) {
  .hierarchical_prior(
    "CHLRM", coefficient_mean, coefficient_covariance,
    coefficient_dispersion_df, coefficient_dispersion_scale,
    residual_degrees_of_freedom, residual_variance_scale_shape,
    residual_variance_scale_rate, mixture_concentration)
}

.inverse_gamma_moments <- function(shape, scale) {
  mean_exists <- shape > 1
  variance_exists <- shape > 2
  list(
    mean = if (mean_exists) scale / (shape - 1) else NA_real_,
    variance = if (variance_exists)
      scale^2 / ((shape - 1)^2 * (shape - 2)) else NA_real_,
    mean_exists = mean_exists,
    variance_exists = variance_exists,
    parameterization = "inverse-gamma shape-scale")
}

.gamma_moments <- function(shape, rate) list(
  mean = shape / rate,
  variance = shape / rate^2,
  mean_exists = TRUE,
  variance_exists = TRUE,
  parameterization = "gamma shape-rate")

.dirichlet_moments <- function(concentration) {
  total <- sum(concentration)
  mean <- concentration / total
  covariance <- -tcrossprod(concentration) /
    (total^2 * (total + 1))
  diag(covariance) <- concentration * (total - concentration) /
    (total^2 * (total + 1))
  list(mean = mean, covariance = covariance,
       mean_exists = TRUE, covariance_exists = TRUE,
       parameterization = "Dirichlet concentration")
}

#' @rdname model-specifications
#' @export
prior_moments <- function(prior) {
  if (!inherits(prior, "vb_prior_spec"))
    stop("prior must be created by lrm_prior(), hlrm_prior() or chlrm_prior()")
  coefficient_moments <- list(
    mean = prior$coefficient_mean,
    covariance = prior$coefficient_covariance,
    mean_exists = TRUE,
    covariance_exists = TRUE,
    parameterization = "multivariate normal covariance")
  if (prior$model == "LRM") {
    residual <- .inverse_gamma_moments(
      prior$residual_degrees_of_freedom / 2,
      prior$residual_degrees_of_freedom *
        prior$residual_variance_scale / 2)
    return(structure(list(
      schema_version = "prior-moments-1", model = "LRM",
      parameters = list(
        coefficients = coefficient_moments,
        residual_variance = residual)),
      class = c("vb_prior_moments", "list")))
  }

  p <- length(prior$coefficient_mean)
  covariance_mean_exists <- prior$coefficient_dispersion_df > p + 1
  covariance_mean <- if (covariance_mean_exists)
    prior$coefficient_dispersion_scale /
      (prior$coefficient_dispersion_df - p - 1) else
        matrix(NA_real_, p, p)
  population_covariance <- list(
    mean = covariance_mean,
    mean_exists = covariance_mean_exists,
    parameterization = "inverse-Wishart degrees-of-freedom and scale",
    note = "Elementwise variances are not summarized by this contract")
  variance_scale <- .gamma_moments(
    prior$residual_variance_scale_shape,
    prior$residual_variance_scale_rate)

  nu <- prior$residual_degrees_of_freedom
  local_variance_mean_exists <- nu > 2
  local_variance_variance_exists <- nu > 4
  expected_xi <- variance_scale$mean
  expected_xi_squared <- variance_scale$variance + expected_xi^2
  local_variance_mean <- if (local_variance_mean_exists)
    nu * expected_xi / (nu - 2) else NA_real_
  local_variance_second <- if (local_variance_variance_exists)
    nu^2 * expected_xi_squared / ((nu - 2) * (nu - 4)) else NA_real_
  local_variance <- list(
    mean = local_variance_mean,
    variance = if (local_variance_variance_exists)
      local_variance_second - local_variance_mean^2 else NA_real_,
    mean_exists = local_variance_mean_exists,
    variance_exists = local_variance_variance_exists,
    parameterization = paste0(
      "marginal under inverse-gamma residual variance conditional on ",
      "gamma variance_scale"))
  local_coefficients <- list(
    mean = prior$coefficient_mean,
    covariance = if (covariance_mean_exists)
      prior$coefficient_covariance + covariance_mean else
        matrix(NA_real_, p, p),
    mean_exists = TRUE,
    covariance_exists = covariance_mean_exists,
    parameterization = "marginal normal hierarchy")

  parameters <- list(
    population_coefficients = coefficient_moments,
    population_covariance = population_covariance,
    variance_scale = variance_scale)
  local_name <- if (prior$model == "HLRM")
    "group_coefficients" else "component_coefficients"
  variance_name <- if (prior$model == "HLRM")
    "group_variance" else "component_variance"
  parameters[[local_name]] <- local_coefficients
  parameters[[variance_name]] <- local_variance
  if (prior$model == "CHLRM") parameters$mixture_weights <-
    .dirichlet_moments(prior$mixture_concentration)
  structure(list(
    schema_version = "prior-moments-1", model = prior$model,
    parameters = parameters), class = c("vb_prior_moments", "list"))
}

