# Explicit, training-only affine scaling shared by every inference engine.

data_scaling <- function(
    predictors = "all", response = TRUE,
    center_predictors = TRUE, scale_predictors = TRUE,
    center_response = TRUE, scale_response = TRUE,
    intercept_column = 1L, zero_scale_tolerance = sqrt(.Machine$double.eps)) {
  if (!(identical(predictors, "all") || is.numeric(predictors) ||
        is.character(predictors)))
    stop("predictors must be 'all', column indices or column names")
  if (!identical(predictors, "all") &&
      (anyDuplicated(predictors) || anyNA(predictors)))
    stop("predictors must be unique and non-missing")
  logical_fields <- list(
    response = response, center_predictors = center_predictors,
    scale_predictors = scale_predictors, center_response = center_response,
    scale_response = scale_response)
  logical_fields <- Map(.logical_scalar, logical_fields, names(logical_fields))
  if (is.null(intercept_column)) {
    if (center_predictors || (response && center_response))
      stop("Centering requires an explicit intercept_column")
  } else if (!(is.numeric(intercept_column) || is.character(intercept_column)) ||
             length(intercept_column) != 1L || is.na(intercept_column)) {
    stop("intercept_column must be NULL, one index or one column name")
  }
  zero_scale_tolerance <- .positive_finite_scalar(
    zero_scale_tolerance, "zero_scale_tolerance")
  .new_vb_spec(c(list(
    schema_version = "data-scaling-spec-1", predictors = predictors,
    intercept_column = intercept_column,
    zero_scale_tolerance = zero_scale_tolerance), logical_fields),
    "vb_data_scaling_spec")
}

.scaling_data_parts <- function(data) {
  if (!is.list(data) || is.null(names(data)) || anyDuplicated(names(data)) ||
      !all(c("X", "y") %in% names(data)))
    stop("data must be a uniquely named list containing X and y")
  grouped <- is.list(data$X) && !is.data.frame(data$X)
  if (grouped) {
    info <- validate_grouped_data(data$X, data$y)
    X <- lapply(data$X, as.matrix)
    y <- lapply(data$y, as.numeric)
    list(grouped = TRUE, X = X, y = y, p = info$p,
         group_sizes = info$nj)
  } else {
    X <- as.matrix(data$X)
    y <- as.numeric(data$y)
    if (!is.numeric(X) || !nrow(X) || !ncol(X) || nrow(X) != length(y) ||
        any(!is.finite(X)) || any(!is.finite(y)))
      stop("Ungrouped data require a finite nonempty X matrix and matching y")
    storage.mode(X) <- "double"
    list(grouped = FALSE, X = list(X), y = list(y), p = ncol(X),
         group_sizes = nrow(X))
  }
}

.resolve_scaling_column <- function(column, names, p, label, allow_null = FALSE) {
  if (is.null(column) && allow_null) return(NULL)
  if (is.character(column)) {
    if (is.null(names)) stop(label, " uses names but X has no column names")
    index <- match(column, names)
    if (anyNA(index)) stop(label, " contains unknown column names")
  } else {
    if (!is.numeric(column) || any(!is.finite(column)) ||
        any(column != floor(column)) || any(column < 1L | column > p))
      stop(label, " contains invalid column indices")
    index <- as.integer(column)
  }
  index
}

.coefficient_affine_map <- function(predictor_center, predictor_scale,
                                    response_center, response_scale,
                                    intercept_column) {
  p <- length(predictor_center)
  restore_matrix <- diag(response_scale / predictor_scale, p)
  restore_offset <- numeric(p)
  to_model_matrix <- diag(predictor_scale / response_scale, p)
  to_model_offset <- numeric(p)
  if (!is.null(intercept_column)) {
    others <- setdiff(seq_len(p), intercept_column)
    restore_offset[intercept_column] <- response_center
    restore_matrix[intercept_column, others] <-
      -predictor_center[others] * response_scale / predictor_scale[others]
    to_model_offset[intercept_column] <- -response_center / response_scale
    to_model_matrix[intercept_column, others] <-
      predictor_center[others] / response_scale
  }
  list(
    restore_offset = restore_offset,
    restore_matrix = restore_matrix,
    to_model_offset = to_model_offset,
    to_model_matrix = to_model_matrix)
}

fit_scaling <- function(scaling, data) {
  if (!inherits(scaling, "vb_data_scaling_spec"))
    stop("scaling must be created by data_scaling()")
  parts <- .scaling_data_parts(data)
  reference_names <- colnames(parts$X[[1L]])
  for (index in seq_along(parts$X))
    if (!identical(colnames(parts$X[[index]]), reference_names))
      stop("All design matrices must use identical column names and order")
  intercept <- .resolve_scaling_column(
    scaling$intercept_column, reference_names, parts$p,
    "intercept_column", allow_null = TRUE)
  if (!is.null(intercept) && any(vapply(parts$X, function(X)
      any(abs(X[, intercept] - 1) > scaling$zero_scale_tolerance), logical(1L))))
    stop("intercept_column must contain ones in all training rows")
  predictors <- if (identical(scaling$predictors, "all"))
    setdiff(seq_len(parts$p), intercept) else .resolve_scaling_column(
      scaling$predictors, reference_names, parts$p, "predictors")
  if (length(intersect(predictors, intercept)))
    stop("predictors cannot include intercept_column")

  pooled_X <- do.call(rbind, parts$X)
  pooled_y <- unlist(parts$y, use.names = FALSE)
  predictor_center <- numeric(parts$p)
  predictor_scale <- rep(1, parts$p)
  if (length(predictors)) {
    raw_center <- colMeans(pooled_X[, predictors, drop = FALSE])
    raw_scale <- apply(pooled_X[, predictors, drop = FALSE], 2L, stats::sd)
    if (any(!is.finite(raw_scale)) ||
        any(raw_scale <= scaling$zero_scale_tolerance))
      stop("Selected predictors must have nonzero finite training standard deviations")
    if (scaling$center_predictors) predictor_center[predictors] <- raw_center
    if (scaling$scale_predictors) predictor_scale[predictors] <- raw_scale
  }
  response_center <- if (scaling$response && scaling$center_response)
    mean(pooled_y) else 0
  response_scale <- if (scaling$response && scaling$scale_response)
    stats::sd(pooled_y) else 1
  if (!is.finite(response_scale) ||
      response_scale <= scaling$zero_scale_tolerance)
    stop("Response must have a nonzero finite training standard deviation")
  affine <- .coefficient_affine_map(
    predictor_center, predictor_scale, response_center, response_scale,
    intercept)
  structure(c(list(
    schema_version = "data-scaling-map-1",
    specification = scaling,
    coefficient_count = parts$p,
    column_names = reference_names,
    intercept_column = intercept,
    scaled_predictor_columns = predictors,
    unscaled_predictor_columns = setdiff(seq_len(parts$p), c(intercept, predictors)),
    predictor_center = predictor_center,
    predictor_scale = predictor_scale,
    response_center = response_center,
    response_scale = response_scale,
    training = list(
      observation_count = length(pooled_y),
      grouped = parts$grouped,
      group_sizes = as.integer(parts$group_sizes))), affine),
    class = c("vb_data_scaling_map", "list"))
}

.validate_scaling_X <- function(X, map) {
  X <- as.matrix(X)
  if (!is.numeric(X) || ncol(X) != map$coefficient_count ||
      any(!is.finite(X)))
    stop("X does not match the fitted scaling map")
  if (!is.null(map$column_names) && !identical(colnames(X), map$column_names))
    stop("X column names or order differ from the fitted scaling map")
  if (!is.null(map$intercept_column) && any(
      abs(X[, map$intercept_column] - 1) >
        map$specification$zero_scale_tolerance))
    stop("The fitted intercept column must remain equal to one")
  transformed <- sweep(X, 2L, map$predictor_center, `-`)
  transformed <- sweep(transformed, 2L, map$predictor_scale, `/`)
  if (!is.null(map$intercept_column))
    transformed[, map$intercept_column] <- 1
  transformed
}

transform_data <- function(data, map) {
  if (!inherits(map, "vb_data_scaling_map"))
    stop("map must be created by fit_scaling()")
  parts <- .scaling_data_parts(data)
  if (parts$p != map$coefficient_count)
    stop("data coefficient count differs from the fitted scaling map")
  X <- lapply(parts$X, .validate_scaling_X, map = map)
  y <- lapply(parts$y, function(value)
    (value - map$response_center) / map$response_scale)
  out <- data
  out$X <- if (parts$grouped) X else X[[1L]]
  out$y <- if (parts$grouped) y else y[[1L]]
  out
}

.transform_coefficient_vector <- function(value, map, direction) {
  vector <- is.null(dim(value))
  value <- if (vector) matrix(as.numeric(value), nrow = 1L) else as.matrix(value)
  if (ncol(value) != map$coefficient_count || any(!is.finite(value)))
    stop("Coefficient dimensions differ from the fitted scaling map")
  if (direction == "to_model") {
    out <- sweep(value %*% t(map$to_model_matrix), 2L,
                 map$to_model_offset, `+`)
  } else {
    out <- sweep(value %*% t(map$restore_matrix), 2L,
                 map$restore_offset, `+`)
  }
  if (vector) drop(out) else out
}

transform_coefficients <- function(value, map) {
  .transform_coefficient_vector(value, map, "to_model")
}

restore_coefficients <- function(value, map) {
  .transform_coefficient_vector(value, map, "to_input")
}

.transform_covariance <- function(value, matrix) {
  value <- as.matrix(value)
  if (!identical(dim(value), dim(matrix)) || any(!is.finite(value)))
    stop("Covariance dimensions differ from the fitted scaling map")
  symmetrize_numeric(matrix %*% value %*% t(matrix),
                     "transformed covariance")
}

transform_prior <- function(prior, map) {
  if (!inherits(prior, "vb_prior_spec"))
    stop("prior must be created by a prior constructor")
  if (!inherits(map, "vb_data_scaling_map"))
    stop("map must be created by fit_scaling()")
  if (length(prior$coefficient_mean) != map$coefficient_count)
    stop("Prior coefficient dimension differs from the fitted scaling map")
  mean <- transform_coefficients(prior$coefficient_mean, map)
  covariance <- .transform_covariance(
    prior$coefficient_covariance, map$to_model_matrix)
  if (prior$model == "LRM") return(lrm_prior(
    mean, covariance, prior$residual_degrees_of_freedom,
    prior$residual_variance_scale / map$response_scale^2))
  arguments <- list(
    coefficient_mean = mean,
    coefficient_covariance = covariance,
    coefficient_dispersion_df = prior$coefficient_dispersion_df,
    coefficient_dispersion_scale = .transform_covariance(
      prior$coefficient_dispersion_scale, map$to_model_matrix),
    residual_degrees_of_freedom = prior$residual_degrees_of_freedom,
    residual_variance_scale_shape = prior$residual_variance_scale_shape,
    residual_variance_scale_rate =
      prior$residual_variance_scale_rate * map$response_scale^2)
  if (prior$model == "CHLRM") {
    arguments$mixture_concentration <- prior$mixture_concentration
    return(do.call(chlrm_prior, arguments))
  }
  do.call(hlrm_prior, arguments)
}

restore_covariance <- function(value, map) {
  .transform_covariance(value, map$restore_matrix)
}

restore_residual_variance <- function(value, map) {
  value * map$response_scale^2
}

transform_residual_variance <- function(value, map) {
  value / map$response_scale^2
}

scaling_map <- function(fit) {
  if (!inherits(fit, "vb_model_fit"))
    stop("fit must be produced by fit_model()")
  map <- attr(fit, "scaling_map", exact = TRUE)
  if (is.null(map)) stop("This fit did not use data scaling")
  map
}

.transform_coefficient_array <- function(value, map, direction) {
  dimensions <- dim(value)
  if (is.null(dimensions) || length(dimensions) < 2L ||
      tail(dimensions, 1L) != map$coefficient_count ||
      !is.numeric(value) || any(!is.finite(value)))
    stop("Coefficient array must use its last dimension for coefficients")
  dim_names <- dimnames(value)
  rows <- matrix(value, ncol = map$coefficient_count)
  transformed <- .transform_coefficient_vector(rows, map, direction)
  dim(transformed) <- dimensions
  dimnames(transformed) <- dim_names
  transformed
}

.transform_covariance_collection <- function(value, matrix) {
  p <- nrow(matrix)
  if (is.list(value))
    return(lapply(value, .transform_covariance, matrix = matrix))
  dimensions <- dim(value)
  if (identical(dimensions, c(p, p)))
    return(.transform_covariance(value, matrix))
  if (is.null(dimensions) || length(dimensions) != 3L ||
      !identical(as.integer(dimensions[2:3]), c(p, p)) ||
      !is.numeric(value) || any(!is.finite(value)))
    stop("Covariance collection must be a matrix, list, or [draw,p,p] array")
  out <- array(NA_real_, dimensions, dimnames = dimnames(value))
  for (index in seq_len(dimensions[1L]))
    out[index, , ] <- .transform_covariance(value[index, , ], matrix)
  out
}

# Explicit user values are expressed in the data's input units. Initializers
# derived from the prior/data need no conversion because those inputs are
# already transformed before the backend sees them.
transform_initialization <- function(initialization, map, model, method) {
  if (!inherits(initialization, "vb_initialization_spec"))
    stop("initialization must be created by an initialization constructor")
  if (!inherits(map, "vb_data_scaling_map"))
    stop("map must be created by fit_scaling()")
  if (initialization$strategy != "values") return(initialization)
  model <- match.arg(model, c("LRM", "HLRM", "CHLRM"))
  method <- match.arg(method, c("MCMC", "CAVI", "SVI"))
  values <- unserialize(serialize(initialization$values, NULL))

  coefficient_vectors <- intersect(names(values), c(
    "coefficients", "coefficient_mean", "population_coefficients",
    "population_coefficient_mean"))
  for (name in coefficient_vectors)
    values[[name]] <- transform_coefficients(values[[name]], map)
  coefficient_arrays <- intersect(names(values), c(
    "group_coefficients", "component_coefficients",
    "group_coefficient_means", "component_coefficient_means"))
  for (name in coefficient_arrays)
    values[[name]] <- .transform_coefficient_array(
      values[[name]], map, "to_model")

  covariance_fields <- intersect(names(values), c(
    "population_covariance", "coefficient_covariance",
    "population_coefficient_covariance", "population_covariance_scale",
    "group_coefficient_covariances", "component_coefficient_covariances"))
  for (name in covariance_fields)
    values[[name]] <- .transform_covariance_collection(
      values[[name]], map$to_model_matrix)

  variance_fields <- intersect(names(values), c(
    "residual_variance", "variance_scale", "group_residual_variances",
    "component_residual_variances"))
  for (name in variance_fields)
    values[[name]] <- transform_residual_variance(values[[name]], map)
  inverse_gamma_scale_fields <- intersect(names(values), c(
    "residual_variance_rate", "group_residual_variance_rates",
    "component_residual_variance_rates"))
  for (name in inverse_gamma_scale_fields)
    values[[name]] <- values[[name]] / map$response_scale^2
  if ("variance_scale_rate" %in% names(values))
    values$variance_scale_rate <-
      values$variance_scale_rate * map$response_scale^2
  initialization_values(values)
}

restore_fit_parameters <- function(parameters, model, method, map) {
  if (!is.list(parameters) ||
      !identical(names(parameters), parameter_names_for_model(model)))
    stop("parameters do not match the selected model")
  method <- match.arg(method, c("MCMC", "CAVI", "SVI"))
  restored <- unserialize(serialize(parameters, NULL))
  coefficient_names <- if (model == "LRM") "coefficients" else c(
    "population_coefficients",
    if (model == "HLRM") "group_coefficients" else "component_coefficients")
  covariance_name <- if (model == "LRM") NULL else "population_covariance"
  residual_names <- switch(model,
    LRM = "residual_variance",
    HLRM = c("variance_scale", "group_variances"),
    CHLRM = c("variance_scale", "component_variances"))

  if (method == "MCMC") {
    for (name in coefficient_names)
      restored[[name]]$draws <- .transform_coefficient_array(
        restored[[name]]$draws, map, "to_input")
    if (!is.null(covariance_name))
      restored[[covariance_name]]$draws <- .transform_covariance_collection(
        restored[[covariance_name]]$draws, map$restore_matrix)
    for (name in residual_names)
      restored[[name]]$draws <- restore_residual_variance(
        restored[[name]]$draws, map)
    return(restored)
  }

  for (name in coefficient_names) {
    global <- name %in% c("coefficients", "population_coefficients")
    if (global) {
      original_dimensions <- dim(restored[[name]]$mean)
      restored[[name]]$mean <- restore_coefficients(
        as.numeric(restored[[name]]$mean), map)
      if (!is.null(original_dimensions))
        dim(restored[[name]]$mean) <- original_dimensions
    } else {
      restored[[name]]$mean <- .transform_coefficient_array(
        restored[[name]]$mean, map, "to_input")
    }
    restored[[name]]$covariance <- .transform_covariance_collection(
      restored[[name]]$covariance, map$restore_matrix)
  }
  if (!is.null(covariance_name))
    restored[[covariance_name]]$scale <- .transform_covariance(
      restored[[covariance_name]]$scale, map$restore_matrix)
  for (name in setdiff(residual_names, "variance_scale"))
    restored[[name]]$scale <- restore_residual_variance(
      restored[[name]]$scale, map)
  if ("variance_scale" %in% residual_names)
    restored$variance_scale$rate <-
      restored$variance_scale$rate / map$response_scale^2
  restored
}
