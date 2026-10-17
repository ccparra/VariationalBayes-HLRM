# Typed data-generating processes. One call generates one dataset; replication,
# scenario grids and paper-specific coefficient geometry live outside R/.

.simulation_spec <- function(fields, class) {
  structure(fields, class = c(class, "vb_simulation_spec", "vb_spec", "list"))
}

.validate_dgp_model <- function(model) {
  if (!inherits(model, "vb_model_spec"))
    stop("model must be created by lrm_model(), hlrm_model() or chlrm_model()")
  model
}

.named_parameter_list <- function(parameters) {
  if (!is.list(parameters) || !length(parameters) || is.null(names(parameters)) ||
      any(!nzchar(names(parameters))) || anyDuplicated(names(parameters)) ||
      any(vapply(parameters, is.null, logical(1L))))
    stop("parameters must be a nonempty named list without NULL entries")
  parameters
}

conditional_dgp <- function(model, parameters) {
  model <- .validate_dgp_model(model)
  parameters <- .named_parameter_list(parameters)
  required <- switch(model$model,
    LRM = c("coefficients", "residual_variance"),
    HLRM = c("group_coefficients", "group_variances"),
    CHLRM = c("component_coefficients", "component_variances",
              "mixture_weights"))
  optional <- if (model$model == "CHLRM") "allocations" else character()
  missing <- setdiff(required, names(parameters))
  unknown <- setdiff(names(parameters), c(required, optional))
  if (length(missing)) stop("Missing conditional parameters: ",
                            paste(missing, collapse = ", "))
  if (length(unknown)) stop("Unknown conditional parameters: ",
                            paste(unknown, collapse = ", "))
  .simulation_spec(list(
    schema_version = "dgp-spec-1", type = "conditional",
    model = model, parameters = parameters), "vb_conditional_dgp")
}

hierarchical_dgp <- function(model, parameters) {
  model <- .validate_dgp_model(model)
  if (model$model == "LRM")
    stop("hierarchical_dgp() applies only to HLRM and CHLRM")
  parameters <- .named_parameter_list(parameters)
  required <- c("population_coefficients", "population_covariance",
                "variance_scale", "residual_degrees_of_freedom")
  if (model$model == "CHLRM") required <- c(required, "mixture_weights")
  missing <- setdiff(required, names(parameters))
  unknown <- setdiff(names(parameters), required)
  if (length(missing)) stop("Missing hierarchical parameters: ",
                            paste(missing, collapse = ", "))
  if (length(unknown)) stop("Unknown hierarchical parameters: ",
                            paste(unknown, collapse = ", "))
  .simulation_spec(list(
    schema_version = "dgp-spec-1", type = "hierarchical",
    model = model, parameters = parameters), "vb_hierarchical_dgp")
}

prior_predictive_dgp <- function(model, prior) {
  model <- .validate_dgp_model(model)
  if (!inherits(prior, "vb_prior_spec") || !identical(prior$model, model$model))
    stop("prior must be a model-matched object from a prior constructor")
  .simulation_spec(list(
    schema_version = "dgp-spec-1", type = "prior_predictive",
    model = model, prior = prior), "vb_prior_predictive_dgp")
}

simulation_design <- function(
    distribution = c("independent_normal", "multivariate_normal", "mixed",
                     "grid", "fixed"),
    include_intercept = TRUE, mean = 0, sd = 1, covariance = NULL,
    binary_probabilities = numeric(), bounds = NULL, jitter_sd = 0,
    binary_position = c("last", "first"),
    fixed_train = NULL, fixed_test = NULL, column_names = NULL) {
  distribution <- match.arg(distribution)
  include_intercept <- .logical_scalar(include_intercept, "include_intercept")
  fields <- list(schema_version = "simulation-design-1",
                 distribution = distribution,
                 include_intercept = include_intercept)
  if (!is.null(column_names)) {
    if (!is.character(column_names) || !length(column_names) ||
        any(!nzchar(column_names)) || anyDuplicated(column_names))
      stop("column_names must contain unique nonempty names")
    fields$column_names <- column_names
  }
  if (distribution == "fixed") {
    if (is.null(fixed_train)) stop("fixed_train is required for a fixed design")
    fields$fixed_train <- fixed_train
    if (!is.null(fixed_test)) fields$fixed_test <- fixed_test
    return(.simulation_spec(fields, "vb_simulation_design"))
  }
  if (!is.numeric(mean) || !length(mean) || any(!is.finite(mean)))
    stop("mean must contain finite numeric values")
  fields$mean <- as.numeric(mean)
  if (distribution == "multivariate_normal") {
    if (is.null(covariance) || !is.matrix(covariance) ||
        !is.numeric(covariance) ||
        nrow(covariance) != ncol(covariance))
      stop("covariance must be a square numeric matrix")
    fields$covariance <- if (nrow(covariance)) .spd_matrix(
      covariance, nrow(covariance), "covariance") else
        matrix(numeric(), 0L, 0L)
  } else if (distribution %in% c("independent_normal", "mixed")) {
    if (distribution == "mixed") {
      fields$binary_position <- match.arg(binary_position)
      if (!is.numeric(binary_probabilities) ||
          any(!is.finite(binary_probabilities)) ||
          any(binary_probabilities <= 0 | binary_probabilities >= 1))
        stop("binary_probabilities must lie strictly between zero and one")
      fields$binary_probabilities <- as.numeric(binary_probabilities)
      if (!is.null(covariance)) {
        if (!is.matrix(covariance) || !is.numeric(covariance) ||
            nrow(covariance) != ncol(covariance) || !nrow(covariance))
          stop("covariance must be a nonempty square numeric matrix")
        fields$covariance <- .spd_matrix(
          covariance, nrow(covariance), "covariance")
      }
    }
    if (is.null(fields$covariance)) {
      if (!is.numeric(sd) || !length(sd) || any(!is.finite(sd)) || any(sd <= 0))
        stop("sd must contain positive finite values")
      fields$sd <- as.numeric(sd)
    }
  } else {
    bounds <- as.matrix(bounds)
    if (!is.numeric(bounds) || ncol(bounds) != 2L ||
        any(!is.finite(bounds)) || any(bounds[, 1L] >= bounds[, 2L]))
      stop("bounds must contain one increasing finite interval per covariate")
    if (!is.numeric(jitter_sd) || !length(jitter_sd) ||
        any(!is.finite(jitter_sd)) || any(jitter_sd < 0))
      stop("jitter_sd must contain nonnegative finite values")
    fields$bounds <- bounds
    fields$jitter_sd <- as.numeric(jitter_sd)
  }
  .simulation_spec(fields, "vb_simulation_design")
}

simulation_sample <- function(sample_size = NULL, group_count = NULL,
                              observations_per_group = NULL,
                              test_sample_size = NULL,
                              test_observations_per_group = NULL,
                              test_group_mode = c("same", "new"),
                              test_group_count = NULL) {
  lrm_mode <- !is.null(sample_size)
  grouped_mode <- !is.null(group_count) || !is.null(observations_per_group)
  if (lrm_mode == grouped_mode)
    stop("Specify either sample_size or group_count plus observations_per_group")
  fields <- list(schema_version = "simulation-sample-1")
  if (lrm_mode) {
    fields$sample_size <- validate_integer_scalar(sample_size, "sample_size")
    if (!is.null(test_sample_size))
      fields$test_sample_size <- validate_integer_scalar(
        test_sample_size, "test_sample_size")
  } else {
    fields$group_count <- validate_integer_scalar(group_count, "group_count")
    sizes <- observations_per_group
    if (length(sizes) == 1L) sizes <- rep(sizes, fields$group_count)
    fields$observations_per_group <- validate_integer_vector(
      sizes, "observations_per_group", fields$group_count)
    if (!is.null(test_observations_per_group)) {
      fields$test_group_mode <- match.arg(test_group_mode)
      count <- if (fields$test_group_mode == "same") fields$group_count else
        validate_integer_scalar(
          if (is.null(test_group_count)) fields$group_count else test_group_count,
          "test_group_count")
      test_sizes <- test_observations_per_group
      if (length(test_sizes) == 1L) test_sizes <- rep(test_sizes, count)
      fields$test_group_count <- count
      fields$test_observations_per_group <- validate_integer_vector(
        test_sizes, "test_observations_per_group", count)
    }
  }
  .simulation_spec(fields, "vb_simulation_sample")
}

.simulation_observation_count <- function(sample) {
  train <- if (!is.null(sample$sample_size)) sample$sample_size else
    sum(sample$observations_per_group)
  test <- if (!is.null(sample$test_sample_size)) sample$test_sample_size else
    if (!is.null(sample$test_observations_per_group))
      sum(sample$test_observations_per_group) else 0
  c(train = as.double(train), test = as.double(test))
}

.check_simulation_budget <- function(sample, coefficient_count,
                                     retain_observation_truth,
                                     max_output_mb) {
  if (!is.numeric(max_output_mb) || length(max_output_mb) != 1L ||
      !is.finite(max_output_mb) || max_output_mb <= 0)
    stop("max_output_mb must be a positive finite scalar")
  observations <- sum(.simulation_observation_count(sample))
  # X, y and group indices plus the conditional mean and innovation used while
  # generating y. Retaining observation truth keeps the latter two in output.
  transient_vectors <- 2
  retained_vectors <- if (retain_observation_truth) 2 else 0
  required_doubles <- observations *
    (coefficient_count + 2 + transient_vectors + retained_vectors)
  required_mb <- 8 * required_doubles / 1024^2
  if (required_mb > max_output_mb)
    stop(sprintf(
      paste0("Simulation arrays require approximately %.1f MB, exceeding ",
             "max_output_mb = %.1f; reduce the sample or increase the ",
             "explicit budget"), required_mb, max_output_mb))
  invisible(required_mb)
}

