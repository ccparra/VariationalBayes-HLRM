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

.expand_design_parameter <- function(value, q, name, positive = FALSE) {
  if (q == 0L) return(numeric())
  if (length(value) == 1L) value <- rep(value, q)
  if (!is.numeric(value) || length(value) != q || any(!is.finite(value)) ||
      (positive && any(value <= 0)))
    stop(name, " must have one ", if (positive) "positive " else "",
         "finite value per non-binary covariate")
  as.numeric(value)
}

.design_matrix <- function(design, n, coefficient_count, split, group = NULL) {
  if (design$distribution == "fixed") {
    source <- if (split == "train") design$fixed_train else design$fixed_test
    if (is.null(source)) stop("A fixed test design requires fixed_test")
    value <- if (is.list(source) && !is.data.frame(source)) source[[group]] else source
    value <- as.matrix(value)
    if (!is.numeric(value) || nrow(value) != n ||
        ncol(value) != coefficient_count || any(!is.finite(value)))
      stop("Fixed design dimensions must match the sample and coefficient count")
    if (!is.null(design$column_names)) {
      if (length(design$column_names) != ncol(value))
        stop("column_names length must equal the final design column count")
      colnames(value) <- design$column_names
    }
    return(value)
  }
  q <- coefficient_count - as.integer(design$include_intercept)
  if (q < 0L) stop("coefficient_count is incompatible with include_intercept")
  if (design$distribution == "multivariate_normal") {
    if (nrow(design$covariance) != q)
      stop("covariance dimension must match the model covariate count")
    location <- .expand_design_parameter(design$mean, q, "mean")
    covariates <- if (q) mvtnorm::rmvnorm(
      n, location, design$covariance, method = "chol") else matrix(numeric(), n, 0L)
  } else if (design$distribution == "independent_normal") {
    location <- .expand_design_parameter(design$mean, q, "mean")
    scale <- .expand_design_parameter(design$sd, q, "sd", TRUE)
    covariates <- if (q) matrix(stats::rnorm(n * q), n, q) else
      matrix(numeric(), n, 0L)
    if (q) covariates <- sweep(sweep(covariates, 2L, scale, `*`),
                               2L, location, `+`)
  } else if (design$distribution == "mixed") {
    binary_count <- length(design$binary_probabilities)
    numeric_count <- q - binary_count
    if (numeric_count < 0L)
      stop("binary_probabilities exceed the model covariate count")
    binary <- if (binary_count) vapply(
      design$binary_probabilities,
      function(probability) stats::rbinom(n, 1L, probability), numeric(n)) else
        matrix(numeric(), n, 0L)
    numeric_mean <- .expand_design_parameter(design$mean, numeric_count, "mean")
    if (!is.null(design$covariance)) {
      if (nrow(design$covariance) != numeric_count)
        stop("covariance dimension must match the continuous covariate count")
      continuous <- if (numeric_count) mvtnorm::rmvnorm(
        n, numeric_mean, design$covariance, method = "chol") else
          matrix(numeric(), n, 0L)
    } else {
      numeric_sd <- .expand_design_parameter(design$sd, numeric_count, "sd", TRUE)
      continuous <- if (numeric_count) matrix(
        stats::rnorm(n * numeric_count), n, numeric_count) else
          matrix(numeric(), n, 0L)
      if (numeric_count) continuous <- sweep(
        sweep(continuous, 2L, numeric_sd, `*`), 2L, numeric_mean, `+`)
    }
    covariates <- if (design$binary_position == "first")
      cbind(binary, continuous) else cbind(continuous, binary)
  } else {
    if (nrow(design$bounds) != q)
      stop("bounds must have one row per model covariate")
    if (!q) covariates <- matrix(numeric(), n, 0L) else {
      base <- floor(n^(1 / q)); sizes <- rep(max(1L, base), q)
      while (prod(sizes) < n) for (index in seq_len(q)) {
        if (prod(sizes) >= n) break
        sizes[index] <- sizes[index] + 1L
      }
      axes <- lapply(seq_len(q), function(index) seq(
        design$bounds[index, 1L], design$bounds[index, 2L],
        length.out = sizes[index]))
      covariates <- as.matrix(do.call(expand.grid, axes))[seq_len(n), , drop = FALSE]
      jitter <- .expand_design_parameter(
        design$jitter_sd, q, "jitter_sd")
      for (index in seq_len(q)) covariates[, index] <- pmin(pmax(
        covariates[, index] + stats::rnorm(n, sd = jitter[index]),
        design$bounds[index, 1L]), design$bounds[index, 2L])
    }
  }
  X <- if (design$include_intercept) cbind(intercept = 1, covariates) else covariates
  if (!is.null(design$column_names)) {
    if (length(design$column_names) != ncol(X))
      stop("column_names length must equal the final design column count")
    colnames(X) <- design$column_names
  } else {
    dimnames(X) <- NULL
  }
  X
}

.validate_probability_vector <- function(value, count, name) {
  value <- as.numeric(value)
  if (length(value) != count || any(!is.finite(value)) || any(value < 0) ||
      abs(sum(value) - 1) > 1e-10)
    stop(name, " must be a nonnegative vector summing to one")
  value
}

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

new_simulation_result <- function(model, train, truth, configuration,
                                  conditioning, parameter_source,
                                  design_source, test = NULL) {
  model <- match.arg(model, c("LRM", "HLRM", "CHLRM"))
  validate_simulation_split(train, model, "train")
  if (!is.null(test)) {
    validate_simulation_split(test, model, "test")
    train_p <- if (model == "LRM") ncol(train$X) else ncol(train$X[[1L]])
    test_p <- if (model == "LRM") ncol(test$X) else ncol(test$X[[1L]])
    if (train_p != test_p)
      stop("train and test must have the same design column count")
  }
  if (!is.list(truth) || !length(truth) || is.null(names(truth)) ||
      any(!nzchar(names(truth))) || anyDuplicated(names(truth)))
    stop("truth must be a nonempty uniquely named list")
  blocks <- list(truth = truth, configuration = configuration,
                 conditioning = conditioning,
                 parameter_source = parameter_source,
                 design_source = design_source)
  if (any(!vapply(blocks, is.list, logical(1L))) ||
      any(vapply(blocks, function(value)
        is.null(names(value)) || any(!nzchar(names(value))), logical(1L))))
    stop("simulation result metadata blocks must be lists")
  out <- list(schema_version = "simulation-result-2", model = model,
              train = train)
  if (!is.null(test)) out$test <- test
  out <- c(out, blocks)
  structure(out, class = c("vb_simulation_result", "list"))
}

.truth_source <- function(simulation, quantity) {
  source <- simulation$parameter_source
  if (grepl("conditional_mean|residuals", quantity)) return("generated")
  if (grepl("allocations|mixture_proportions", quantity))
    return(if (is.null(source$allocations)) "generated" else source$allocations)
  if (grepl("^(test_)?group_", quantity))
    return(if (is.null(source$local)) "generated" else source$local)
  if (grepl("^component_", quantity))
    return(if (is.null(source$components)) "generated" else source$components)
  source$global
}

.flatten_truth_value <- function(quantity, value, source) {
  if (is.list(value) && !is.data.frame(value)) {
    rows <- lapply(seq_along(value), function(group) {
      nested <- .flatten_truth_value(quantity, value[[group]], source)
      nested$index_3 <- nested$index_2
      nested$index_2 <- nested$index_1
      nested$index_1 <- group
      nested
    })
    return(do.call(rbind, rows))
  }
  if (!is.numeric(value) && !is.integer(value))
    stop("truth values must be numeric, integer or lists of numeric values")
  dimensions <- dim(value)
  if (is.null(dimensions)) {
    indices <- cbind(seq_along(value), NA_integer_, NA_integer_)
  } else {
    if (length(dimensions) > 3L)
      stop("tidy_truth() supports truth arrays with at most three dimensions")
    indices <- arrayInd(seq_along(value), dimensions)
    indices <- cbind(indices, matrix(
      NA_integer_, nrow(indices), 3L - ncol(indices)))
  }
  data.frame(
    quantity = quantity,
    index_1 = indices[, 1L], index_2 = indices[, 2L],
    index_3 = indices[, 3L], value = as.numeric(value), source = source,
    stringsAsFactors = FALSE)
}

tidy_truth <- function(simulation, include_observation_level = FALSE) {
  if (!inherits(simulation, "vb_simulation_result") ||
      !identical(simulation$schema_version, "simulation-result-2"))
    stop("simulation must be a simulation-result-2 object")
  include_observation_level <- .logical_scalar(
    include_observation_level, "include_observation_level")
  quantities <- names(simulation$truth)
  if (!include_observation_level) quantities <- quantities[
    !grepl("conditional_mean|residuals", quantities)]
  rows <- lapply(quantities, function(quantity) .flatten_truth_value(
    quantity, simulation$truth[[quantity]],
    .truth_source(simulation, quantity)))
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  structure(out, class = c("vb_tidy_truth", "data.frame"))
}

simulate_data <- function(generator, design, sample, seed = NULL,
                          retain_observation_truth = FALSE,
                          max_output_mb = 256) {
  if (!inherits(generator, "vb_simulation_spec") ||
      !inherits(design, "vb_simulation_design") ||
      !inherits(sample, "vb_simulation_sample"))
    stop("generator, design and sample must be created by simulation constructors")
  validate_rng_seed(seed)
  retain_observation_truth <- .logical_scalar(
    retain_observation_truth, "retain_observation_truth")
  model <- generator$model
  grouped <- model$model != "LRM"
  if (grouped && is.null(sample$group_count))
    stop("HLRM/CHLRM require a grouped simulation_sample()")
  if (!grouped && is.null(sample$sample_size))
    stop("LRM requires an ungrouped simulation_sample()")
  p <- .dgp_dimension(generator)
  estimated_output_mb <- .check_simulation_budget(
    sample, p, retain_observation_truth, max_output_mb)
  scope <- inference_rng_scope(seed)
  on.exit(scope$restore(), add = TRUE)

  parameter_source <- list()
  if (generator$type == "prior_predictive") {
    parameters <- .draw_prior_parameters(model, generator$prior)
    parameter_source$global <- "prior"
  } else {
    parameters <- generator$parameters
    parameter_source$global <- "conditioned"
  }

  has_test <- if (grouped) !is.null(sample$test_observations_per_group) else
    !is.null(sample$test_sample_size)
  truth <- list()
  if (model$model == "LRM") {
    if (generator$type == "hierarchical") stop("LRM has no hierarchical DGP")
    parameters <- .validate_conditional_parameters(
      conditional_dgp(model, parameters), sample, p)
    train <- .simulate_lrm_split(
      design, sample$sample_size, p, parameters, "train")
    truth <- parameters
    if (retain_observation_truth) truth <- c(truth, list(
      train_conditional_mean = train$conditional_mean,
      train_residuals = train$residuals))
    test <- NULL
    if (has_test) {
      test_generated <- .simulate_lrm_split(
        design, sample$test_sample_size, p, parameters, "test")
      test <- test_generated$data
      if (retain_observation_truth) {
        truth$test_conditional_mean <- test_generated$conditional_mean
        truth$test_residuals <- test_generated$residuals
      }
    }
    train <- train$data
  } else if (model$model == "HLRM") {
    if (generator$type == "conditional") {
      parameters <- .validate_conditional_parameters(generator, sample, p)
      locals <- parameters
      parameter_source$local <- "conditioned"
      population <- list()
    } else {
      population <- .validate_population_parameters(parameters, model, p)
      locals <- .draw_group_truth(population, sample$group_count)
      parameter_source$local <- "hierarchy"
    }
    train_generated <- .simulate_grouped_split(
      design, sample$observations_per_group, p, locals$group_coefficients,
      locals$group_variances, "train")
    truth <- c(population, locals)
    if (retain_observation_truth) truth <- c(truth, list(
      train_conditional_mean = train_generated$conditional_mean,
      train_residuals = train_generated$residuals))
    train <- train_generated$data; test <- NULL
    if (has_test) {
      test_locals <- if (sample$test_group_mode == "same") locals else
        .draw_group_truth(population, sample$test_group_count)
      test_generated <- .simulate_grouped_split(
        design, sample$test_observations_per_group, p,
        test_locals$group_coefficients, test_locals$group_variances, "test")
      test <- test_generated$data
      if (sample$test_group_mode == "new") {
        truth$test_group_coefficients <- test_locals$group_coefficients
        truth$test_group_variances <- test_locals$group_variances
      }
      if (retain_observation_truth) {
        truth$test_conditional_mean <- test_generated$conditional_mean
        truth$test_residuals <- test_generated$residuals
      }
    }
  } else {
    if (generator$type == "conditional") {
      components <- .validate_conditional_parameters(generator, sample, p)
      population <- list()
      parameter_source$components <- "conditioned"
    } else {
      population <- .validate_population_parameters(parameters, model, p)
      components <- c(.draw_component_truth(population, model),
                      list(mixture_weights = population$mixture_weights))
      parameter_source$components <- "hierarchy"
    }
    K <- model$component_count
    allocations <- components$allocations
    if (is.null(allocations)) allocations <- sample.int(
      K, sample$group_count, replace = TRUE, prob = components$mixture_weights)
    local_coefficients <- components$component_coefficients[allocations, , drop = FALSE]
    local_variances <- components$component_variances[allocations]
    train_generated <- .simulate_grouped_split(
      design, sample$observations_per_group, p, local_coefficients,
      local_variances, "train")
    component_truth <- components[c(
      "component_coefficients", "component_variances")]
    if (!length(population)) component_truth$mixture_weights <-
      components$mixture_weights
    parameter_source$allocations <- if (!is.null(components$allocations))
      "conditioned" else "mixture_distribution"
    truth <- c(population, component_truth,
      list(allocations = allocations,
           realized_mixture_proportions = tabulate(allocations, K) / length(allocations)))
    if (retain_observation_truth) truth <- c(truth, list(
      train_conditional_mean = train_generated$conditional_mean,
      train_residuals = train_generated$residuals))
    train <- train_generated$data; test <- NULL
    if (has_test) {
      test_allocations <- if (sample$test_group_mode == "same") allocations else
        sample.int(K, sample$test_group_count, replace = TRUE,
                   prob = components$mixture_weights)
      test_generated <- .simulate_grouped_split(
        design, sample$test_observations_per_group, p,
        components$component_coefficients[test_allocations, , drop = FALSE],
        components$component_variances[test_allocations], "test")
      test <- test_generated$data
      if (sample$test_group_mode == "new") {
        truth$test_allocations <- test_allocations
        truth$test_realized_mixture_proportions <-
          tabulate(test_allocations, K) / length(test_allocations)
      }
      if (retain_observation_truth) {
        truth$test_conditional_mean <- test_generated$conditional_mean
        truth$test_residuals <- test_generated$residuals
      }
    }
  }

  conditioning <- list(
    dgp = generator$type,
    global_parameters = parameter_source$global,
    test_units = if (!has_test) "not_generated" else if (!grouped) "new_observations" else
      paste0(sample$test_group_mode, "_groups"))
  design_source <- list(
    distribution = design$distribution,
    random = design$distribution != "fixed",
    include_intercept = design$include_intercept)
  new_simulation_result(
    model$model, train, truth,
    configuration = list(
      generator = generator, design = design, sample = sample, seed = seed,
      retain_observation_truth = retain_observation_truth,
      max_output_mb = max_output_mb,
      estimated_peak_memory_mb = estimated_output_mb),
    conditioning = conditioning,
    parameter_source = parameter_source,
    design_source = design_source, test = test)
}
