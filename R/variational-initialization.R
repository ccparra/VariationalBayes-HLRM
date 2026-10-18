.variational_initialization_payload <- function(initialization_set) {
  payload <- attr(initialization_set, "serialized_states", exact = TRUE)
  if (!is.list(payload) || length(payload) != initialization_set$starts ||
      any(!vapply(payload, is.raw, logical(1L))))
    stop("initialization_set does not contain valid serialized states")
  lapply(payload, unserialize)
}

.validate_variational_initialization_set <- function(
    initialization_set, model, data, starts, scaling_map) {
  if (!inherits(initialization_set, "vb_variational_initialization_set") ||
      !identical(initialization_set$schema_version,
                 "variational-initialization-set-1"))
    stop("initialization must be a valid variational initialization set")
  if (!identical(initialization_set$model, model$model) ||
      !identical(initialization_set$component_count,
                 model$component_count))
    stop("initialization_set belongs to a different model specification")
  expected_dimensions <- .vi_data_dimensions(
    model$model, data$y, data$X,
    if (model$model == "CHLRM") model$component_count else NULL)
  if (!identical(initialization_set$dimensions, expected_dimensions))
    stop("initialization_set is incompatible with the data dimensions")
  if (!identical(initialization_set$starts, as.integer(starts)))
    stop("control starts must match initialization_set starts")

  stored_map <- initialization_set$scaling_map
  if (xor(is.null(stored_map), is.null(scaling_map)))
    stop("initialization_set and fit must use the same scaling convention")
  if (!is.null(stored_map) && !isTRUE(all.equal(
      stored_map, scaling_map, check.attributes = TRUE)))
    stop("initialization_set was created with a different scaling map")

  states <- .variational_initialization_payload(initialization_set)
  structure(list(
    schema_version = "vi-initial-states-1",
    model = initialization_set$model,
    dimensions = initialization_set$dimensions,
    states = states,
    metadata = initialization_set$metadata,
    initialization = initialization_set$initialization$strategy
  ), class = "vi_initial_states")
}

#' @rdname model-specifications
#' @export
variational_initialization_set <- function(
    model, data, prior,
    initialization = initialization_from_prior(),
    starts = 1L, scaling = NULL, seed) {
  if (!inherits(model, "vb_model_spec"))
    stop("model must be created by lrm_model(), hlrm_model() or chlrm_model()")
  if (!inherits(initialization, "vb_initialization_spec") ||
      inherits(initialization, "vb_variational_initialization_set"))
    stop("initialization must be a base initialization specification")
  if (!is.null(scaling) && !inherits(scaling, "vb_data_scaling_spec"))
    stop("scaling must be NULL or created by data_scaling()")
  starts <- validate_integer_scalar(starts, "starts", 1L)
  if (missing(seed) || is.null(seed))
    stop("seed is required to construct a reproducible initialization set")
  validate_rng_seed(seed)

  input_data <- .validate_fit_data(model, data)
  .validate_prior_for_model(prior, model, input_data$dimensions)
  scaling_map <- if (is.null(scaling)) NULL else fit_scaling(scaling, data)
  working_data <- if (is.null(scaling_map)) data else
    transform_data(data, scaling_map)
  working_prior <- if (is.null(scaling_map)) prior else
    transform_prior(prior, scaling_map)
  working_initialization <- if (is.null(scaling_map)) initialization else
    transform_initialization(initialization, scaling_map, model$model, "CAVI")
  validated <- .validate_fit_data(model, working_data)
  .validate_prior_for_model(working_prior, model, validated$dimensions)

  if (initialization$strategy == "kmeans" && model$model != "CHLRM")
    stop("initialization_kmeans() applies only to CHLRM")
  if (initialization$strategy == "ols" && model$model == "CHLRM")
    stop("CHLRM requires initialization_kmeans() rather than an unclustered OLS state")

  lrm_ols <- NULL
  grouped_ols <- NULL
  if (working_initialization$strategy %in% c("ols", "kmeans")) {
    if (model$model == "LRM")
      lrm_ols <- .lrm_ols_initialization(
        validated, working_initialization$jitter)
    else
      grouped_ols <- .group_ols_initialization(
        validated, working_initialization$jitter)
  }
  native_prior <- .prior_to_native(
    working_prior, model, validated$dimensions, grouped_ols)
  control <- cavi_control(
    max_iterations = 1L, evaluation_draws = 2L, starts = 1L)
  native_control <- .control_to_native(control, seed)
  if (model$model == "CHLRM")
    native_control$initialization <- if (
      working_initialization$strategy == "kmeans") "kmeans_ols" else
      "random_responsibilities"
  initialization_seeds <- derive_independent_seeds(starts, seed)
  native <- .typed_vi_initial_states(
    model, validated, native_prior, native_control,
    working_initialization, "CAVI", starts, initialization_seeds,
    lrm_ols, grouped_ols)

  out <- structure(list(
    schema_version = "variational-initialization-set-1",
    model = model$model,
    component_count = model$component_count,
    dimensions = native$dimensions,
    starts = starts,
    initialization = initialization,
    scaling = scaling,
    scaling_map = scaling_map,
    metadata = transform(
      native$metadata,
      coordinate_scale = if (is.null(scaling_map)) "input" else
        "standardized_backend")
  ), class = c(
    "vb_variational_initialization_set", "vb_initialization_spec",
    "vb_spec", "list"))
  attr(out, "serialized_states") <- lapply(native$states, serialize, NULL)
  out
}
