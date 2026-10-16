# Conservative, dimension-based resource estimates for typed model fits.
#
# These estimates are guards against obviously unsafe requests, not promises
# about resident-set size. R object overhead, BLAS workspaces and allocator
# behaviour vary by platform, so the formulas deliberately include headroom.

.numeric_mb <- function(value_count, bytes_per_value = 8) {
  as.numeric(value_count) * bytes_per_value / 1024^2
}

.mcmc_values_per_draw <- function(model, dimensions, component_count = NULL) {
  p <- dimensions$p
  m <- if (is.null(dimensions$m)) 1L else dimensions$m
  switch(model,
    LRM = p + 3,
    HLRM = m * p + m + p + p^2 + 3,
    CHLRM = m * p + 2 * m + component_count * p +
      2 * component_count + p + p^2 + 3)
}

.mcmc_diagnostic_values <- function(model, dimensions, component_count = NULL,
                                    scope = "global") {
  p <- dimensions$p
  m <- if (is.null(dimensions$m)) 1L else dimensions$m
  switch(model,
    LRM = p + 1,
    HLRM = p + p * (p + 1) / 2 + 1 +
      if (scope == "all") m * (p + 1) else 0,
    CHLRM = component_count * (p + 1) + max(component_count - 1, 0) +
      p + p * (p + 1) / 2 + 1 +
      if (scope == "all") m * (p + 1) else 0)
}

.variational_state_values <- function(model, dimensions,
                                      component_count = NULL) {
  p <- dimensions$p
  m <- if (is.null(dimensions$m)) 1L else dimensions$m
  switch(model,
    LRM = p + p^2 + 2,
    HLRM = m * p + m * p^2 + p + 3 * p^2 + 2 * m + 3,
    CHLRM = m * component_count + component_count +
      component_count * p + component_count * p^2 + p + 3 * p^2 +
      2 * component_count + 3)
}

.resource_dimensions <- function(model, validated) {
  dimensions <- validated$dimensions
  observations <- if (model$model == "LRM") dimensions$n else
    sum(lengths(validated$y))
  list(
    observations = as.integer(observations),
    groups = if (model$model == "LRM") NA_integer_ else
      as.integer(dimensions$m),
    coefficients = as.integer(dimensions$p),
    components = if (model$model == "CHLRM")
      as.integer(model$component_count) else NA_integer_)
}

.estimate_fit_resources_validated <- function(model, validated, control,
                                              input_data, scaling = NULL) {
  family <- model$model
  dimensions <- validated$dimensions
  shape <- .resource_dimensions(model, validated)
  input_data_mb <- as.numeric(utils::object.size(input_data)) / 1024^2
  # The backend retains the input and sufficient statistics. Grouped models
  # additionally retain one p x p cross-product per group.
  sufficient_values <- if (family == "LRM") {
    dimensions$p^2 + dimensions$p
  } else {
    dimensions$m * (dimensions$p^2 + dimensions$p + 2)
  }
  sufficient_statistics_mb <- .numeric_mb(sufficient_values)
  workers <- if (control$method == "MCMC") control$chains else control$starts
  scaling_copy_mb <- if (is.null(scaling)) 0 else input_data_mb

  if (control$method == "MCMC") {
    values_per_draw <- .mcmc_values_per_draw(
      family, dimensions, model$component_count)
    retained_draws_mb <- .numeric_mb(
      control$post_warmup_draws * workers * values_per_draw)
    diagnostic_scope <- if (workers > 1L) control$diagnostic_scope else "global"
    diagnostic_mb <- if (workers > 1L) .numeric_mb(
      control$post_warmup_draws * workers *
        .mcmc_diagnostic_values(
          family, dimensions, model$component_count, diagnostic_scope)) else 0
    # During multichain finalization, aligned/combined draws coexist with the
    # individual chains. A single chain needs only a smaller result wrapper.
    result_copy_mb <- retained_draws_mb * if (workers > 1L) 1 else 0.25
    state_mb <- .numeric_mb(workers * .mcmc_values_per_draw(
      family, dimensions, model$component_count))
    history_mb <- .numeric_mb(
      max(1, ceiling(control$post_warmup_draws /
        if (workers > 1L) control$stopping$check_every else
          control$post_warmup_draws)) * 8)
    fit_storage_mb <- retained_draws_mb
  } else {
    state_values <- .variational_state_values(
      family, dimensions, model$component_count)
    state_mb <- .numeric_mb(state_values * workers)
    # Multistart retains initial states, fitted backends and result snapshots.
    # A single start still has a mutable state and an immutable result.
    result_copy_mb <- state_mb * if (workers > 1L) 3 else 1.5
    diagnostic_mb <- 0
    checks <- ceiling(control$max_iterations /
      control$stopping$check_every) + 1L
    history_mb <- .numeric_mb(checks * workers * 8)
    fit_storage_mb <- state_mb + history_mb
  }

  backend_data_mb <- workers * (input_data_mb + sufficient_statistics_mb)
  raw_peak_mb <- input_data_mb + scaling_copy_mb + backend_data_mb +
    fit_storage_mb + diagnostic_mb + result_copy_mb + state_mb
  estimated_peak_mb <- 1.35 * raw_peak_mb

  # This is not allocated by fit_model(). It quantifies a plausible later
  # request retaining posterior parameter draws, pointwise log likelihood and
  # y_rep simultaneously, so users can budget evaluation separately.
  evaluation_parameter_values <- if (control$method == "MCMC") {
    .mcmc_values_per_draw(family, dimensions, model$component_count)
  } else {
    .variational_state_values(family, dimensions, model$component_count)
  }
  evaluation_storage_mb <- .numeric_mb(control$evaluation_draws *
    (evaluation_parameter_values + 2 * shape$observations))

  components <- data.frame(
    component = c(
      "input_data", "scaling_copy", "backend_data_and_statistics",
      "retained_fit_state", "diagnostics", "result_copies"),
    estimated_mb = c(
      input_data_mb, scaling_copy_mb, backend_data_mb, fit_storage_mb,
      diagnostic_mb, result_copy_mb),
    stringsAsFactors = FALSE)
  budget <- control$max_memory_mb
  structure(list(
    schema_version = "fit-resource-estimate-1",
    model = family,
    method = control$method,
    dimensions = shape,
    workers = as.integer(workers),
    components = components,
    estimated_fit_storage_mb = as.numeric(fit_storage_mb),
    estimated_peak_memory_mb = as.numeric(estimated_peak_mb),
    estimated_evaluation_storage_mb = as.numeric(evaluation_storage_mb),
    max_memory_mb = as.numeric(budget),
    within_budget = estimated_peak_mb <= budget,
    assumptions = c(
      "Sequential chains/starts; no parallel worker duplication.",
      "Peak estimate includes 35% headroom for R object overhead.",
      "Evaluation storage is reported separately and is not allocated by fit_model().",
      "The estimate is a preflight guard, not a measured operating-system memory bound.")),
    class = c("vb_fit_resource_estimate", "list"))
}

#' Estimate resources required by a typed model fit
#'
#' @param model A model specification from `lrm_model()`, `hlrm_model()` or
#'   `chlrm_model()`.
#' @param data A list containing `X` and `y` in the representation required by
#'   the model.
#' @param control An inference control created by `mcmc_control()`,
#'   `cavi_control()` or `svi_control()`.
#' @param scaling Optional preprocessing specification from `data_scaling()`.
#'
#' @return A `vb_fit_resource_estimate` containing a conservative peak-memory
#'   estimate, retained-fit storage, potential downstream evaluation storage
#'   and the declared memory budget.
estimate_fit_resources <- function(model, data, control, scaling = NULL) {
  if (!inherits(model, "vb_model_spec"))
    stop("model must be created by lrm_model(), hlrm_model() or chlrm_model()")
  if (!inherits(control, "vb_inference_control"))
    stop("control must be created by mcmc_control(), cavi_control() or svi_control()")
  if (!is.null(scaling) && !inherits(scaling, "vb_data_scaling_spec"))
    stop("scaling must be NULL or created by data_scaling()")
  validated <- .validate_fit_data(model, data)
  .estimate_fit_resources_validated(model, validated, control, data, scaling)
}

#' @export
print.vb_fit_resource_estimate <- function(x, ...) {
  cat(sprintf(
    "%s %s resource estimate: %.1f MB peak (budget %.1f MB; %s)\n",
    x$model, x$method, x$estimated_peak_memory_mb, x$max_memory_mb,
    if (isTRUE(x$within_budget)) "within budget" else "exceeds budget"))
  cat(sprintf(
    "Retained fit: %.1f MB; potential downstream evaluation: %.1f MB\n",
    x$estimated_fit_storage_mb, x$estimated_evaluation_storage_mb))
  invisible(x)
}
