.reference_for_public_fit <- function(reference, fit) {
  if (is.null(reference)) return(NULL)
  out <- unserialize(serialize(reference, NULL))
  map <- .fit_scaling_map(fit)
  out$metadata <- list(
    coordinate_scale = if (is.null(map)) "input" else "standardized_backend",
    scaling_map = map)
  out
}

.validate_reference_for_public_fit <- function(reference, fit) {
  .require_chlrm_fit(fit)
  specification <- fit$specifications$model
  coefficient_count <- specification$coefficient_count
  if (is.null(coefficient_count))
    coefficient_count <- .fit_training_stats(fit)$p
  .validate_chlrm_label_reference(
    reference, specification$component_count,
    coefficient_count)
  metadata <- reference$metadata
  if (!is.list(metadata) || is.null(metadata$coordinate_scale) ||
      !metadata$coordinate_scale %in% c("input", "standardized_backend"))
    stop(paste(
      "reference must carry coordinate metadata from",
      "fit_orchestration()"))
  map <- .fit_scaling_map(fit)
  expected_scale <- if (is.null(map)) "input" else "standardized_backend"
  if (!identical(metadata$coordinate_scale, expected_scale))
    stop("reference and fit use different coordinate scales")
  if (!is.null(map) && !isTRUE(all.equal(
      metadata$scaling_map, map, check.attributes = TRUE)))
    stop("reference and fit use different scaling maps")
  invisible(reference)
}

#' @rdname chlrm-accessors
#' @export
relabel_fit <- function(fit, reference) {
  .validate_reference_for_public_fit(reference, fit)
  reference <- unserialize(serialize(reference, NULL))
  backend <- .fit_backend(fit)
  native <- if (fit$method == "MCMC") backend$chains else backend$results
  aligned <- align_to_reference_chlrm(native, reference)

  out <- fit
  out$parameters <- canonical_fit_parameters("CHLRM", fit$method, aligned)
  out$latent_variables <- canonical_fit_latent_variables(
    "CHLRM", fit$method, aligned)
  map <- .fit_scaling_map(fit)
  if (!is.null(map))
    out$parameters <- restore_fit_parameters(
      out$parameters, "CHLRM", fit$method, map)
  out$metadata$label_reference <- list(
    applied = TRUE,
    algorithm = reference$algorithm,
    coordinate_scale = reference$metadata$coordinate_scale)

  reporting_backend <- backend$clone(deep = TRUE)
  reporting_backend$label_reference <- reference
  if (fit$method == "MCMC")
    reporting_backend$relabeled_chains <- aligned else
    reporting_backend$aligned_results <- aligned
  attr(out, "backend") <- reporting_backend
  attr(out, "label_reference") <- reference
  out
}
