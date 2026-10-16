# Semantic accessors for hierarchical layers. They expose only quantities that
# exist in the fitted family and avoid borrowing ambiguous fixed/random-effect
# names from a different model parameterization.

.require_hierarchical_fit <- function(fit, models = c("HLRM", "CHLRM"),
                                      operation = "This operation") {
  if (!inherits(fit, "vb_model_fit"))
    stop("fit must be produced by fit_model()")
  if (!fit$model %in% models)
    stop(operation, " applies only to ", paste(models, collapse = " or "),
         " fits")
  invisible(fit)
}

.hierarchical_summary <- function(fit, ...) {
  .require_hierarchical_fit(fit)
  if (fit$model == "CHLRM")
    summarize_posterior(fit, include_incidence = FALSE, ...) else
    summarize_posterior(fit, ...)
}

population_parameters <- function(fit, ...) {
  UseMethod("population_parameters")
}

population_parameters.vb_model_fit <- function(fit, ...) {
  summary <- .hierarchical_summary(fit, ...)
  summary$parameters[c(
    "population_coefficients", "population_covariance", "variance_scale")]
}

group_parameters <- function(fit, ...) {
  UseMethod("group_parameters")
}

group_parameters.vb_model_fit <- function(fit, ...) {
  .require_hierarchical_fit(fit, "HLRM", "group_parameters()")
  summary <- summarize_posterior(fit, ...)
  summary$parameters[c("group_coefficients", "group_variances")]
}

group_deviations <- function(fit, ...) {
  UseMethod("group_deviations")
}

group_deviations.vb_model_fit <- function(fit, ...) {
  .require_hierarchical_fit(fit, "HLRM", "group_deviations()")
  summary <- summarize_posterior(fit, ...)
  group_mean <- as.matrix(summary$parameters$group_coefficients$mean)
  population_mean <- as.numeric(
    summary$parameters$population_coefficients$mean)
  deviations <- sweep(group_mean, 2L, population_mean, `-`)
  names <- .coefficient_names(fit)
  colnames(deviations) <- names
  rownames(deviations) <- paste0("group_", seq_len(nrow(deviations)))
  attr(deviations, "estimand") <-
    "posterior mean of group coefficient minus population coefficient"
  deviations
}
