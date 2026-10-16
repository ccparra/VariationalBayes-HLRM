# Idiomatic methods for fitted models. These are deliberately thin views over
# the versioned package contracts: they never refit a model or expose the R6
# backend. Quantities that differ structurally by model are returned only where
# they apply.

.fit_training_stats <- function(fit) {
  backend <- .fit_backend(fit)
  field <- .backend_state_field(fit$model, fit$method)
  state <- backend[[field]]
  if (!is.list(state) || !is.list(state$stats))
    stop("The fitted backend does not retain valid training statistics")
  state$stats
}

.coefficient_names <- function(fit) {
  stats <- .fit_training_stats(fit)
  X <- if (fit$model == "LRM") stats$X else stats$X[[1L]]
  names <- colnames(X)
  if (is.null(names)) paste0("coefficient_", seq_len(ncol(X))) else names
}

.name_coefficient_vector <- function(value, names) {
  value <- as.numeric(value)
  names(value) <- names
  value
}

.name_coefficient_matrix <- function(value, column_names, row_prefix) {
  value <- as.matrix(value)
  colnames(value) <- column_names
  if (is.null(rownames(value)))
    rownames(value) <- paste0(row_prefix, seq_len(nrow(value)))
  value
}

.flatten_observation_blocks <- function(values) {
  if (!is.list(values) || !length(values))
    stop("Expected one numeric block per observed group")
  out <- unlist(values, use.names = FALSE)
  names(out) <- unlist(Map(function(value, group)
    paste0("group_", group, ".observation_", seq_along(value)),
    values, seq_along(values)), use.names = FALSE)
  out
}

print.vb_model_fit <- function(x, ...) {
  convergence <- x$diagnostics
  cat("Bayesian ", x$model, " fit\n", sep = "")
  cat("  Method: ", x$method, "\n", sep = "")
  cat("  Observations: ", nobs(x), "\n", sep = "")
  if (!is.null(convergence$status))
    cat("  Convergence: ", convergence$status, " (",
        convergence$reason, ")\n", sep = "")
  cat("  Elapsed: ", format(round(x$timing$total_seconds, 3), nsmall = 3),
      " seconds\n", sep = "")
  invisible(x)
}

summary.vb_model_fit <- function(object, ...) {
  summarize_posterior(object, ...)
}

coef.vb_model_fit <- function(object, ...) {
  posterior <- summarize_posterior(object, ...)
  coefficient_names <- .coefficient_names(object)
  if (object$model == "LRM")
    return(.name_coefficient_vector(
      posterior$parameters$coefficients$mean, coefficient_names))

  out <- list(population_coefficients = .name_coefficient_vector(
    posterior$parameters$population_coefficients$mean, coefficient_names))
  if (object$model == "HLRM") {
    out$group_coefficients <- .name_coefficient_matrix(
      posterior$parameters$group_coefficients$mean,
      coefficient_names, "group_")
  } else {
    out$component_coefficients <- .name_coefficient_matrix(
      posterior$parameters$component_coefficients$mean,
      coefficient_names, "component_")
  }
  out
}

fitted.vb_model_fit <- function(object, ...) {
  details <- evaluate_in_sample_fit(object, return_details = TRUE)
  if (object$model == "LRM") {
    value <- as.numeric(details$fitted_values[[1L]])
    names(value) <- paste0("observation_", seq_along(value))
    return(value)
  }
  .flatten_observation_blocks(details$fitted_values)
}

residuals.vb_model_fit <- function(object, ...) {
  details <- evaluate_in_sample_fit(object, return_details = TRUE)
  if (object$model == "LRM") {
    value <- as.numeric(details$residuals[[1L]])
    names(value) <- paste0("observation_", seq_along(value))
    return(value)
  }
  .flatten_observation_blocks(details$residuals)
}

nobs.vb_model_fit <- function(object, ...) {
  stats <- .fit_training_stats(object)
  if (object$model == "LRM") as.integer(stats$n) else
    as.integer(sum(stats$nj))
}
