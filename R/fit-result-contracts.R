new_fit_result <- function(model, method, representation, parameters,
                           diagnostics, total_seconds, cumulative_seconds,
                           evaluation_iterations = NULL,
                           update_seconds = NULL, objective_seconds = NULL,
                           metadata = list()) {
  if (!model %in% c("LRM", "HLRM", "CHLRM")) stop("unknown model")
  if (!method %in% c("MCMC", "CAVI", "SVI")) stop("unknown inference method")
  if (!representation %in% c("posterior draws", "variational parameters"))
    stop("unknown posterior representation")
  if (!is.list(parameters) || !length(parameters))
    stop("native parameters must be a nonempty list")
  if (!is.list(diagnostics)) stop("diagnostics must be a list")
  if (!is.numeric(total_seconds) || length(total_seconds) != 1L ||
      !is.finite(total_seconds) || total_seconds < 0)
    stop("total_seconds must be a nonnegative finite scalar")
  native_parameters <- parameters
  parameters <- canonical_fit_parameters(model, method, native_parameters)
  timing <- list(total_seconds = total_seconds,
                 cumulative_seconds = cumulative_seconds)
  if (!is.null(evaluation_iterations))
    timing$evaluation_iterations <- evaluation_iterations
  if (!is.null(update_seconds)) timing$update_seconds <- update_seconds
  if (!is.null(objective_seconds)) timing$objective_seconds <- objective_seconds
  out <- list(
    schema_version = "fit-result-1",
    model = model,
    method = method,
    status = "fitted",
    representation = representation,
    parameters = parameters,
    objective = canonical_fit_objective(method, native_parameters),
    diagnostics = diagnostics,
    timing = timing,
    metadata = metadata
  )
  latent_variables <- canonical_fit_latent_variables(
    model, method, native_parameters)
  if (!is.null(latent_variables)) out$latent_variables <- latent_variables
  structure(out, class = c("vb_fit_result", "list"))
}

new_convergence_result <- function(method, status, reason, work_completed,
                                   checks_completed, metrics = list(),
                                   thresholds = list(), history = data.frame(),
                                   work_unit) {
  if (!method %in% c("MCMC", "CAVI", "SVI"))
    stop("unknown inference method")
  if (!status %in% c("converged", "max_work", "failed", "interrupted"))
    stop("unknown convergence status")
  if (!is.character(reason) || length(reason) != 1L || !nzchar(reason))
    stop("reason must be a nonempty string")
  work_completed <- validate_integer_scalar(
    work_completed, "work_completed", 0L)
  checks_completed <- validate_integer_scalar(
    checks_completed, "checks_completed", 0L)
  if (!is.list(metrics) || !is.list(thresholds) || !is.data.frame(history))
    stop("metrics and thresholds must be lists and history a data frame")
  if (!is.character(work_unit) || length(work_unit) != 1L || !nzchar(work_unit))
    stop("work_unit must be a nonempty string")
  structure(list(
    schema_version = "convergence-1", method = method, status = status,
    reason = reason, work_completed = work_completed,
    work_unit = work_unit, checks_completed = checks_completed,
    metrics = metrics, thresholds = thresholds, history = history
  ), class = c("vb_convergence", "list"))
}

