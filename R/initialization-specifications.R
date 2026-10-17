.initialization_spec <- function(strategy, values = NULL, jitter = NULL) {
  fields <- list(schema_version = "initialization-spec-1",
                 strategy = strategy)
  if (!is.null(values)) fields$values <- values
  if (!is.null(jitter)) fields$jitter <- jitter
  .new_vb_spec(fields, "vb_initialization_spec")
}

#' @rdname model-specifications
#' @export
initialization_from_prior <- function() .initialization_spec("from_prior")

#' @rdname model-specifications
#' @export
initialization_ols <- function(jitter = 1e-6) {
  jitter <- .positive_scalar(jitter, "jitter")
  .initialization_spec("ols", jitter = jitter)
}

#' @rdname model-specifications
#' @export
initialization_kmeans <- function(jitter = 1e-6) {
  jitter <- .positive_scalar(jitter, "jitter")
  .initialization_spec("kmeans", jitter = jitter)
}

#' @rdname model-specifications
#' @export
initialization_values <- function(values) {
  if (!is.list(values) || !length(values) || is.null(names(values)) ||
      any(!nzchar(names(values))) || anyDuplicated(names(values)) ||
      any(vapply(values, is.null, logical(1L))))
    stop("values must be a nonempty named list without NULL entries")
  .initialization_spec("values", values = values)
}

