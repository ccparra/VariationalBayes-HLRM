# Dependency-free long representation of posterior summaries for tables,
# plotting, and simulation studies.

.parameter_indices <- function(value) {
  dimensions <- dim(value)
  if (is.null(dimensions))
    return(cbind(index_1 = seq_along(value), index_2 = NA_integer_))
  if (length(dimensions) > 2L)
    stop("tidy_parameters() supports summary parameters with at most two dimensions")
  indices <- arrayInd(seq_along(value), dimensions)
  if (ncol(indices) == 1L)
    indices <- cbind(indices, NA_integer_)
  colnames(indices) <- c("index_1", "index_2")
  indices
}

.tidy_parameter_block <- function(name, parameter, model, method) {
  required <- c("mean", "sd", "lower", "upper", "mean_exists",
                "contains_zero", "level", "distribution")
  if (!is.list(parameter) || !all(required %in% names(parameter)))
    stop("Every parameter must be a canonical posterior summary")
  indices <- .parameter_indices(parameter$mean)
  length_required <- length(parameter$mean)
  fields <- lapply(parameter[c(
    "mean", "sd", "lower", "upper", "mean_exists", "contains_zero")],
    as.vector)
  if (any(lengths(fields) != length_required))
    stop("Posterior summary fields have incompatible dimensions")
  data.frame(
    model = model,
    method = method,
    parameter = name,
    index_1 = indices[, "index_1"],
    index_2 = indices[, "index_2"],
    estimate = as.numeric(fields$mean),
    sd = as.numeric(fields$sd),
    lower = as.numeric(fields$lower),
    upper = as.numeric(fields$upper),
    mean_exists = as.logical(fields$mean_exists),
    contains_zero = as.logical(fields$contains_zero),
    distribution = parameter$distribution,
    interval_level = parameter$level,
    stringsAsFactors = FALSE)
}

#' @rdname tidy-parameters
#' @export
tidy_parameters <- function(x, ...) {
  UseMethod("tidy_parameters")
}

#' @rdname tidy-parameters
#' @export
tidy_parameters.vb_model_fit <- function(x, ...) {
  tidy_parameters(summarize_posterior(x, ...))
}

#' @rdname tidy-parameters
#' @export
tidy_parameters.vb_posterior_summary <- function(x, ...) {
  if (!identical(x$schema_version, "posterior-summary-1") ||
      !is.list(x$parameters) || !length(x$parameters))
    stop("x must be a posterior-summary-1 object")
  rows <- Map(
    .tidy_parameter_block, names(x$parameters), x$parameters,
    MoreArgs = list(model = x$model, method = x$method))
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  structure(out, class = c("vb_tidy_parameters", "data.frame"))
}
