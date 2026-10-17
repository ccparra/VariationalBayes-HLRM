# CHLRM-only public accessors. LRM and HLRM fits fail explicitly instead of
# carrying latent-variable slots filled with NULL.

.require_chlrm_fit <- function(fit) {
  if (!inherits(fit, "vb_model_fit"))
    stop("fit must be produced by fit_model()")
  if (!identical(fit$model, "CHLRM"))
    stop("This operation applies only to CHLRM fits")
  invisible(fit)
}

.chlrm_posterior_summary <- function(fit, include_incidence = TRUE, ...) {
  .require_chlrm_fit(fit)
  summarize_posterior(
    fit, include_incidence = include_incidence, ...)
}

#' @rdname chlrm-accessors
#' @export
components <- function(fit, ...) {
  UseMethod("components")
}

#' @rdname chlrm-accessors
#' @export
components.vb_model_fit <- function(fit, ...) {
  summary <- .chlrm_posterior_summary(fit, include_incidence = FALSE, ...)
  summary$parameters[c(
    "component_coefficients", "component_variances", "mixture_weights")]
}

#' @rdname chlrm-accessors
#' @export
responsibilities <- function(fit, ...) {
  UseMethod("responsibilities")
}

#' @rdname chlrm-accessors
#' @export
responsibilities.vb_model_fit <- function(fit, ...) {
  summary <- .chlrm_posterior_summary(fit, include_incidence = FALSE, ...)
  probabilities <- summary$allocations$probabilities
  rownames(probabilities) <- paste0("group_", seq_len(nrow(probabilities)))
  colnames(probabilities) <- paste0(
    "component_", seq_len(ncol(probabilities)))
  probabilities
}

#' @rdname chlrm-accessors
#' @export
classify_groups <- function(fit, ...) {
  UseMethod("classify_groups")
}

#' @rdname chlrm-accessors
#' @export
classify_groups.vb_model_fit <- function(fit, ...) {
  probabilities <- responsibilities(fit, ...)
  labels <- max.col(probabilities, ties.method = "first")
  selected <- probabilities[cbind(seq_len(nrow(probabilities)), labels)]
  entropy <- -rowSums(ifelse(
    probabilities > 0, probabilities * log(probabilities), 0))
  data.frame(
    group = seq_len(nrow(probabilities)),
    component = labels,
    assignment_probability = selected,
    uncertainty = 1 - selected,
    entropy = entropy,
    stringsAsFactors = FALSE)
}

#' @rdname chlrm-accessors
#' @export
classification_uncertainty <- function(fit, ...) {
  classify_groups(fit, ...)[c(
    "group", "component", "assignment_probability", "uncertainty", "entropy")]
}

#' @rdname chlrm-accessors
#' @export
co_clustering <- function(fit, ...) {
  UseMethod("co_clustering")
}

#' @rdname chlrm-accessors
#' @export
co_clustering.vb_model_fit <- function(fit, ...) {
  summary <- .chlrm_posterior_summary(fit, include_incidence = TRUE, ...)
  matrix <- summary$allocations$co_clustering
  if (is.null(matrix))
    stop("The selected CHLRM summary did not compute co-clustering")
  labels <- paste0("group_", seq_len(nrow(matrix)))
  dimnames(matrix) <- list(labels, labels)
  matrix
}
