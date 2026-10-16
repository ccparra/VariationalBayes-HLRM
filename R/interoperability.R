# Explicit orientation adapters for the wider Bayesian R ecosystem. Internal
# contracts keep predictive units in rows and draws in columns; bayesplot and
# loo conventionally consume draws in rows. The transpose is therefore named,
# validated, and tested rather than left to callers.

as_yrep <- function(x, ...) {
  UseMethod("as_yrep")
}

as_yrep.vb_prediction_result <- function(x, ...) {
  if (!identical(x$target, "response"))
    stop("as_yrep() requires response predictions, not latent-mean predictions")
  if (is.null(x$draws) || !is.list(x$draws) || !length(x$draws))
    stop("Prediction draws were not retained; call posterior_predict(..., return_draws = TRUE)")
  draw_counts <- vapply(x$draws, ncol, integer(1L))
  if (length(unique(draw_counts)) != 1L)
    stop("All prediction blocks must contain the same number of draws")
  unit_draw <- do.call(rbind, x$draws)
  if (!is.numeric(unit_draw) || any(!is.finite(unit_draw)))
    stop("Prediction draws must be finite numeric matrices")
  out <- t(unit_draw)
  if (nrow(x$summary) == ncol(out)) {
    colnames(out) <- paste0(
      "group_", x$summary$group, ".observation_", x$summary$row)
  }
  attr(out, "model") <- x$model
  attr(out, "method") <- x$method
  attr(out, "target") <- x$target
  out
}

as_yrep.vb_posterior_predictive_check <- function(x, ...) {
  if (is.null(x$y_rep))
    stop("Replicated responses were not retained; set return_yrep = TRUE")
  values <- as.matrix(x$y_rep)
  if (!is.numeric(values) || !nrow(values) || !ncol(values) ||
      any(!is.finite(values)))
    stop("Stored replicated responses must be a finite matrix")
  out <- t(values)
  colnames(out) <- paste0("observation_", seq_len(ncol(out)))
  attr(out, "model") <- x$model
  attr(out, "method") <- x$method
  attr(out, "target") <- x$predictive_target
  out
}

as_log_lik_matrix <- function(x, ...) {
  UseMethod("as_log_lik_matrix")
}

as_log_lik_matrix.vb_log_likelihood <- function(x, ...) {
  values <- as.matrix(x$values)
  out <- t(values)
  identifiers <- apply(x$predictive_units, 1L, function(row) paste(
    paste0(names(x$predictive_units), "_", row), collapse = "."))
  colnames(out) <- make.unique(identifiers)
  attr(out, "predictive_units") <- x$predictive_units
  attr(out, "estimand") <- log_likelihood_estimand(x)
  out
}
