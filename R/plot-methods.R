# Dependency-light diagnostic plots. Each method returns the plotted data
# invisibly so figures can be reproduced or rebuilt with another graphics
# system without extracting private fields.

.plot_call <- function(fun, defaults, dots) {
  arguments <- utils::modifyList(defaults, dots)
  do.call(fun, arguments)
}

.parameter_plot_labels <- function(data) {
  index <- ifelse(
    is.na(data$index_2), as.character(data$index_1),
    paste(data$index_1, data$index_2, sep = ","))
  paste0(data$parameter, "[", index, "]")
}

#' Plot package results
#'
#' Base-R diagnostic views of fitted models, simulations, predictions,
#' posterior predictive checks, pointwise log likelihoods, and orchestration
#' audits. Plotting never refits a model or changes a result object. Every
#' method invisibly returns the data used in the figure.
#'
#' @param x A supported `VariationalBayesHLRM` result object.
#' @param type Plot type. Fits accept `"convergence"`, `"coefficients"`, or
#'   `"responsibilities"`; orchestration audits currently accept
#'   `"elapsed"`.
#' @param parameter Optional parameter-block names for a coefficient plot.
#'   For posterior draws, these names select the distributions to display.
#' @param metric Optional convergence-history column to display.
#' @param max_panels Maximum number of scalar posterior densities in one plot.
#' @param observed Optional observed response vector to overlay on a prediction
#'   plot.
#' @param predictor Predictor column name or index for a simulation plot. By
#'   default the first nonconstant column is used.
#' @param ... Additional graphical arguments passed to the underlying base-R
#'   plotting function.
#' @return The plotted data, invisibly.
#' @examples
#' simulation <- simulate_data(
#'   conditional_dgp(lrm_model(2), list(
#'     coefficients = c(1, 2), residual_variance = 1)),
#'   simulation_design("grid", bounds = matrix(c(-1, 1), nrow = 1)),
#'   simulation_sample(sample_size = 20), seed = 31
#' )
#' plot(simulation)
#' @name package-plots
NULL

#' @rdname package-plots
#' @export
plot.vb_model_fit <- function(
    x, type = c("convergence", "coefficients", "responsibilities"),
    parameter = NULL, metric = NULL, ...) {
  type <- match.arg(type)
  dots <- list(...)

  if (type == "convergence") {
    history <- x$diagnostics$history
    if (!is.data.frame(history) || !nrow(history))
      stop("This fit has no convergence history to plot")
    step_name <- if ("iteration" %in% names(history)) "iteration" else
      if ("stored_samples" %in% names(history)) "stored_samples" else
        names(history)[1L]
    candidates <- if (x$method == "MCMC") c(
      "max_rhat", "min_ess_bulk", "min_ess_tail",
      "max_relative_mcse_mean", "max_relative_mcse_sd") else
      c("smoothed_objective", "objective", "objective_relative_change",
        "parameter_relative_change")
    candidates <- candidates[candidates %in% names(history)]
    if (is.null(metric)) {
      metric <- candidates[vapply(candidates, function(name)
        any(is.finite(history[[name]])), logical(1L))][1L]
    }
    if (length(metric) != 1L || is.na(metric) ||
        !metric %in% names(history) || !is.numeric(history[[metric]]))
      stop("metric must identify a numeric convergence-history column")
    data <- data.frame(
      step = history[[step_name]], value = history[[metric]],
      metric = metric, stringsAsFactors = FALSE)
    .plot_call(graphics::plot, list(
      x = data$step, y = data$value, type = "l",
      xlab = step_name, ylab = metric,
      main = paste(x$model, x$method, "convergence")), dots)
    return(invisible(data))
  }

  if (type == "coefficients") {
    data <- as.data.frame(tidy_parameters(x))
    if (is.null(parameter)) {
      selected <- grepl("coefficients$", data$parameter)
      if (any(selected)) data <- data[selected, , drop = FALSE]
    } else {
      if (!is.character(parameter) || !length(parameter) || anyNA(parameter))
        stop("parameter must contain one or more parameter-block names")
      unknown <- setdiff(parameter, unique(data$parameter))
      if (length(unknown))
        stop("Unknown parameter block: ", paste(unknown, collapse = ", "))
      data <- data[data$parameter %in% parameter, , drop = FALSE]
    }
    keep <- is.finite(data$estimate) & is.finite(data$lower) &
      is.finite(data$upper)
    data <- data[keep, , drop = FALSE]
    if (!nrow(data)) stop("No finite parameter intervals are available to plot")
    data$label <- .parameter_plot_labels(data)
    position <- rev(seq_len(nrow(data)))
    .plot_call(graphics::plot, list(
      x = data$estimate, y = position,
      xlim = range(c(data$lower, data$upper)), ylim = c(.5, nrow(data) + .5),
      yaxt = "n", ylab = "", xlab = "posterior estimate",
      main = paste(x$model, x$method, "coefficient intervals"), pch = 19), dots)
    graphics::segments(data$lower, position, data$upper, position)
    graphics::axis(2, at = position, labels = data$label, las = 2)
    graphics::abline(v = 0, lty = 3, col = "grey60")
    return(invisible(data))
  }

  .require_chlrm_fit(x)
  data <- responsibilities(x)
  palette <- grDevices::hcl.colors(32, "Blues 3")
  .plot_call(graphics::image, list(
    x = seq_len(ncol(data)), y = seq_len(nrow(data)), z = t(data),
    col = palette, xlab = "component", ylab = "group",
    main = paste(x$method, "allocation probabilities"), axes = FALSE), dots)
  graphics::axis(1, at = seq_len(ncol(data)), labels = colnames(data))
  graphics::axis(2, at = seq_len(nrow(data)), labels = rownames(data), las = 2)
  invisible(data)
}

#' @rdname package-plots
#' @export
plot.vb_prediction_result <- function(x, observed = NULL, ...) {
  data <- x$summary
  required <- c("mean", "lower", "upper")
  if (!is.data.frame(data) || any(!required %in% names(data)))
    stop("x does not contain a valid posterior prediction summary")
  if (!is.null(observed) &&
      (!is.numeric(observed) || length(observed) != nrow(data) ||
       any(!is.finite(observed))))
    stop("observed must be a finite vector matching the prediction rows")
  position <- seq_len(nrow(data))
  dots <- list(...)
  .plot_call(graphics::plot, list(
    x = position, y = data$mean, ylim = range(c(data$lower, data$upper,
      if (is.null(observed)) numeric() else observed)),
    xlab = "prediction row", ylab = "response", pch = 19,
    main = paste(x$model, x$method, x$target, "prediction")), dots)
  graphics::segments(position, data$lower, position, data$upper)
  if (!is.null(observed)) {
    graphics::points(position, observed, pch = 4, col = "firebrick")
  }
  invisible(data)
}

#' @rdname package-plots
#' @export
plot.vb_posterior_draws <- function(
    x, parameter = NULL, max_panels = 9L, ...) {
  if (is.null(parameter)) parameter <- names(x$parameters)[1L]
  if (!is.character(parameter) || !length(parameter) || anyNA(parameter))
    stop("parameter must contain one or more posterior parameter names")
  unknown <- setdiff(parameter, names(x$parameters))
  if (length(unknown))
    stop("Unknown posterior parameter: ", paste(unknown, collapse = ", "))
  max_panels <- validate_integer_scalar(max_panels, "max_panels", 1L)
  draw_count <- x$metadata$draw_count
  matrices <- Map(
    .flatten_draw_parameter, x$parameters[parameter], parameter,
    MoreArgs = list(draw_count = draw_count))
  values <- do.call(cbind, matrices)
  if (ncol(values) > max_panels)
    stop("Selected parameters exceed max_panels; request fewer blocks or increase the explicit limit")

  previous <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(previous), add = TRUE)
  graphics::par(
    mfrow = grDevices::n2mfrow(ncol(values)),
    mar = c(3, 3, 2, 1), mgp = c(1.8, .6, 0))
  dots <- list(...)
  for (column in seq_len(ncol(values))) {
    estimate <- stats::density(values[, column])
    .plot_call(graphics::plot, list(
      x = estimate, xlab = "value", ylab = "density",
      main = colnames(values)[column]), dots)
  }
  data <- data.frame(
    draw = rep.int(seq_len(draw_count), ncol(values)),
    parameter = rep(colnames(values), each = draw_count),
    value = as.vector(values), stringsAsFactors = FALSE)
  invisible(data)
}

#' @rdname package-plots
#' @export
plot.vb_posterior_predictive_check <- function(x, ...) {
  replicated <- as.data.frame(x$discrepancies$repdata)
  observed <- x$discrepancies$t0
  if (!ncol(replicated) || length(observed) != ncol(replicated))
    stop("x does not contain valid posterior predictive discrepancies")
  dots <- list(...)
  .plot_call(graphics::boxplot, list(
    x = replicated, xlab = "discrepancy", ylab = "replicated value",
    main = paste(x$model, x$method, "posterior predictive check")), dots)
  graphics::points(seq_along(observed), observed, pch = 4,
                   col = "firebrick", lwd = 2)
  invisible(list(
    replicated = replicated, observed = observed,
    posterior_predictive_p_value = x$discrepancies$ppp))
}

#' @rdname package-plots
#' @export
plot.vb_log_likelihood <- function(x, ...) {
  values <- as.matrix(x)
  data <- cbind(
    x$predictive_units,
    mean = rowMeans(values),
    lower = apply(values, 1L, stats::quantile, probs = .025),
    upper = apply(values, 1L, stats::quantile, probs = .975))
  position <- seq_len(nrow(data))
  dots <- list(...)
  .plot_call(graphics::plot, list(
    x = position, y = data$mean,
    ylim = range(c(data$lower, data$upper)),
    xlab = x$unit, ylab = "log likelihood", pch = 19,
    main = paste(x$model, x$conditioning, "log likelihood")), dots)
  graphics::segments(position, data$lower, position, data$upper)
  invisible(data)
}

#' @rdname package-plots
#' @export
plot.vb_fit_orchestration <- function(x, type = "elapsed", ...) {
  type <- match.arg(type, "elapsed")
  data <- if (identical(x$type, "multistart")) x$summary else
    x$chain_metadata
  label <- if (identical(x$type, "multistart")) data$start else data$chain
  color <- if ("selected" %in% names(data))
    ifelse(data$selected, "steelblue", "grey70") else "steelblue"
  dots <- list(...)
  .plot_call(graphics::barplot, list(
    height = data$elapsed_seconds, names.arg = label, col = color,
    xlab = if (identical(x$type, "multistart")) "start" else "chain",
    ylab = "elapsed seconds",
    main = paste(x$model, x$method, x$type, "timing")), dots)
  invisible(data)
}

.simulation_plot_data <- function(x, predictor) {
  grouped <- is.list(x$train$X) && !is.data.frame(x$train$X)
  X <- if (grouped) do.call(rbind, x$train$X) else as.matrix(x$train$X)
  response <- if (grouped) unlist(x$train$y, use.names = FALSE) else
    as.numeric(x$train$y)
  if (is.null(predictor)) {
    variable <- which(apply(X, 2L, function(value)
      stats::sd(value) > sqrt(.Machine$double.eps)))[1L]
    if (is.na(variable)) variable <- 1L
  } else if (is.character(predictor)) {
    if (length(predictor) != 1L || is.null(colnames(X)) ||
        is.na(variable <- match(predictor, colnames(X))))
      stop("predictor must name one design-matrix column")
  } else {
    variable <- validate_integer_scalar(predictor, "predictor", 1L)
    if (variable > ncol(X)) stop("predictor exceeds the design-matrix columns")
  }
  group <- if (grouped)
    rep.int(seq_along(x$train$X), lengths(x$train$y)) else
      rep.int(1L, nrow(X))
  component <- if (x$model == "CHLRM" &&
      !is.null(x$truth$allocations)) x$truth$allocations[group] else
        rep.int(NA_integer_, nrow(X))
  data.frame(
    predictor = X[, variable], response = response,
    group = group, component = component,
    stringsAsFactors = FALSE)
}

#' @rdname package-plots
#' @export
plot.vb_simulation_result <- function(x, predictor = NULL, ...) {
  data <- .simulation_plot_data(x, predictor)
  color_index <- if (x$model == "CHLRM" && any(!is.na(data$component)))
    data$component else if (x$model == "HLRM") data$group else rep.int(1L, nrow(data))
  colors <- grDevices::hcl.colors(max(color_index, na.rm = TRUE), "Dark 3")
  dots <- list(...)
  .plot_call(graphics::plot, list(
    x = data$predictor, y = data$response, col = colors[color_index],
    pch = 19, xlab = "predictor", ylab = "response",
    main = paste(x$model, "simulated data")), dots)
  invisible(data)
}
