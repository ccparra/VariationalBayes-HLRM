# Reusable log-likelihood contract for predictive evaluation. Gibbs kernels
# keep their scalar conditional likelihoods; this module materializes the
# pointwise matrix needed by WAIC, CV and predictive scoring.

new_log_likelihood_result <- function(values, model, unit, conditioning,
                                      predictive_units,
                                      integrated_over = character(),
                                      metadata = list()) {
  model <- match.arg(model, c("LRM", "HLRM", "CHLRM"))
  unit <- match.arg(unit, c("observation", "group"))
  conditioning <- match.arg(conditioning,
    c("parameter_conditional", "group_marginal"))
  values <- as.matrix(values)
  if (!is.numeric(values) || !nrow(values) || !ncol(values) ||
      any(!is.finite(values)))
    stop("values must be a nonempty finite [predictive unit, draw] matrix")
  if (!is.data.frame(predictive_units) ||
      nrow(predictive_units) != nrow(values) || !ncol(predictive_units))
    stop("predictive_units must identify every row of values")
  if (!is.character(integrated_over) || anyNA(integrated_over) ||
      any(!nzchar(integrated_over)))
    stop("integrated_over must contain nonempty names")
  if (!is.list(metadata)) stop("metadata must be a list")
  structure(list(
    schema_version = "log-likelihood-1",
    model = model,
    unit = unit,
    conditioning = conditioning,
    values = values,
    predictive_units = predictive_units,
    integrated_over = unique(integrated_over),
    metadata = metadata
  ), class = c("vb_log_likelihood", "list"))
}

#' @rdname interoperability
#' @export
as.matrix.vb_log_likelihood <- function(x, ...) x$values

log_likelihood_estimand <- function(x) {
  if (!inherits(x, "vb_log_likelihood"))
    stop("x must be a log-likelihood-1 object")
  x[c("schema_version", "model", "unit", "conditioning",
      "integrated_over", "metadata")]
}

.observation_ids <- function(group_sizes) {
  data.frame(
    group = rep.int(seq_along(group_sizes), group_sizes),
    observation = sequence(group_sizes),
    stringsAsFactors = FALSE)
}

.group_ids <- function(labels) {
  data.frame(group = labels, stringsAsFactors = FALSE)
}

.validate_lrm_draws <- function(draws, coefficient_count) {
  if (!is.list(draws) || !is.matrix(draws$Beta) ||
      !is.numeric(draws$Beta) || nrow(draws$Beta) < 1L ||
      ncol(draws$Beta) != coefficient_count ||
      any(!is.finite(draws$Beta)))
    stop("LRM draws$Beta must be finite [S,p]")
  S <- nrow(draws$Beta)
  if (!is.numeric(draws$sigma2) || length(draws$sigma2) != S ||
      any(!is.finite(draws$sigma2)) || any(draws$sigma2 <= 0))
    stop("LRM draws$sigma2 must be a positive finite S-vector")
  S
}

.log_likelihood_lrm <- function(y, X, draws, unit, groups,
                                max_output_mb) {
  X <- as.matrix(X); y <- as.numeric(y)
  if (!is.numeric(X) || !nrow(X) || nrow(X) != length(y) ||
      any(!is.finite(X)) || any(!is.finite(y)))
    stop("LRM X and y must be finite with compatible dimensions")
  S <- .validate_lrm_draws(draws, ncol(X))
  if (unit == "observation") {
    validate_evaluation_size(max(2L, S), max_output_mb,
                             as.double(length(y)) * S)
    values <- vapply(seq_len(S), function(s)
      stats::dnorm(y, drop(X %*% draws$Beta[s, ]),
                   sqrt(draws$sigma2[s]), log = TRUE), numeric(length(y)))
    values <- matrix(values, nrow = length(y), ncol = S)
    return(new_log_likelihood_result(
      values, "LRM", unit, "parameter_conditional",
      .observation_ids(length(y)),
      metadata = list(conditioned_on = c("coefficients", "residual_variance"))))
  }
  if (is.null(groups) || !is.atomic(groups) || length(groups) != length(y) ||
      anyNA(groups) || (is.numeric(groups) && any(!is.finite(groups))))
    stop("groups must give one nonmissing label per observation for group units")
  labels <- unique(groups)
  group_index <- match(groups, labels)
  validate_evaluation_size(max(2L, S), max_output_mb,
                           as.double(length(y) + length(labels)) * S)
  values <- vapply(seq_len(S), function(s) {
    ll <- stats::dnorm(y, drop(X %*% draws$Beta[s, ]),
                       sqrt(draws$sigma2[s]), log = TRUE)
    drop(rowsum(ll, group_index, reorder = FALSE))
  }, numeric(length(labels)))
  values <- matrix(values, nrow = length(labels), ncol = S)
  new_log_likelihood_result(
    values, "LRM", unit, "parameter_conditional", .group_ids(labels),
    metadata = list(conditioned_on = c("coefficients", "residual_variance"),
                    aggregation = "joint likelihood within supplied groups"))
}

.log_likelihood_hlrm_conditional <- function(y, X, draws, unit,
                                              max_output_mb) {
  info <- validate_grouped_data(X, y)
  S <- validate_gaussian_draws(draws, info$m, info$p)
  required <- as.double(if (unit == "observation") sum(info$nj) else info$m) * S
  validate_evaluation_size(S, max_output_mb, required)
  if (unit == "observation") {
    values <- matrix(NA_real_, sum(info$nj), S)
    start <- cumsum(c(1L, head(info$nj, -1L)))
    for (j in seq_len(info$m)) for (s in seq_len(S)) {
      rows <- start[j]:(start[j] + info$nj[j] - 1L)
      values[rows, s] <- stats::dnorm(
        y[[j]], drop(X[[j]] %*% draws$Betaj[s, j, ]),
        sqrt(draws$sigma2j[s, j]), log = TRUE)
    }
    units <- .observation_ids(info$nj)
  } else {
    values <- matrix(NA_real_, info$m, S)
    for (j in seq_len(info$m)) for (s in seq_len(S))
      values[j, s] <- sum(stats::dnorm(
        y[[j]], drop(X[[j]] %*% draws$Betaj[s, j, ]),
        sqrt(draws$sigma2j[s, j]), log = TRUE))
    units <- .group_ids(seq_len(info$m))
  }
  new_log_likelihood_result(
    values, "HLRM", unit, "parameter_conditional", units,
    metadata = list(conditioned_on = c("group_coefficients", "group_variances")))
}

#' Pointwise log likelihood under an explicit predictive estimand
#'
#' @return A `log-likelihood-1` object retaining pointwise values and estimand.
#' @rdname predictive-evaluation
#' @export
log_likelihood <- function(model, y, X, draws, unit, conditioning,
                           groups = NULL, nu0 = NULL,
                           marginal_control = list(), max_output_mb = 256) {
  model <- match.arg(model, c("LRM", "HLRM", "CHLRM"))
  unit <- match.arg(unit, c("observation", "group"))
  conditioning <- match.arg(conditioning,
    c("parameter_conditional", "group_marginal"))
  if (!is.list(marginal_control) || (length(marginal_control) &&
      (is.null(names(marginal_control)) || anyDuplicated(names(marginal_control)) ||
       any(!names(marginal_control) %in%
             c("tolerance", "max_nodes", "max_evaluations")))))
    stop("Unknown or duplicated marginal_control fields")

  if (model == "LRM") {
    if (conditioning != "parameter_conditional")
      stop("LRM has no group latent quantity to marginalize")
    return(.log_likelihood_lrm(y, X, draws, unit, groups, max_output_mb))
  }
  if (conditioning == "parameter_conditional") {
    if (model != "HLRM")
      stop("CHLRM predictive scoring uses group_marginal allocation likelihood")
    return(.log_likelihood_hlrm_conditional(y, X, draws, unit,
                                             max_output_mb))
  }
  if (unit != "group")
    stop("group_marginal likelihood requires unit = 'group'")
  info <- validate_grouped_data(X, y)
  if (model == "CHLRM") {
    values <- loglik_group_marginal_chlrm(y, X, draws, max_output_mb)
    return(new_log_likelihood_result(
      values, model, unit, conditioning, .group_ids(seq_len(info$m)),
      integrated_over = "component_allocation",
      metadata = list(marginalization = "finite mixture sum")))
  }
  if (is.null(nu0)) stop("nu0 is required for HLRM group marginalization")
  values <- do.call(loglik_group_marginal_hlrm,
    c(list(y_list = y, X_list = X, draws = draws, nu0 = nu0,
           max_output_mb = max_output_mb), marginal_control))
  new_log_likelihood_result(
    values, model, unit, conditioning, .group_ids(seq_len(info$m)),
    integrated_over = c("group_coefficients", "group_variances"),
    metadata = list(marginalization = "Gaussian analytic plus inverse-gamma quadrature",
                    quadrature = attr(values, "quadrature")))
}
