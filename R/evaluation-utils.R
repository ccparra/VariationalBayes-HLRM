# Evaluation primitives shared by all model families. No model fitting or I/O.
validate_evaluation_size <- function(S, max_output_mb, required_doubles = 0) {
  if (!is.numeric(S) || length(S) != 1L || !is.finite(S) || S < 2 ||
      S != floor(S) || S > .Machine$integer.max) stop("S must be an integer >= 2")
  if (!is.numeric(max_output_mb) || length(max_output_mb) != 1L ||
      !is.finite(max_output_mb) || max_output_mb <= 0)
    stop("max_output_mb must be a positive finite scalar")
  if (8 * required_doubles / 1024^2 > max_output_mb)
    stop("Evaluation arrays exceed max_output_mb; reduce S or increase the explicit budget")
  invisible(TRUE)
}

evaluation_rng_streams <- function(seed = NULL) {
  validate_rng_seed(seed)
  if (is.null(seed)) return(list(draws = NULL, replicates = NULL))
  seeds <- derive_independent_seeds(2L, seed = seed)
  list(draws = seeds[1L], replicates = seeds[2L])
}

# A PPC discrepancy specification is deliberately independent of a fitted
# model.  The standard preset reproduces the five checks used in the paper:
# T1 global mean, T2 global SD, T3 variance of group means, T4 mean within-
# group variance and T5 MSE against a fixed observed-data OLS reference.
#' @rdname predictive-evaluation
#' @export
ppc_discrepancies <- function(statistics = "standard") {
  if (!is.character(statistics) || !length(statistics) || anyNA(statistics))
    stop("statistics must be a nonempty character vector")
  aliases <- c(
    T1 = "mean", T2 = "sd",
    T3 = "between_group_mean_variance",
    T4 = "mean_within_group_variance", T5 = "ols_mse")
  statistics <- unname(ifelse(statistics %in% names(aliases),
                              aliases[statistics], statistics))
  allowed <- c("standard", "distribution", "mean", "sd",
               "between_group_mean_variance",
               "mean_within_group_variance", "ols_mse")
  unknown <- setdiff(statistics, allowed)
  if (length(unknown))
    stop("Unknown PPC discrepancy: ", paste(unknown, collapse = ", "))
  if ("standard" %in% statistics && length(statistics) > 1L)
    stop("'standard' cannot be combined with named discrepancies")
  if ("distribution" %in% statistics && length(statistics) > 1L)
    stop("'distribution' cannot be combined with named discrepancies")
  structure(list(statistics = unique(statistics)),
            class = c("vb_ppc_discrepancies", "list"))
}

.resolve_ppc_statistics <- function(specification, y_list, X_list) {
  if (is.null(specification)) specification <- ppc_discrepancies()
  if (!inherits(specification, "vb_ppc_discrepancies"))
    stop("discrepancies must be produced by ppc_discrepancies()")
  statistics <- specification$statistics
  if (identical(statistics, "distribution")) statistics <- c("mean", "sd")
  if (identical(statistics, "standard")) {
    statistics <- c("mean", "sd")
    if (length(y_list) > 1L) {
      statistics <- c(statistics, "between_group_mean_variance")
      if (all(lengths(y_list) >= 2L))
        statistics <- c(statistics, "mean_within_group_variance")
    }
    if (!is.null(X_list)) statistics <- c(statistics, "ols_mse")
  }
  grouped <- c("between_group_mean_variance",
               "mean_within_group_variance")
  if (length(y_list) < 2L && any(statistics %in% grouped))
    stop("between/within-group discrepancies require at least two groups")
  if ("mean_within_group_variance" %in% statistics &&
      any(lengths(y_list) < 2L))
    stop("mean_within_group_variance requires at least two observations per group")
  if ("ols_mse" %in% statistics && is.null(X_list))
    stop("ols_mse requires X_list")
  statistics
}

.ols_fitted_reference <- function(y_list, X_list) {
  if (!is.list(X_list) || length(X_list) != length(y_list))
    stop("X_list must match y_list when ols_mse is requested")
  Map(function(y, X) {
    X <- as.matrix(X)
    if (!is.numeric(X) || nrow(X) != length(y) || !ncol(X) ||
        any(!is.finite(X)))
      stop("X_list must contain finite design matrices matching y_list")
    qr.fitted(qr(X, tol = 1e-10), y)
  }, y_list, X_list)
}

.ppc_statistic_function <- function(statistics, y_list, X_list) {
  ols_fitted <- if ("ols_mse" %in% statistics)
    unlist(.ols_fitted_reference(y_list, X_list), use.names = FALSE) else NULL
  function(values) {
    flat <- unlist(values, use.names = FALSE)
    result <- numeric(length(statistics)); names(result) <- statistics
    if ("mean" %in% statistics) result["mean"] <- mean(flat)
    if ("sd" %in% statistics) result["sd"] <- stats::sd(flat)
    if ("between_group_mean_variance" %in% statistics)
      result["between_group_mean_variance"] <-
        stats::var(vapply(values, mean, numeric(1L)))
    if ("mean_within_group_variance" %in% statistics)
      result["mean_within_group_variance"] <-
        mean(vapply(values, stats::var, numeric(1L)))
    if ("ols_mse" %in% statistics)
      result["ols_mse"] <- mean((flat - ols_fitted)^2)
    result
  }
}

# simulate_draw(s) returns a list with one response vector per observed group.
# Retains only discrepancy summaries unless replicates are requested explicitly.
pp_check_stream <- function(y_list, S, simulate_draw, target, seed = NULL,
                            return_yrep = FALSE, max_output_mb = 256,
                            X_list = NULL, discrepancies = NULL) {
  if (!is.list(y_list) || !length(y_list) || any(lengths(y_list) < 1L) ||
      !all(vapply(y_list, function(y) is.numeric(y) && all(is.finite(y)), logical(1))))
    stop("y_list must contain finite nonempty numeric response vectors")
  if (!is.logical(return_yrep) || length(return_yrep) != 1L || is.na(return_yrep))
    stop("return_yrep must be TRUE or FALSE")
  n <- sum(lengths(y_list))
  statistics <- .resolve_ppc_statistics(discrepancies, y_list, X_list)
  summarize <- .ppc_statistic_function(statistics, y_list, X_list)
  validate_evaluation_size(
    S, max_output_mb,
    S * (length(statistics) + if (return_yrep) as.double(n) else 0))
  scope <- inference_rng_scope(seed)
  on.exit(scope$restore(), add = TRUE)
  observed <- summarize(y_list)
  repdata <- matrix(NA_real_, S, length(observed),
                    dimnames = list(NULL, names(observed)))
  Yrep <- if (return_yrep) matrix(NA_real_, n, S) else NULL
  for (s in seq_len(S)) {
    yr <- simulate_draw(s)
    if (!is.list(yr) || !identical(lengths(yr), lengths(y_list)) ||
        !all(vapply(yr, function(y) is.numeric(y) && all(is.finite(y)), logical(1))))
      stop("simulate_draw returned invalid responses at draw ", s)
    repdata[s, ] <- summarize(yr)
    if (return_yrep) Yrep[, s] <- unlist(yr, use.names = FALSE)
  }
  ppp <- colMeans(sweep(repdata, 2, observed, `>=`))
  out <- list(repdata = as.data.frame(repdata), t0 = observed, ppp = ppp,
              mcse_iid = sqrt(ppp * (1 - ppp) / S),
              target = target,
              statistic_definition = c(
                mean = "T1: global response mean",
                sd = "T2: global response standard deviation",
                between_group_mean_variance =
                  "T3: sample variance of observed-group response means",
                mean_within_group_variance =
                  "T4: mean observed-group sample response variance",
                ols_mse =
                  "T5: MSE against the fixed observed-data groupwise OLS fit"
              )[statistics])
  if (return_yrep) out$y_rep <- Yrep
  out
}

# Reuse response draws produced by any model$posterior_predict() call. This is
# useful for plots and custom validation sets; posterior_predictive_check()
# keeps a streaming simulator so large runs do not materialise y_rep.
#' @rdname predictive-evaluation
#' @export
pp_check_from_prediction <- function(y, prediction, X = NULL,
                                     discrepancies = ppc_discrepancies("distribution"),
                                     return_yrep = FALSE,
                                     max_output_mb = 256) {
  if (!is.list(prediction) ||
      !identical(prediction$schema_version, "posterior-prediction-2") ||
      !identical(prediction$target, "response") ||
      is.null(prediction$draws) || !is.list(prediction$draws)) {
    stop("prediction must be a response posterior-prediction-2 object with stored draws")
  }
  y_list <- if (is.list(y)) y else list(as.numeric(y))
  if (length(y_list) != length(prediction$draws) ||
      !identical(lengths(y_list), vapply(prediction$draws, nrow, integer(1))))
    stop("y must match the groups and rows in prediction$draws")
  S <- unique(vapply(prediction$draws, ncol, integer(1)))
  if (length(S) != 1L) stop("All predictive draw matrices must have the same number of columns")
  simulate_draw <- function(s) lapply(prediction$draws, function(x) x[, s])
  X_list <- if (is.null(X)) NULL else if (is.list(X)) X else list(as.matrix(X))
  out <- pp_check_stream(y_list, S, simulate_draw,
                         paste0("external_", prediction$target, "_draws"),
                         return_yrep = return_yrep, max_output_mb = max_output_mb,
                         X_list = X_list, discrepancies = discrepancies)
  out$prediction_metadata <- prediction[c(
    "schema_version", "model", "method", "target", "interval", "metadata")]
  out
}

summarize_waic_pointwise <- function(pointwise) {
  rownames(pointwise) <- NULL
  lppd <- pointwise$lppd; penalty <- pointwise$p_waic; elpd <- lppd - penalty
  if (!length(lppd) || any(!is.finite(c(lppd, penalty, elpd))))
    stop("WAIC pointwise values must be finite and nonempty")
  pointwise$elpd_waic <- elpd
  list(WAIC = -2 * sum(elpd), lppd = sum(lppd), p_waic = sum(penalty),
       se = if (length(elpd) > 1L) 2 * sqrt(length(elpd) * var(elpd)) else NA_real_,
       pointwise = pointwise,
       warning = if (any(penalty > .4))
         "Some predictive units have Var(log likelihood) > 0.4; inspect PSIS-LOO." else NULL)
}
mean_squared_error <- function(observed, predicted) {
  mean((observed - predicted)^2)
}

coefficient_of_determination <- function(observed, predicted) {
  observed <- as.numeric(observed)
  predicted <- as.numeric(predicted)
  residual_sum_squares <- sum((observed - predicted)^2)
  total_sum_squares <- sum((observed - mean(observed))^2)
  if (total_sum_squares == 0) return(NA_real_)
  1 - residual_sum_squares / total_sum_squares
}
