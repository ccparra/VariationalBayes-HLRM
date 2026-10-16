# Metrics that summarize repeated experiments against known truth or an MCMC
# reference. They belong to the paper suite, not to the fitted-model API.

.metric_numeric_vector <- function(value, name, minimum_length = 1L) {
  value <- as.numeric(value)
  if (length(value) < minimum_length || any(!is.finite(value)))
    stop(name, " must be a finite numeric vector")
  value
}

.choose_two <- function(value) value * (value - 1) / 2

.mode_integer <- function(value) {
  counts <- table(as.integer(value))
  as.integer(names(counts)[which.max(counts)])
}

adjusted_rand_index <- function(truth, estimate) {
  truth <- as.character(truth)
  estimate <- as.character(estimate)
  if (!length(truth) || length(truth) != length(estimate) ||
      anyNA(truth) || anyNA(estimate))
    stop("truth and estimate must be nonempty, complete vectors of equal length")
  if (length(truth) == 1L) return(1)
  contingency <- table(truth, estimate)
  pair_agreement <- sum(.choose_two(contingency))
  truth_pairs <- sum(.choose_two(rowSums(contingency)))
  estimate_pairs <- sum(.choose_two(colSums(contingency)))
  total_pairs <- .choose_two(length(truth))
  expected <- truth_pairs * estimate_pairs / total_pairs
  maximum <- (truth_pairs + estimate_pairs) / 2
  denominator <- maximum - expected
  if (abs(denominator) <= .Machine$double.eps)
    return(as.numeric(identical(
      outer(truth, truth, `==`), outer(estimate, estimate, `==`))))
  (pair_agreement - expected) / denominator
}

co_clustering_error <- function(truth, co_clustering) {
  truth <- as.character(truth)
  probability <- as.matrix(co_clustering)
  m <- length(truth)
  if (m < 2L || anyNA(truth))
    stop("truth must contain at least two complete labels")
  if (!identical(dim(probability), as.integer(c(m, m))) ||
      any(!is.finite(probability)) || any(probability < 0) ||
      any(probability > 1))
    stop("co_clustering must be a finite m by m probability matrix")
  if (!isTRUE(all.equal(probability, t(probability), tolerance = 1e-12)))
    stop("co_clustering must be symmetric")
  if (any(abs(diag(probability) - 1) > 1e-12))
    stop("co_clustering must have a unit diagonal")
  actual <- outer(truth, truth, `==`) * 1
  2 * sum((probability[upper.tri(probability)] -
             actual[upper.tri(actual)])^2) / (m * (m - 1))
}

normalized_classification_entropy <- function(probabilities) {
  probabilities <- as.matrix(probabilities)
  if (!nrow(probabilities) || !ncol(probabilities) ||
      any(!is.finite(probabilities)) || any(probabilities < 0) ||
      any(probabilities > 1) ||
      any(abs(rowSums(probabilities) - 1) > 1e-10))
    stop("probabilities must be a finite row-stochastic matrix")
  if (ncol(probabilities) == 1L) return(0)
  terms <- ifelse(probabilities == 0, 0,
                  probabilities * log(probabilities))
  -sum(terms) / (nrow(probabilities) * log(ncol(probabilities)))
}

scalar_recovery_metrics <- function(estimates, truth, lower = NULL,
                                    upper = NULL) {
  estimates <- .metric_numeric_vector(estimates, "estimates")
  if (!is.numeric(truth) || length(truth) != 1L || !is.finite(truth))
    stop("truth must be one finite scalar")
  error <- estimates - truth
  out <- list(
    replicates = length(estimates),
    bias = mean(error),
    rmse = sqrt(mean(error^2)),
    bias_mcse = if (length(error) > 1L)
      stats::sd(error) / sqrt(length(error)) else NA_real_)
  if (xor(is.null(lower), is.null(upper)))
    stop("lower and upper must be supplied together")
  if (!is.null(lower)) {
    lower <- .metric_numeric_vector(lower, "lower")
    upper <- .metric_numeric_vector(upper, "upper")
    if (length(lower) != length(estimates) ||
        length(upper) != length(estimates) || any(lower > upper))
      stop("interval bounds must match estimates and satisfy lower <= upper")
    covered <- lower <= truth & truth <= upper
    out$coverage <- mean(covered)
    out$coverage_mcse <- sqrt(out$coverage * (1 - out$coverage) /
                                length(covered))
    out$mean_width <- mean(upper - lower)
  }
  out
}

block_interval_coverage <- function(truth, lower, upper) {
  truth <- .metric_numeric_vector(truth, "truth")
  lower <- as.matrix(lower); upper <- as.matrix(upper)
  if (!identical(dim(lower), dim(upper)) || ncol(lower) != length(truth) ||
      !nrow(lower) || any(!is.finite(lower)) || any(!is.finite(upper)) ||
      any(lower > upper))
    stop("lower and upper must be finite replicate by parameter matrices")
  covered <- sweep(lower, 2L, truth, `<=`) &
    sweep(upper, 2L, truth, `>=`)
  per_replicate <- rowMeans(covered)
  width_per_replicate <- rowMeans(upper - lower)
  list(
    replicates = nrow(lower),
    parameters = ncol(lower),
    coverage = mean(per_replicate),
    coverage_mcse = if (nrow(lower) > 1L)
      stats::sd(per_replicate) / sqrt(nrow(lower)) else NA_real_,
    mean_width = mean(width_per_replicate),
    per_replicate = data.frame(
      replicate = seq_len(nrow(lower)), coverage = per_replicate,
      mean_width = width_per_replicate))
}

posterior_fidelity <- function(reference_draws, approximate_draws) {
  reference <- as.matrix(reference_draws)
  approximate <- as.matrix(approximate_draws)
  if (!nrow(reference) || !nrow(approximate) ||
      ncol(reference) != ncol(approximate) || !ncol(reference) ||
      any(!is.finite(reference)) || any(!is.finite(approximate)))
    stop("draw matrices must be finite and have the same variables in columns")
  if (!is.null(colnames(reference)) && !is.null(colnames(approximate)) &&
      !identical(colnames(reference), colnames(approximate)))
    stop("draw matrices must use the same variable order")
  dimension <- ncol(reference)
  reference_sd <- apply(reference, 2L, stats::sd)
  approximate_sd <- apply(approximate, 2L, stats::sd)
  if (any(!is.finite(reference_sd)) || any(reference_sd <= 0) ||
      any(!is.finite(approximate_sd)))
    stop("reference variables must have positive finite standard deviations")
  reference_correlation <- stats::cor(reference)
  approximate_correlation <- stats::cor(approximate)
  difference <- approximate_correlation - reference_correlation
  off_diagonal <- row(difference) != col(difference)
  list(
    dimension = dimension,
    mean_distance = sqrt(sum((colMeans(approximate) -
                                colMeans(reference))^2) / dimension),
    mean_relative_sd_difference = mean(
      (approximate_sd - reference_sd) / reference_sd),
    mean_sd_ratio = mean(approximate_sd / reference_sd),
    correlation_frobenius = sqrt(sum(difference^2)) / dimension,
    maximum_off_diagonal_correlation_error = if (dimension == 1L) 0 else
      max(abs(difference[off_diagonal])),
    reference_correlation = reference_correlation,
    approximate_correlation = approximate_correlation)
}

predictive_regret <- function(mlpd, selected_index) {
  mlpd <- .metric_numeric_vector(mlpd, "mlpd")
  if (!is.numeric(selected_index) || length(selected_index) != 1L ||
      !is.finite(selected_index) || selected_index != floor(selected_index) ||
      selected_index < 1L || selected_index > length(mlpd))
    stop("selected_index must identify one MLPD entry")
  max(mlpd) - mlpd[as.integer(selected_index)]
}

selection_recovery <- function(selected, truth, candidate_values = NULL,
                               regret = NULL) {
  selected <- .metric_numeric_vector(selected, "selected")
  if (!is.numeric(truth) || length(truth) != 1L || !is.finite(truth))
    stop("truth must be one finite value")
  if (is.null(candidate_values)) candidate_values <- sort(unique(selected))
  candidate_values <- .metric_numeric_vector(candidate_values, "candidate_values")
  if (any(!selected %in% candidate_values))
    stop("all selections must belong to candidate_values")
  distribution <- data.frame(
    selected = candidate_values,
    count = vapply(candidate_values, function(value) sum(selected == value),
                   integer(1L)),
    probability = vapply(candidate_values, function(value)
      mean(selected == value), numeric(1L)))
  out <- list(
    replicates = length(selected),
    exact_recovery = mean(selected == truth),
    exact_recovery_mcse = sqrt(mean(selected == truth) *
      (1 - mean(selected == truth)) / length(selected)),
    mean_absolute_error = mean(abs(selected - truth)),
    distribution = distribution)
  if (!is.null(regret)) {
    regret <- .metric_numeric_vector(regret, "regret")
    if (length(regret) != length(selected) || any(regret < 0))
      stop("regret must be nonnegative with one value per selection")
    out$mean_predictive_regret <- mean(regret)
    out$predictive_regret_mcse <- if (length(regret) > 1L)
      stats::sd(regret) / sqrt(length(regret)) else NA_real_
  }
  out
}

match_components_to_truth <- function(
    estimated_coefficients, estimated_variances,
    true_coefficients, true_variances, coordinate_scale = NULL) {
  estimated_coefficients <- as.matrix(estimated_coefficients)
  true_coefficients <- as.matrix(true_coefficients)
  estimated_variances <- .metric_numeric_vector(
    estimated_variances, "estimated_variances")
  true_variances <- .metric_numeric_vector(true_variances, "true_variances")
  if (!identical(dim(estimated_coefficients), dim(true_coefficients)) ||
      nrow(estimated_coefficients) != length(estimated_variances) ||
      length(estimated_variances) != length(true_variances) ||
      any(!is.finite(estimated_coefficients)) ||
      any(!is.finite(true_coefficients)) ||
      any(estimated_variances <= 0) || any(true_variances <= 0))
    stop("estimated and true component parameters must have compatible finite dimensions")
  estimated_location <- cbind(estimated_coefficients,
                              log(estimated_variances))
  true_location <- cbind(true_coefficients, log(true_variances))
  if (is.null(coordinate_scale)) {
    coordinate_scale <- apply(
      rbind(estimated_location, true_location), 2L, stats::sd)
    coordinate_scale[!is.finite(coordinate_scale) |
                       coordinate_scale <= sqrt(.Machine$double.eps)] <- 1
  } else {
    coordinate_scale <- .metric_numeric_vector(
      coordinate_scale, "coordinate_scale")
    if (length(coordinate_scale) != ncol(estimated_location) ||
        any(coordinate_scale <= 0))
      stop("coordinate_scale must contain one positive value per location coordinate")
  }
  K <- nrow(true_location)
  cost <- matrix(NA_real_, K, K)
  for (truth_index in seq_len(K)) for (estimated_index in seq_len(K)) {
    difference <- (estimated_location[estimated_index, ] -
                     true_location[truth_index, ]) / coordinate_scale
    cost[truth_index, estimated_index] <- sum(difference^2)
  }
  assignment <- as.integer(clue::solve_LSAP(cost))
  structure(list(
    schema_version = "paper-component-alignment-1",
    estimated_index_for_truth = assignment,
    coordinate_scale = coordinate_scale,
    total_cost = sum(cost[cbind(seq_len(K), assignment)]),
    cost = cost), class = c("paper_component_alignment", "list"))
}

.chlrm_parameter_matrix <- function(parameters, name, value = "estimate") {
  rows <- parameters[parameters$parameter == name, , drop = FALSE]
  if (!nrow(rows) || !value %in% names(rows))
    stop("parameter table does not contain ", name, " with field ", value)
  if (all(is.na(rows$index_2))) {
    output <- rows[[value]][order(rows$index_1)]
    return(as.numeric(output))
  }
  K <- max(rows$index_1); p <- max(rows$index_2)
  output <- matrix(NA_real_, K, p)
  output[cbind(rows$index_1, rows$index_2)] <- rows[[value]]
  if (any(!is.finite(output)))
    stop("parameter table has an incomplete component matrix")
  output
}

aligned_chlrm_recovery <- function(parameters, truth,
                                   coordinate_scale = NULL) {
  if (!inherits(parameters, "vb_tidy_parameters") ||
      !identical(unique(parameters$model), "CHLRM"))
    stop("parameters must be CHLRM output from tidy_parameters()")
  required_truth <- c("component_coefficients", "component_variances",
                      "mixture_weights")
  if (!is.list(truth) || !all(required_truth %in% names(truth)))
    stop("truth must contain CHLRM component parameters")
  estimated_coefficients <- .chlrm_parameter_matrix(
    parameters, "component_coefficients")
  estimated_variances <- .chlrm_parameter_matrix(
    parameters, "component_variances")
  if (any(!is.finite(estimated_variances) | estimated_variances <= 0)) {
    variance_rows <- parameters[
      parameters$parameter == "component_variances", , drop = FALSE]
    order <- order(variance_rows$index_1)
    variance_rows <- variance_rows[order, , drop = FALSE]
    fallback <- exp((log(variance_rows$lower) + log(variance_rows$upper)) / 2)
    replace <- !is.finite(estimated_variances) | estimated_variances <= 0
    estimated_variances[replace] <- fallback[replace]
  }
  alignment <- match_components_to_truth(
    estimated_coefficients, estimated_variances,
    truth$component_coefficients, truth$component_variances,
    coordinate_scale)
  permutation <- alignment$estimated_index_for_truth
  blocks <- c("component_coefficients", "component_variances",
              "mixture_weights")
  rows <- lapply(blocks, function(name) {
    block <- parameters[parameters$parameter == name, , drop = FALSE]
    aligned <- do.call(rbind, lapply(seq_along(permutation), function(k) {
      selected <- block[block$index_1 == permutation[k], , drop = FALSE]
      selected$estimated_component <- selected$index_1
      selected$true_component <- k
      selected$index_1 <- k
      selected
    }))
    truth_value <- switch(name,
      component_coefficients = truth$component_coefficients[
        cbind(aligned$true_component, aligned$index_2)],
      component_variances = truth$component_variances[aligned$true_component],
      mixture_weights = truth$mixture_weights[aligned$true_component])
    aligned$truth <- as.numeric(truth_value)
    aligned$error <- aligned$estimate - aligned$truth
    aligned$covered <- aligned$lower <= aligned$truth &
      aligned$truth <= aligned$upper
    aligned$interval_width <- aligned$upper - aligned$lower
    aligned
  })
  aligned <- do.call(rbind, rows)
  rownames(aligned) <- NULL
  finite_mean <- function(value)
    if (any(is.finite(value))) mean(value[is.finite(value)]) else NA_real_
  summaries <- do.call(rbind, lapply(split(aligned, aligned$parameter), function(x)
    data.frame(
      parameter = x$parameter[1L], scalar_parameters = nrow(x),
      finite_posterior_means = sum(is.finite(x$estimate)),
      mean_error = finite_mean(x$error),
      mean_absolute_error = finite_mean(abs(x$error)),
      root_mean_squared_error = sqrt(finite_mean(x$error^2)),
      coverage = mean(x$covered), mean_width = mean(x$interval_width),
      mean_posterior_sd = finite_mean(x$sd), stringsAsFactors = FALSE)))
  rownames(summaries) <- NULL
  structure(list(
    schema_version = "paper-chlrm-recovery-1",
    alignment = alignment, pointwise = aligned, summary = summaries),
    class = c("paper_chlrm_recovery", "list"))
}

clustering_stability <- function(orchestration) {
  if (!inherits(orchestration, "vb_fit_orchestration") ||
      !identical(orchestration$type, "multistart") ||
      !identical(orchestration$model, "CHLRM") ||
      !is.data.frame(orchestration$start_classifications))
    stop("orchestration must be a CHLRM multistart audit object")
  classifications <- orchestration$start_classifications
  starts <- sort(unique(classifications$start))
  if (length(starts) < 2L)
    return(list(starts = length(starts), pairs = 0L, stability = NA_real_,
                pairwise = data.frame()))
  labels <- lapply(starts, function(start) {
    values <- classifications[classifications$start == start, ]
    values <- values[order(values$group), ]
    values$component
  })
  if (length(unique(lengths(labels))) != 1L)
    stop("completed starts contain incompatible group classifications")
  pairs <- utils::combn(seq_along(starts), 2L)
  values <- vapply(seq_len(ncol(pairs)), function(index)
    adjusted_rand_index(labels[[pairs[1L, index]]],
                        labels[[pairs[2L, index]]]), numeric(1L))
  pairwise <- data.frame(
    start_1 = starts[pairs[1L, ]], start_2 = starts[pairs[2L, ]],
    adjusted_rand_index = values)
  list(starts = length(starts), pairs = length(values),
       stability = mean(values), pairwise = pairwise)
}
