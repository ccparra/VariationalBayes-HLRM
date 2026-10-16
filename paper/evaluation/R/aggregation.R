# Aggregation from immutable per-task artifacts. No function in this file fits
# a model or changes an estimand.

e1_task_metric_row <- function(result) {
  if (!inherits(result, "paper_task_result") ||
      !identical(result$schema_version, "paper-task-result-4"))
    stop("result must be a paper-task-result-4 object")
  held_out <- result$evaluation$held_out$statistics
  clustering <- result$allocations$metrics
  recovery <- result$recovery$summary
  recovery_value <- function(parameter, field) {
    row <- recovery[recovery$parameter == parameter, , drop = FALSE]
    if (nrow(row) != 1L) stop("Task recovery summary is incomplete")
    row[[field]]
  }
  data.frame(
    task_id = result$task$task_id,
    scenario = result$task$scenario,
    replicate = result$task$replicate,
    method = result$task$method,
    convergence_status = result$fit$diagnostics$status,
    fit_seconds = result$fit$timing$total_seconds,
    ari = clustering$adjusted_rand_index,
    cpe = clustering$co_clustering_error,
    normalized_entropy = clustering$normalized_entropy,
    coefficient_rmse_per_dataset = recovery_value(
      "component_coefficients", "root_mean_squared_error"),
    coefficient_coverage_per_dataset = recovery_value(
      "component_coefficients", "coverage"),
    variance_rmse_per_dataset = recovery_value(
      "component_variances", "root_mean_squared_error"),
    variance_coverage_per_dataset = recovery_value(
      "component_variances", "coverage"),
    weight_rmse_per_dataset = recovery_value(
      "mixture_weights", "root_mean_squared_error"),
    mlpd = held_out$mean_log_predictive_density,
    test_rmse = held_out$rmse,
    test_mae = held_out$mae,
    predictive_coverage = held_out$empirical_coverage,
    mean_interval_width = held_out$mean_interval_width,
    stringsAsFactors = FALSE)
}

collect_e1_metrics <- function(manifest, output_dir) {
  if (!is.data.frame(manifest) ||
      !all(c("task_id", "path") %in% names(manifest)))
    stop("manifest must contain task_id and path")
  paths <- file.path(output_dir, manifest$path)
  if (any(!file.exists(paths)))
    stop("manifest references missing task artifacts")
  rows <- lapply(paths, function(path) e1_task_metric_row(readRDS(path)))
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out[match(manifest$task_id, out$task_id), , drop = FALSE]
}

aggregate_e1_metrics <- function(metrics) {
  required <- c("scenario", "method", "ari", "cpe", "normalized_entropy",
                "coefficient_rmse_per_dataset",
                "coefficient_coverage_per_dataset",
                "variance_rmse_per_dataset",
                "variance_coverage_per_dataset", "weight_rmse_per_dataset",
                "mlpd", "test_rmse", "test_mae", "predictive_coverage",
                "mean_interval_width", "fit_seconds")
  if (!is.data.frame(metrics) || !all(required %in% names(metrics)) ||
      !nrow(metrics))
    stop("metrics must contain the columns produced by collect_e1_metrics()")
  numeric_metrics <- setdiff(required, c("scenario", "method"))
  groups <- split(seq_len(nrow(metrics)),
                  interaction(metrics$scenario, metrics$method, drop = TRUE))
  rows <- lapply(groups, function(index) {
    group <- metrics[index, , drop = FALSE]
    values <- c(
      list(scenario = group$scenario[1L], method = group$method[1L],
           replicates = nrow(group)),
      unlist(lapply(numeric_metrics, function(name) {
        value <- group[[name]]
        stats::setNames(c(
          mean(value),
          if (length(value) > 1L) stats::sd(value) / sqrt(length(value)) else
            NA_real_), c(paste0(name, "_mean"), paste0(name, "_mcse")))
      }), use.names = TRUE))
    as.data.frame(as.list(values), stringsAsFactors = FALSE,
                  check.names = FALSE)
  })
  out <- do.call(rbind, rows)
  numeric_columns <- setdiff(names(out), c("scenario", "method"))
  out[numeric_columns] <- lapply(out[numeric_columns], as.numeric)
  rownames(out) <- NULL
  out[order(out$scenario,
            match(out$method, c("MCMC", "CAVI", "SVI"))), , drop = FALSE]
}

e1_parameter_recovery_rows <- function(result) {
  if (!inherits(result, "paper_task_result") ||
      !identical(result$schema_version, "paper-task-result-4"))
    stop("result must be a paper-task-result-4 object")
  rows <- result$recovery$pointwise
  required <- c("parameter", "index_1", "index_2", "truth", "estimate",
                "sd", "lower", "upper", "error", "covered",
                "interval_width")
  if (!is.data.frame(rows) || !all(required %in% names(rows)))
    stop("result does not contain aligned pointwise recovery")
  rows <- rows[required]
  rows$task_id <- result$task$task_id
  rows$scenario <- result$task$scenario
  rows$replicate <- result$task$replicate
  rows$method <- result$task$method
  rows$absolute_error <- abs(rows$error)
  rows$squared_error <- rows$error^2
  provenance <- result$provenance
  if (inherits(provenance, "paper_run_provenance")) {
    rows$run_fingerprint <- provenance$run_fingerprint
    rows$configuration_fingerprint <- provenance$configuration_fingerprint
    rows$source_fingerprint <- provenance$source_fingerprint
    rows$suite_source_fingerprint <- provenance$suite_source_fingerprint
    rows$output_fingerprint <- provenance$output_fingerprint
  }
  leading <- c("task_id", "scenario", "replicate", "method")
  rows[c(leading, setdiff(names(rows), leading))]
}

collect_e1_parameter_recovery <- function(manifest, output_dir) {
  if (!is.data.frame(manifest) ||
      !all(c("task_id", "path") %in% names(manifest)))
    stop("manifest must contain task_id and path")
  paths <- file.path(output_dir, manifest$path)
  if (any(!file.exists(paths)))
    stop("manifest references missing task artifacts")
  rows <- lapply(paths, function(path)
    e1_parameter_recovery_rows(readRDS(path)))
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

.paper_mean_mcse <- function(value) {
  value <- value[is.finite(value)]
  if (!length(value)) return(c(mean = NA_real_, mcse = NA_real_))
  c(mean = mean(value), mcse = if (length(value) > 1L)
    stats::sd(value) / sqrt(length(value)) else NA_real_)
}

aggregate_e1_parameter_recovery <- function(recovery) {
  required <- c("scenario", "method", "replicate", "parameter", "index_1",
                "index_2", "truth", "estimate", "sd", "error", "covered",
                "interval_width")
  if (!is.data.frame(recovery) || !nrow(recovery) ||
      !all(required %in% names(recovery)))
    stop("recovery must contain aligned pointwise parameter rows")
  key <- paste(
    recovery$scenario, recovery$method, recovery$parameter,
    ifelse(is.na(recovery$index_1), "", recovery$index_1),
    ifelse(is.na(recovery$index_2), "", recovery$index_2), sep = "\r")
  groups <- split(seq_len(nrow(recovery)), key)
  rows <- lapply(groups, function(index) {
    group <- recovery[index, , drop = FALSE]
    truth <- unique(group$truth[is.finite(group$truth)])
    if (length(truth) != 1L)
      stop("Each aligned parameter coordinate must have one finite truth")
    error <- .paper_mean_mcse(group$error)
    absolute_error <- .paper_mean_mcse(abs(group$error))
    squared_error <- .paper_mean_mcse(group$error^2)
    coverage <- .paper_mean_mcse(as.numeric(group$covered))
    width <- .paper_mean_mcse(group$interval_width)
    posterior_sd <- .paper_mean_mcse(group$sd)
    data.frame(
      scenario = group$scenario[1L], method = group$method[1L],
      parameter = group$parameter[1L], index_1 = group$index_1[1L],
      index_2 = group$index_2[1L], truth = truth,
      replicates = length(unique(group$replicate)),
      finite_posterior_means = sum(is.finite(group$estimate)),
      bias = error["mean"], bias_mcse = error["mcse"],
      mean_absolute_error = absolute_error["mean"],
      mean_absolute_error_mcse = absolute_error["mcse"],
      rmse = sqrt(squared_error["mean"]),
      coverage = coverage["mean"], coverage_mcse = coverage["mcse"],
      mean_interval_width = width["mean"],
      mean_interval_width_mcse = width["mcse"],
      mean_posterior_sd = posterior_sd["mean"],
      mean_posterior_sd_mcse = posterior_sd["mcse"],
      run_fingerprint = if ("run_fingerprint" %in% names(group))
        unique(group$run_fingerprint)[1L] else NA_character_,
      stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out[order(out$scenario, match(out$method, c("MCMC", "CAVI", "SVI")),
            out$parameter, out$index_1, out$index_2,
            na.last = TRUE), , drop = FALSE]
}

write_e1_metric_outputs <- function(manifest, output_dir) {
  raw <- collect_e1_metrics(manifest, output_dir)
  summary <- aggregate_e1_metrics(raw)
  parameter_raw <- collect_e1_parameter_recovery(manifest, output_dir)
  parameter_summary <- aggregate_e1_parameter_recovery(parameter_raw)
  utils::write.csv(raw, file.path(output_dir, "metrics-raw.csv"),
                   row.names = FALSE)
  utils::write.csv(summary, file.path(output_dir, "metrics-summary.csv"),
                   row.names = FALSE)
  utils::write.csv(
    parameter_raw, file.path(output_dir, "parameter-recovery-raw.csv"),
    row.names = FALSE)
  utils::write.csv(
    parameter_summary,
    file.path(output_dir, "parameter-recovery-summary.csv"), row.names = FALSE)
  invisible(list(
    raw = raw, summary = summary,
    parameter_raw = parameter_raw, parameter_summary = parameter_summary))
}
