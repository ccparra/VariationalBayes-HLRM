# Resource review and execution helpers for the E3 pilot stages.

.e3_required_screen_columns <- c(
  "task_id", "method", "minibatch_size", "exponent", "offset", "status",
  "final_elbo", "mlpd", "ari", "cpe", "runtime",
  "iterations", "convergence_status", "result_file_mb",
  "worker_wall_seconds", "measured_peak_rss_mb", "rss_measurement",
  "output_fingerprint", "git_commit", "git_dirty")

review_e3_screen <- function(results, factorial_task_count = 109L,
                             max_memory_mb = 2048,
                             rss_safety_factor = 2,
                             first_eligible_svi_stop = 1200L) {
  if (!is.data.frame(results) || !nrow(results))
    stop("results must be a nonempty E3 screen data frame")
  missing <- setdiff(.e3_required_screen_columns, names(results))
  if (length(missing))
    stop("E3 screen results are missing: ", paste(missing, collapse = ", "))
  factorial_task_count <- .paper_positive_integer(
    factorial_task_count, "factorial_task_count")
  if (!is.numeric(max_memory_mb) || length(max_memory_mb) != 1L ||
      !is.finite(max_memory_mb) || max_memory_mb <= 0)
    stop("max_memory_mb must be positive and finite")
  if (!is.numeric(rss_safety_factor) || length(rss_safety_factor) != 1L ||
      !is.finite(rss_safety_factor) || rss_safety_factor < 1)
    stop("rss_safety_factor must be finite and at least one")
  first_eligible_svi_stop <- .paper_positive_integer(
    first_eligible_svi_stop, "first_eligible_svi_stop")

  cavi <- results[results$method == "CAVI", , drop = FALSE]
  svi <- results[results$method == "SVI", , drop = FALSE]
  unique_value <- function(value) {
    candidate <- unique(value)
    if (length(candidate) == 1L) candidate[[1L]] else NA
  }
  all_complete <- all(results$status %in% c("complete", "reused"))
  resources_measured <- all(
    is.finite(results$worker_wall_seconds) &
      is.finite(results$measured_peak_rss_mb) &
      nzchar(results$rss_measurement))
  max_rss <- if (resources_measured) max(results$measured_peak_rss_mb) else NA_real_
  memory_safe <- resources_measured &&
    rss_safety_factor * max_rss <= max_memory_mb
  clean_provenance <- all(!is.na(results$git_dirty) & !results$git_dirty) &&
    length(unique(results$git_commit)) == 1L &&
    length(unique(results$output_fingerprint)) == 1L
  all_converged <- all(results$convergence_status == "converged")
  shape_valid <- nrow(results) == 10L && nrow(cavi) == 1L && nrow(svi) == 9L &&
    anyDuplicated(results$task_id) == 0L

  mean_svi_wall <- if (nrow(svi) && all(is.finite(svi$worker_wall_seconds)))
    mean(svi$worker_wall_seconds) else NA_real_
  projected_wall <- if (nrow(cavi) == 1L && is.finite(mean_svi_wall))
    cavi$worker_wall_seconds[[1L]] +
      (factorial_task_count - 1L) * mean_svi_wall else NA_real_
  mean_svi_storage <- if (nrow(svi) && all(is.finite(svi$result_file_mb)))
    mean(svi$result_file_mb) else NA_real_
  projected_storage <- if (nrow(cavi) == 1L && is.finite(mean_svi_storage))
    cavi$result_file_mb[[1L]] +
      (factorial_task_count - 1L) * mean_svi_storage else NA_real_
  earliest_fraction <- if (nrow(svi))
    mean(svi$iterations == first_eligible_svi_stop) else NA_real_

  gate <- data.frame(
    requirement = c(
      "screen_shape", "all_tasks_complete", "all_tasks_converged",
      "resource_measurements", "rss_safety_margin", "clean_provenance"),
    status = c(
      if (shape_valid) "pass" else "fail",
      if (all_complete) "pass" else "fail",
      if (all_converged) "pass" else "fail",
      if (resources_measured) "pass" else "fail",
      if (memory_safe) "pass" else "fail",
      if (clean_provenance) "pass" else "fail"),
    detail = c(
      sprintf("%d tasks: %d CAVI and %d SVI", nrow(results), nrow(cavi), nrow(svi)),
      sprintf("%d of %d complete/reused", sum(results$status %in% c("complete", "reused")), nrow(results)),
      sprintf("%d of %d diagnosed converged", sum(results$convergence_status == "converged"), nrow(results)),
      sprintf("%d of %d with finite wall time and RSS", sum(
        is.finite(results$worker_wall_seconds) &
          is.finite(results$measured_peak_rss_mb)), nrow(results)),
      if (resources_measured) sprintf(
        "%.1f x %.1f MiB = %.1f MiB <= %.1f MiB",
        rss_safety_factor, max_rss, rss_safety_factor * max_rss,
        max_memory_mb) else "RSS unavailable",
      if (clean_provenance) paste(
        unique_value(results$git_commit),
        unique_value(results$output_fingerprint)) else
          "dirty, missing, or inconsistent provenance"),
    stringsAsFactors = FALSE)
  ready <- all(gate$status == "pass")
  summary <- data.frame(
    schema_version = "paper-e3-screen-review-1",
    screen_task_count = nrow(results),
    factorial_task_count = factorial_task_count,
    all_tasks_converged = all_converged,
    svi_first_eligible_stop_fraction = earliest_fraction,
    observed_total_wall_seconds = sum(results$worker_wall_seconds),
    observed_max_rss_mb = max_rss,
    rss_safety_factor = rss_safety_factor,
    memory_budget_mb = max_memory_mb,
    projected_factorial_wall_seconds = projected_wall,
    projected_factorial_storage_mb = projected_storage,
    screen_output_fingerprint = unique_value(results$output_fingerprint),
    screen_git_commit = unique_value(results$git_commit),
    ready_for_factorial = ready,
    stringsAsFactors = FALSE)
  ranking <- svi[order(-svi$mlpd, -svi$ari, svi$cpe), c(
    "task_id", "minibatch_size", "exponent", "offset", "final_elbo",
    "mlpd", "ari", "cpe", "runtime", "iterations",
    "convergence_status"), drop = FALSE]
  ranking$predictive_rank <- seq_len(nrow(ranking))
  structure(list(
    schema_version = "paper-e3-screen-review-1",
    summary = summary,
    gate = gate,
    ranking = ranking), class = c("paper_e3_screen_review", "list"))
}

write_e3_screen_review <- function(screen_dir, factorial_task_count = 109L,
                                   max_memory_mb = 2048,
                                   rss_safety_factor = 2) {
  results_path <- file.path(screen_dir, "fit-results.csv")
  if (!file.exists(results_path))
    stop("screen_dir does not contain fit-results.csv")
  results <- utils::read.csv(results_path, stringsAsFactors = FALSE,
                             check.names = FALSE)
  review <- review_e3_screen(
    results, factorial_task_count = factorial_task_count,
    max_memory_mb = max_memory_mb, rss_safety_factor = rss_safety_factor)
  utils::write.csv(review$summary, file.path(screen_dir, "screen-review.csv"),
                   row.names = FALSE)
  utils::write.csv(review$gate, file.path(screen_dir, "screen-gate.csv"),
                   row.names = FALSE)
  utils::write.csv(review$ranking, file.path(screen_dir, "screen-ranking.csv"),
                   row.names = FALSE)
  saveRDS(review, file.path(screen_dir, "screen-review.rds"), version = 3)
  invisible(review)
}

.paper_peak_rss_mb <- function(output, platform = Sys.info()[["sysname"]]) {
  if (identical(platform, "Darwin")) {
    line <- output[grepl("maximum resident set size", output, fixed = TRUE)]
    if (!length(line)) return(NA_real_)
    bytes <- suppressWarnings(as.numeric(sub(
      "^\\s*([0-9]+).*$", "\\1", line[[1L]])))
    return(bytes / 1024^2)
  }
  line <- output[grepl("Maximum resident set size", output, fixed = TRUE)]
  if (!length(line)) return(NA_real_)
  kilobytes <- suppressWarnings(as.numeric(sub(
    "^.*:\\s*([0-9]+).*$", "\\1", line[[1L]])))
  kilobytes / 1024
}

