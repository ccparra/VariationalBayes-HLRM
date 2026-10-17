# Execute one frozen E3 pilot stage sequentially, one task per fresh R process.
# The factorial stage is unavailable until a passing screen review has been
# explicitly accepted and embedded in its frozen configuration.
.cli_value <- function(arguments, prefix, required = TRUE) {
  values <- arguments[startsWith(arguments, prefix)]
  if (!length(values)) {
    if (required) stop("Missing required argument ", prefix)
    return(NULL)
  }
  sub(prefix, "", values[[1L]], fixed = TRUE)
}

.atomic_write_csv <- function(value, path) {
  temporary <- tempfile(paste0(basename(path), "-"), tmpdir = dirname(path),
                        fileext = ".tmp")
  on.exit(if (file.exists(temporary)) unlink(temporary), add = TRUE)
  utils::write.csv(value, temporary, row.names = FALSE)
  if (!file.rename(temporary, path))
    stop("Could not atomically replace ", path)
  invisible(path)
}

arguments <- commandArgs(trailingOnly = TRUE)
if (!"--execute" %in% arguments)
  stop(paste(
    "E3 pilot not started. Review/freeze its preflight, then pass --execute",
    "and --accept-configuration-fingerprint=<md5:...>."))
stage <- match.arg(.cli_value(arguments, "--stage="), c("screen", "factorial"))
accepted <- .cli_value(arguments, "--accept-configuration-fingerprint=")
accepted_screen <- .cli_value(
  arguments, "--accept-screen-output-fingerprint=", required = FALSE)
root <- normalizePath(".", mustWork = TRUE)
source(file.path(root, "paper", "evaluation", "R", "project.R"))
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1")
preflight_dir <- file.path(
  root, "paper", "evaluation", "output", "preflight", "E3",
  paste0("pilot-", stage))
frozen_path <- file.path(preflight_dir, "frozen-configuration.rds")
if (!file.exists(frozen_path))
  stop("No frozen E3 ", stage, " configuration exists")

for (file in c("scenarios.R", "metrics.R", "aggregation.R", "provenance.R",
               "tasks.R", "e3.R", "e3-execution.R"))
  source(file.path(root, "paper", "evaluation", "R", file))
frozen <- readRDS(frozen_path)
validate_e3_frozen_configuration(frozen, stage)
if (!identical(accepted, frozen$configuration_fingerprint))
  stop("Explicitly accepted fingerprint does not match the frozen E3 stage")
if (identical(stage, "factorial")) {
  if (!is.character(accepted_screen) || length(accepted_screen) != 1L ||
      !identical(accepted_screen,
                 frozen$upstream_evidence$screen_output_fingerprint))
    stop("Explicitly accepted screen fingerprint does not match factorial evidence")
} else if (!is.null(accepted_screen)) {
  stop("--accept-screen-output-fingerprint applies only to factorial")
}

summary_path <- file.path(preflight_dir, "budget-summary.csv")
budget <- utils::read.csv(summary_path, stringsAsFactors = FALSE)
output_dir <- file.path(
  root, "paper", "evaluation", "output", "pilot", "E3", stage)
dir.create(file.path(output_dir, "tasks"), recursive = TRUE, showWarnings = FALSE)
disk <- system2("df", c("-Pk", shQuote(output_dir)), stdout = TRUE, stderr = TRUE)
disk_fields <- strsplit(trimws(tail(disk, 1L)), "[[:space:]]+")[[1L]]
available_mb <- suppressWarnings(as.numeric(disk_fields[[4L]]) / 1024)
projected_storage_mb <- if (identical(stage, "factorial"))
  frozen$upstream_evidence$projected_factorial_storage_mb else
    budget$total_estimated_fit_storage_mb
required_mb <- max(1024, 2 * projected_storage_mb)
if (!is.finite(available_mb) || available_mb < required_mb)
  stop(sprintf(
    "Insufficient disk: %.1f MiB available, %.1f MiB safety requirement",
    available_mb, required_mb))

validation_library <- tempfile("vb-hlrm-paper-e3-pilot-library-")
on.exit(unlink(validation_library, recursive = TRUE, force = TRUE), add = TRUE)
install_paper_package(root, validation_library)

tasks <- frozen$configuration$tasks
provenance <- paper_run_provenance(
  "E3", paste0("pilot-", stage), frozen$configuration)
worker <- file.path(root, "paper", "evaluation", "run-e3-task.R")
time_binary <- "/usr/bin/time"
if (!file.exists(time_binary))
  stop("/usr/bin/time is required for operating-system peak RSS measurement")
time_flag <- if (identical(Sys.info()[["sysname"]], "Darwin")) "-l" else "-v"
worker_environment <- c(
  "OMP_NUM_THREADS=1", "OPENBLAS_NUM_THREADS=1", "MKL_NUM_THREADS=1")
rows <- vector("list", nrow(tasks))
for (index in seq_len(nrow(tasks))) {
  task <- tasks[index, , drop = FALSE]
  cat(sprintf("[%d/%d] %s\n", index, nrow(tasks), task$task_id))
  flush.console()
  result_path <- file.path(output_dir, "tasks", paste0(task$task_id, ".rds"))
  reused <- FALSE
  if (file.exists(result_path)) {
    result <- readRDS(result_path)
    reused <- .paper_result_compatible(
      result, "paper_e3_task_result", "paper-e3-task-result-3",
      task$task_id, provenance) &&
      identical(result$frozen_configuration$configuration_fingerprint,
                frozen$configuration_fingerprint)
  }
  worker_wall_seconds <- NA_real_
  peak_rss_mb <- NA_real_
  if (!reused) {
    started <- proc.time()[["elapsed"]]
    worker_output <- system2(
      time_binary,
      c(time_flag, file.path(R.home("bin"), "Rscript"), shQuote(worker),
        "--execute", paste0("--root=", shQuote(root)),
        paste0("--library=", shQuote(validation_library)),
        paste0("--frozen=", shQuote(frozen_path)),
        paste0("--task-id=", shQuote(task$task_id)),
        paste0("--output=", shQuote(result_path))),
      stdout = TRUE, stderr = TRUE, env = worker_environment)
    worker_wall_seconds <- proc.time()[["elapsed"]] - started
    status <- attr(worker_output, "status")
    writeLines(worker_output, file.path(
      output_dir, "tasks", paste0(task$task_id, ".log")))
    if (!is.null(status) && status != 0L)
      stop("E3 ", stage, " stopped after failed task ", task$task_id,
           "; inspect its log. Completed tasks remain reusable.")
    peak_rss_mb <- .paper_peak_rss_mb(worker_output)
    if (!is.finite(peak_rss_mb))
      stop("Could not parse peak RSS for task ", task$task_id)
    result <- readRDS(result_path)
    result$execution_measurement <- list(
      schema_version = "paper-task-resource-measurement-1",
      worker_wall_seconds = worker_wall_seconds,
      measured_peak_rss_mb = peak_rss_mb,
      rss_measurement = if (identical(Sys.info()[["sysname"]], "Darwin"))
        "usr-bin-time-darwin-bytes" else "usr-bin-time-linux-kbytes")
    .atomic_save_rds(result, result_path)
  } else if (!is.null(result$execution_measurement)) {
    worker_wall_seconds <- result$execution_measurement$worker_wall_seconds
    peak_rss_mb <- result$execution_measurement$measured_peak_rss_mb
  }
  row <- .e3_result_row(
    result, result_path, if (reused) "reused" else "complete")
  row$worker_wall_seconds <- worker_wall_seconds
  row$measured_peak_rss_mb <- peak_rss_mb
  row$rss_measurement <- if (!is.null(result$execution_measurement))
    result$execution_measurement$rss_measurement else NA_character_
  rows[[index]] <- row
  partial <- .paper_annotate_manifest(do.call(rbind, rows[seq_len(index)]),
                                      provenance)
  .atomic_write_csv(partial, file.path(output_dir, "manifest.csv"))
}
manifest <- .paper_annotate_manifest(do.call(rbind, rows), provenance)
comparison <- .paper_annotate_manifest(
  compare_e3_to_reference(manifest), provenance)
.atomic_write_csv(manifest, file.path(output_dir, "fit-results.csv"))
.atomic_write_csv(comparison, file.path(output_dir, "svi-comparison.csv"))
write_paper_provenance(provenance, output_dir)
cat("E3 pilot ", stage, " completed sequentially. These are ",
    if (identical(stage, "screen")) "profiling" else "pilot",
    " results, not confirmatory evidence.\n", sep = "")

