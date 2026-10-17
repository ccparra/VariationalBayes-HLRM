# Internal worker for one frozen E3 pilot task. The orchestrator invokes this
# in a separate process so operating-system peak RSS can be measured.
.cli_value <- function(arguments, prefix, required = TRUE) {
  values <- arguments[startsWith(arguments, prefix)]
  if (!length(values)) {
    if (required) stop("Missing required argument ", prefix)
    return(NULL)
  }
  sub(prefix, "", values[[1L]], fixed = TRUE)
}

arguments <- commandArgs(trailingOnly = TRUE)
if (!"--execute" %in% arguments)
  stop("This internal worker requires the explicit --execute flag")
root <- normalizePath(.cli_value(arguments, "--root="), mustWork = TRUE)
library_path <- normalizePath(.cli_value(arguments, "--library="), mustWork = TRUE)
frozen_path <- normalizePath(.cli_value(arguments, "--frozen="), mustWork = TRUE)
task_id <- .cli_value(arguments, "--task-id=")
output_path <- .cli_value(arguments, "--output=")
.libPaths(c(library_path, .libPaths()))
suppressPackageStartupMessages(library(VariationalBayesHLRM))
source(file.path(root, "paper", "evaluation", "R", "project.R"))
for (file in c("scenarios.R", "metrics.R", "aggregation.R", "provenance.R",
               "tasks.R", "e3.R", "e3-execution.R"))
  source(file.path(root, "paper", "evaluation", "R", file))

frozen <- readRDS(frozen_path)
stage <- frozen$stage
validate_e3_frozen_configuration(frozen, stage)
tasks <- frozen$configuration$tasks
task <- tasks[tasks$task_id == task_id, , drop = FALSE]
if (nrow(task) != 1L)
  stop("task-id must identify exactly one task in the frozen E3 stage")
provenance <- paper_run_provenance(
  "E3", paste0("pilot-", stage), frozen$configuration)
result <- fit_e3_task(
  task, pilot_authorized = TRUE,
  max_memory_mb = unique(frozen$configuration$controls$max_memory_mb))
result$provenance <- provenance
result$frozen_configuration <- list(
  schema_version = frozen$schema_version,
  stage = stage,
  configuration_fingerprint = frozen$configuration_fingerprint,
  source_fingerprint = frozen$source_fingerprint,
  git_commit = frozen$git_commit,
  upstream_evidence = frozen$upstream_evidence)
dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
.atomic_save_rds(result, output_path)
cat("completed_task_id=", task_id, "\n", sep = "")
cat("result_path=", normalizePath(output_path), "\n", sep = "")

