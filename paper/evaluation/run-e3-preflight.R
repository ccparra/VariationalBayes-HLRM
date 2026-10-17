# Materialize E3 controls and conservative resource estimates without fitting.
# This runner never executes CAVI or SVI. Run from the repository root.
root <- normalizePath(".", mustWork = TRUE)
source(file.path(root, "paper", "evaluation", "R", "project.R"))
validation_library <- tempfile("vb-hlrm-paper-e3-preflight-library-")
on.exit(unlink(validation_library, recursive = TRUE, force = TRUE), add = TRUE)
install_paper_package(root, validation_library)
suppressPackageStartupMessages(library(VariationalBayesHLRM))
for (file in c("scenarios.R", "metrics.R", "aggregation.R", "provenance.R",
               "tasks.R", "e3.R", "e3-execution.R"))
  source(file.path(root, "paper", "evaluation", "R", file))

arguments <- commandArgs(trailingOnly = TRUE)
profile_argument <- arguments[grepl("^--profile=", arguments)]
stage_argument <- arguments[grepl("^--stage=", arguments)]
memory_argument <- arguments[grepl("^--max-memory-mb=", arguments)]
profile <- if (length(profile_argument)) sub("^--profile=", "", profile_argument[[1L]]) else
  "pilot"
stage <- if (length(stage_argument)) sub("^--stage=", "", stage_argument[[1L]]) else
  "screen"
max_memory_mb <- if (length(memory_argument))
  as.numeric(sub("^--max-memory-mb=", "", memory_argument[[1L]])) else 2048
if (!is.finite(max_memory_mb) || max_memory_mb <= 0)
  stop("--max-memory-mb must be positive and finite")

output_dir <- file.path(root, "paper", "evaluation", "output", "preflight",
                        "E3", paste(profile, stage, sep = "-"))
screen_review <- if (identical(profile, "pilot") &&
                     identical(stage, "factorial")) {
  write_e3_screen_review(file.path(
    root, "paper", "evaluation", "output", "pilot", "E3", "screen"))
} else NULL
result <- prepare_e3_preflight(
  output_dir, profile = profile, stage = stage,
  max_memory_mb = max_memory_mb, screen_review = screen_review)
print(result$summary, row.names = FALSE)
cat("\nExecution gate (preflight never runs model fits):\n")
print(result$gate[c("requirement", "status", "detail")], row.names = FALSE)
