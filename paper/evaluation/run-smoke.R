# Run from the repository root after validating the installed package boundary.
library(testthat)

root <- normalizePath(".", mustWork = TRUE)
Sys.setenv(VB_HLRM_ROOT = root)

# Reuse the installation boundary established by the validation runner in this
# fresh process by installing into another temporary library.
validation_library <- tempfile("vb-hlrm-paper-smoke-library-")
dir.create(validation_library, recursive = TRUE, showWarnings = FALSE)
on.exit(unlink(validation_library, recursive = TRUE, force = TRUE), add = TRUE)
install_log <- system2(
  file.path(R.home("bin"), "R"),
  c("CMD", "INSTALL", "--no-multiarch",
    paste0("--library=", shQuote(validation_library)), shQuote(root)),
  stdout = TRUE, stderr = TRUE)
install_status <- attr(install_log, "status")
if (!is.null(install_status) && install_status != 0L)
  stop("Temporary package installation failed:\n", paste(install_log, collapse = "\n"))
.libPaths(c(validation_library, .libPaths()))
suppressPackageStartupMessages(library(VariationalBayesHLRM))
source(file.path(root, "paper", "evaluation", "R", "scenarios.R"))
source(file.path(root, "paper", "evaluation", "R", "metrics.R"))
source(file.path(root, "paper", "evaluation", "R", "aggregation.R"))
source(file.path(root, "paper", "evaluation", "R", "provenance.R"))
source(file.path(root, "paper", "evaluation", "R", "tasks.R"))

output_dir <- file.path(root, "paper", "evaluation", "output", "smoke")
resume <- !"--fresh" %in% commandArgs(trailingOnly = TRUE)
manifest <- run_e1_smoke(output_dir, resume = resume)
print(manifest[c("task_id", "status", "execution_seconds", "fit_seconds",
                 "estimated_peak_memory_mb", "result_file_mb")],
      row.names = FALSE)
cat("\nResource summary:\n")
print(summarize_e1_resources(manifest), row.names = FALSE)
