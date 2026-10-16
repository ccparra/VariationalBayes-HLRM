# Run the bounded E3 integration profile from the repository root.
root <- normalizePath(".", mustWork = TRUE)
validation_library <- tempfile("vb-hlrm-paper-e3-library-")
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
for (file in c("scenarios.R", "metrics.R", "aggregation.R", "provenance.R",
               "tasks.R", "e3.R"))
  source(file.path(root, "paper", "evaluation", "R", file))
output_dir <- file.path(root, "paper", "evaluation", "output", "smoke", "E3")
resume <- !"--fresh" %in% commandArgs(trailingOnly = TRUE)
result <- run_e3_smoke(output_dir, resume)
print(result$fits, row.names = FALSE)
cat("\nSVI relative to CAVI:\n")
print(result$comparison, row.names = FALSE)
