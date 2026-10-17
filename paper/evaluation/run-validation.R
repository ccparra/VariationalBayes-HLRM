# Lightweight infrastructure checks only; this does not execute confirmatory
# paper simulations. Run from the repository root.
library(testthat)

root <- normalizePath(".", mustWork = TRUE)
Sys.setenv(VB_HLRM_RESEARCH_ROOT = root, VB_HLRM_ROOT = root)
source(file.path(root, "paper", "evaluation", "R", "project.R"))
output_dir <- file.path(root, "paper", "evaluation", "output")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Validate the same boundary used by a paper replication: install the package
# into an isolated library and expose no source-tree internals to the suite.
validation_library <- tempfile("vb-hlrm-paper-library-")
on.exit(unlink(validation_library, recursive = TRUE, force = TRUE), add = TRUE)
install_paper_package(root, validation_library)
Sys.setenv(VB_HLRM_USE_INSTALLED = "true")

started <- Sys.time()
results <- test_dir(
  file.path(root, "paper", "evaluation", "tests", "testthat"),
  reporter = "summary")
finished <- Sys.time()
details <- as.data.frame(results)
utils::write.csv(details[c("file", "context", "test", "nb", "failed",
                           "skipped", "error", "warning", "real")],
                 file.path(output_dir, "results.csv"), row.names = FALSE)
summary <- c(
  "Package boundary: temporary installation",
  paste("Started:", format(started, tz = "UTC", usetz = TRUE)),
  paste("Finished:", format(finished, tz = "UTC", usetz = TRUE)),
  sprintf("Elapsed seconds: %.3f",
          as.numeric(difftime(finished, started, units = "secs"))),
  paste("Test blocks:", nrow(details)),
  paste("Failed expectations:", sum(details$failed, na.rm = TRUE)),
  paste("Errors:", sum(details$error, na.rm = TRUE)),
  paste("Warnings:", sum(details$warning, na.rm = TRUE)),
  paste("Skipped:", sum(details$skipped, na.rm = TRUE)))
writeLines(summary, file.path(output_dir, "latest-summary.txt"))
writeLines(capture.output(utils::sessionInfo()),
           file.path(output_dir, "session-info.txt"))
cat(paste(summary, collapse = "\n"), "\n")
