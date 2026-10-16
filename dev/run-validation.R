# Run from the repository root. Generated files live in tests/output/, which is
# ignored by git but remains available for local inspection.
library(testthat)

root <- normalizePath(".", mustWork = TRUE)
Sys.setenv(VB_HLRM_ROOT = root)
output_dir <- file.path(root, "tests", "output")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

source(file.path(root, "dev", "check-consumer-boundary.R"))
boundary <- check_consumer_boundary(root)

started <- Sys.time()
results <- test_dir(file.path(root, "tests", "testthat"), reporter = "summary")
finished <- Sys.time()
details <- as.data.frame(results)
csv_details <- details
list_columns <- vapply(csv_details, is.list, logical(1L))
csv_details[list_columns] <- lapply(csv_details[list_columns], function(column)
  vapply(column, function(value) {
    if (!length(value)) return("")
    paste(vapply(value, function(item) {
      message <- tryCatch(conditionMessage(item), error = function(e) "")
      paste(c(class(item)[1L], message), collapse = ": ")
    }, character(1L)), collapse = " | ")
  }, character(1L)))
utils::write.csv(csv_details, file.path(output_dir, "results.csv"), row.names = FALSE)

count <- function(name) if (name %in% names(details)) sum(details[[name]], na.rm = TRUE) else NA
summary <- c(
  paste("Started:", format(started, tz = "UTC", usetz = TRUE)),
  paste("Finished:", format(finished, tz = "UTC", usetz = TRUE)),
  sprintf("Elapsed seconds: %.3f", as.numeric(difftime(finished, started, units = "secs"))),
  paste("Consumer R files checked:", boundary$files_checked),
  paste("Consumer boundary violations:", boundary$violations),
  paste("Test blocks:", nrow(details)),
  paste("Failed expectations:", count("failed")),
  paste("Errors:", count("error")),
  paste("Warnings:", count("warning")),
  paste("Skipped:", count("skipped"))
)
writeLines(summary, file.path(output_dir, "latest-summary.txt"))
writeLines(capture.output(utils::sessionInfo()), file.path(output_dir, "session-info.txt"))

cat(paste(summary, collapse = "\n"), "\n")
