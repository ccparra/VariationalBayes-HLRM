# Freeze a reviewed E3 pilot candidate. Factorial freezing additionally
# requires an accepted, passing screen review. This records no scientific
# result and refuses dirty/uncommitted source trees.
.cli_value <- function(arguments, prefix, required = TRUE) {
  values <- arguments[startsWith(arguments, prefix)]
  if (!length(values)) {
    if (required) stop("Missing required argument ", prefix)
    return(NULL)
  }
  sub(prefix, "", values[[1L]], fixed = TRUE)
}

arguments <- commandArgs(trailingOnly = TRUE)
accepted <- .cli_value(arguments, "--accept-configuration-fingerprint=")
stage <- match.arg(.cli_value(arguments, "--stage="), c("screen", "factorial"))
accepted_screen <- .cli_value(
  arguments, "--accept-screen-output-fingerprint=", required = FALSE)
root <- normalizePath(".", mustWork = TRUE)
for (file in c("scenarios.R", "metrics.R", "aggregation.R", "provenance.R",
               "tasks.R", "e3.R", "e3-execution.R"))
  source(file.path(root, "paper", "evaluation", "R", file))
preflight_dir <- file.path(
  root, "paper", "evaluation", "output", "preflight", "E3",
  paste0("pilot-", stage))
if (identical(stage, "factorial")) {
  screen_dir <- file.path(
    root, "paper", "evaluation", "output", "pilot", "E3", "screen")
  screen_review <- write_e3_screen_review(screen_dir)
} else {
  screen_review <- NULL
}
frozen <- freeze_e3_configuration(
  preflight_dir, accepted, screen_review = screen_review,
  accepted_screen_output_fingerprint = accepted_screen)
cat("Frozen E3 ", stage, " configuration\n", sep = "")
cat("configuration_fingerprint=", frozen$configuration_fingerprint, "\n", sep = "")
cat("source_fingerprint=", frozen$source_fingerprint, "\n", sep = "")
if (!is.null(frozen$upstream_evidence))
  cat("screen_output_fingerprint=",
      frozen$upstream_evidence$screen_output_fingerprint, "\n", sep = "")
cat("git_commit=", frozen$git_commit, "\n", sep = "")
