# Freeze a reviewed E3 screen candidate. This records no scientific result and
# refuses dirty/uncommitted source trees.
.cli_value <- function(arguments, prefix) {
  values <- arguments[startsWith(arguments, prefix)]
  if (!length(values)) stop("Missing required argument ", prefix)
  sub(prefix, "", values[[1L]], fixed = TRUE)
}

arguments <- commandArgs(trailingOnly = TRUE)
accepted <- .cli_value(arguments, "--accept-configuration-fingerprint=")
root <- normalizePath(".", mustWork = TRUE)
for (file in c("scenarios.R", "metrics.R", "aggregation.R", "provenance.R",
               "tasks.R", "e3.R"))
  source(file.path(root, "paper", "evaluation", "R", file))
preflight_dir <- file.path(
  root, "paper", "evaluation", "output", "preflight", "E3", "pilot-screen")
frozen <- freeze_e3_configuration(preflight_dir, accepted)
cat("Frozen E3 screen configuration\n")
cat("configuration_fingerprint=", frozen$configuration_fingerprint, "\n", sep = "")
cat("source_fingerprint=", frozen$source_fingerprint, "\n", sep = "")
cat("git_commit=", frozen$git_commit, "\n", sep = "")
