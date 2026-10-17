# Run the bounded E2 integration profile from the repository root.
root <- normalizePath(".", mustWork = TRUE)
source(file.path(root, "paper", "evaluation", "R", "project.R"))
validation_library <- tempfile("vb-hlrm-paper-e2-library-")
on.exit(unlink(validation_library, recursive = TRUE, force = TRUE), add = TRUE)
install_paper_package(root, validation_library)
suppressPackageStartupMessages(library(VariationalBayesHLRM))
for (file in c("scenarios.R", "metrics.R", "aggregation.R", "provenance.R",
               "tasks.R", "e2.R"))
  source(file.path(root, "paper", "evaluation", "R", file))

output_dir <- file.path(root, "paper", "evaluation", "output", "smoke", "E2")
resume <- !"--fresh" %in% commandArgs(trailingOnly = TRUE)
result <- run_e2_smoke(output_dir, resume)
print(result$fits, row.names = FALSE)
cat("\nSelections:\n")
print(result$selection, row.names = FALSE)
