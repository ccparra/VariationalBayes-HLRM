# Build reporting-only smoke figures from completed CSV artifacts. This script
# never loads the inference package and cannot fit a model.
root <- normalizePath(".", mustWork = TRUE)
source(file.path(root, "paper", "evaluation", "R", "provenance.R"))
source(file.path(root, "paper", "evaluation", "R", "reporting.R"))
input_dir <- file.path(root, "paper", "evaluation", "output", "smoke")
manifest <- generate_smoke_figures(input_dir)
print(manifest, row.names = FALSE)
