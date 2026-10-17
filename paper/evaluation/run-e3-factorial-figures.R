# Generate pilot-factorial figures from completed E3 CSV artifacts. This
# script never fits a model and marks every output as non-confirmatory.
root <- normalizePath(".", mustWork = TRUE)
source(file.path(root, "paper", "evaluation", "R", "provenance.R"))
source(file.path(root, "paper", "evaluation", "R", "reporting.R"))
input_dir <- file.path(
  root, "paper", "evaluation", "output", "pilot", "E3", "factorial")
manifest <- generate_e3_factorial_figures(input_dir)
print(manifest, row.names = FALSE)
