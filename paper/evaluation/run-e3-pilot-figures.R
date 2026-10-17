# Generate reporting-only figures from the completed E3 screen artifacts. This
# script never fits a model and labels every product as preliminary.
root <- normalizePath(".", mustWork = TRUE)
source(file.path(root, "paper", "evaluation", "R", "provenance.R"))
source(file.path(root, "paper", "evaluation", "R", "reporting.R"))
input_dir <- file.path(
  root, "paper", "evaluation", "output", "pilot", "E3", "screen")
manifest <- generate_e3_screen_figures(input_dir)
print(manifest, row.names = FALSE)
