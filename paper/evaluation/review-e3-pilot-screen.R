# Convert the completed E3 screen measurements into a reproducible execution
# gate for the factorial pilot. This reviews resources and provenance; it does
# not select a scientifically preferred SVI configuration.

root <- normalizePath(".", mustWork = TRUE)
for (file in c("scenarios.R", "metrics.R", "aggregation.R", "provenance.R",
               "tasks.R", "e3.R", "e3-execution.R"))
  source(file.path(root, "paper", "evaluation", "R", file))

screen_dir <- file.path(
  root, "paper", "evaluation", "output", "pilot", "E3", "screen")
review <- write_e3_screen_review(screen_dir)
print(review$summary, row.names = FALSE)
print(review$gate, row.names = FALSE)
if (!isTRUE(review$summary$ready_for_factorial))
  stop("E3 screen review did not pass; factorial remains blocked")
cat("E3 screen resource/provenance review passed. Statistical conclusions\n",
    "still require the factorial pilot and confirmatory runs.\n", sep = "")

