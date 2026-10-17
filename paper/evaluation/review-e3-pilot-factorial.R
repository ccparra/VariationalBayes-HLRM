# Aggregate the completed 4 x 3 x 3 E3 factorial across its three pilot seeds.
# Ranking is exploratory and cannot by itself authorize confirmatory execution.
root <- normalizePath(".", mustWork = TRUE)
for (file in c("scenarios.R", "metrics.R", "aggregation.R", "provenance.R",
               "tasks.R", "e3.R", "e3-execution.R"))
  source(file.path(root, "paper", "evaluation", "R", file))

factorial_dir <- file.path(
  root, "paper", "evaluation", "output", "pilot", "E3", "factorial")
review <- write_e3_factorial_review(factorial_dir)
print(review$review, row.names = FALSE)
cat("\nTop exploratory cells by mean MLPD, then stability and clustering:\n")
print(utils::head(review$ranking[c(
  "pilot_rank", "minibatch_size", "exponent", "offset", "mlpd_mean",
  "mlpd_sd", "ari_mean", "ari_sd", "zero_ari_rate")], 10L),
  row.names = FALSE)
if (!isTRUE(review$review$execution_valid))
  stop("E3 factorial execution review failed")
cat("E3 factorial execution is complete and internally valid. It remains\n",
    "pilot evidence; confirmatory execution is still blocked.\n", sep = "")
