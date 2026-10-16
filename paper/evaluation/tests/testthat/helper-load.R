root <- normalizePath(Sys.getenv("VB_HLRM_ROOT"), mustWork = TRUE)

# The default runner installs the package in a temporary library. Developers
# may still run individual tests against the source tree, but source loading
# exposes only exported symbols.
if (identical(Sys.getenv("VB_HLRM_USE_INSTALLED"), "true")) {
  suppressPackageStartupMessages(library(
    "VariationalBayesHLRM", character.only = TRUE))
} else {
  pkgload::load_all(
    root, quiet = TRUE, export_all = FALSE, helpers = FALSE,
    attach_testthat = FALSE)
}
source(file.path(root, "paper", "evaluation", "R", "sbc.R"))
source(file.path(root, "paper", "evaluation", "R", "metrics.R"))
source(file.path(root, "paper", "evaluation", "R", "scenarios.R"))
source(file.path(root, "paper", "evaluation", "R", "aggregation.R"))
source(file.path(root, "paper", "evaluation", "R", "provenance.R"))
source(file.path(root, "paper", "evaluation", "R", "tasks.R"))
source(file.path(root, "paper", "evaluation", "R", "e2.R"))
source(file.path(root, "paper", "evaluation", "R", "e3.R"))
source(file.path(root, "paper", "evaluation", "R", "e4.R"))
source(file.path(root, "paper", "evaluation", "R", "e5.R"))

.test_rng_state <- function() {
  if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
    get(".Random.seed", envir = .GlobalEnv, inherits = FALSE) else NULL
}
