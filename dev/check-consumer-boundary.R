# Static guard for the one-way dependency from scientific consumers to the
# installed public package. This check does not inspect thesis/v0, which is a
# deliberately frozen historical artifact.

check_consumer_boundary <- function(root = normalizePath(".", mustWork = TRUE)) {
  consumer_roots <- file.path(root, c("paper/evaluation", "thesis/final/code"))
  files <- unlist(lapply(consumer_roots, function(path)
    list.files(path, pattern = "\\.[Rr]$", recursive = TRUE,
               full.names = TRUE)), use.names = FALSE)
  forbidden <- c(
    "source\\([^\\n]*R/(LRM|HLRM|CHLRM)",
    "VariationalBayesHLRM:::",
    "\\bhyperparams_(lrm|hlrm|chlrm)\\b",
    "\\brun_mcmc_chains\\b",
    "\\brun_vi_multistart\\b",
    "\\b(mcmc|vi)_(lrm|hlrm|chlrm)\\$new\\b"
  )
  violations <- do.call(rbind, lapply(files, function(path) {
    lines <- readLines(path, warn = FALSE)
    hit <- unique(unlist(lapply(forbidden, function(pattern)
      grep(pattern, lines, perl = TRUE)), use.names = FALSE))
    if (!length(hit)) return(NULL)
    data.frame(
      file = substring(path, nchar(root) + 2L),
      line = hit,
      text = trimws(lines[hit]),
      stringsAsFactors = FALSE)
  }))
  if (!is.null(violations) && nrow(violations)) {
    message(paste(capture.output(print(violations, row.names = FALSE)),
                  collapse = "\n"))
    stop("Active consumers depend on private package implementation")
  }
  invisible(data.frame(files_checked = length(files), violations = 0L))
}

if (sys.nframe() == 0L) {
  result <- check_consumer_boundary()
  cat(sprintf("Consumer boundary: %d R files checked; no violations.\n",
              result$files_checked))
}
