# Project-boundary helpers for the research repository. During the repository
# split the evaluation suite may live either beside the package source or in
# the historical monorepo. Callers can always override discovery explicitly
# with VB_HLRM_PACKAGE_SOURCE.

paper_package_source <- function(
    research_root = Sys.getenv("VB_HLRM_RESEARCH_ROOT", unset = "."),
    explicit = Sys.getenv("VB_HLRM_PACKAGE_SOURCE", unset = "")) {
  research_root <- normalizePath(research_root, mustWork = TRUE)
  candidates <- unique(c(
    explicit,
    research_root,
    file.path(dirname(research_root), "VariationalBayes-HLRM")))
  candidates <- candidates[nzchar(candidates)]
  valid <- vapply(candidates, function(path) {
    description <- file.path(path, "DESCRIPTION")
    if (!file.exists(description)) return(FALSE)
    package <- tryCatch(
      unname(read.dcf(description, fields = "Package")[1L, 1L]),
      error = function(condition) NA_character_)
    identical(package, "VariationalBayesHLRM")
  }, logical(1L))
  if (!any(valid)) stop(paste(
    "Cannot locate the VariationalBayesHLRM package source.",
    "Set VB_HLRM_PACKAGE_SOURCE to its checkout."))
  normalizePath(candidates[which(valid)[1L]], mustWork = TRUE)
}

install_paper_package <- function(research_root, library_path) {
  package_source <- paper_package_source(research_root)
  dir.create(library_path, recursive = TRUE, showWarnings = FALSE)
  install_log <- system2(
    file.path(R.home("bin"), "R"),
    c("CMD", "INSTALL", "--no-multiarch",
      paste0("--library=", shQuote(library_path)), shQuote(package_source)),
    stdout = TRUE, stderr = TRUE)
  install_status <- attr(install_log, "status")
  if (!is.null(install_status) && install_status != 0L)
    stop("Temporary package installation failed:\n",
         paste(install_log, collapse = "\n"))
  .libPaths(c(library_path, .libPaths()))
  Sys.setenv(
    VB_HLRM_RESEARCH_ROOT = normalizePath(research_root, mustWork = TRUE),
    VB_HLRM_PACKAGE_SOURCE = package_source)
  invisible(package_source)
}
