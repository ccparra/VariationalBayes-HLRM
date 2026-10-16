# Reproducibility metadata for the paper suite. Fingerprints deliberately cover
# both the declared configuration and the exact package/suite source files, so
# a cached task is not reused after an uncommitted implementation change.

.paper_canonicalize <- function(value) {
  if (is.environment(value) || is.function(value) || typeof(value) == "externalptr")
    stop("configuration cannot contain environments, functions, or external pointers")
  if (is.data.frame(value)) {
    columns <- lapply(value, .paper_canonicalize)
    names(columns) <- names(value)
    return(list(type = "data.frame", rows = nrow(value), columns = columns))
  }
  if (is.factor(value)) value <- as.character(value)
  if (is.list(value)) {
    value <- lapply(value, .paper_canonicalize)
    if (!is.null(names(value)) && all(nzchar(names(value))))
      value <- value[order(names(value))]
    return(value)
  }
  attributes(value) <- NULL
  value
}

.paper_fingerprint <- function(value) {
  payload <- serialize(.paper_canonicalize(value), NULL, version = 3)
  path <- tempfile("paper-fingerprint-", fileext = ".bin")
  on.exit(unlink(path), add = TRUE)
  connection <- file(path, open = "wb")
  writeBin(payload, connection)
  close(connection)
  paste0("md5:", unname(tools::md5sum(path)))
}

.paper_source_fingerprint <- function(root, experiment, scope = c("task", "suite")) {
  scope <- match.arg(scope)
  root <- normalizePath(root, mustWork = TRUE)
  package_files <- c(
    file.path(root, c("DESCRIPTION", "NAMESPACE")),
    list.files(file.path(root, "R"), pattern = "\\.[Rr]$",
               recursive = TRUE, full.names = TRUE))
  evaluation_root <- file.path(root, "paper", "evaluation")
  if (scope == "suite") {
    evaluation_files <- c(
      list.files(file.path(evaluation_root, "R"), pattern = "\\.[Rr]$",
                 full.names = TRUE),
      list.files(evaluation_root, pattern = "\\.[Rr]$", full.names = TRUE))
  } else {
    common <- file.path(evaluation_root, "R", c(
      "scenarios.R", "metrics.R", "provenance.R", "tasks.R"))
    experiment_file <- if (experiment %in% paste0("E", 2:5))
      file.path(evaluation_root, "R", paste0(tolower(experiment), ".R")) else
      character()
    task_worker <- if (identical(experiment, "E3"))
      file.path(evaluation_root, "run-e3-screen-task.R") else character()
    evaluation_files <- c(common, experiment_file, task_worker)
  }
  files <- c(package_files, evaluation_files)
  files <- sort(unique(files[file.exists(files)]))
  relative <- substring(files, nchar(root) + 2L)
  hashes <- unname(tools::md5sum(files))
  .paper_fingerprint(as.list(stats::setNames(hashes, relative)))
}

.paper_git_value <- function(root, arguments) {
  value <- tryCatch(
    suppressWarnings(system2("git", c("-C", root, arguments),
                             stdout = TRUE, stderr = FALSE)),
    error = function(condition) character())
  status <- attr(value, "status")
  if (!length(value) || (!is.null(status) && status != 0L)) NA_character_ else
    paste(value, collapse = "\n")
}

paper_run_provenance <- function(experiment, profile, configuration,
                                 root = Sys.getenv("VB_HLRM_ROOT", unset = ".")) {
  if (!is.character(experiment) || length(experiment) != 1L || !nzchar(experiment))
    stop("experiment must be one nonempty string")
  if (!is.character(profile) || length(profile) != 1L || !nzchar(profile))
    stop("profile must be one nonempty string")
  root <- normalizePath(root, mustWork = TRUE)
  configuration_fingerprint <- .paper_fingerprint(configuration)
  source_fingerprint <- .paper_source_fingerprint(
    root, experiment, scope = "task")
  suite_source_fingerprint <- .paper_source_fingerprint(
    root, experiment, scope = "suite")
  package_version <- tryCatch(
    as.character(utils::packageVersion("VariationalBayesHLRM")),
    error = function(condition) {
      description <- read.dcf(file.path(root, "DESCRIPTION"))
      unname(description[1L, "Version"])
    })
  git_commit <- .paper_git_value(root, c("rev-parse", "HEAD"))
  git_status <- .paper_git_value(root, c("status", "--porcelain"))
  git_dirty <- if (is.na(git_status)) NA else nzchar(git_status)
  run_fingerprint <- .paper_fingerprint(list(
    experiment = experiment, profile = profile,
    configuration = configuration_fingerprint,
    source = source_fingerprint, package_version = package_version))
  output_fingerprint <- .paper_fingerprint(list(
    run = run_fingerprint, suite_source = suite_source_fingerprint))
  structure(list(
    schema_version = "paper-run-provenance-1",
    experiment = experiment, profile = profile,
    configuration_fingerprint = configuration_fingerprint,
    source_fingerprint = source_fingerprint,
    suite_source_fingerprint = suite_source_fingerprint,
    run_fingerprint = run_fingerprint,
    output_fingerprint = output_fingerprint,
    package_version = package_version,
    git_commit = git_commit, git_dirty = git_dirty,
    r_version = R.version.string, platform = R.version$platform,
    generated_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    threads = list(
      omp = Sys.getenv("OMP_NUM_THREADS", unset = ""),
      openblas = Sys.getenv("OPENBLAS_NUM_THREADS", unset = ""),
      mkl = Sys.getenv("MKL_NUM_THREADS", unset = ""))),
    class = c("paper_run_provenance", "list"))
}

.paper_result_compatible <- function(result, result_class, schema_version,
                                     task_id, provenance) {
  inherits(result, result_class) &&
    identical(result$schema_version, schema_version) &&
    identical(result$task$task_id, task_id) &&
    inherits(result$provenance, "paper_run_provenance") &&
    identical(result$provenance$run_fingerprint,
              provenance$run_fingerprint)
}

.paper_annotate_manifest <- function(manifest, provenance) {
  manifest$configuration_fingerprint <- provenance$configuration_fingerprint
  manifest$source_fingerprint <- provenance$source_fingerprint
  manifest$suite_source_fingerprint <- provenance$suite_source_fingerprint
  manifest$run_fingerprint <- provenance$run_fingerprint
  manifest$output_fingerprint <- provenance$output_fingerprint
  manifest$package_version <- provenance$package_version
  manifest$git_commit <- provenance$git_commit
  manifest$git_dirty <- provenance$git_dirty
  manifest
}

write_paper_provenance <- function(provenance, output_dir) {
  if (!inherits(provenance, "paper_run_provenance"))
    stop("provenance must be produced by paper_run_provenance()")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  row <- data.frame(
    schema_version = provenance$schema_version,
    experiment = provenance$experiment, profile = provenance$profile,
    configuration_fingerprint = provenance$configuration_fingerprint,
    source_fingerprint = provenance$source_fingerprint,
    suite_source_fingerprint = provenance$suite_source_fingerprint,
    run_fingerprint = provenance$run_fingerprint,
    output_fingerprint = provenance$output_fingerprint,
    package_version = provenance$package_version,
    git_commit = provenance$git_commit, git_dirty = provenance$git_dirty,
    r_version = provenance$r_version, platform = provenance$platform,
    generated_at_utc = provenance$generated_at_utc,
    omp_threads = provenance$threads$omp,
    openblas_threads = provenance$threads$openblas,
    mkl_threads = provenance$threads$mkl,
    stringsAsFactors = FALSE)
  utils::write.csv(row, file.path(output_dir, "provenance.csv"), row.names = FALSE)
  saveRDS(provenance, file.path(output_dir, "provenance.rds"), version = 3)
  invisible(row)
}
