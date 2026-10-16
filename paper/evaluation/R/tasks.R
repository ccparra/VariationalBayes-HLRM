# Declarative task expansion and lightweight execution for the paper suite.
# Confirmatory settings are never inferred from package defaults.

evaluation_profile <- function(profile = c("smoke", "pilot", "full")) {
  profile <- match.arg(profile)
  catalog <- e1_scenario_catalog()
  switch(profile,
    smoke = list(
      name = "smoke", scenarios = "S1", replicates = 1:2,
      methods = c("MCMC", "CAVI", "SVI"), confirmatory = FALSE),
    pilot = list(
      name = "pilot", scenarios = catalog$scenario, replicates = 1:5,
      methods = c("MCMC", "CAVI", "SVI"), confirmatory = FALSE),
    full = list(
      name = "full", scenarios = catalog$scenario, replicates = 1:200,
      methods = c("MCMC", "CAVI", "SVI"), confirmatory = TRUE))
}

expand_e1_tasks <- function(profile = c("smoke", "pilot", "full"),
                            base_seed = 210000L) {
  specification <- evaluation_profile(match.arg(profile))
  grid <- expand.grid(
    scenario = specification$scenarios,
    replicate = specification$replicates,
    method = specification$methods,
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  grid <- grid[order(match(grid$scenario, e1_scenario_catalog()$scenario),
                     grid$replicate,
                     match(grid$method, specification$methods)), ]
  rownames(grid) <- NULL
  dataset_seed <- mapply(
    .e1_task_seed, grid$scenario, grid$replicate,
    MoreArgs = list(base_seed = base_seed))
  method_offset <- match(grid$method, specification$methods) * 10000L
  fit_seed <- dataset_seed + method_offset
  if (any(fit_seed > .Machine$integer.max))
    stop("Task seeds exceed the valid R seed range")
  data.frame(
    schema_version = "paper-task-manifest-1",
    profile = specification$name,
    experiment = "E1",
    scenario = grid$scenario,
    replicate = as.integer(grid$replicate),
    method = grid$method,
    dataset_seed = as.integer(dataset_seed),
    fit_seed = as.integer(fit_seed),
    task_id = sprintf("E1-%s-r%03d-%s", grid$scenario, grid$replicate,
                      tolower(grid$method)),
    confirmatory = specification$confirmatory,
    stringsAsFactors = FALSE)
}

.smoke_control <- function(method, group_count) {
  switch(method,
    MCMC = mcmc_control(
      warmup = 2, post_warmup_draws = 8, chains = 2,
      evaluation_draws = 4,
      stopping = mcmc_stopping(
        min_post_warmup_draws = 4, check_every = 4,
        rhat_max = 100, ess_bulk_min = .1, ess_tail_min = .1,
        relative_mcse_mean_max = 100, relative_mcse_sd_max = 100),
      max_memory_mb = 256),
    CAVI = cavi_control(
      max_iterations = 4, evaluation_draws = 4,
      stopping = cavi_stopping(
        objective_tolerance = 1e20, parameter_tolerance = 1e20,
        min_iterations = 2, patience = 1),
      starts = 2, require_converged = FALSE, max_memory_mb = 256),
    SVI = svi_control(
      max_iterations = 4, minibatch_size = min(10L, group_count),
      learning_rate = power_schedule(offset = 19, exponent = .7),
      evaluation_draws = 4,
      stopping = svi_stopping(
        objective_tolerance = 1e20, parameter_tolerance = 1e20,
        min_iterations = 2, patience = 1, check_every = 1,
        smoothing_window = 1),
      starts = 2, require_converged = FALSE, max_memory_mb = 256),
    stop("Unknown inference method"))
}

fit_e1_task <- function(task) {
  if (!is.data.frame(task) || nrow(task) != 1L ||
      !all(c("profile", "scenario", "replicate", "method", "dataset_seed",
             "fit_seed", "task_id") %in% names(task)))
    stop("task must be one row from expand_e1_tasks()")
  if (!identical(task$profile, "smoke"))
    stop("Only smoke tasks may run before pilot controls and budgets are frozen")
  scenario <- scenario_e1(
    task$scenario, task$replicate,
    base_seed = task$dataset_seed -
      100000L * match(task$scenario, e1_scenario_catalog()$scenario) -
      task$replicate)
  if (!identical(scenario$seed, task$dataset_seed))
    stop("Task and scenario dataset seeds are inconsistent")
  simulation <- simulate_data(
    scenario$generator, scenario$design, scenario$sample,
    seed = task$dataset_seed)
  fit <- fit_model(
    scenario$model, simulation$train, scenario$prior,
    method = tolower(task$method),
    initialization = initialization_from_prior(),
    control = .smoke_control(task$method, scenario$configuration$group_count),
    seed = task$fit_seed, verbose = FALSE)
  summary <- summarize_posterior(fit)
  parameter_table <- tidy_parameters(summary)
  in_sample <- evaluate_in_sample_fit(fit)
  predictive <- predictive_criteria(
    fit, n_draws = 4, seed = task$fit_seed + 1L,
    return_log_likelihood = TRUE)
  ppc <- posterior_predictive_check(
    fit, n_draws = 4, seed = task$fit_seed + 2L,
    return_yrep = FALSE)
  held_out <- predictive_performance(
    fit, simulation$test, predictive_target = "observed_groups",
    predictive_unit = "observation",
    group_ids = seq_along(simulation$test$y), n_draws = 4,
    seed = task$fit_seed + 3L)
  classification <- classify_groups(fit)
  co_clustering_matrix <- co_clustering(fit)
  probabilities <- responsibilities(fit)
  clustering_metrics <- list(
    adjusted_rand_index = adjusted_rand_index(
      simulation$truth$allocations, classification$component),
    co_clustering_error = co_clustering_error(
      simulation$truth$allocations, co_clustering_matrix),
    normalized_entropy = normalized_classification_entropy(probabilities))
  recovery <- aligned_chlrm_recovery(parameter_table, simulation$truth)
  structure(list(
    schema_version = "paper-task-result-4",
    task = task,
    configuration = scenario$configuration,
    truth = simulation$truth,
    fit = list(
      model = fit$model, method = fit$method, status = fit$status,
      diagnostics = fit$diagnostics, timing = fit$timing,
      metadata = fit$metadata),
    parameters = parameter_table,
    recovery = recovery,
    allocations = list(
      classification = classification,
      co_clustering = co_clustering_matrix,
      metrics = clustering_metrics),
    evaluation = list(
      in_sample = in_sample, predictive = predictive, ppc = ppc,
      held_out = held_out)),
    class = c("paper_task_result", "list"))
}

.atomic_save_rds <- function(object, path) {
  directory <- dirname(path)
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(".paper-task-", tmpdir = directory, fileext = ".rds")
  on.exit(if (file.exists(temporary)) unlink(temporary), add = TRUE)
  saveRDS(object, temporary, version = 3)
  backup <- NULL
  if (file.exists(path)) {
    backup <- tempfile(".paper-task-backup-", tmpdir = directory,
                       fileext = ".rds")
    if (!file.rename(path, backup))
      stop("Could not stage the previous task result for replacement: ", path)
    on.exit({
      if (!file.exists(path) && file.exists(backup)) file.rename(backup, path)
      if (file.exists(path) && file.exists(backup)) unlink(backup)
    }, add = TRUE)
  }
  if (!file.rename(temporary, path))
    stop("Could not atomically publish task result: ", path)
  invisible(normalizePath(path, mustWork = TRUE))
}

.paper_task_record <- function(result, path, execution_status,
                               execution_seconds) {
  estimate <- result$fit$metadata$resource_estimate
  convergence <- result$fit$diagnostics
  data.frame(
    task_id = result$task$task_id,
    status = execution_status,
    execution_seconds = as.numeric(execution_seconds),
    fit_seconds = as.numeric(result$fit$timing$total_seconds),
    convergence_status = convergence$status,
    estimated_peak_memory_mb = if (is.null(estimate)) NA_real_ else
      estimate$estimated_peak_memory_mb,
    memory_budget_mb = if (is.null(estimate)) NA_real_ else
      estimate$max_memory_mb,
    estimated_evaluation_storage_mb = if (is.null(estimate)) NA_real_ else
      estimate$estimated_evaluation_storage_mb,
    result_file_mb = file.info(path)$size / 1024^2,
    path = file.path("tasks", basename(path)),
    stringsAsFactors = FALSE)
}

summarize_e1_resources <- function(manifest) {
  required <- c(
    "method", "fit_seconds", "estimated_peak_memory_mb",
    "estimated_evaluation_storage_mb", "result_file_mb")
  if (!is.data.frame(manifest) || !all(required %in% names(manifest)))
    stop("manifest does not contain the resource columns produced by run_e1_smoke()")
  groups <- split(manifest, manifest$method)
  rows <- lapply(groups, function(group) data.frame(
    method = group$method[1L],
    tasks = nrow(group),
    mean_fit_seconds = mean(group$fit_seconds),
    max_fit_seconds = max(group$fit_seconds),
    max_estimated_peak_memory_mb = max(group$estimated_peak_memory_mb),
    max_estimated_evaluation_storage_mb =
      max(group$estimated_evaluation_storage_mb),
    total_result_file_mb = sum(group$result_file_mb),
    stringsAsFactors = FALSE))
  summary <- do.call(rbind, rows)
  rownames(summary) <- NULL
  summary[order(match(summary$method, c("MCMC", "CAVI", "SVI"))), ]
}

run_e1_smoke <- function(output_dir, resume = TRUE) {
  if (!is.character(output_dir) || length(output_dir) != 1L ||
      !nzchar(output_dir)) stop("output_dir must be one nonempty path")
  if (!is.logical(resume) || length(resume) != 1L || is.na(resume))
    stop("resume must be TRUE or FALSE")
  tasks <- expand_e1_tasks("smoke")
  provenance <- paper_run_provenance(
    "E1", "smoke",
    list(profile = evaluation_profile("smoke"),
         scenario_catalog = e1_scenario_catalog(), tasks = tasks))
  task_dir <- file.path(output_dir, "tasks")
  dir.create(task_dir, recursive = TRUE, showWarnings = FALSE)
  records <- vector("list", nrow(tasks))
  for (index in seq_len(nrow(tasks))) {
    task <- tasks[index, , drop = FALSE]
    path <- file.path(task_dir, paste0(task$task_id, ".rds"))
    started <- proc.time()[["elapsed"]]
    status <- "complete"
    if (resume && file.exists(path)) {
      result <- readRDS(path)
      valid <- .paper_result_compatible(
        result, "paper_task_result", "paper-task-result-4",
        task$task_id, provenance)
      if (valid) {
        status <- "reused"
      } else {
        result <- fit_e1_task(task)
        result$provenance <- provenance
        .atomic_save_rds(result, path)
        status <- "regenerated_incompatible"
      }
    } else {
      result <- fit_e1_task(task)
      result$provenance <- provenance
      .atomic_save_rds(result, path)
    }
    records[[index]] <- .paper_task_record(
      result, path, status, proc.time()[["elapsed"]] - started)
  }
  manifest <- merge(tasks, do.call(rbind, records), by = "task_id",
                    sort = FALSE)
  manifest <- manifest[match(tasks$task_id, manifest$task_id), ]
  manifest <- .paper_annotate_manifest(manifest, provenance)
  write.csv(manifest, file.path(output_dir, "manifest.csv"), row.names = FALSE)
  write_paper_provenance(provenance, output_dir)
  write.csv(summarize_e1_resources(manifest),
            file.path(output_dir, "resource-summary.csv"), row.names = FALSE)
  write_e1_metric_outputs(manifest, output_dir)
  manifest
}
