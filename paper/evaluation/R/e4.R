# E4 paired scalability tasks over nested prefixes of one maximum dataset.

expand_e4_tasks <- function(profile = c("smoke", "pilot", "full"),
                            base_seed = 510000L) {
  profile <- .paper_profile(profile)
  protocol <- e4_protocol(profile)
  grid <- expand.grid(
    replicate = protocol$execution$replicates,
    group_count = protocol$execution$nested_group_counts,
    method = c("CAVI", "SVI"), KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE)
  if (profile == "smoke") grid <- rbind(
    data.frame(replicate = 1L,
               group_count = min(protocol$execution$nested_group_counts),
               method = "MCMC"), grid)
  grid <- grid[order(grid$replicate, grid$group_count,
                     match(grid$method, c("MCMC", "CAVI", "SVI"))), ]
  rownames(grid) <- NULL
  dataset_seed <- as.integer(base_seed + grid$replicate)
  fit_seed <- as.integer(dataset_seed + grid$group_count * 100L +
    match(grid$method, c("MCMC", "CAVI", "SVI")) * 10000L)
  data.frame(
    schema_version = "paper-task-manifest-1", profile = profile,
    experiment = "E4", replicate = as.integer(grid$replicate),
    group_count = as.integer(grid$group_count), method = grid$method,
    dataset_seed = dataset_seed, fit_seed = fit_seed,
    task_id = sprintf("E4-r%03d-m%04d-%s", grid$replicate,
                      grid$group_count, tolower(grid$method)),
    confirmatory = identical(profile, "full"), stringsAsFactors = FALSE)
}

.subset_grouped_simulation <- function(simulation, group_count) {
  keep <- seq_len(group_count)
  list(
    train = list(X = simulation$train$X[keep], y = simulation$train$y[keep]),
    test = list(X = simulation$test$X[keep], y = simulation$test$y[keep]),
    allocations = simulation$truth$allocations[keep])
}

.e4_smoke_control <- function(method, group_count) {
  if (method != "SVI") return(.smoke_control(method, group_count))
  svi_control(
    max_iterations = 4L, minibatch_size = min(100L, group_count),
    learning_rate = power_schedule(10, .7), evaluation_draws = 4L,
    stopping = svi_stopping(
      objective_tolerance = 1e20, parameter_tolerance = 1e20,
      min_iterations = 2L, patience = 1L, check_every = 1L,
      smoothing_window = 1L),
    starts = 2L, require_converged = FALSE, max_memory_mb = 256)
}

fit_e4_task <- function(task) {
  if (!is.data.frame(task) || nrow(task) != 1L ||
      !all(c("profile", "replicate", "group_count", "method",
             "dataset_seed", "fit_seed", "task_id") %in% names(task)))
    stop("task must be one row from expand_e4_tasks()")
  if (!identical(task$profile, "smoke"))
    stop("Only E4 smoke tasks may run before pilot settings are frozen")
  scenario <- scenario_e4(task$profile, task$replicate,
                          base_seed = task$dataset_seed - task$replicate)
  maximum <- simulate_data(
    scenario$generator, scenario$design, scenario$sample,
    seed = task$dataset_seed)
  data <- .subset_grouped_simulation(maximum, task$group_count)
  fit <- fit_model(
    scenario$model, data$train,
    review_chlrm_prior(scenario$model$coefficient_count,
                       scenario$model$component_count),
    method = tolower(task$method), initialization = initialization_from_prior(),
    control = .e4_smoke_control(task$method, task$group_count),
    seed = task$fit_seed, verbose = FALSE)
  held_out <- predictive_performance(
    fit, data$test, predictive_target = "observed_groups",
    predictive_unit = "observation", group_ids = seq_len(task$group_count),
    n_draws = 4, seed = task$fit_seed + 1L)
  classification <- classify_groups(fit)
  structure(list(
    schema_version = "paper-e4-task-result-2", task = task,
    observations = sum(lengths(data$train$y)),
    fit = list(diagnostics = fit$diagnostics, timing = fit$timing,
               metadata = fit$metadata), held_out = held_out,
    clustering = list(
      adjusted_rand_index = adjusted_rand_index(
        data$allocations, classification$component),
      co_clustering_error = co_clustering_error(
        data$allocations, co_clustering(fit)))),
    class = c("paper_e4_task_result", "list"))
}

.e4_result_row <- function(result, path, status) data.frame(
  task_id = result$task$task_id, replicate = result$task$replicate,
  group_count = result$task$group_count, observations = result$observations,
  method = result$task$method, runtime = result$fit$timing$total_seconds,
  mlpd = result$held_out$statistics$mean_log_predictive_density,
  ari = result$clustering$adjusted_rand_index,
  cpe = result$clustering$co_clustering_error,
  convergence_status = result$fit$diagnostics$status,
  estimated_peak_memory_mb =
    result$fit$metadata$resource_estimate$estimated_peak_memory_mb,
  result_file_mb = file.info(path)$size / 1024^2,
  status = status, path = file.path("tasks", basename(path)),
  stringsAsFactors = FALSE)

e4_scaling_slopes <- function(results) {
  groups <- split(results[results$method %in% c("CAVI", "SVI"), ],
                  results$method[results$method %in% c("CAVI", "SVI")])
  rows <- lapply(groups, function(values) {
    if (length(unique(values$observations)) < 2L)
      return(data.frame(method = values$method[1L], slope = NA_real_,
                        sizes = nrow(values)))
    model <- stats::lm(log(runtime) ~ log(observations), data = values)
    data.frame(method = values$method[1L],
               slope = unname(stats::coef(model)[2L]), sizes = nrow(values))
  })
  do.call(rbind, rows)
}

run_e4_smoke <- function(output_dir, resume = TRUE) {
  tasks <- expand_e4_tasks("smoke")
  provenance <- paper_run_provenance(
    "E4", "smoke", list(protocol = e4_protocol("smoke"), tasks = tasks))
  task_dir <- file.path(output_dir, "tasks")
  dir.create(task_dir, recursive = TRUE, showWarnings = FALSE)
  rows <- vector("list", nrow(tasks))
  for (i in seq_len(nrow(tasks))) {
    task <- tasks[i, , drop = FALSE]
    path <- file.path(task_dir, paste0(task$task_id, ".rds"))
    reused <- FALSE
    if (resume && file.exists(path)) {
      result <- readRDS(path)
      reused <- .paper_result_compatible(
        result, "paper_e4_task_result", "paper-e4-task-result-2",
        task$task_id, provenance)
    }
    if (!reused) {
      result <- fit_e4_task(task)
      result$provenance <- provenance
      .atomic_save_rds(result, path)
    }
    rows[[i]] <- .e4_result_row(
      result, path, if (reused) "reused" else "complete")
  }
  raw <- .paper_annotate_manifest(do.call(rbind, rows), provenance)
  slopes <- e4_scaling_slopes(raw)
  slopes <- .paper_annotate_manifest(slopes, provenance)
  utils::write.csv(raw, file.path(output_dir, "manifest.csv"), row.names = FALSE)
  utils::write.csv(raw, file.path(output_dir, "fit-results.csv"), row.names = FALSE)
  utils::write.csv(slopes, file.path(output_dir, "scaling-slopes.csv"),
                   row.names = FALSE)
  write_paper_provenance(provenance, output_dir)
  invisible(list(fits = raw, slopes = slopes))
}
