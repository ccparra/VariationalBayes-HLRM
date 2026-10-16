# E2 identification and K-selection tasks. Only the smoke profile is executable
# until pilot controls and resource budgets are frozen.

expand_e2_tasks <- function(profile = c("smoke", "pilot", "full"),
                            base_seed = 310000L) {
  profile <- .paper_profile(profile)
  protocol <- e2_protocol(profile)
  fixed <- expand.grid(
    minimum_distance = protocol$data_generating$minimum_distance,
    replicate = protocol$execution$replicates,
    method = c("CAVI", "SVI"),
    component_count = protocol$execution$candidate_components,
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  fixed$task_type <- "fixed_k"
  overfit <- expand.grid(
    minimum_distance = protocol$data_generating$minimum_distance,
    replicate = protocol$execution$replicates,
    method = "MCMC",
    component_count = protocol$execution$overfitted_mcmc_components,
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  overfit$task_type <- "overfitted_mcmc"
  grid <- rbind(fixed, overfit)
  grid <- grid[order(grid$minimum_distance, grid$replicate,
                     match(grid$method, c("MCMC", "CAVI", "SVI")),
                     grid$component_count), ]
  rownames(grid) <- NULL
  distance_index <- match(
    grid$minimum_distance, protocol$data_generating$minimum_distance)
  dataset_seed <- as.integer(base_seed + 100000L * distance_index +
                               grid$replicate)
  fit_seed <- as.integer(dataset_seed +
    match(grid$method, c("MCMC", "CAVI", "SVI")) * 10000L +
    grid$component_count * 100L)
  data.frame(
    schema_version = "paper-task-manifest-1", profile = profile,
    experiment = "E2", task_type = grid$task_type,
    minimum_distance = grid$minimum_distance,
    replicate = as.integer(grid$replicate), method = grid$method,
    component_count = as.integer(grid$component_count),
    dataset_seed = dataset_seed, fit_seed = fit_seed,
    task_id = sprintf(
      "E2-d%s-r%03d-%s-K%02d",
      sub("\\.", "p", format(grid$minimum_distance, trim = TRUE)),
      grid$replicate, tolower(grid$method), grid$component_count),
    confirmatory = identical(profile, "full"), stringsAsFactors = FALSE)
}

.e2_smoke_control <- function(method, group_count) {
  control <- .smoke_control(method, group_count)
  if (method %in% c("CAVI", "SVI")) {
    control$starts <- 2L
    control$require_converged <- FALSE
  }
  control
}

fit_e2_task <- function(task) {
  if (!is.data.frame(task) || nrow(task) != 1L ||
      !all(c("profile", "minimum_distance", "replicate", "method",
             "component_count", "dataset_seed", "fit_seed", "task_id") %in%
           names(task)))
    stop("task must be one row from expand_e2_tasks()")
  if (!identical(task$profile, "smoke"))
    stop("Only E2 smoke tasks may run before pilot settings are frozen")
  scenario <- scenario_e2(
    task$minimum_distance, task$replicate, task$profile,
    base_seed = task$dataset_seed -
      100000L * match(task$minimum_distance, c(1.5, .5)) - task$replicate)
  simulation <- simulate_data(
    scenario$generator, scenario$design, scenario$sample,
    seed = task$dataset_seed)
  candidate_model <- chlrm_model(task$component_count,
                                 scenario$model$coefficient_count)
  fit <- fit_model(
    candidate_model, simulation$train,
    review_chlrm_prior(candidate_model$coefficient_count,
                       candidate_model$component_count),
    method = tolower(task$method),
    initialization = if (task$method == "MCMC") initialization_kmeans() else
      initialization_from_prior(),
    control = .e2_smoke_control(task$method,
                                length(simulation$train$y)),
    seed = task$fit_seed, verbose = FALSE)
  held_out <- predictive_performance(
    fit, simulation$test, predictive_target = "observed_groups",
    predictive_unit = "observation",
    group_ids = seq_along(simulation$test$y), n_draws = 4,
    seed = task$fit_seed + 1L)
  classification <- classify_groups(fit)
  clustering <- list(
    adjusted_rand_index = adjusted_rand_index(
      simulation$truth$allocations, classification$component),
    co_clustering_error = co_clustering_error(
      simulation$truth$allocations, co_clustering(fit)),
    normalized_entropy = normalized_classification_entropy(
      responsibilities(fit)))
  result <- list(
    schema_version = "paper-e2-task-result-2", task = task,
    truth_component_count = scenario$model$component_count,
    fit = list(method = fit$method, status = fit$status,
               diagnostics = fit$diagnostics, timing = fit$timing,
               metadata = fit$metadata),
    held_out = held_out, clustering = clustering)
  if (task$method %in% c("CAVI", "SVI")) {
    result$selection_statistic <- list(
      final_elbo = tail(fit$objective$elbo[
        is.finite(fit$objective$elbo)], 1L),
      mlpd = held_out$statistics$mean_log_predictive_density)
  } else {
    draws <- posterior_draws(
      fit, n_draws = 4, seed = task$fit_seed + 2L, include_latent = TRUE)
    occupied <- apply(
      draws$latent_variables$group_allocations, 1L,
      function(labels) length(unique(labels)))
    result$selection_statistic <- list(
      occupied_components = occupied,
      posterior_mode_occupied_components = .mode_integer(occupied))
  }
  structure(result, class = c("paper_e2_task_result", "list"))
}

.e2_result_row <- function(result) {
  if (!inherits(result, "paper_e2_task_result"))
    stop("result must be produced by fit_e2_task()")
  variational <- result$task$method %in% c("CAVI", "SVI")
  data.frame(
    task_id = result$task$task_id,
    minimum_distance = result$task$minimum_distance,
    replicate = result$task$replicate, method = result$task$method,
    component_count = result$task$component_count,
    final_elbo = if (variational)
      result$selection_statistic$final_elbo else NA_real_,
    mlpd = result$held_out$statistics$mean_log_predictive_density,
    posterior_mode_occupied_components = if (!variational)
      result$selection_statistic$posterior_mode_occupied_components else NA_integer_,
    ari = result$clustering$adjusted_rand_index,
    cpe = result$clustering$co_clustering_error,
    normalized_entropy = result$clustering$normalized_entropy,
    fit_seconds = result$fit$timing$total_seconds,
    estimated_peak_memory_mb =
      result$fit$metadata$resource_estimate$estimated_peak_memory_mb,
    convergence_status = result$fit$diagnostics$status,
    stringsAsFactors = FALSE)
}

summarize_e2_selection <- function(fit_results, truth = 3L) {
  if (!is.data.frame(fit_results) || !nrow(fit_results))
    stop("fit_results must be a nonempty data frame")
  keys <- unique(fit_results[c("minimum_distance", "replicate")])
  rows <- list(); index <- 0L
  for (i in seq_len(nrow(keys))) {
    same_dataset <- fit_results$minimum_distance == keys$minimum_distance[i] &
      fit_results$replicate == keys$replicate[i]
    for (method in c("CAVI", "SVI")) {
      values <- fit_results[same_dataset & fit_results$method == method, ]
      if (!nrow(values)) next
      for (criterion in c("ELBO", "MLPD")) {
        score <- if (criterion == "ELBO") values$final_elbo else values$mlpd
        selected <- which.max(score)
        index <- index + 1L
        rows[[index]] <- data.frame(
          minimum_distance = keys$minimum_distance[i],
          replicate = keys$replicate[i], method = method,
          criterion = criterion,
          selected_components = values$component_count[selected],
          correct = values$component_count[selected] == truth,
          absolute_error = abs(values$component_count[selected] - truth),
          predictive_regret = predictive_regret(values$mlpd, selected),
          stringsAsFactors = FALSE)
      }
    }
    mcmc <- fit_results[same_dataset & fit_results$method == "MCMC", ]
    if (nrow(mcmc)) {
      index <- index + 1L
      selected <- mcmc$posterior_mode_occupied_components[1L]
      rows[[index]] <- data.frame(
        minimum_distance = keys$minimum_distance[i],
        replicate = keys$replicate[i], method = "MCMC",
        criterion = "occupied_components_mode",
        selected_components = selected, correct = selected == truth,
        absolute_error = abs(selected - truth), predictive_regret = NA_real_,
        stringsAsFactors = FALSE)
    }
  }
  do.call(rbind, rows)
}

run_e2_smoke <- function(output_dir, resume = TRUE) {
  tasks <- expand_e2_tasks("smoke")
  provenance <- paper_run_provenance(
    "E2", "smoke", list(protocol = e2_protocol("smoke"), tasks = tasks))
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
        result, "paper_e2_task_result", "paper-e2-task-result-2",
        task$task_id, provenance)
    }
    if (!reused) {
      result <- fit_e2_task(task)
      result$provenance <- provenance
      .atomic_save_rds(result, path)
    }
    row <- .e2_result_row(result)
    row$status <- if (reused) "reused" else "complete"
    row$path <- file.path("tasks", basename(path))
    row$result_file_mb <- file.info(path)$size / 1024^2
    rows[[i]] <- row
  }
  raw <- .paper_annotate_manifest(do.call(rbind, rows), provenance)
  selection <- summarize_e2_selection(raw)
  selection <- .paper_annotate_manifest(selection, provenance)
  manifest_columns <- c(
    "task_id", "minimum_distance", "replicate", "method",
    "component_count", "status", "convergence_status", "fit_seconds",
    "estimated_peak_memory_mb", "result_file_mb", "path",
    "configuration_fingerprint", "source_fingerprint",
    "suite_source_fingerprint", "run_fingerprint", "output_fingerprint",
    "package_version", "git_commit", "git_dirty")
  utils::write.csv(raw[manifest_columns],
                   file.path(output_dir, "manifest.csv"), row.names = FALSE)
  utils::write.csv(raw, file.path(output_dir, "fit-results.csv"), row.names = FALSE)
  utils::write.csv(selection, file.path(output_dir, "selection-results.csv"),
                   row.names = FALSE)
  write_paper_provenance(provenance, output_dir)
  invisible(list(fits = raw, selection = selection))
}
