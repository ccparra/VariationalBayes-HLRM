# E5 Farms reanalysis tasks. The smoke profile exercises full-data fitting,
# within-farm CV, multistart stability, WAIC, and overfitted MCMC occupancy.

expand_e5_tasks <- function(profile = c("smoke", "pilot", "full"),
                            fit_seed_base = 610000L) {
  profile <- .paper_profile(profile)
  protocol <- e5_protocol(profile)
  fixed <- expand.grid(
    task_type = "fixed_k", method = c("MCMC", "CAVI", "SVI"),
    component_count = protocol$execution$candidate_components,
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  overfit <- data.frame(
    task_type = "overfitted_mcmc", method = "MCMC",
    component_count = protocol$execution$overfitted_mcmc_components)
  grid <- rbind(fixed, overfit)
  grid <- grid[order(grid$task_type, grid$component_count,
                     match(grid$method, c("MCMC", "CAVI", "SVI"))), ]
  rownames(grid) <- NULL
  fit_seed <- as.integer(fit_seed_base + seq_len(nrow(grid)))
  data.frame(
    schema_version = "paper-task-manifest-1", profile = profile,
    experiment = "E5", task_type = grid$task_type, method = grid$method,
    component_count = as.integer(grid$component_count), fit_seed = fit_seed,
    task_id = sprintf("E5-%s-K%02d-%s",
      ifelse(grid$task_type == "fixed_k", "fixed", "overfit"),
      grid$component_count, tolower(grid$method)),
    confirmatory = identical(profile, "full"), stringsAsFactors = FALSE)
}

.e5_control <- function(method, group_count) .smoke_control(method, group_count)

fit_e5_task <- function(task, data, cv_plan = NULL) {
  if (!is.data.frame(task) || nrow(task) != 1L ||
      !all(c("profile", "task_type", "method", "component_count",
             "fit_seed", "task_id") %in% names(task)))
    stop("task must be one row from expand_e5_tasks()")
  if (!identical(task$profile, "smoke"))
    stop("Only E5 smoke tasks may run before pilot settings are frozen")
  model <- chlrm_model(task$component_count, ncol(data$X[[1L]]))
  prior <- review_chlrm_prior(model$coefficient_count, model$component_count)
  control <- .e5_control(task$method, length(data$y))
  fit <- fit_model(
    model, data, prior, tolower(task$method), initialization_kmeans(), control,
    scaling = data_scaling(), seed = task$fit_seed, verbose = FALSE)
  result <- list(
    schema_version = "paper-e5-task-result-3", task = task,
    fit = list(diagnostics = fit$diagnostics, timing = fit$timing,
               metadata = fit$metadata))
  if (task$task_type == "overfitted_mcmc") {
    draws <- posterior_draws(
      fit, n_draws = 4, seed = task$fit_seed + 1L, include_latent = TRUE)
    occupied <- apply(draws$latent_variables$group_allocations, 1L,
                      function(labels) length(unique(labels)))
    result$occupied_components <- occupied
    result$posterior_mode_occupied_components <- .mode_integer(occupied)
  } else {
    criteria <- predictive_criteria(
      fit, n_draws = 4, seed = task$fit_seed + 1L,
      return_log_likelihood = FALSE)
    if (is.null(cv_plan))
      cv_plan <- cross_validation_plan(
        data, "within_group_kfold", folds = 5L, seed = 610001L)
    cv <- cross_validate_model(
      model, data, prior, tolower(task$method), initialization_kmeans(),
      control, cv_plan, scaling = data_scaling(), n_draws = 4,
      seed = task$fit_seed + 2L, max_fits = 5L,
      compute_performance = TRUE, retain_fits = FALSE, verbose = FALSE)
    result$waic_group_marginal <-
      criteria$criteria$waic_group_marginal$WAIC
    result$cross_validation <- cv
    if (task$method %in% c("CAVI", "SVI")) {
      orchestration <- fit_orchestration(fit)
      elbo <- orchestration$summary$final_elbo
      result$elbo <- list(
        maximum = max(elbo, na.rm = TRUE), median = stats::median(elbo, na.rm = TRUE),
        iqr = stats::IQR(elbo, na.rm = TRUE))
      result$clustering_stability <- clustering_stability(orchestration)
    }
  }
  structure(result, class = c("paper_e5_task_result", "list"))
}

.e5_result_row <- function(result, path, status) {
  fixed <- result$task$task_type == "fixed_k"
  data.frame(
    task_id = result$task$task_id, task_type = result$task$task_type,
    method = result$task$method,
    component_count = result$task$component_count,
    fit_seconds = result$fit$timing$total_seconds,
    convergence_status = result$fit$diagnostics$status,
    waic_group_marginal = if (fixed) result$waic_group_marginal else NA_real_,
    cv_mlpd = if (fixed)
      result$cross_validation$statistics$mean_log_predictive_density else NA_real_,
    cv_mse = if (fixed) result$cross_validation$statistics$rmse^2 else NA_real_,
    elbo_maximum = if (fixed && !is.null(result$elbo))
      result$elbo$maximum else NA_real_,
    elbo_median = if (fixed && !is.null(result$elbo))
      result$elbo$median else NA_real_,
    elbo_iqr = if (fixed && !is.null(result$elbo))
      result$elbo$iqr else NA_real_,
    clustering_stability = if (fixed && !is.null(result$clustering_stability))
      result$clustering_stability$stability else NA_real_,
    posterior_mode_occupied_components = if (!fixed)
      result$posterior_mode_occupied_components else NA_integer_,
    estimated_peak_memory_mb =
      result$fit$metadata$resource_estimate$estimated_peak_memory_mb,
    result_file_mb = file.info(path)$size / 1024^2,
    status = status, path = file.path("tasks", basename(path)),
    stringsAsFactors = FALSE)
}

run_e5_smoke <- function(output_dir, data, resume = TRUE) {
  tasks <- expand_e5_tasks("smoke")
  provenance <- paper_run_provenance(
    "E5", "smoke",
    list(protocol = e5_protocol("smoke"), tasks = tasks, data = data))
  plan <- cross_validation_plan(
    data, "within_group_kfold", folds = 5L, seed = 610001L)
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
        result, "paper_e5_task_result", "paper-e5-task-result-3",
        task$task_id, provenance)
    }
    if (!reused) {
      result <- fit_e5_task(task, data, plan)
      result$provenance <- provenance
      .atomic_save_rds(result, path)
    }
    rows[[i]] <- .e5_result_row(
      result, path, if (reused) "reused" else "complete")
  }
  raw <- .paper_annotate_manifest(do.call(rbind, rows), provenance)
  utils::write.csv(raw, file.path(output_dir, "manifest.csv"), row.names = FALSE)
  utils::write.csv(raw, file.path(output_dir, "fit-results.csv"), row.names = FALSE)
  utils::write.csv(raw[raw$task_type == "fixed_k", ],
                   file.path(output_dir, "K-comparison.csv"), row.names = FALSE)
  utils::write.csv(raw[raw$task_type == "overfitted_mcmc", ],
                   file.path(output_dir, "occupied-components.csv"), row.names = FALSE)
  write_paper_provenance(provenance, output_dir)
  invisible(raw)
}
