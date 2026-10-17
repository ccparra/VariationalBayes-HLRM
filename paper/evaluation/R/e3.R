# E3 SVI-sensitivity tasks. A single CAVI fit is the dataset-level reference.

expand_e3_tasks <- function(profile = c("smoke", "pilot", "full"),
                            dataset_seed = 410001L,
                            fit_seed_base = 420000L) {
  profile <- .paper_profile(profile)
  protocol <- e3_protocol(profile)
  svi <- merge(
    protocol$execution$grid,
    data.frame(seed_index = protocol$execution$seeds), all = TRUE)
  svi$method <- "SVI"
  reference <- data.frame(
    minibatch_size = NA_integer_, exponent = NA_real_, offset = NA_real_,
    seed_index = 0L, method = "CAVI")
  grid <- rbind(reference, svi)
  grid$fit_seed <- as.integer(fit_seed_base + seq_len(nrow(grid)))
  data.frame(
    schema_version = "paper-task-manifest-1", profile = profile,
    experiment = "E3", method = grid$method,
    minibatch_size = as.integer(grid$minibatch_size),
    exponent = grid$exponent, offset = grid$offset,
    seed_index = as.integer(grid$seed_index),
    dataset_seed = as.integer(dataset_seed), fit_seed = grid$fit_seed,
    task_id = ifelse(
      grid$method == "CAVI", "E3-cavi-reference",
      sprintf("E3-svi-b%03d-x%s-t%03d-s%02d",
              grid$minibatch_size,
              sub("\\.", "p", format(grid$exponent, trim = TRUE)),
              as.integer(grid$offset), grid$seed_index)),
    confirmatory = identical(profile, "full"), stringsAsFactors = FALSE)
}

e3_inference_control <- function(task, group_count,
                                 max_memory_mb = NULL) {
  if (!is.data.frame(task) || nrow(task) != 1L ||
      !all(c("profile", "method", "minibatch_size", "exponent", "offset") %in%
           names(task)))
    stop("task must be one row from expand_e3_tasks()")
  profile <- .paper_profile(task$profile)
  method <- match.arg(task$method, c("CAVI", "SVI"))
  if (!is.numeric(group_count) || length(group_count) != 1L ||
      !is.finite(group_count) || group_count < 1L || group_count != floor(group_count))
    stop("group_count must be a positive integer")
  group_count <- as.integer(group_count)

  if (profile == "smoke") {
    if (method == "CAVI") return(.smoke_control("CAVI", group_count))
    return(svi_control(
      max_iterations = 4L, minibatch_size = task$minibatch_size,
      learning_rate = power_schedule(task$offset, task$exponent),
      evaluation_draws = 4L,
      stopping = svi_stopping(
        objective_tolerance = 1e20, parameter_tolerance = 1e20,
        min_iterations = 2L, patience = 1L, check_every = 1L,
        smoothing_window = 1L),
      starts = 2L, require_converged = FALSE, max_memory_mb = 256))
  }

  evaluation_draws <- if (profile == "pilot") 250L else 5000L
  memory_budget <- if (is.null(max_memory_mb)) 2048 else max_memory_mb
  if (method == "CAVI") {
    return(cavi_control(
      max_iterations = 10000L, evaluation_draws = evaluation_draws,
      stopping = cavi_stopping(
        objective_tolerance = 1e-8, parameter_tolerance = NULL,
        min_iterations = 1L, patience = 5L),
      starts = 5L, require_converged = TRUE,
      max_memory_mb = memory_budget))
  }
  if (task$minibatch_size > group_count)
    stop("E3 minibatch_size cannot exceed group_count")
  # In E3 each of the 20 declared seeds is one independent SVI run. Setting
  # starts = 5 here would silently turn the specified 720 runs into 3600 and
  # would discard the per-seed failures that the experiment must estimate.
  svi_control(
    max_iterations = 50000L, minibatch_size = task$minibatch_size,
    learning_rate = power_schedule(task$offset, task$exponent),
    evaluation_draws = evaluation_draws,
    stopping = svi_stopping(
      objective_tolerance = 1e-4, parameter_tolerance = NULL,
      min_iterations = 1000L, patience = 5L, check_every = 50L,
      smoothing_window = 5L),
    starts = 1L, max_memory_mb = memory_budget)
}

expand_e3_pilot_stage <- function(stage = c("screen", "factorial")) {
  stage <- match.arg(stage)
  tasks <- expand_e3_tasks("pilot")
  if (stage == "factorial") return(tasks)
  svi <- tasks[tasks$method == "SVI" & tasks$seed_index == 1L, , drop = FALSE]
  corner <- svi$minibatch_size %in% c(25L, 200L) &
    svi$exponent %in% c(.6, .8) & svi$offset %in% c(1, 100)
  centre <- svi$minibatch_size == 100L & svi$exponent == .7 &
    svi$offset == 10
  rbind(tasks[tasks$method == "CAVI", , drop = FALSE],
        svi[corner | centre, , drop = FALSE])
}

.e3_control_row <- function(task, control, estimate) {
  stopping <- control$stopping
  data.frame(
    task_id = task$task_id, profile = task$profile, method = task$method,
    minibatch_size = task$minibatch_size,
    exponent = task$exponent, offset = task$offset,
    seed_index = task$seed_index, dataset_seed = task$dataset_seed,
    fit_seed = task$fit_seed,
    max_iterations = control$max_iterations,
    starts = control$starts, evaluation_draws = control$evaluation_draws,
    objective_tolerance = stopping$objective_tolerance,
    parameter_tolerance = if (is.null(stopping$parameter_tolerance)) NA_real_ else
      stopping$parameter_tolerance,
    min_iterations = stopping$min_iterations,
    patience = stopping$patience, check_every = stopping$check_every,
    smoothing_window = stopping$smoothing_window,
    max_memory_mb = control$max_memory_mb,
    estimated_peak_memory_mb = estimate$estimated_peak_memory_mb,
    estimated_fit_storage_mb = estimate$estimated_fit_storage_mb,
    estimated_evaluation_storage_mb = estimate$estimated_evaluation_storage_mb,
    within_memory_budget = estimate$within_budget,
    stringsAsFactors = FALSE)
}

prepare_e3_preflight <- function(output_dir, profile = c("pilot", "full"),
                                 stage = c("screen", "factorial"),
                                 max_memory_mb = 2048,
                                 screen_review = NULL) {
  profile <- match.arg(profile)
  stage <- match.arg(stage)
  if (profile == "full" && stage == "screen")
    stop("screen is a pilot-only stage")
  pilot_factorial <- identical(profile, "pilot") &&
    identical(stage, "factorial")
  if (pilot_factorial &&
      (!inherits(screen_review, "paper_e3_screen_review") ||
       !isTRUE(screen_review$summary$ready_for_factorial)))
    stop("A passing E3 screen review is required for factorial preflight")
  if (!pilot_factorial && !is.null(screen_review))
    stop("screen_review applies only to the pilot factorial preflight")
  tasks <- if (profile == "pilot") expand_e3_pilot_stage(stage) else
    expand_e3_tasks("full")
  scenario <- scenario_e3(profile, tasks$dataset_seed[[1L]])
  simulation <- simulate_data(
    scenario$generator, scenario$design, scenario$sample,
    seed = tasks$dataset_seed[[1L]])
  group_count <- length(simulation$train$y)
  rows <- lapply(seq_len(nrow(tasks)), function(index) {
    task <- tasks[index, , drop = FALSE]
    control <- e3_inference_control(task, group_count, max_memory_mb)
    estimate <- estimate_fit_resources(scenario$model, simulation$train, control)
    .e3_control_row(task, control, estimate)
  })
  plan <- do.call(rbind, rows)
  configuration <- list(
    protocol = e3_protocol(profile), stage = stage, tasks = tasks,
    controls = plan[c(
      "task_id", "max_iterations", "starts", "evaluation_draws",
      "objective_tolerance", "parameter_tolerance", "min_iterations",
      "patience", "check_every", "smoothing_window", "max_memory_mb")])
  provenance <- paper_run_provenance("E3", paste0(profile, "-preflight"),
                                     configuration)
  plan <- .paper_annotate_manifest(plan, provenance)
  data_size_mb <- as.numeric(utils::object.size(simulation$train)) / 1024^2
  summary <- data.frame(
    schema_version = "paper-e3-preflight-1", profile = profile, stage = stage,
    task_count = nrow(plan), cavi_tasks = sum(plan$method == "CAVI"),
    svi_tasks = sum(plan$method == "SVI"), group_count = group_count,
    observation_count = sum(lengths(simulation$train$y)),
    coefficient_count = scenario$model$coefficient_count,
    component_count = scenario$model$component_count,
    data_object_mb = data_size_mb,
    maximum_estimated_peak_memory_mb = max(plan$estimated_peak_memory_mb),
    total_estimated_fit_storage_mb = sum(plan$estimated_fit_storage_mb),
    maximum_estimated_evaluation_storage_mb =
      max(plan$estimated_evaluation_storage_mb),
    tasks_within_memory_budget = sum(plan$within_memory_budget),
    sequential_execution_only = TRUE,
    runtime_estimate_status = if (pilot_factorial)
      "projected_from_empirical_screen" else "requires empirical screen",
    screen_output_fingerprint = if (pilot_factorial)
      screen_review$summary$screen_output_fingerprint else NA_character_,
    observed_screen_max_rss_mb = if (pilot_factorial)
      screen_review$summary$observed_max_rss_mb else NA_real_,
    projected_wall_seconds = if (pilot_factorial)
      screen_review$summary$projected_factorial_wall_seconds else NA_real_,
    projected_storage_mb = if (pilot_factorial)
      screen_review$summary$projected_factorial_storage_mb else NA_real_,
    scientific_result = FALSE,
    stringsAsFactors = FALSE)
  summary <- .paper_annotate_manifest(summary, provenance)
  gate <- data.frame(
    requirement = c(
      "task_plan_created", "deterministic_fingerprints",
      "per_fit_memory_preflight", "sequential_execution",
      "empirical_wall_time", "empirical_peak_rss", "disk_capacity",
      "clean_git_commit", "configuration_frozen", "explicit_execution"),
    status = c(
      "pass", "pass",
      if (all(plan$within_memory_budget)) "pass" else "fail", "pass",
      if (pilot_factorial) "pass" else "pending",
      if (pilot_factorial) "pass" else "pending", "pending",
      if (isFALSE(provenance$git_dirty)) "pass" else "pending",
      "pending", "blocked"),
    detail = c(
      paste(nrow(plan), "tasks materialized"),
      provenance$run_fingerprint,
      paste("maximum dimension estimate",
            sprintf("%.1f MB", max(plan$estimated_peak_memory_mb))),
      "max_concurrent_fits = 1",
      if (pilot_factorial) sprintf(
        "screen projects %.1f seconds for %d tasks",
        screen_review$summary$projected_factorial_wall_seconds, nrow(plan)) else
        "measure with the screen stage before factorial execution",
      if (pilot_factorial) sprintf(
        "screen maximum %.1f MiB; %.1fx margin passes %.1f MiB budget",
        screen_review$summary$observed_max_rss_mb,
        screen_review$summary$rss_safety_factor,
        screen_review$summary$memory_budget_mb) else
        "measure operating-system RSS; the package estimate is not RSS",
      "verify free space after measuring screen artifacts",
      if (isFALSE(provenance$git_dirty)) provenance$git_commit else
        "working tree is dirty",
      if (identical(stage, "screen"))
        "freeze the reviewed preflight fingerprint" else
        "freeze only with an accepted passing screen review",
      paste("run-e3-pilot.R additionally requires --execute and both accepted",
            "fingerprints for factorial; full remains blocked")),
    stringsAsFactors = FALSE)
  gate$ready_to_execute <- FALSE
  gate <- .paper_annotate_manifest(gate, provenance)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(plan, file.path(output_dir, "task-plan.csv"), row.names = FALSE)
  utils::write.csv(summary, file.path(output_dir, "budget-summary.csv"),
                   row.names = FALSE)
  utils::write.csv(gate, file.path(output_dir, "execution-gate.csv"),
                   row.names = FALSE)
  saveRDS(configuration, file.path(output_dir, "configuration-candidate.rds"),
          version = 3)
  write_paper_provenance(provenance, output_dir)
  invisible(list(plan = plan, summary = summary, gate = gate,
                 provenance = provenance))
}

freeze_e3_configuration <- function(preflight_dir,
                                    accepted_configuration_fingerprint,
                                    screen_review = NULL,
                                    accepted_screen_output_fingerprint = NULL) {
  if (!is.character(accepted_configuration_fingerprint) ||
      length(accepted_configuration_fingerprint) != 1L ||
      !nzchar(accepted_configuration_fingerprint))
    stop("accepted_configuration_fingerprint must be one nonempty string")
  configuration_path <- file.path(preflight_dir, "configuration-candidate.rds")
  provenance_path <- file.path(preflight_dir, "provenance.rds")
  if (!file.exists(configuration_path) || !file.exists(provenance_path))
    stop("preflight_dir must contain configuration-candidate.rds and provenance.rds")
  configuration <- readRDS(configuration_path)
  recorded <- readRDS(provenance_path)
  if (identical(configuration$stage, "factorial")) {
    if (!inherits(screen_review, "paper_e3_screen_review") ||
        !isTRUE(screen_review$summary$ready_for_factorial))
      stop("A passing E3 screen review is required before freezing factorial")
    if (!is.character(accepted_screen_output_fingerprint) ||
        length(accepted_screen_output_fingerprint) != 1L ||
        !nzchar(accepted_screen_output_fingerprint) ||
        !identical(screen_review$summary$screen_output_fingerprint,
                   accepted_screen_output_fingerprint))
      stop("Accepted screen output fingerprint does not match its review")
    upstream_evidence <- as.list(screen_review$summary[1L, c(
      "schema_version", "screen_output_fingerprint", "screen_git_commit",
      "observed_max_rss_mb", "rss_safety_factor", "memory_budget_mb",
      "projected_factorial_wall_seconds",
      "projected_factorial_storage_mb", "ready_for_factorial")])
  } else {
    if (!is.null(screen_review) || !is.null(accepted_screen_output_fingerprint))
      stop("Screen evidence applies only to the factorial stage")
    upstream_evidence <- NULL
  }
  current <- paper_run_provenance("E3", recorded$profile, configuration)
  if (!identical(recorded$configuration_fingerprint,
                 accepted_configuration_fingerprint) ||
      !identical(current$configuration_fingerprint,
                 accepted_configuration_fingerprint))
    stop("Accepted fingerprint does not match the current E3 configuration")
  if (!identical(recorded$source_fingerprint, current$source_fingerprint))
    stop("E3 calculation sources changed after preflight; regenerate it")
  if (!isFALSE(recorded$git_dirty) ||
      !identical(recorded$git_commit, current$git_commit))
    stop("Regenerate the E3 preflight from the current clean Git commit")
  if (!isFALSE(current$git_dirty))
    stop("E3 configuration can only be frozen from a clean Git commit")
  frozen <- structure(list(
    schema_version = "paper-e3-frozen-configuration-2",
    experiment = "E3", profile = configuration$protocol$profile,
    stage = configuration$stage,
    configuration = configuration,
    upstream_evidence = upstream_evidence,
    configuration_fingerprint = current$configuration_fingerprint,
    source_fingerprint = current$source_fingerprint,
    git_commit = current$git_commit,
    frozen_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)),
    class = c("paper_e3_frozen_configuration", "list"))
  saveRDS(frozen, file.path(preflight_dir, "frozen-configuration.rds"),
          version = 3)
  utils::write.csv(data.frame(
    schema_version = frozen$schema_version, experiment = frozen$experiment,
    profile = frozen$profile, stage = frozen$stage,
    configuration_fingerprint = frozen$configuration_fingerprint,
    source_fingerprint = frozen$source_fingerprint,
    screen_output_fingerprint = if (is.null(frozen$upstream_evidence)) NA_character_ else
      frozen$upstream_evidence$screen_output_fingerprint,
    git_commit = frozen$git_commit, frozen_at_utc = frozen$frozen_at_utc,
    stringsAsFactors = FALSE),
    file.path(preflight_dir, "frozen-configuration.csv"), row.names = FALSE)
  invisible(frozen)
}

validate_e3_frozen_configuration <- function(frozen, stage = "screen") {
  if (!inherits(frozen, "paper_e3_frozen_configuration") ||
      !identical(frozen$schema_version, "paper-e3-frozen-configuration-2"))
    stop("frozen must be produced by freeze_e3_configuration()")
  if (!identical(frozen$profile, "pilot") || !identical(frozen$stage, stage))
    stop("Frozen E3 configuration does not match the requested pilot stage")
  current <- paper_run_provenance(
    "E3", paste0(frozen$profile, "-preflight"), frozen$configuration)
  if (!identical(current$configuration_fingerprint,
                 frozen$configuration_fingerprint) ||
      !identical(current$source_fingerprint, frozen$source_fingerprint) ||
      !identical(current$git_commit, frozen$git_commit) ||
      !isFALSE(current$git_dirty))
    stop("Frozen E3 configuration no longer matches the clean source commit")
  if (identical(stage, "factorial") &&
      (is.null(frozen$upstream_evidence) ||
       !isTRUE(frozen$upstream_evidence$ready_for_factorial) ||
       !is.character(frozen$upstream_evidence$screen_output_fingerprint) ||
       !nzchar(frozen$upstream_evidence$screen_output_fingerprint)))
    stop("Frozen factorial lacks accepted E3 screen evidence")
  invisible(current)
}

fit_e3_task <- function(task, pilot_authorized = FALSE,
                        max_memory_mb = NULL) {
  if (!is.data.frame(task) || nrow(task) != 1L ||
      !all(c("profile", "method", "dataset_seed", "fit_seed", "task_id") %in%
           names(task)))
    stop("task must be one row from expand_e3_tasks()")
  if (identical(task$profile, "full"))
    stop("E3 full tasks remain blocked until the pilot is accepted")
  if (!identical(task$profile, "smoke") && !isTRUE(pilot_authorized))
    stop("E3 pilot tasks require a validated frozen configuration")
  scenario <- scenario_e3(task$profile, task$dataset_seed)
  simulation <- simulate_data(
    scenario$generator, scenario$design, scenario$sample,
    seed = task$dataset_seed)
  fit <- fit_model(
    scenario$model, simulation$train,
    review_chlrm_prior(scenario$model$coefficient_count,
                       scenario$model$component_count),
    method = tolower(task$method), initialization = initialization_from_prior(),
    control = e3_inference_control(
      task, length(simulation$train$y), max_memory_mb = max_memory_mb),
    seed = task$fit_seed, verbose = FALSE)
  held_out <- predictive_performance(
    fit, simulation$test, predictive_target = "observed_groups",
    predictive_unit = "observation",
    group_ids = seq_along(simulation$test$y),
    n_draws = if (identical(task$profile, "smoke")) 4L else 250L,
    seed = task$fit_seed + 1L)
  classification <- classify_groups(fit)
  structure(list(
    schema_version = "paper-e3-task-result-3", task = task,
    fit = list(diagnostics = fit$diagnostics, timing = fit$timing,
               metadata = fit$metadata),
    final_elbo = tail(fit$objective$elbo[is.finite(fit$objective$elbo)], 1L),
    held_out = held_out,
    clustering = list(
      adjusted_rand_index = adjusted_rand_index(
        simulation$truth$allocations, classification$component),
      co_clustering_error = co_clustering_error(
        simulation$truth$allocations, co_clustering(fit)))),
    class = c("paper_e3_task_result", "list"))
}

.e3_result_row <- function(result, path, status) data.frame(
  task_id = result$task$task_id, method = result$task$method,
  minibatch_size = result$task$minibatch_size,
  exponent = result$task$exponent, offset = result$task$offset,
  seed_index = result$task$seed_index,
  final_elbo = result$final_elbo,
  mlpd = result$held_out$statistics$mean_log_predictive_density,
  ari = result$clustering$adjusted_rand_index,
  cpe = result$clustering$co_clustering_error,
  runtime = result$fit$timing$total_seconds,
  iterations = result$fit$diagnostics$work_completed,
  convergence_status = result$fit$diagnostics$status,
  estimated_peak_memory_mb =
    result$fit$metadata$resource_estimate$estimated_peak_memory_mb,
  result_file_mb = file.info(path)$size / 1024^2,
  status = status, path = file.path("tasks", basename(path)),
  stringsAsFactors = FALSE)

compare_e3_to_reference <- function(results) {
  reference <- results[results$method == "CAVI", , drop = FALSE]
  if (nrow(reference) != 1L)
    stop("E3 results require exactly one CAVI reference")
  svi <- results[results$method == "SVI", , drop = FALSE]
  svi$delta_elbo <- reference$final_elbo - svi$final_elbo
  svi$delta_mlpd <- svi$mlpd - reference$mlpd
  svi$delta_ari <- svi$ari - reference$ari
  svi
}

run_e3_smoke <- function(output_dir, resume = TRUE) {
  tasks <- expand_e3_tasks("smoke")
  provenance <- paper_run_provenance(
    "E3", "smoke", list(protocol = e3_protocol("smoke"), tasks = tasks))
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
        result, "paper_e3_task_result", "paper-e3-task-result-3",
        task$task_id, provenance)
    }
    if (!reused) {
      result <- fit_e3_task(task)
      result$provenance <- provenance
      .atomic_save_rds(result, path)
    }
    rows[[i]] <- .e3_result_row(
      result, path, if (reused) "reused" else "complete")
  }
  raw <- .paper_annotate_manifest(do.call(rbind, rows), provenance)
  comparison <- compare_e3_to_reference(raw)
  comparison <- .paper_annotate_manifest(comparison, provenance)
  utils::write.csv(raw, file.path(output_dir, "manifest.csv"), row.names = FALSE)
  utils::write.csv(raw, file.path(output_dir, "fit-results.csv"), row.names = FALSE)
  utils::write.csv(comparison, file.path(output_dir, "svi-comparison.csv"),
                   row.names = FALSE)
  write_paper_provenance(provenance, output_dir)
  invisible(list(fits = raw, comparison = comparison))
}
