simulate_data <- function(generator, design, sample, seed = NULL,
                          retain_observation_truth = FALSE,
                          max_output_mb = 256) {
  if (!inherits(generator, "vb_simulation_spec") ||
      !inherits(design, "vb_simulation_design") ||
      !inherits(sample, "vb_simulation_sample"))
    stop("generator, design and sample must be created by simulation constructors")
  validate_rng_seed(seed)
  retain_observation_truth <- .logical_scalar(
    retain_observation_truth, "retain_observation_truth")
  model <- generator$model
  grouped <- model$model != "LRM"
  if (grouped && is.null(sample$group_count))
    stop("HLRM/CHLRM require a grouped simulation_sample()")
  if (!grouped && is.null(sample$sample_size))
    stop("LRM requires an ungrouped simulation_sample()")
  p <- .dgp_dimension(generator)
  estimated_output_mb <- .check_simulation_budget(
    sample, p, retain_observation_truth, max_output_mb)
  scope <- inference_rng_scope(seed)
  on.exit(scope$restore(), add = TRUE)

  parameter_source <- list()
  if (generator$type == "prior_predictive") {
    parameters <- .draw_prior_parameters(model, generator$prior)
    parameter_source$global <- "prior"
  } else {
    parameters <- generator$parameters
    parameter_source$global <- "conditioned"
  }

  has_test <- if (grouped) !is.null(sample$test_observations_per_group) else
    !is.null(sample$test_sample_size)
  truth <- list()
  if (model$model == "LRM") {
    if (generator$type == "hierarchical") stop("LRM has no hierarchical DGP")
    parameters <- .validate_conditional_parameters(
      conditional_dgp(model, parameters), sample, p)
    train <- .simulate_lrm_split(
      design, sample$sample_size, p, parameters, "train")
    truth <- parameters
    if (retain_observation_truth) truth <- c(truth, list(
      train_conditional_mean = train$conditional_mean,
      train_residuals = train$residuals))
    test <- NULL
    if (has_test) {
      test_generated <- .simulate_lrm_split(
        design, sample$test_sample_size, p, parameters, "test")
      test <- test_generated$data
      if (retain_observation_truth) {
        truth$test_conditional_mean <- test_generated$conditional_mean
        truth$test_residuals <- test_generated$residuals
      }
    }
    train <- train$data
  } else if (model$model == "HLRM") {
    if (generator$type == "conditional") {
      parameters <- .validate_conditional_parameters(generator, sample, p)
      locals <- parameters
      parameter_source$local <- "conditioned"
      population <- list()
    } else {
      population <- .validate_population_parameters(parameters, model, p)
      locals <- .draw_group_truth(population, sample$group_count)
      parameter_source$local <- "hierarchy"
    }
    train_generated <- .simulate_grouped_split(
      design, sample$observations_per_group, p, locals$group_coefficients,
      locals$group_variances, "train")
    truth <- c(population, locals)
    if (retain_observation_truth) truth <- c(truth, list(
      train_conditional_mean = train_generated$conditional_mean,
      train_residuals = train_generated$residuals))
    train <- train_generated$data; test <- NULL
    if (has_test) {
      test_locals <- if (sample$test_group_mode == "same") locals else
        .draw_group_truth(population, sample$test_group_count)
      test_generated <- .simulate_grouped_split(
        design, sample$test_observations_per_group, p,
        test_locals$group_coefficients, test_locals$group_variances, "test")
      test <- test_generated$data
      if (sample$test_group_mode == "new") {
        truth$test_group_coefficients <- test_locals$group_coefficients
        truth$test_group_variances <- test_locals$group_variances
      }
      if (retain_observation_truth) {
        truth$test_conditional_mean <- test_generated$conditional_mean
        truth$test_residuals <- test_generated$residuals
      }
    }
  } else {
    if (generator$type == "conditional") {
      components <- .validate_conditional_parameters(generator, sample, p)
      population <- list()
      parameter_source$components <- "conditioned"
    } else {
      population <- .validate_population_parameters(parameters, model, p)
      components <- c(.draw_component_truth(population, model),
                      list(mixture_weights = population$mixture_weights))
      parameter_source$components <- "hierarchy"
    }
    K <- model$component_count
    allocations <- components$allocations
    if (is.null(allocations)) allocations <- sample.int(
      K, sample$group_count, replace = TRUE, prob = components$mixture_weights)
    local_coefficients <- components$component_coefficients[allocations, , drop = FALSE]
    local_variances <- components$component_variances[allocations]
    train_generated <- .simulate_grouped_split(
      design, sample$observations_per_group, p, local_coefficients,
      local_variances, "train")
    component_truth <- components[c(
      "component_coefficients", "component_variances")]
    if (!length(population)) component_truth$mixture_weights <-
      components$mixture_weights
    parameter_source$allocations <- if (!is.null(components$allocations))
      "conditioned" else "mixture_distribution"
    truth <- c(population, component_truth,
      list(allocations = allocations,
           realized_mixture_proportions = tabulate(allocations, K) / length(allocations)))
    if (retain_observation_truth) truth <- c(truth, list(
      train_conditional_mean = train_generated$conditional_mean,
      train_residuals = train_generated$residuals))
    train <- train_generated$data; test <- NULL
    if (has_test) {
      test_allocations <- if (sample$test_group_mode == "same") allocations else
        sample.int(K, sample$test_group_count, replace = TRUE,
                   prob = components$mixture_weights)
      test_generated <- .simulate_grouped_split(
        design, sample$test_observations_per_group, p,
        components$component_coefficients[test_allocations, , drop = FALSE],
        components$component_variances[test_allocations], "test")
      test <- test_generated$data
      if (sample$test_group_mode == "new") {
        truth$test_allocations <- test_allocations
        truth$test_realized_mixture_proportions <-
          tabulate(test_allocations, K) / length(test_allocations)
      }
      if (retain_observation_truth) {
        truth$test_conditional_mean <- test_generated$conditional_mean
        truth$test_residuals <- test_generated$residuals
      }
    }
  }

  conditioning <- list(
    dgp = generator$type,
    global_parameters = parameter_source$global,
    test_units = if (!has_test) "not_generated" else if (!grouped) "new_observations" else
      paste0(sample$test_group_mode, "_groups"))
  design_source <- list(
    distribution = design$distribution,
    random = design$distribution != "fixed",
    include_intercept = design$include_intercept)
  new_simulation_result(
    model$model, train, truth,
    configuration = list(
      generator = generator, design = design, sample = sample, seed = seed,
      retain_observation_truth = retain_observation_truth,
      max_output_mb = max_output_mb,
      estimated_peak_memory_mb = estimated_output_mb),
    conditioning = conditioning,
    parameter_source = parameter_source,
    design_source = design_source, test = test)
}
