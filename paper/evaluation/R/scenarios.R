# Paper-specific recipes. These functions configure the package API but are not
# part of the reusable model implementation.

.paper_positive_integer <- function(value, name) {
  if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
      value < 1 || value != floor(value))
    stop(name, " must be a positive integer")
  as.integer(value)
}

.paper_positive_scalar <- function(value, name) {
  if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
      value <= 0)
    stop(name, " must be a positive finite scalar")
  as.numeric(value)
}

make_review_component_coefficients <- function(component_count,
                                               coefficient_count,
                                               minimum_distance) {
  component_count <- .paper_positive_integer(
    component_count, "component_count")
  coefficient_count <- .paper_positive_integer(
    coefficient_count, "coefficient_count")
  if (component_count < 2L || coefficient_count < 5L)
    stop("The review geometry requires component_count >= 2 and coefficient_count >= 5")
  minimum_distance <- .paper_positive_scalar(
    minimum_distance, "minimum_distance")
  angle <- 2 * pi * (seq_len(component_count) - 1L) / component_count
  coefficients <- cbind(
    0, cos(angle), sin(angle), cos(2 * angle) / sqrt(2),
    sin(2 * angle) / sqrt(2))
  coefficients <- coefficients *
    (minimum_distance / min(stats::dist(coefficients)))
  if (coefficient_count > 5L) coefficients <- cbind(
    coefficients, matrix(0, component_count, coefficient_count - 5L))
  unname(coefficients)
}

review_chlrm_prior <- function(coefficient_count, component_count) {
  coefficient_count <- .paper_positive_integer(
    coefficient_count, "coefficient_count")
  component_count <- .paper_positive_integer(
    component_count, "component_count")
  chlrm_prior(
    coefficient_mean = rep(0, coefficient_count),
    coefficient_covariance = 100 * diag(coefficient_count),
    coefficient_dispersion_df = coefficient_count + 2,
    coefficient_dispersion_scale = diag(coefficient_count),
    residual_degrees_of_freedom = 1,
    residual_variance_scale_shape = 1,
    residual_variance_scale_rate = 1,
    mixture_concentration = rep(1 / component_count, component_count))
}

e1_scenario_catalog <- function() {
  data.frame(
    scenario = c("S1", "S2", "S3", "S4", "S5", "S6",
                 "S7a", "S7b", "S7c"),
    description = c(
      "baseline", "weak_separation", "strong_predictor_dependence",
      "unbalanced_groups", "weak_separation_unbalanced_rare_component",
      "high_dimension", "increasing_components_3",
      "increasing_components_6", "increasing_components_10"),
    component_count = c(3L, 3L, 3L, 3L, 3L, 3L, 3L, 6L, 10L),
    group_count = c(50L, 50L, 50L, 50L, 80L, 100L, 90L, 180L, 300L),
    group_size_mode = c("fixed", "fixed", "fixed", "discrete_uniform",
                        "discrete_uniform", "fixed", "fixed", "fixed",
                        "fixed"),
    observations_per_group = c(20L, 20L, 20L, NA, NA, 20L, 20L, 20L, 20L),
    coefficient_count = c(5L, 5L, 5L, 5L, 5L, 20L, 10L, 10L, 10L),
    predictor_correlation = c(0, 0, .8, 0, 0, .5, 0, 0, 0),
    minimum_distance = c(1.5, .5, 1.5, 1.5, .75, 1, 1.5, 1.5, 1.5),
    weight_scheme = c(
      "0.4_0.3_0.3", "0.4_0.3_0.3", "0.4_0.3_0.3",
      "0.4_0.3_0.3", "0.75_0.20_0.05", "0.4_0.3_0.3",
      "uniform", "uniform", "uniform"),
    stringsAsFactors = FALSE)
}

.paper_seed_scope <- function(seed) {
  existed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  previous <- if (existed)
    get(".Random.seed", envir = .GlobalEnv, inherits = FALSE) else NULL
  set.seed(seed)
  function() {
    if (existed) {
      assign(".Random.seed", previous, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
    invisible(NULL)
  }
}

.e1_task_seed <- function(scenario, replicate, base_seed) {
  catalog <- e1_scenario_catalog()
  scenario_index <- match(scenario, catalog$scenario)
  if (is.na(scenario_index))
    stop("scenario must be one of: ", paste(catalog$scenario, collapse = ", "))
  replicate <- .paper_positive_integer(replicate, "replicate")
  base_seed <- .paper_positive_integer(base_seed, "base_seed")
  seed <- base_seed + 100000L * scenario_index + replicate
  if (seed > .Machine$integer.max)
    stop("base_seed, scenario and replicate produce an invalid R seed")
  as.integer(seed)
}

.ar1_covariance <- function(dimension, correlation) {
  if (dimension < 1L) return(matrix(numeric(), 0L, 0L))
  outer(seq_len(dimension), seq_len(dimension), function(r, s)
    correlation^abs(r - s))
}

.e1_weights <- function(scheme, component_count) {
  switch(scheme,
    "0.4_0.3_0.3" = c(.4, .3, .3),
    "0.75_0.20_0.05" = c(.75, .20, .05),
    uniform = rep(1 / component_count, component_count),
    stop("Unknown E1 weight scheme"))
}

scenario_e1 <- function(scenario, replicate = 1L, base_seed = 210000L) {
  catalog <- e1_scenario_catalog()
  row <- catalog[catalog$scenario == scenario, , drop = FALSE]
  if (nrow(row) != 1L)
    stop("scenario must be one of: ", paste(catalog$scenario, collapse = ", "))
  seed <- .e1_task_seed(scenario, replicate, base_seed)
  K <- row$component_count
  p <- row$coefficient_count
  q <- p - 1L
  covariance <- .ar1_covariance(q, row$predictor_correlation)
  group_sizes <- if (row$group_size_mode == "fixed") {
    rep.int(as.integer(row$observations_per_group), row$group_count)
  } else {
    restore <- .paper_seed_scope(seed)
    on.exit(restore(), add = TRUE)
    sample.int(46L, row$group_count, replace = TRUE) + 4L
  }
  weights <- .e1_weights(row$weight_scheme, K)
  coefficients <- make_review_component_coefficients(
    K, p, row$minimum_distance)
  model <- chlrm_model(K, p)
  configuration <- as.list(row[1L, ])
  configuration$replicate <- .paper_positive_integer(replicate, "replicate")
  configuration$seed <- seed
  configuration$test_observations_per_group <- 20L
  configuration$residual_variance <- 1
  structure(list(
    schema_version = "paper-e1-scenario-1",
    experiment = "E1",
    scenario = scenario,
    replicate = configuration$replicate,
    seed = seed,
    model = model,
    generator = conditional_dgp(model, list(
      component_coefficients = coefficients,
      component_variances = rep(1, K),
      mixture_weights = weights)),
    design = simulation_design(
      "multivariate_normal", include_intercept = TRUE,
      mean = rep(0, q), covariance = covariance),
    sample = simulation_sample(
      group_count = row$group_count,
      observations_per_group = group_sizes,
      test_observations_per_group = 20L,
      test_group_mode = "same"),
    prior = review_chlrm_prior(p, K),
    configuration = configuration),
    class = c("paper_e1_scenario", "list"))
}

.paper_profile <- function(profile) match.arg(profile, c("smoke", "pilot", "full"))

e2_protocol <- function(profile = c("smoke", "pilot", "full")) {
  profile <- .paper_profile(profile)
  execution <- switch(profile,
    smoke = list(replicates = 1L, candidate_components = c(2L, 3L),
                 group_count = 12L, observations_per_group = 4L,
                 overfitted_mcmc_components = 4L),
    pilot = list(replicates = 1:5, candidate_components = 1:8,
                 group_count = 100L, observations_per_group = 20L,
                 overfitted_mcmc_components = 8L),
    full = list(replicates = 1:100, candidate_components = 1:8,
                group_count = 100L, observations_per_group = 20L,
                overfitted_mcmc_components = 8L))
  structure(list(
    schema_version = "paper-e2-protocol-1", experiment = "E2",
    profile = profile,
    data_generating = list(
      component_count = 3L, group_count = 100L,
      observations_per_group = 20L, coefficient_count = 5L,
      minimum_distance = c(1.5, .5)),
    inference = list(
      candidate_components = 1:8,
      cavi_starts = 10L, svi_starts = 10L,
      overfitted_mcmc_components = 8L),
    selection = list(
      variational = c("maximum_elbo", "maximum_mlpd"),
      mcmc = "posterior_mode_occupied_components"),
    target_replicates = 100L,
    execution = execution),
    class = c("paper_experiment_protocol", "list"))
}

scenario_e2 <- function(minimum_distance, replicate = 1L,
                        profile = c("smoke", "pilot", "full"),
                        base_seed = 310000L) {
  profile <- .paper_profile(profile)
  protocol <- e2_protocol(profile)
  if (!is.numeric(minimum_distance) || length(minimum_distance) != 1L ||
      !minimum_distance %in% protocol$data_generating$minimum_distance)
    stop("minimum_distance must be 1.5 or 0.5")
  replicate <- .paper_positive_integer(replicate, "replicate")
  distance_index <- match(minimum_distance,
                          protocol$data_generating$minimum_distance)
  seed <- as.integer(base_seed + 100000L * distance_index + replicate)
  K <- protocol$data_generating$component_count
  p <- protocol$data_generating$coefficient_count
  model <- chlrm_model(K, p)
  structure(list(
    schema_version = "paper-e2-scenario-1", experiment = "E2",
    profile = profile, replicate = replicate, seed = seed,
    minimum_distance = minimum_distance,
    model = model,
    generator = conditional_dgp(model, list(
      component_coefficients = make_review_component_coefficients(
        K, p, minimum_distance),
      component_variances = rep(1, K),
      mixture_weights = c(.4, .3, .3))),
    design = simulation_design(
      "multivariate_normal", include_intercept = TRUE,
      mean = rep(0, p - 1L), covariance = diag(p - 1L)),
    sample = simulation_sample(
      group_count = protocol$execution$group_count,
      observations_per_group = protocol$execution$observations_per_group,
      test_observations_per_group = protocol$execution$observations_per_group,
      test_group_mode = "same"),
    protocol = protocol), class = c("paper_e2_scenario", "list"))
}

e3_protocol <- function(profile = c("smoke", "pilot", "full")) {
  profile <- .paper_profile(profile)
  full_grid <- expand.grid(
    minibatch_size = c(25L, 50L, 100L, 200L),
    exponent = c(.6, .7, .8), offset = c(1, 10, 100),
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  execution <- switch(profile,
    smoke = list(group_count = 30L, seeds = 1L,
                 grid = full_grid[full_grid$minibatch_size == 25L &
                                    full_grid$exponent == .7 &
                                    full_grid$offset == 10, , drop = FALSE]),
    pilot = list(group_count = 1000L, seeds = 1:3,
                 grid = full_grid),
    full = list(group_count = 1000L, seeds = 1:20,
                grid = full_grid))
  structure(list(
    schema_version = "paper-e3-protocol-1", experiment = "E3",
    profile = profile,
    data_generating = list(
      component_count = 3L, group_count = 1000L,
      group_size_support = 10:30, coefficient_count = 10L,
      predictor_correlation = .8, minimum_distance = .75),
    full_grid = full_grid, target_seeds = 1:20,
    reference_method = "CAVI",
    outcomes = c("final_elbo", "mlpd", "ari", "cpe", "runtime",
                 "iterations", "converged"),
    execution = execution),
    class = c("paper_experiment_protocol", "list"))
}

scenario_e3 <- function(profile = c("smoke", "pilot", "full"),
                        seed = 410001L) {
  profile <- .paper_profile(profile)
  protocol <- e3_protocol(profile)
  group_count <- protocol$execution$group_count
  restore <- .paper_seed_scope(seed)
  on.exit(restore(), add = TRUE)
  group_sizes <- sample(protocol$data_generating$group_size_support,
                        group_count, replace = TRUE)
  K <- protocol$data_generating$component_count
  p <- protocol$data_generating$coefficient_count
  correlation <- protocol$data_generating$predictor_correlation
  model <- chlrm_model(K, p)
  structure(list(
    schema_version = "paper-e3-scenario-1", experiment = "E3",
    profile = profile, seed = as.integer(seed), model = model,
    generator = conditional_dgp(model, list(
      component_coefficients = make_review_component_coefficients(
        K, p, protocol$data_generating$minimum_distance),
      component_variances = rep(1, K),
      mixture_weights = rep(1 / K, K))),
    design = simulation_design(
      "multivariate_normal", include_intercept = TRUE,
      mean = rep(0, p - 1L),
      covariance = .ar1_covariance(p - 1L, correlation)),
    sample = simulation_sample(
      group_count = group_count, observations_per_group = group_sizes,
      test_observations_per_group = group_sizes, test_group_mode = "same"),
    protocol = protocol), class = c("paper_e3_scenario", "list"))
}

e4_protocol <- function(profile = c("smoke", "pilot", "full")) {
  profile <- .paper_profile(profile)
  execution <- switch(profile,
    smoke = list(replicates = 1L, nested_group_counts = c(10L, 20L)),
    pilot = list(replicates = 1:5,
                 nested_group_counts = c(50L, 200L, 1000L, 5000L)),
    full = list(replicates = 1:20,
                nested_group_counts = c(50L, 200L, 1000L, 5000L)))
  structure(list(
    schema_version = "paper-e4-protocol-1", experiment = "E4",
    profile = profile,
    data_generating = list(
      component_count = 3L, maximum_group_count = 5000L,
      group_size_support = 10:30, coefficient_count = 10L,
      minimum_distance = 1),
    nested_group_counts = c(50L, 200L, 1000L, 5000L),
    target_replicates = 20L,
    methods = c("CAVI", "SVI", "MCMC_when_diagnostics_pass"),
    svi_minibatch_rule = "min(100, group_count)",
    outcomes = c("runtime", "mlpd", "ari", "cpe"),
    scaling_model = "log(runtime) = intercept + slope * log(observations)",
    execution = execution),
    class = c("paper_experiment_protocol", "list"))
}

scenario_e4 <- function(profile = c("smoke", "pilot", "full"),
                        replicate = 1L, base_seed = 510000L) {
  profile <- .paper_profile(profile)
  replicate <- .paper_positive_integer(replicate, "replicate")
  protocol <- e4_protocol(profile)
  maximum_groups <- max(protocol$execution$nested_group_counts)
  seed <- as.integer(base_seed + replicate)
  restore <- .paper_seed_scope(seed)
  on.exit(restore(), add = TRUE)
  group_sizes <- sample(protocol$data_generating$group_size_support,
                        maximum_groups, replace = TRUE)
  K <- protocol$data_generating$component_count
  p <- protocol$data_generating$coefficient_count
  model <- chlrm_model(K, p)
  structure(list(
    schema_version = "paper-e4-scenario-1", experiment = "E4",
    profile = profile, replicate = replicate, seed = seed, model = model,
    generator = conditional_dgp(model, list(
      component_coefficients = make_review_component_coefficients(
        K, p, protocol$data_generating$minimum_distance),
      component_variances = rep(1, K),
      mixture_weights = rep(1 / K, K))),
    design = simulation_design(
      "independent_normal", include_intercept = TRUE),
    sample = simulation_sample(
      group_count = maximum_groups, observations_per_group = group_sizes,
      test_observations_per_group = group_sizes, test_group_mode = "same"),
    protocol = protocol), class = c("paper_e4_scenario", "list"))
}

e5_protocol <- function(profile = c("smoke", "pilot", "full")) {
  profile <- .paper_profile(profile)
  execution <- switch(profile,
    smoke = list(candidate_components = 1:2, evaluation_draws = 4L,
                 overfitted_mcmc_components = 4L),
    pilot = list(candidate_components = 1:4, evaluation_draws = 250L,
                 overfitted_mcmc_components = 8L),
    full = list(candidate_components = 1:12, evaluation_draws = 5000L,
                overfitted_mcmc_components = 23L))
  structure(list(
    schema_version = "paper-e5-protocol-1", experiment = "E5",
    profile = profile, dataset = "Farms",
    preprocessing = list(
      standardize_predictor = TRUE, standardize_response = TRUE,
      learn_within_training_fold = TRUE,
      report_parameters_on_original_scale = TRUE),
    selection = list(
      initial_components = 1:12,
      extend_to = 18L, extend_when_boundary_selected = TRUE,
      mcmc_chains = 4L, cavi_starts = 20L, svi_starts = 20L,
      criteria = c("elbo", "within_group_cv_mlpd", "clustering_stability")),
    cross_validation = list(
      folds = 5L, scheme = "one_observation_per_group",
      predictive_target = "observed_groups"),
    overfitted_mcmc_components = 23L,
    waic_unit = "group", waic_conditioning = "group_marginal",
    waic_draws = 5000L, posterior_predictive_draws = 5000L,
    execution = execution),
    class = c("paper_experiment_protocol", "list"))
}
