# Test-only fixtures built on the public simulation API. They keep individual
# numerical tests concise without reintroducing production simulator aliases.

.test_lrm_simulation <- function(sample_size, coefficients, residual_variance,
                                 covariate_range, seed = NULL) {
  q <- length(coefficients) - 1L
  bounds <- matrix(covariate_range, q, 2L, byrow = TRUE)
  simulate_data(
    conditional_dgp(lrm_model(length(coefficients)), list(
      coefficients = coefficients, residual_variance = residual_variance)),
    simulation_design("grid", bounds = bounds,
                      jitter_sd = .01 * (bounds[, 2L] - bounds[, 1L])),
    simulation_sample(sample_size = sample_size), seed = seed)
}

.test_hlrm_simulation <- function(group_sizes, group_coefficients,
                                  group_variances,
                                  binary_probabilities = numeric(),
                                  covariate_mean = 0, covariate_sd = 1,
                                  seed = NULL) {
  group_coefficients <- as.matrix(group_coefficients)
  simulate_data(
    conditional_dgp(hlrm_model(ncol(group_coefficients)), list(
      group_coefficients = group_coefficients,
      group_variances = group_variances)),
    simulation_design(
      if (length(binary_probabilities)) "mixed" else "independent_normal",
      mean = covariate_mean, sd = covariate_sd,
      binary_probabilities = binary_probabilities, binary_position = "first"),
    simulation_sample(group_count = nrow(group_coefficients),
                      observations_per_group = group_sizes), seed = seed)
}

.test_hlrm_hierarchy <- function(group_sizes, population_coefficients,
                                 population_covariance, variance_scale,
                                 residual_degrees_of_freedom = 4,
                                 covariate_mean = 0, covariate_sd = 1,
                                 seed = NULL) {
  simulate_data(
    hierarchical_dgp(hlrm_model(length(population_coefficients)), list(
      population_coefficients = population_coefficients,
      population_covariance = population_covariance,
      variance_scale = variance_scale,
      residual_degrees_of_freedom = residual_degrees_of_freedom)),
    simulation_design("independent_normal", mean = covariate_mean,
                      sd = covariate_sd),
    simulation_sample(group_count = length(group_sizes),
                      observations_per_group = group_sizes), seed = seed)
}

.test_component_coefficients <- function(component_count, coefficient_count,
                                         separation = 3) {
  if (component_count == 1L) return(matrix(0, 1L, coefficient_count))
  angles <- 2 * pi * (seq_len(component_count) - 1L) / component_count
  basis <- cbind(cos(angles), sin(angles))
  basis <- basis * separation / min(stats::dist(basis))
  out <- matrix(0, component_count, coefficient_count)
  usable <- min(2L, max(0L, coefficient_count - 1L))
  if (usable) out[, 1L + seq_len(usable)] <- basis[, seq_len(usable), drop = FALSE]
  if (coefficient_count == 1L)
    out[, 1L] <- seq(-separation / 2, separation / 2,
                     length.out = component_count)
  out
}

.test_chlrm_simulation <- function(
    group_count, group_sizes, component_count, coefficient_count,
    covariate_covariance = diag(coefficient_count - 1L),
    mixture_weights = rep(1 / component_count, component_count),
    component_variances = 1, test_group_sizes = 20L,
    component_coefficients = NULL, allocations = NULL, seed = NULL) {
  if (is.null(component_coefficients)) component_coefficients <-
    .test_component_coefficients(component_count, coefficient_count)
  simulate_data(
    conditional_dgp(chlrm_model(component_count, coefficient_count), list(
      component_coefficients = component_coefficients,
      component_variances = component_variances,
      mixture_weights = mixture_weights,
      allocations = if (is.null(allocations)) NULL else allocations)[
        c("component_coefficients", "component_variances", "mixture_weights",
          if (!is.null(allocations)) "allocations")]),
    simulation_design("multivariate_normal",
                      covariance = covariate_covariance),
    simulation_sample(
      group_count = group_count, observations_per_group = group_sizes,
      test_observations_per_group = test_group_sizes), seed = seed)
}

.test_chlrm_hierarchy <- function(
    group_sizes, component_count, population_coefficients,
    population_covariance, variance_scale,
    residual_degrees_of_freedom = 4,
    mixture_weights = rep(1 / component_count, component_count),
    covariate_mean = 0, covariate_sd = 1, seed = NULL) {
  simulate_data(
    hierarchical_dgp(chlrm_model(
      component_count, length(population_coefficients)), list(
        population_coefficients = population_coefficients,
        population_covariance = population_covariance,
        variance_scale = variance_scale,
        residual_degrees_of_freedom = residual_degrees_of_freedom,
        mixture_weights = mixture_weights)),
    simulation_design("independent_normal", mean = covariate_mean,
                      sd = covariate_sd),
    simulation_sample(group_count = length(group_sizes),
                      observations_per_group = group_sizes), seed = seed)
}

.test_chlrm_native_prior <- function(p, K) {
  validate_hlrm_hyperparams(list(
    Mu0 = rep(0, p), Lambda0 = 100 * diag(p), n0 = p + 2,
    S0 = diag(p), nu0 = 1, a0 = 1, b0 = 1,
    alpha0 = rep(1 / K, K), type = "test_fixture"),
    p, m = 1L, model = "CHLRM", K = K)
}
