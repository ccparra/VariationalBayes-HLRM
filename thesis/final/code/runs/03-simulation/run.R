library(VariationalBayesHLRM)
source("thesis/final/code/R/reporting.R")

out <- ensure_dir(file.path(
  repo_root(), "thesis", "final", "code", "runs", "03-simulation", "output"))
profile <- thesis_execution_profile()

align_chlrm <- function(hat, component_coefficients) {
  K <- nrow(component_coefficients)
  cost <- outer(seq_len(K), seq_len(K), Vectorize(function(i, j)
    sum((hat$parameters$component_coefficients$mean[i, ] -
         component_coefficients[j, ])^2)))
  estimated_to_true <- as.integer(clue::solve_LSAP(cost))
  ord <- match(seq_len(K), estimated_to_true)
  list(beta = hat$parameters$component_coefficients$mean[ord, , drop = FALSE],
       sigma2 = hat$parameters$component_variances$mean[ord],
       z = estimated_to_true[hat$allocations$labels])
}

# HLRM: toda la jerarquía tiene verdad conocida.
dh_simulation <- simulate_data(
  hierarchical_dgp(hlrm_model(2), list(
    population_coefficients = c(1, .6),
    population_covariance = matrix(c(2, .2, .2, .4), 2),
    variance_scale = 1, residual_degrees_of_freedom = 6)),
  simulation_design("independent_normal"),
  simulation_sample(group_count = 40L, observations_per_group = 30L),
  seed = 1301)
dh <- dh_simulation$train
truth_h <- dh_simulation$truth
prior_h <- thesis_hierarchical_prior(
  "HLRM", coefficient_mean = c(0, 0),
  coefficient_covariance = 100 * diag(2),
  coefficient_dispersion_scale = diag(2))
controls_h <- if (profile == "smoke")
  thesis_smoke_controls(length(dh$y)) else list(
    evaluation_draws = 100L,
    cavi = cavi_control(
      max_iterations = 2500, evaluation_draws = 100,
      starts = 3, require_converged = FALSE),
    svi = svi_control(
      max_iterations = 1000, minibatch_size = 10,
      learning_rate = power_schedule(.05^(-1 / .7) - 1, .7),
      evaluation_draws = 100,
      stopping = svi_stopping(
        min_iterations = 200, check_every = 20, smoothing_window = 5),
      starts = 3, require_converged = FALSE))
fit_vh <- fit_model(
  hlrm_model(2), dh, prior_h, "cavi", initialization_ols(),
  controls_h$cavi,
  seed = 1302, verbose = FALSE)
hat_vh <- summarize_posterior(fit_vh)
fit_eval_vh <- evaluate_in_sample_fit(fit_vh)
criteria_vh <- predictive_criteria(
  fit_vh, n_draws = controls_h$evaluation_draws, seed = 1304,
  return_log_likelihood = TRUE)
fit_sh <- fit_model(
  hlrm_model(2), dh, prior_h, "svi", initialization_ols(),
  controls_h$svi,
  seed = 1303, verbose = FALSE)
hat_sh <- summarize_posterior(fit_sh)
fit_eval_sh <- evaluate_in_sample_fit(fit_sh)
criteria_sh <- predictive_criteria(
  fit_sh, n_draws = controls_h$evaluation_draws, seed = 1305,
  return_log_likelihood = TRUE)

metrics_h <- rbind(
  metric_row("simulation_hlrm", "CAVI", fit_vh, fit_eval_vh, criteria_vh),
  metric_row("simulation_hlrm", "SVI", fit_sh, fit_eval_sh, criteria_sh))
metrics_h$beta_local_RMSE <- c(rmse(hat_vh$parameters$group_coefficients$mean, truth_h$group_coefficients),
                               rmse(hat_sh$parameters$group_coefficients$mean, truth_h$group_coefficients))
metrics_h$variance_local_RMSE <- c(rmse(hat_vh$parameters$group_variances$mean, truth_h$group_variances),
                                   rmse(hat_sh$parameters$group_variances$mean, truth_h$group_variances))
metrics_h$beta_global_RMSE <- c(rmse(hat_vh$parameters$population_coefficients$mean, truth_h$population_coefficients),
                                rmse(hat_sh$parameters$population_coefficients$mean, truth_h$population_coefficients))
metrics_h$beta_component_RMSE <- NA_real_
metrics_h$variance_component_RMSE <- NA_real_
metrics_h$allocation_ARI <- NA_real_

# CHLRM: componentes, pesos y partición con verdad conocida.
dc_simulation <- simulate_data(
  hierarchical_dgp(chlrm_model(3L, 2L), list(
    population_coefficients = c(0, .3),
    population_covariance = matrix(c(12, 0, 0, .5), 2),
    variance_scale = 1, residual_degrees_of_freedom = 8,
    mixture_weights = c(.3, .4, .3))),
  simulation_design("independent_normal"),
  simulation_sample(group_count = 90L, observations_per_group = 30L),
  seed = 1310)
dc <- dc_simulation$train
truth_c <- dc_simulation$truth
prior_c <- thesis_hierarchical_prior(
  "CHLRM", coefficient_mean = c(0, 0),
  coefficient_covariance = 100 * diag(2),
  coefficient_dispersion_scale = diag(2), component_count = 3)
controls_c <- if (profile == "smoke")
  thesis_smoke_controls(length(dc$y)) else list(
    evaluation_draws = 100L,
    cavi = cavi_control(
      max_iterations = 3000, evaluation_draws = 100,
      starts = 5, require_converged = FALSE),
    svi = svi_control(
      max_iterations = 1200, minibatch_size = 20,
      learning_rate = power_schedule(.05^(-1 / .7) - 1, .7),
      evaluation_draws = 100,
      stopping = svi_stopping(
        min_iterations = 200, check_every = 20, smoothing_window = 5),
      starts = 5, require_converged = FALSE))
fit_vc <- fit_model(
  chlrm_model(3, 2), dc, prior_c, "cavi", initialization_kmeans(),
  controls_c$cavi,
  seed = 1311, verbose = FALSE)
hat_vc <- summarize_posterior(fit_vc)
fit_eval_vc <- evaluate_in_sample_fit(fit_vc)
criteria_vc <- predictive_criteria(
  fit_vc, n_draws = controls_c$evaluation_draws, seed = 1313,
  return_log_likelihood = TRUE)
fit_sc <- fit_model(
  chlrm_model(3, 2), dc, prior_c, "svi", initialization_kmeans(),
  controls_c$svi,
  seed = 1312, verbose = FALSE)
hat_sc <- summarize_posterior(fit_sc)
fit_eval_sc <- evaluate_in_sample_fit(fit_sc)
criteria_sc <- predictive_criteria(
  fit_sc, n_draws = controls_c$evaluation_draws, seed = 1314,
  return_log_likelihood = TRUE)
av <- align_chlrm(hat_vc, truth_c$component_coefficients)
as <- align_chlrm(hat_sc, truth_c$component_coefficients)

metrics_c <- rbind(
  metric_row("simulation_chlrm", "CAVI", fit_vc, fit_eval_vc, criteria_vc),
  metric_row("simulation_chlrm", "SVI", fit_sc, fit_eval_sc, criteria_sc))
metrics_c$beta_component_RMSE <- c(rmse(av$beta, truth_c$component_coefficients),
                                   rmse(as$beta, truth_c$component_coefficients))
metrics_c$variance_component_RMSE <- c(rmse(av$sigma2, truth_c$component_variances),
                                       rmse(as$sigma2, truth_c$component_variances))
metrics_c$allocation_ARI <- c(ari(av$z, truth_c$allocations), ari(as$z, truth_c$allocations))
metrics_c$beta_local_RMSE <- NA_real_
metrics_c$variance_local_RMSE <- NA_real_
metrics_c$beta_global_RMSE <- NA_real_

metrics_c <- metrics_c[names(metrics_h)]
metrics <- rbind(metrics_h, metrics_c)
metrics$profile <- profile
write.csv(metrics, file.path(out, "metrics.csv"), row.names = FALSE)
write_parameter_summary(hat_vh, "CAVI", file.path(out, "hlrm-parameters-cavi.csv"))
write_parameter_summary(hat_sh, "SVI", file.path(out, "hlrm-parameters-svi.csv"))
write_parameter_summary(hat_vc, "CAVI", file.path(out, "chlrm-parameters-cavi.csv"))
write_parameter_summary(hat_sc, "SVI", file.path(out, "chlrm-parameters-svi.csv"))

png(file.path(out, "recovery.png"), width = 1800, height = 900, res = 160)
par(mfrow = c(1, 3), mar = c(4.3, 4.3, 3, 1))
plot(truth_h$group_coefficients, hat_vh$parameters$group_coefficients$mean, pch = 16, col = "#0072B2",
     xlab = "Verdad", ylab = "Estimación", main = "HLRM: coeficientes locales"); abline(0, 1, lty = 2)
plot(truth_c$component_coefficients, av$beta, pch = 16, col = "#0072B2",
     xlab = "Verdad", ylab = "Estimación", main = "CHLRM: coeficientes CAVI"); abline(0, 1, lty = 2)
plot(truth_c$component_coefficients, as$beta, pch = 17, col = "#D55E00",
     xlab = "Verdad", ylab = "Estimación", main = "CHLRM: coeficientes SVI"); abline(0, 1, lty = 2)
dev.off()

plot_psm(co_clustering(fit_vc), file.path(out, "chlrm-coclustering-cavi.png"),
         "Simulación CHLRM: coclustering CAVI")
plot_psm(hat_sc$allocations$co_clustering, file.path(out, "chlrm-coclustering-svi.png"),
         "Simulación CHLRM: coclustering SVI")
saveRDS(list(profile = profile, metrics = metrics,
             truth = list(HLRM = truth_h, CHLRM = truth_c)),
        file.path(out, "results-light.rds"))
saveRDS(list(
  HLRM = list(
    CAVI = list(in_sample = fit_eval_vh, predictive = criteria_vh),
    SVI = list(in_sample = fit_eval_sh, predictive = criteria_sh)),
  CHLRM = list(
    CAVI = list(in_sample = fit_eval_vc, predictive = criteria_vc),
    SVI = list(in_sample = fit_eval_sc, predictive = criteria_sc))),
  file.path(out, "evaluation-with-log-likelihood.rds"), compress = "xz")
