library(VariationalBayesHLRM)
source("thesis/final/code/R/reporting.R")

out <- ensure_dir(file.path(
  repo_root(), "thesis", "final", "code", "runs", "05-select-k-chlrm", "output"))
farms <- read.delim(file.path(repo_root(), "data", "farm_trees.txt"))
d <- grouped_regression_data(farms, "size", "N", "farm")
profile <- thesis_execution_profile()
if (profile == "smoke") {
  keep <- seq_len(min(4L, length(d$y)))
  d$X <- d$X[keep]
  d$y <- d$y[keep]
  d$mapping <- d$mapping[keep, , drop = FALSE]
}
Ks <- if (profile == "smoke") 1:2 else 1:8
scaling <- data_scaling()
logo_plan <- cross_validation_plan(
  d, scheme = "leave_one_group_out", seed = 5001)

prior_for_k <- function(K) thesis_hierarchical_prior(
  "CHLRM", coefficient_mean = c(100, 0),
  coefficient_covariance = diag(c(50^2, 5^2)),
  coefficient_dispersion_scale = diag(c(20^2, 1^2)),
  component_count = K)

svi_control_for_k <- function(K, starts, epochs, group_count = length(d$y)) {
  if (profile == "smoke") return(svi_control(
    max_iterations = 4, minibatch_size = min(4L, group_count),
    learning_rate = power_schedule(19, .7), evaluation_draws = 4,
    stopping = svi_stopping(
      objective_tolerance = 1e20, parameter_tolerance = 1e20,
      min_iterations = 2, patience = 1, check_every = 1,
      smoothing_window = 1),
    starts = starts,
    require_converged = if (starts > 1L) FALSE else NULL,
    max_memory_mb = 256))
  groups <- group_count
  batch <- max(K, ceiling(.5 * groups))
  iterations <- ceiling(epochs * groups / batch)
  svi_control(
    max_iterations = iterations,
    minibatch_size = batch,
    learning_rate = power_schedule(.05^(-1 / .7) - 1, .7),
    evaluation_draws = 100,
    stopping = svi_stopping(
      objective_tolerance = 1e-4,
      min_iterations = min(100L, iterations), patience = 5,
      check_every = max(1L, ceiling(groups / batch)), smoothing_window = 5),
    starts = starts, require_converged = FALSE, max_memory_mb = 2048)
}

full_rows <- vector("list", length(Ks))
full_evaluations <- vector("list", length(Ks))
best_fits <- vector("list", length(Ks)); names(best_fits) <- as.character(Ks)
logo_results <- vector("list", length(Ks)); names(logo_results) <- as.character(Ks)

for (K in Ks) {
  model <- chlrm_model(K, 2)
  prior <- prior_for_k(K)
  fit <- fit_model(
    model, d, prior, "svi", initialization_kmeans(),
    svi_control_for_k(K, starts = if (profile == "smoke") 2 else 10,
                      epochs = if (profile == "smoke") 2 else 300),
    scaling = scaling, seed = 5000L + K, verbose = FALSE)
  orchestration <- fit_orchestration(fit)
  start_summary <- orchestration$summary
  classification <- classify_groups(fit)
  probabilities <- responsibilities(fit)
  effective_sizes <- colSums(probabilities)
  fit_description <- evaluate_in_sample_fit(fit)
  criteria <- predictive_criteria(
    fit, n_draws = if (profile == "smoke") 4 else 150,
    seed = 7000L + K,
    return_log_likelihood = TRUE)

  full_rows[[K]] <- data.frame(
    K = K,
    selected_start = orchestration$selected_start,
    ELBO_best = max(start_summary$final_elbo, na.rm = TRUE),
    ELBO_mean = mean(start_summary$final_elbo, na.rm = TRUE),
    ELBO_sd = if (nrow(start_summary) > 1L)
      stats::sd(start_summary$final_elbo, na.rm = TRUE) else 0,
    converged_starts = sum(start_summary$converged),
    WAIC_group_marginal = criteria$criteria$waic_group_marginal$WAIC,
    WAIC_group_marginal_se = criteria$criteria$waic_group_marginal$se,
    MSE = fit_description$statistics$mse,
    R2 = fit_description$statistics$r2,
    effective_components = sum(effective_sizes >= 1),
    min_effective_size = min(effective_sizes),
    normalized_entropy = if (K == 1L) 0 else
      mean(classification$entropy) / log(K),
    time_total_starts = sum(start_summary$elapsed_seconds, na.rm = TRUE),
    profile = profile, stringsAsFactors = FALSE)
  full_evaluations[[as.character(K)]] <- list(
    in_sample = fit_description, predictive = criteria,
    starts = start_summary)
  best_fits[[as.character(K)]] <- fit

  logo_started <- proc.time()[["elapsed"]]
  logo <- cross_validate_model(
    model, d, prior, "svi", initialization_kmeans(),
    svi_control_for_k(K, starts = if (profile == "smoke") 1 else 3,
                      epochs = if (profile == "smoke") 2 else 200,
                      group_count = length(d$y) - 1L), logo_plan,
    scaling = scaling, n_draws = if (profile == "smoke") 4 else 250,
    seed = 900000L + K,
    max_fits = logo_plan$fold_count, retain_fits = FALSE, verbose = FALSE)
  logo$metadata$total_seconds <- proc.time()[["elapsed"]] - logo_started
  logo_results[[as.character(K)]] <- logo
  message(sprintf(
    "K=%d: ELBO=%.2f WAIC=%.2f LOGO-ELPD=%.2f effective=%d",
    K, full_rows[[K]]$ELBO_best,
    full_rows[[K]]$WAIC_group_marginal,
    logo$statistics$elpd, full_rows[[K]]$effective_components))
}

full <- do.call(rbind, full_rows)
write.csv(full, file.path(out, "full-data-results.csv"), row.names = FALSE)

logo_pointwise <- do.call(rbind, lapply(names(logo_results), function(name) {
  values <- logo_results[[name]]$pointwise
  values$K <- as.integer(name)
  values
}))
write.csv(logo_pointwise, file.path(out, "logo-pointwise.csv"), row.names = FALSE)

logo_summary <- do.call(rbind, lapply(names(logo_results), function(name) {
  result <- logo_results[[name]]
  values <- result$pointwise$log_predictive_density
  data.frame(
    K = as.integer(name),
    LOGO_ELPD = result$statistics$elpd,
    LOGO_SE = sqrt(length(values) * stats::var(values)),
    mean_group_log_score = result$statistics$mean_log_predictive_density,
    total_fit_seconds = result$metadata$total_seconds)
}))
best_K_predictive <- logo_summary$K[which.max(logo_summary$LOGO_ELPD)]
best_scores <- logo_pointwise$log_predictive_density[
  logo_pointwise$K == best_K_predictive]
logo_summary$delta_from_best <- NA_real_
logo_summary$SE_delta <- NA_real_
for (index in seq_len(nrow(logo_summary))) {
  scores <- logo_pointwise$log_predictive_density[
    logo_pointwise$K == logo_summary$K[index]]
  difference <- scores - best_scores
  logo_summary$delta_from_best[index] <- sum(difference)
  logo_summary$SE_delta[index] <- sqrt(length(difference) * stats::var(difference))
}
eligible <- logo_summary$delta_from_best >= -logo_summary$SE_delta
selected_K <- min(logo_summary$K[eligible])
logo_summary$within_one_se <- eligible
logo_summary$selected_smallest_one_se <- logo_summary$K == selected_K
write.csv(logo_summary, file.path(out, "logo-summary.csv"), row.names = FALSE)

selection <- merge(full, logo_summary, by = "K")
write.csv(selection, file.path(out, "K-selection-summary.csv"), row.names = FALSE)
write.csv(selection[selection$K == selected_K, ],
          file.path(out, "recommended-K.csv"), row.names = FALSE)
saveRDS(list(profile = profile, full_data = full_evaluations,
             logo = logo_results),
        file.path(out, "evaluation-with-log-likelihood.rds"), compress = "xz")

png(file.path(out, "K-selection.png"), width = 1700, height = 1000, res = 160)
par(mfrow = c(1, 3), mar = c(4.2, 4.5, 2.8, 1))
plot(full$K, full$WAIC_group_marginal, type = "b", pch = 19, col = "#0072B2",
     xlab = "K", ylab = "WAIC marginal por finca", main = "Ajuste interno")
plot(logo_summary$K, logo_summary$LOGO_ELPD, type = "b", pch = 19,
     col = "#009E73", xlab = "K", ylab = "LOGO-ELPD",
     main = "Predicción de fincas nuevas")
arrows(logo_summary$K, logo_summary$LOGO_ELPD - logo_summary$LOGO_SE,
       logo_summary$K, logo_summary$LOGO_ELPD + logo_summary$LOGO_SE,
       angle = 90, code = 3, length = .05, col = "#009E73")
abline(v = selected_K, lty = 2, col = "#D55E00")
plot(full$K, full$effective_components, type = "b", pch = 19,
     ylim = c(0, max(Ks)), xlab = "K solicitado",
     ylab = "Componentes efectivos", main = "Ocupación posterior")
abline(0, 1, lty = 3, col = "grey60")
dev.off()

selected_psm <- co_clustering(best_fits[[as.character(selected_K)]])
plot_psm(selected_psm, file.path(out, "coclustering-selected-K.png"),
         sprintf("CHLRM: K=%d seleccionado", selected_K))
write.csv(selected_psm, file.path(out, "coclustering-selected-K.csv"),
          row.names = FALSE)

cat(sprintf(
  "\nPredictive-best K: %d; smallest-within-1-SE K: %d\n",
  best_K_predictive, selected_K))
print(selection)
