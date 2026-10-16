library(VariationalBayesHLRM)
source("thesis/final/code/R/reporting.R")

out <- ensure_dir(file.path(
  repo_root(), "thesis", "final", "code", "runs", "02-farms-chlrm", "output"))
farms <- read.delim(file.path(repo_root(), "data", "farm_trees.txt"))
d <- grouped_regression_data(farms, "size", "N", "farm")
K <- 4L
model <- chlrm_model(K, 2)
prior <- thesis_hierarchical_prior(
  "CHLRM", coefficient_mean = c(100, 0),
  coefficient_covariance = diag(c(50^2, 5^2)),
  coefficient_dispersion_scale = diag(c(20^2, 1^2)),
  component_count = K)
scaling <- data_scaling()
profile <- thesis_execution_profile()
controls <- if (profile == "smoke") thesis_smoke_controls(length(d$y)) else
  list(
    evaluation_draws = 150L,
    mcmc = mcmc_control(
      warmup = 800, post_warmup_draws = 1600, chains = 4,
      evaluation_draws = 150,
      stopping = mcmc_stopping(
        min_post_warmup_draws = 500, check_every = 250,
        ess_bulk_min = 400, ess_tail_min = 400),
      max_memory_mb = 2048),
    cavi = cavi_control(
      max_iterations = 4000, evaluation_draws = 150, starts = 5,
      require_converged = FALSE, max_memory_mb = 1024),
    svi = svi_control(
      max_iterations = 1000, minibatch_size = 8,
      learning_rate = power_schedule(
        offset = .05^(-1 / .7) - 1, exponent = .7),
      evaluation_draws = 150,
      stopping = svi_stopping(
        min_iterations = 200, patience = 5, check_every = 20,
        smoothing_window = 5),
      starts = 5, require_converged = FALSE, max_memory_mb = 1024))

mcmc <- suppressWarnings(fit_model(
  model, d, prior, "mcmc", initialization_kmeans(),
  controls$mcmc,
  scaling = scaling, seed = 1201, verbose = FALSE))
cavi <- fit_model(
  model, d, prior, "cavi", initialization_kmeans(),
  controls$cavi,
  scaling = scaling, seed = 1202, verbose = FALSE)
svi <- fit_model(
  model, d, prior, "svi", initialization_kmeans(),
  controls$svi,
  scaling = scaling, seed = 1203, verbose = FALSE)

fits <- list(MCMC = mcmc, CAVI = cavi, SVI = svi)
summaries <- lapply(fits, summarize_posterior)
in_sample <- lapply(fits, evaluate_in_sample_fit)
criteria <- Map(function(fit, seed) predictive_criteria(
  fit, n_draws = controls$evaluation_draws, seed = seed,
  return_log_likelihood = TRUE),
  fits, 1211:1213)
ppc <- Map(function(fit, seed) posterior_predictive_check(
  fit, n_draws = controls$evaluation_draws, seed = seed,
  return_yrep = TRUE),
  fits, 1221:1223)

metrics <- do.call(rbind, lapply(names(fits), function(method)
  metric_row("farms_chlrm", method, fits[[method]],
             in_sample[[method]], criteria[[method]])))
metrics$profile <- profile
write.csv(metrics, file.path(out, "metrics.csv"), row.names = FALSE)
for (method in names(fits))
  write_parameter_summary(
    summaries[[method]], method,
    file.path(out, paste0("parameters-", tolower(method), ".csv")))
plot_elbo_comparison(
  fits[c("CAVI", "SVI")], file.path(out, "elbo.png"),
  "CHLRM fincas: CAVI frente a SVI")
for (method in names(fits)) {
  plot_psm(
    co_clustering(fits[[method]]),
    file.path(out, paste0("coclustering-", tolower(method), ".png")),
    paste(method, "coclustering"))
  plot_predictive_density(
    d$y, ppc[[method]]$y_rep,
    file.path(out, paste0("predictive-", tolower(method), ".png")),
    paste("CHLRM:", method))
}

classifications <- lapply(fits, classify_groups)
assignments <- data.frame(farm = d$mapping$group_label)
for (method in names(classifications)) {
  assignments[[method]] <- classifications[[method]]$component
  assignments[[paste0(method, "_uncertainty")]] <-
    classifications[[method]]$uncertainty
  assignments[[paste0(method, "_entropy")]] <-
    classifications[[method]]$entropy
}
write.csv(assignments, file.path(out, "assignments.csv"), row.names = FALSE)
saveRDS(list(
  profile = profile, metrics = metrics, assignments = assignments,
  parameters = lapply(summaries, `[[`, "parameters"),
  convergence = lapply(fits, `[[`, "diagnostics")),
  file.path(out, "results-light.rds"), compress = "xz")
saveRDS(Map(function(observed, predictive, check)
  list(in_sample = observed, predictive = predictive, ppc = check),
  in_sample, criteria, ppc),
  file.path(out, "evaluation-with-log-likelihood.rds"), compress = "xz")
