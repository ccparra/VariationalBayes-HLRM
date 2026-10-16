library(VariationalBayesHLRM)
source("thesis/final/code/R/reporting.R")

out <- ensure_dir(file.path(
  repo_root(), "thesis", "final", "code", "runs", "01-farms-hlrm", "output"))
farms <- read.delim(file.path(repo_root(), "data", "farm_trees.txt"))
d <- grouped_regression_data(farms, "size", "N", "farm")
model <- hlrm_model(2)
prior <- thesis_hierarchical_prior(
  "HLRM", coefficient_mean = c(100, 0),
  coefficient_covariance = diag(c(50^2, 5^2)),
  coefficient_dispersion_scale = diag(c(20^2, 1^2)))
scaling <- data_scaling()
profile <- thesis_execution_profile()
controls <- if (profile == "smoke") thesis_smoke_controls(length(d$y)) else
  list(
    evaluation_draws = 150L,
    mcmc = mcmc_control(
      warmup = 600, post_warmup_draws = 1200, chains = 4,
      evaluation_draws = 150,
      stopping = mcmc_stopping(
        min_post_warmup_draws = 500, check_every = 250,
        ess_bulk_min = 400, ess_tail_min = 400),
      max_memory_mb = 2048),
    cavi = cavi_control(
      max_iterations = 3000, evaluation_draws = 150, starts = 3,
      require_converged = FALSE, max_memory_mb = 1024),
    svi = svi_control(
      max_iterations = 800, minibatch_size = 8,
      learning_rate = power_schedule(
        offset = .05^(-1 / .7) - 1, exponent = .7),
      evaluation_draws = 150,
      stopping = svi_stopping(
        min_iterations = 200, patience = 5, check_every = 20,
        smoothing_window = 5),
      starts = 3, require_converged = FALSE, max_memory_mb = 1024))

mcmc <- suppressWarnings(fit_model(
  model, d, prior, "mcmc", initialization_ols(),
  controls$mcmc,
  scaling = scaling, seed = 1101, verbose = FALSE))
cavi <- fit_model(
  model, d, prior, "cavi", initialization_ols(),
  controls$cavi,
  scaling = scaling, seed = 1102, verbose = FALSE)
svi <- fit_model(
  model, d, prior, "svi", initialization_ols(),
  controls$svi,
  scaling = scaling, seed = 1103, verbose = FALSE)

fits <- list(MCMC = mcmc, CAVI = cavi, SVI = svi)
summaries <- lapply(fits, summarize_posterior)
in_sample <- lapply(fits, evaluate_in_sample_fit)
criteria <- Map(function(fit, seed) predictive_criteria(
  fit, n_draws = controls$evaluation_draws, seed = seed,
  return_log_likelihood = TRUE),
  fits, 1111:1113)
ppc <- Map(function(fit, seed) posterior_predictive_check(
  fit, n_draws = controls$evaluation_draws, seed = seed,
  return_yrep = TRUE),
  fits, 1121:1123)

metrics <- do.call(rbind, lapply(names(fits), function(method)
  metric_row("farms_hlrm", method, fits[[method]],
             in_sample[[method]], criteria[[method]])))
metrics$profile <- profile
write.csv(metrics, file.path(out, "metrics.csv"), row.names = FALSE)
for (method in names(fits))
  write_parameter_summary(
    summaries[[method]], method,
    file.path(out, paste0("parameters-", tolower(method), ".csv")))

plot_elbo_comparison(
  fits[c("CAVI", "SVI")], file.path(out, "elbo.png"),
  "HLRM fincas: CAVI frente a SVI")
plot_regression_lines(
  farms, summaries, file.path(out, "regressions-by-farm.png"),
  "Regresiones HLRM")
for (method in names(fits))
  plot_predictive_density(
    d$y, ppc[[method]]$y_rep,
    file.path(out, paste0("predictive-", tolower(method), ".png")),
    paste("HLRM:", method))

saveRDS(list(
  profile = profile, metrics = metrics,
  parameters = lapply(summaries, `[[`, "parameters"),
  convergence = lapply(fits, `[[`, "diagnostics")),
  file.path(out, "results-light.rds"), compress = "xz")
saveRDS(Map(function(observed, predictive, check)
  list(in_sample = observed, predictive = predictive, ppc = check),
  in_sample, criteria, ppc),
  file.path(out, "evaluation-with-log-likelihood.rds"), compress = "xz")
