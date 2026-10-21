expected_exports <- c(
  "advance_mcmc_state", "append_mcmc_history", "as_log_lik_matrix", "as_yrep", "cavi_control", "cavi_stopping",
  "chlrm_model", "chlrm_prior", "classification_uncertainty",
  "classify_groups", "co_clustering", "components", "conditional_dgp",
  "create_mcmc_history", "cross_validate_model", "cross_validation_plan", "cross_validation_split",
  "data_scaling", "estimate_fit_resources", "evaluate_in_sample_fit",
  "fit_model", "fit_orchestration", "fit_scaling", "group_deviations",
  "group_parameters", "grouped_regression_data", "hierarchical_dgp",
  "hlrm_model", "hlrm_prior", "initialization_from_prior",
  "initialization_kmeans", "initialization_ols", "initialization_values",
  "log_likelihood", "lrm_model", "lrm_prior", "mcmc_control", "mcmc_reporting_view", "mcmc_restart_state", "mcmc_history_draws", "mcmc_history_restart_state",
  "mcmc_stopping", "population_parameters", "posterior_draws",
  "posterior_predict", "posterior_predictive_check", "power_schedule",
  "pp_check_from_prediction", "ppc_discrepancies", "predictive_criteria",
  "predictive_log_likelihood", "predictive_performance", "predictive_score",
  "prior_moments", "prior_predictive_dgp", "psis_loo",
  "read_mcmc_history", "read_mcmc_checkpoint", "read_mcmc_restart_state", "relabel_fit", "responsibilities", "resume_mcmc",
  "save_mcmc_checkpoint", "save_mcmc_restart_state", "scaling_map", "simulate_data",
  "simulation_design", "simulation_sample", "summarize_posterior",
  "svi_control", "svi_stopping", "tidy_parameters", "tidy_truth",
  "transform_data", "transform_prior", "variational_initialization_set",
  "waic_from_loglik"
)

expected_s3_methods <- c(
  "as_draws.vb_model_fit", "as_draws.vb_posterior_draws",
  "as_log_lik_matrix.vb_log_likelihood",
  "as_yrep.vb_posterior_predictive_check", "as_yrep.vb_prediction_result",
  "as.matrix.vb_log_likelihood", "classify_groups.vb_model_fit",
  "co_clustering.vb_model_fit", "coef.vb_model_fit",
  "components.vb_model_fit", "evaluate_in_sample_fit.vb_model_fit",
  "fitted.vb_model_fit", "group_deviations.vb_model_fit",
  "group_parameters.vb_model_fit", "nobs.vb_model_fit",
  "population_parameters.vb_model_fit", "posterior_draws.vb_model_fit",
  "posterior_predict.vb_model_fit",
  "posterior_predictive_check.vb_model_fit", "predict.vb_model_fit",
  "plot.vb_fit_orchestration", "plot.vb_log_likelihood",
  "plot.vb_model_fit", "plot.vb_posterior_predictive_check",
  "plot.vb_posterior_draws", "plot.vb_prediction_result",
  "plot.vb_simulation_result",
  "predictive_criteria.vb_model_fit", "print.vb_fit_resource_estimate",
  "print.vb_model_fit", "residuals.vb_model_fit",
  "responsibilities.vb_model_fit", "summarize_posterior.vb_model_fit",
  "summary.vb_model_fit", "tidy_parameters.vb_model_fit",
  "tidy_parameters.vb_posterior_summary"
)

test_that("the public namespace matches the frozen pre-1.0 manifest", {
  expect_setequal(
    getNamespaceExports("VariationalBayesHLRM"),
    expected_exports
  )
})

test_that("registered S3 methods match the frozen pre-1.0 manifest", {
  registered <- getNamespaceInfo(
    asNamespace("VariationalBayesHLRM"), "S3methods"
  )
  method_names <- paste(registered[, 1L], registered[, 2L], sep = ".")
  expect_setequal(method_names, expected_s3_methods)
})
