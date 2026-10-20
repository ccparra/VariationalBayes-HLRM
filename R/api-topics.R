#' Model, prior, and initialization specifications
#'
#' Construct validated model and prior specifications and explicit
#' initialization strategies. Priors and initial values remain separate so an
#' OLS start never changes the Bayesian model.
#'
#' @param model A model specification created by `lrm_model()`, `hlrm_model()`,
#'   or `chlrm_model()`.
#' @param coefficient_count Number of regression coefficients, including an
#'   explicit intercept when present.
#' @param component_count Number of CHLRM mixture components.
#' @param coefficient_mean Prior mean of the population coefficient vector.
#' @param coefficient_covariance Prior covariance of the population
#'   coefficient vector.
#' @param residual_degrees_of_freedom Degrees of freedom in the residual
#'   variance hierarchy.
#' @param residual_variance_scale Positive LRM inverse-gamma scale.
#' @param coefficient_dispersion_df Inverse-Wishart degrees of freedom.
#' @param coefficient_dispersion_scale Inverse-Wishart scale matrix.
#' @param residual_variance_scale_shape,residual_variance_scale_rate Shape and
#'   rate of the Gamma prior for the hierarchical variance scale.
#' @param mixture_concentration Positive CHLRM Dirichlet concentrations.
#' @param prior A prior specification.
#' @param jitter Positive regularization used only to construct an OLS or
#'   K-means initialization.
#' @param values Complete semantic initial values for a model and engine.
#' @param data Model data containing `X` and `y`.
#' @param initialization A base initialization specification used to construct
#'   every state in a reusable variational initialization set.
#' @param starts Number of independently generated variational initial states.
#' @param scaling Optional explicit data-scaling specification.
#' @param seed Required seed for reproducible state generation.
#' @return A validated typed specification, or a `prior-moments-1` object from
#'   `prior_moments()`.
#' @examples
#' model <- lrm_model(2)
#' prior <- lrm_prior(c(0, 0), 10 * diag(2), 4, 1)
#' initialization <- initialization_ols()
#' initialization_set <- variational_initialization_set(
#'   model, list(X = diag(2), y = c(0, 1)), prior,
#'   initialization, starts = 2, seed = 10
#' )
#' prior_moments(prior)
#' @name model-specifications
NULL

#' Inference and stopping controls
#'
#' Declare computational budgets, convergence rules, evaluation draw counts,
#' multichain or multistart settings, and an SVI learning-rate schedule.
#'
#' @param warmup Number of MCMC warm-up iterations per chain.
#' @param post_warmup_draws Maximum retained MCMC draws per chain.
#' @param thinning MCMC thinning interval.
#' @param chains Number of independent MCMC chains.
#' @param evaluation_draws Default number of downstream posterior draws.
#' @param stopping A method-specific stopping specification.
#' @param diagnostic_scope MCMC diagnostic scope, `"global"` or `"all"`.
#' @param max_memory_mb Positive preflight memory budget in MiB.
#' @param max_iterations Maximum CAVI or SVI iterations.
#' @param starts Number of independent variational starts.
#' @param require_converged Whether multistart selection requires convergence.
#' @param minibatch_size Number of observations for LRM or groups for
#'   HLRM/CHLRM sampled per SVI iteration.
#' @param learning_rate An SVI schedule such as `power_schedule()`.
#' @param min_post_warmup_draws Minimum retained draws before MCMC may stop.
#' @param rhat_max,ess_bulk_min,ess_tail_min MCMC diagnostic thresholds.
#' @param relative_mcse_mean_max,relative_mcse_sd_max Relative MCSE thresholds.
#' @param objective_tolerance,parameter_tolerance Variational stabilization
#'   tolerances; `NULL` disables the corresponding condition.
#' @param min_iterations Minimum variational iterations before stopping.
#' @param patience Required consecutive successful checks.
#' @param check_every Iterations between comparable objective evaluations.
#' @param smoothing_window Number of objective evaluations used for SVI
#'   smoothing.
#' @param offset,exponent Parameters of the power learning-rate schedule.
#' @return A validated inference, stopping, or learning-rate specification.
#' @examples
#' cavi_control(max_iterations = 1000)
#' svi_control(
#'   max_iterations = 5000, minibatch_size = 10,
#'   learning_rate = power_schedule(offset = 100, exponent = 0.7)
#' )
#' @name inference-controls
NULL

#' Fit models and extract posterior information
#'
#' Fit a typed model, summarize only applicable parameters, extract consistently
#' oriented draws, inspect orchestration provenance, and apply reversible data
#' scaling.
#'
#' @param model A model specification.
#' @param data Model data containing `X` and `y`.
#' @param prior A matching prior specification.
#' @param method One of `"mcmc"`, `"cavi"`, or `"svi"`.
#' @param initialization A compatible initialization specification.
#' @param control A method-specific inference control.
#' @param scaling Optional `data_scaling_spec`.
#' @param seed Optional reproducibility seed.
#' @param verbose Whether to display progress.
#' @param fit A fitted `vb_model_fit`.
#' @param n_draws Number of draws to return; `NULL` uses the fit default.
#' @param include_latent Whether to include CHLRM allocations.
#' @param max_output_mb Allocation budget for materialized draws.
#' @param draw_indices Optional strictly increasing positive indices into the
#'   combined stored MCMC draws, mutually exclusive with `n_draws`. These select
#'   exact draws without random subsampling. Not available for CAVI/SVI.
#' @param parameters Optional unique MCMC parameter block names. Blocks are
#'   returned in canonical model order. `NULL` returns every parameter;
#'   `character()` with `include_latent = TRUE` returns only CHLRM allocations.
#'   Not available for CAVI/SVI.
#' @param predictors Predictor names/indices to scale, or `"all"`.
#' @param center_predictors,scale_predictors Whether to center or scale selected
#'   predictor columns.
#' @param center_response,scale_response Whether to center or scale the response.
#' @param response Whether to scale the response.
#' @param intercept_column Explicit intercept column name or index.
#' @param zero_scale_tolerance Positive tolerance used to reject numerically
#'   constant columns.
#' @param map A learned `data_scaling_map`.
#' @param x A fitted model or posterior-draw object.
#' @param ... Additional method arguments.
#' @return A versioned fit, posterior summary/draws, orchestration audit, or
#'   scaling object, depending on the operation.
#' @details
#' When `fit_model()` receives a reusable variational initialization set and
#' an explicit `seed`, optimization starts from that seed independently of
#' random draws used to construct a temporary backend initialization. This
#' permits common optimization random streams across initialization strategies.
#' Single-start SVI results under this contract differ from versions before
#' 0.2.1, which continued the constructor's random stream. Reproduce historical
#' results with their original package version and source commit.
#'
#' For MCMC, consecutive `draw_indices` requests can extract a complete trace
#' in bounded output blocks. Each result retains the original `stored_draw`,
#' `chain`, and `iteration` in `draw_index`; its local `draw` starts at one.
#' Indices address the combined stored order, not an individual chain's
#' iteration. Explicit indices consume no random draws. Omitting selectors
#' preserves the existing extraction and variational sampling behavior.
#' Parameter subsets use schema `posterior-draws-2` and record
#' `metadata$parameter_selection`; complete parameter results retain schema
#' `posterior-draws-1`. Latent-only objects require `include_latent = TRUE` when
#' converted with `posterior::as_draws()` and have no parameter density plot.
#' The output budget accounts for selected blocks before extraction. It is
#' not a bound on process RSS: the fitted object and its backend must already
#' reside in memory. This API does not stream a checkpoint from disk, change
#' sampler retention, or compute convergence diagnostics independently on chunks.
#' @examples
#' set.seed(1)
#' X <- cbind(intercept = 1, predictor = seq(-1, 1, length.out = 20))
#' data <- list(X = X, y = drop(X %*% c(1, 2)) + rnorm(20, sd = 0.3))
#' prior <- lrm_prior(c(0, 0), 10 * diag(2), 4, 1)
#' fit <- fit_model(
#'   lrm_model(2), data, prior, "cavi", initialization_ols(),
#'   cavi_control(max_iterations = 50), seed = 2, verbose = FALSE
#' )
#' summarize_posterior(fit)
#' posterior_draws(fit, n_draws = 5, seed = 3)
#' @name fitting-and-draws
NULL

#' Build grouped regression data
#'
#' Convert a data frame to the common list-of-matrices representation used by
#' HLRM and CHLRM. No implicit factor-contrast expansion is performed.
#'
#' @param data A nonempty data frame.
#' @param response One response column name or index.
#' @param predictors Numeric predictor column names or indices.
#' @param group One grouping column name or index.
#' @param intercept Whether to prepend an explicit intercept.
#' @return A `vb_grouped_regression_data` object with grouped `X`, `y`, group
#'   mapping, and column provenance.
#' @examples
#' frame <- data.frame(y = 1:6, x = c(-1, 0, 1, -1, 0, 1),
#'                     group = rep(c("a", "b"), each = 3))
#' grouped_regression_data(frame, "y", "x", "group")
#' @name grouped-regression-data
NULL

#' Prediction and Bayesian model evaluation
#'
#' Keep mean prediction, response prediction, descriptive in-sample fit,
#' posterior predictive checking, pointwise likelihood, predictive criteria,
#' held-out scoring, refitted cross-validation, and PSIS-LOO as distinct
#' statistical operations.
#'
#' @param fit A fitted `vb_model_fit`.
#' @param object A fitted `vb_model_fit` for `predict()`.
#' @param data Training data used to construct a cross-validation plan or run.
#' @param prior A prior specification used in cross-validation refits.
#' @param method Inference engine used in cross-validation refits.
#' @param initialization Initialization used in cross-validation refits.
#' @param control Inference control used in cross-validation refits.
#' @param scaling Optional training-only scaling used in each refit.
#' @param newdata A design matrix or grouped design compatible with the fit.
#' @param group_ids Optional observed/new group identifiers.
#' @param n_draws Number of posterior or variational draws.
#' @param level Pointwise credible or predictive interval level.
#' @param seed Optional reproducibility seed.
#' @param return_yrep Whether to retain posterior replicated responses.
#' @param return_log_likelihood Whether to retain the pointwise matrix.
#' @param max_output_mb Allocation budget in MiB.
#' @param discrepancies A discrepancy specification from `ppc_discrepancies()`.
#' @param statistics Named posterior predictive discrepancy statistics.
#' @param prediction A response-prediction result containing `y_rep`.
#' @param model Model name or specification for an advanced likelihood call.
#' @param y,X Response and design in the representation required by `model`.
#' @param draws Posterior or variational draws.
#' @param unit Predictive unit, observation or group.
#' @param conditioning Conditional or marginalized latent-state convention.
#' @param nu0 HLRM residual degrees of freedom used for group marginalization.
#' @param groups Optional group labels.
#' @param predictive_target Explicit observed/new-data target.
#' @param predictive_unit Explicit observation- or group-level scoring unit.
#' @param marginal_control HLRM numerical marginalization controls.
#' @param log_lik A finite pointwise log-likelihood matrix or contract.
#' @param scheme Cross-validation scheme.
#' @param folds Number of cross-validation folds.
#' @param fold Fold index to materialize.
#' @param plan A cross-validation plan.
#' @param compute_performance Whether CV also computes prediction metrics.
#' @param max_fits Explicit cross-validation refit budget.
#' @param retain_fits Whether to retain fitted fold objects.
#' @param prediction_level Interval level for cross-validated predictions.
#' @param require_converged Whether PSIS-LOO requires diagnosed convergence.
#' @param cores Number of cores passed to `loo`.
#' @param verbose Whether to display refit progress.
#' @param ... Additional method arguments.
#' @return A versioned prediction, fit, PPC, likelihood, score, performance,
#'   cross-validation, WAIC, or PSIS-LOO result.
#' @examples
#' set.seed(11)
#' X <- cbind(intercept = 1, predictor = seq(-1, 1, length.out = 20))
#' data <- list(X = X, y = drop(X %*% c(0.5, 1.5)) + rnorm(20, sd = 0.3))
#' fit <- fit_model(
#'   lrm_model(2), data, lrm_prior(c(0, 0), 10 * diag(2), 4, 1),
#'   "cavi", initialization_ols(), cavi_control(max_iterations = 50),
#'   seed = 12, verbose = FALSE
#' )
#' predict(fit, X[1:2, , drop = FALSE], n_draws = 10, seed = 13)
#' posterior_predict(fit, X[1:2, , drop = FALSE], n_draws = 10, seed = 14)
#' evaluate_in_sample_fit(fit)
#' posterior_predictive_check(fit, n_draws = 20, seed = 15)
#' predictive_criteria(fit, n_draws = 20, seed = 16)
#' @name predictive-evaluation
NULL

#' Simulate from the supported Bayesian models
#'
#' Define a conditional, hierarchical, or prior-predictive DGP, a covariate
#' design, and sample structure. Simulation returns data, realized truth, and
#' conditioning metadata without fitting a model.
#'
#' @param model A model specification.
#' @param parameters Fixed conditional or population DGP parameters.
#' @param prior A prior specification for prior-predictive simulation.
#' @param distribution Covariate-design distribution.
#' @param mean,sd Covariate distribution parameters.
#' @param covariance Covariate covariance matrix.
#' @param bounds Bounds for a grid design.
#' @param binary_probabilities Probabilities for binary design columns.
#' @param binary_position Placement of binary columns.
#' @param include_intercept Whether to prepend an intercept.
#' @param fixed_train,fixed_test Optional fixed design matrices.
#' @param column_names Optional design-column names.
#' @param jitter_sd Optional grid jitter standard deviation.
#' @param sample_size LRM training sample size.
#' @param group_count Number of HLRM/CHLRM groups.
#' @param observations_per_group Scalar/vector group size.
#' @param test_sample_size Optional LRM test size.
#' @param test_observations_per_group Optional grouped test size.
#' @param test_group_mode Whether grouped test data reuse training groups or
#'   represent new groups.
#' @param test_group_count Number of new test groups when requested.
#' @param retain_observation_truth Whether to retain observation-level truth.
#' @param generator A typed DGP specification.
#' @param design A `simulation_design`.
#' @param sample A `simulation_sample`.
#' @param seed Optional reproducibility seed.
#' @param max_output_mb Output allocation budget.
#' @param simulation A `vb_simulation_result`.
#' @param include_observation_level Whether `tidy_truth()` includes retained
#'   observation-level truth.
#' @return A typed DGP/design/sample specification, `simulation-result-2`, or
#'   tidy truth table.
#' @examples
#' simulation <- simulate_data(
#'   conditional_dgp(lrm_model(2), list(
#'     coefficients = c(1, 2), residual_variance = 0.5)),
#'   simulation_design("grid", bounds = matrix(c(-1, 1), nrow = 1)),
#'   simulation_sample(sample_size = 20), seed = 21
#' )
#' tidy_truth(simulation)
#' @name simulation-api
NULL

#' Inspect CHLRM components and classifications
#'
#' Extract mixture-specific quantities without exposing internal result lists.
#' LRM and HLRM fits are rejected rather than receiving null fields.
#'
#' @param fit A CHLRM `vb_model_fit`.
#' @param reference A `chlrm_label_reference` obtained from
#'   `fit_orchestration()` for a compatible CHLRM fit. Coordinate-scale
#'   metadata and, when applicable, the fitted scaling map must match.
#' @param ... Arguments forwarded to `summarize_posterior()` where applicable.
#' @return `relabel_fit()` returns a new fit without modifying its input. Other
#'   functions return component summaries, a responsibility matrix, a
#'   classification data frame, uncertainty data frame, or co-clustering
#'   matrix.
#' @name chlrm-accessors
NULL

#' Inspect hierarchical parameter blocks
#'
#' Extract population, observed-group, and group-deviation summaries using the
#' package hierarchy rather than unrelated fixed/random-effect terminology.
#'
#' @param fit An HLRM or CHLRM `vb_model_fit` as applicable.
#' @param ... Arguments forwarded to `summarize_posterior()`.
#' @return Named parameter summaries or a group-by-coefficient deviation matrix.
#' @name hierarchical-accessors
NULL

#' Standard methods for fitted models
#'
#' Inspect a fitted model through conventional R generics without exposing its
#' mutable inference backend.
#'
#' @param x,object A `vb_model_fit`.
#' @param ... Additional arguments forwarded to the versioned operation.
#' @return The fit invisibly, a posterior summary, coefficients, fitted values,
#'   residuals, or the observation count, depending on the generic.
#' @name model-fit-methods
NULL

#' Extract matrices for Bayesian ecosystem tools
#'
#' Convert versioned prediction and likelihood results to the draws-by-unit
#' orientation used by packages such as bayesplot and loo.
#'
#' @param x A response prediction, PPC, or pointwise likelihood object.
#' @param ... Reserved for methods.
#' @return A numeric draws-by-unit matrix with estimand attributes where
#'   applicable.
#' @name interoperability
NULL

#' Checkpoint and resume typed MCMC fits
#'
#' Save plain serializable sampler state, reconstruct it, or continue every
#' chain without repeating warm-up. Continuation returns a new fit.
#'
#' @param fit An MCMC `vb_model_fit`.
#' @param path A nonempty checkpoint path; saving refuses overwrite.
#' @param additional_samples Post-warmup draws to append to every chain.
#' @param seed Optional seed for an explicit new continuation branch.
#' @param verbose Whether to display sampler progress.
#' @return The normalized path invisibly, or a restored/continued fit.
#' @name orchestration
NULL

#' Represent posterior summaries in long form
#'
#' Convert applicable parameter summaries to a dependency-free long data frame
#' for tables, plots, and replicated studies.
#'
#' @param x A fitted model or `posterior-summary-1` object.
#' @param ... For a fit, arguments forwarded to `summarize_posterior()`.
#' @return A `vb_tidy_parameters` data frame.
#' @name tidy-parameters
NULL

#' Estimate resources required by a fit
#'
#' Perform a conservative dimension-based memory preflight. The result is a
#' safety estimate, not a measured RSS value or operating-system limit.
#'
#' @param model A model specification.
#' @param data Data in the representation required by the model.
#' @param control A method-specific inference control.
#' @param scaling Optional data-scaling specification.
#' @param x A `vb_fit_resource_estimate`.
#' @param ... Unused.
#' @return A `fit-resource-estimate-1` object.
#' @name estimate-fit-resources
NULL
