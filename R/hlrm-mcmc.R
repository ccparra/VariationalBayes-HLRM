mcmc_hlrm <- R6::R6Class(
  "mcmc_hlrm",
  public = list(
    hyperparams     = NULL,
    run_settings    = NULL,
    hlrm_mcmc_object = NULL,
    chains          = NULL,
    posterior_summary_result = NULL,
    rng_state       = NULL,
    fit_status      = "initialized",
    fit_metadata    = NULL,

    initialize = function(y,X,hyperparams, run_settings) {
      dimensions <- validate_grouped_data(X, y)
      validate_mcmc_settings(run_settings)
      hyperparams <- validate_hlrm_hyperparams(
        hyperparams, dimensions$p, dimensions$m, "HLRM")
      scope <- inference_rng_scope(seed = run_settings$seed)
      on.exit(scope$restore(), add = TRUE)
      self$run_settings     <- run_settings
      self$hyperparams      <- hyperparams
      self$hlrm_mcmc_object  <- list(hyperparams = hyperparams,
                                     params = list(),
                                     stats = list(
                                       X   = X,
                                       y   = y,
                                       XtX = lapply(X, crossprod),
                                       Xty = Map(crossprod,X,y),
                                       nj = lengths(y),
                                       m = length(y),
                                       p = ncol(X[[1]])
                                     ))
      self$hlrm_mcmc_object  <- initialize_hlrm_mcmc(self$hlrm_mcmc_object)
      self$chains           <- list()
      self$posterior_summary_result <- NULL
      self$fit_metadata     <- NULL
      if (scope$isolated) self$rng_state <- rng_snapshot()
    },

    fit = function(verbose = TRUE, seed = NULL) {
      if (self$fit_status != "initialized")
        stop("This model has already been fitted or failed; construct a new model for a new run")
      validate_mcmc_settings(self$run_settings)
      scope <- inference_rng_scope(seed, self$rng_state)
      on.exit(scope$restore(), add = TRUE)
      self$fit_status <- "running"
      on.exit({ if (self$fit_status == "running") self$fit_status <- "failed" }, add = TRUE)
      B <- n_samples_mcmc(self$run_settings)
      S <- self$run_settings$samples
      m <- self$hlrm_mcmc_object$stats$m
      p <- self$hlrm_mcmc_object$stats$p

      Betaj <- array(NA_real_, dim = c(S, m, p))
      sigma2j  <- array(NA_real_, dim = c(S, m))
      Beta <- matrix(NA_real_, nrow = S, ncol = p)
      Sigma <- array(NA_real_, dim = c(S, p, p))
      xi2 <- numeric(S)
      ll     <- numeric(S)
      lj     <- numeric(S)
      iter_time <- numeric(S)

      t0 <- proc.time()[["elapsed"]]

      if (verbose) pb <- txtProgressBar(min = 0, max = B, width = 50, char = "=", style = 3)
      if (verbose) on.exit(close(pb), add = TRUE)

      for (b in seq_len(B)) {
        self$hlrm_mcmc_object <- sample_group_coefficients_hlrm(self$hlrm_mcmc_object)
        self$hlrm_mcmc_object <- sample_group_variances_hlrm(self$hlrm_mcmc_object)
        self$hlrm_mcmc_object <- sample_xi2_hlrm(self$hlrm_mcmc_object)
        self$hlrm_mcmc_object <- sample_global_coefficients_hlrm(self$hlrm_mcmc_object)
        self$hlrm_mcmc_object <- sample_global_covariance_hlrm(self$hlrm_mcmc_object)

        idx <- thin_chain_mcmc(b, self$run_settings)
        if (!is.null(idx)) {
          Betaj[idx, , ]      <- self$hlrm_mcmc_object$params$Betaj
          sigma2j[idx,]      <- self$hlrm_mcmc_object$params$sigma2j
          Beta[idx, ]        <- self$hlrm_mcmc_object$params$Beta
          Sigma[idx, , ]     <- self$hlrm_mcmc_object$params$Sigma
          xi2[idx]           <- self$hlrm_mcmc_object$params$xi2
          iter_time[idx]   <- proc.time()[["elapsed"]] - t0

          ll[idx] <- loglik_total_hlrm(
            y_list = self$hlrm_mcmc_object$stats$y,
            X_list = self$hlrm_mcmc_object$stats$X,
            Betaj  = self$hlrm_mcmc_object$params$Betaj,
            sigma2j = self$hlrm_mcmc_object$params$sigma2j
          )
          lj[idx] <- logprior_hlrm(self$hlrm_mcmc_object) + ll[idx]
        }

        if (verbose) setTxtProgressBar(pb, b)
      }
      total_time <- proc.time()[["elapsed"]] - t0

      self$chains <- list(
        Betaj     = Betaj,
        sigma2j   = sigma2j,
        Beta      = Beta,
        Sigma     = Sigma,
        xi2       = xi2,
        loglik   = ll,
        logjoint = lj
      )
      if (scope$isolated) self$rng_state <- rng_snapshot()
      self$fit_status <- "fitted"
      self$fit_metadata <- list(total_seconds = total_time,
        warmup_iterations = self$run_settings$burn_in,
        stored_samples = S, extensions = list())

      invisible(new_fit_result(
        "HLRM", "MCMC", "posterior draws", self$chains,
        diagnostics = list(convergence = "not assessed for a single chain"),
        total_seconds = total_time, cumulative_seconds = iter_time,
        metadata = self$fit_metadata))
    },
    summarize_posterior = function() {
      if (self$fit_status != "fitted") stop("Call fit() successfully before summarize_posterior()")
      parameters <- new_parameter_slots("HLRM")
      parameters$population_coefficients <- summarize_draws_parameter(
        self$chains$Beta, margin = 2)
      parameters$population_covariance <- summarize_draws_parameter(
        self$chains$Sigma, margin = c(2, 3))
      parameters$variance_scale <- summarize_draws_parameter(self$chains$xi2)
      parameters$group_coefficients <- summarize_draws_parameter(
        self$chains$Betaj, margin = c(2, 3))
      parameters$group_variances <- summarize_draws_parameter(
        self$chains$sigma2j, margin = 2)
      self$posterior_summary_result <- new_posterior_summary_result(
        "HLRM", "MCMC", parameters)
      invisible(self$posterior_summary_result)
    },

    evaluate_in_sample_fit = function(return_details = FALSE) {
      if (self$fit_status != "fitted")
        stop("Call fit() successfully before evaluate_in_sample_fit()")
      stats <- self$hlrm_mcmc_object$stats
      point_beta <- matrix(apply(self$chains$Betaj, c(2, 3), mean),
                           stats$m, stats$p)
      invisible(evaluate_in_sample_fit_gaussian(
        stats$y, stats$X, point_beta, "HLRM", "MCMC", return_details))
    },

    posterior_predictive_check = function(n_draws = NULL, seed = NULL,
                                           return_yrep = FALSE,
                                           max_output_mb = 256,
                                           discrepancies = ppc_discrepancies()) {
      if (self$fit_status != "fitted")
        stop("Call fit() successfully before posterior_predictive_check()")
      stats <- self$hlrm_mcmc_object$stats
      seeds <- evaluation_rng_streams(seed)
      prepared <- prepare_gaussian_mcmc_evaluation(
        "HLRM", self$chains, stats, n_draws, self$run_settings$samples_gof,
        seeds$draws, max_output_mb)
      result <- posterior_predictive_check_gaussian(
        stats$y, stats$X, prepared$draws, "HLRM", "MCMC",
        seeds$replicates,
        return_yrep, max_output_mb, discrepancies)
      result$metadata$draw_indices <- prepared$indices
      invisible(result)
    },

    predictive_criteria = function(n_draws = NULL, seed = NULL,
                                   predictive_target = c("observed_groups", "new_groups"),
                                   return_log_likelihood = FALSE,
                                   max_output_mb = 256,
                                   marginal_control = list()) {
      if (self$fit_status != "fitted")
        stop("Call fit() successfully before predictive_criteria()")
      predictive_target <- match.arg(predictive_target)
      stats <- self$hlrm_mcmc_object$stats
      prepared <- prepare_gaussian_mcmc_evaluation(
        "HLRM", self$chains, stats, n_draws, self$run_settings$samples_gof,
        seed, max_output_mb)
      if (predictive_target == "observed_groups") {
        result <- predictive_criteria_gaussian(
          stats$y, stats$X, prepared$draws, "HLRM", "MCMC",
          predictive_target = predictive_target,
          return_log_likelihood = return_log_likelihood,
          max_output_mb = max_output_mb)
      } else {
        global <- list(
          Beta = self$chains$Beta[prepared$indices, , drop = FALSE],
          Sigma = self$chains$Sigma[prepared$indices, , , drop = FALSE],
          xi2 = self$chains$xi2[prepared$indices])
        result <- predictive_criteria_hlrm_marginal(
          stats$y, stats$X, global, self$hlrm_mcmc_object$hyperparams$nu0,
          method = "MCMC", control = marginal_control,
          return_log_likelihood = return_log_likelihood,
          max_output_mb = max_output_mb)
      }
      result$metadata$draw_indices <- prepared$indices
      invisible(result)
    },

    predict = function(X_new, group_ids = NULL, n_draws = NULL, level = 0.95,
                       seed = NULL, return_draws = FALSE, max_output_mb = 256) {
      run_hlrm_prediction_mcmc(self$chains, X_new, group_ids,
        self$hlrm_mcmc_object$hyperparams$nu0, n_draws,
        self$run_settings$samples_gof, seed, level, return_draws,
        max_output_mb, include_noise = FALSE)
    },

    posterior_predict = function(X_new, group_ids = NULL, n_draws = NULL,
                                 level = 0.95, seed = NULL,
                                 return_draws = TRUE, max_output_mb = 256) {
      run_hlrm_prediction_mcmc(self$chains, X_new, group_ids,
        self$hlrm_mcmc_object$hyperparams$nu0, n_draws,
        self$run_settings$samples_gof, seed, level, return_draws,
        max_output_mb, include_noise = TRUE)
    }
  )
)
