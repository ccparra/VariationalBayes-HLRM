mcmc_lrm <- R6::R6Class(
  "mcmc_lrm",
  public = list(
    hyperparams     = NULL,
    run_settings    = NULL,
    lrm_mcmc_object = NULL,
    chains          = NULL,
    posterior_summary_result = NULL,
    rng_state       = NULL,
    fit_status      = "initialized",
    fit_metadata    = NULL,

    initialize = function(y, X, hyperparams, run_settings) {
      validate_mcmc_settings(run_settings)
      X <- as.matrix(X); y <- as.numeric(y)
      if (!nrow(X) || !ncol(X) || nrow(X) != length(y) ||
          any(!is.finite(X)) || any(!is.finite(y)))
        stop("LRM requires finite y and X with compatible nonzero dimensions")
      hyperparams <- validate_lrm_hyperparams(hyperparams, ncol(X))
      scope <- inference_rng_scope(seed = run_settings$seed)
      on.exit(scope$restore(), add = TRUE)
      self$run_settings     <- run_settings
      self$hyperparams      <- hyperparams
      self$lrm_mcmc_object  <- list(hyperparams = hyperparams,
                                    params = list(),
                                    stats = list(
                                      X   = as.matrix(X),
                                      y   = drop(y),
                                      XtX = crossprod(X),
                                      Xty = crossprod(X, y),
                                      n = length(y),
                                      p = ncol(X)
                                    ))
      self$lrm_mcmc_object  <- initialize_lrm_mcmc(self$lrm_mcmc_object)
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
      p <- self$lrm_mcmc_object$stats$p


      Beta      <- matrix(NA_real_, nrow = S, ncol = p)
      sigma2    <- numeric(S)
      ll        <- numeric(S)
      lj        <- numeric(S)
      iter_time <- numeric(S)
      t0 <- proc.time()[["elapsed"]]

      if (verbose) pb <- txtProgressBar(min = 0, max = B, width = 50, char = "=", style = 3)
      if (verbose) on.exit(close(pb), add = TRUE)

      for (b in seq_len(B)) {
        self$lrm_mcmc_object <- sample_coefficients_lrm(self$lrm_mcmc_object)
        self$lrm_mcmc_object <- sample_residual_variance_lrm(self$lrm_mcmc_object)

        idx <- thin_chain_mcmc(b, self$run_settings)
        if (!is.null(idx)) {
          Beta[idx, ]      <- self$lrm_mcmc_object$params$Beta
          sigma2[idx]      <- self$lrm_mcmc_object$params$sigma2
          iter_time[idx]   <- proc.time()[["elapsed"]] - t0

          ll[idx] <- loglik_total_lrm(
            y     = self$lrm_mcmc_object$stats$y,
            X     = self$lrm_mcmc_object$stats$X,
            beta  = self$lrm_mcmc_object$params$Beta,
            sigma2= self$lrm_mcmc_object$params$sigma2
          )
          lj[idx] <- logprior_lrm(self$lrm_mcmc_object) + ll[idx]
        }

        if (verbose) setTxtProgressBar(pb, b)
      }
      total_time <- proc.time()[["elapsed"]] - t0

      self$chains <- list(
        Beta     = Beta,
        sigma2   = sigma2,
        loglik   = ll,
        logjoint = lj
      )
      if (scope$isolated) self$rng_state <- rng_snapshot()
      self$fit_status <- "fitted"
      self$fit_metadata <- list(total_seconds = total_time,
        warmup_iterations = self$run_settings$burn_in,
        stored_samples = S, extensions = list())

      invisible(new_fit_result(
        "LRM", "MCMC", "posterior draws", self$chains,
        diagnostics = list(convergence = "not assessed for a single chain"),
        total_seconds = total_time, cumulative_seconds = iter_time,
        metadata = self$fit_metadata))
    },

    summarize_posterior = function() {
      if (self$fit_status != "fitted") stop("Call fit() successfully before summarize_posterior()")
      parameters <- new_parameter_slots("LRM")
      parameters$coefficients <- summarize_draws_parameter(
        self$chains$Beta, margin = 2)
      parameters$residual_variance <- summarize_draws_parameter(self$chains$sigma2)
      self$posterior_summary_result <- new_posterior_summary_result(
        "LRM", "MCMC", parameters)
      invisible(self$posterior_summary_result)
    },

    evaluate_in_sample_fit = function(return_details = FALSE) {
      if (self$fit_status != "fitted")
        stop("Call fit() successfully before evaluate_in_sample_fit()")
      stats <- self$lrm_mcmc_object$stats
      invisible(evaluate_in_sample_fit_gaussian(
        list(stats$y), list(stats$X),
        matrix(colMeans(self$chains$Beta), 1L),
        "LRM", "MCMC", return_details))
    },

    posterior_predictive_check = function(n_draws = NULL, seed = NULL,
                                           return_yrep = FALSE,
                                           max_output_mb = 256,
                                           discrepancies = ppc_discrepancies()) {
      if (self$fit_status != "fitted")
        stop("Call fit() successfully before posterior_predictive_check()")
      stats <- self$lrm_mcmc_object$stats
      seeds <- evaluation_rng_streams(seed)
      prepared <- prepare_gaussian_mcmc_evaluation(
        "LRM", self$chains, stats, n_draws, self$run_settings$samples_gof,
        seeds$draws, max_output_mb)
      result <- posterior_predictive_check_gaussian(
        list(stats$y), list(stats$X), prepared$draws, "LRM", "MCMC",
        seeds$replicates, return_yrep, max_output_mb, discrepancies)
      result$metadata$draw_indices <- prepared$indices
      invisible(result)
    },

    predictive_criteria = function(n_draws = NULL, seed = NULL,
                                   predictive_target = "observed_data",
                                   predictive_unit = c("observation", "group"),
                                   return_log_likelihood = FALSE,
                                   max_output_mb = 256,
                                   groups = NULL) {
      if (self$fit_status != "fitted")
        stop("Call fit() successfully before predictive_criteria()")
      if (!identical(predictive_target, "observed_data"))
        stop("LRM predictive_criteria() targets observed_data")
      predictive_unit <- match.arg(predictive_unit)
      stats <- self$lrm_mcmc_object$stats
      prepared <- prepare_gaussian_mcmc_evaluation(
        "LRM", self$chains, stats, n_draws, self$run_settings$samples_gof,
        seed, max_output_mb)
      result <- predictive_criteria_gaussian(
        list(stats$y), list(stats$X), prepared$draws, "LRM", "MCMC",
        predictive_target = predictive_target,
        predictive_unit = predictive_unit,
        return_log_likelihood = return_log_likelihood,
        max_output_mb = max_output_mb, groups = groups)
      result$metadata$draw_indices <- prepared$indices
      invisible(result)
    },

    predict = function(X_new, n_draws = NULL, level = 0.95, seed = NULL,
                       return_draws = FALSE, max_output_mb = 256) {
      run_lrm_prediction_mcmc(self$chains, X_new, n_draws,
        self$run_settings$samples_gof, seed, level, return_draws,
        max_output_mb, include_noise = FALSE)
    },

    posterior_predict = function(X_new, n_draws = NULL, level = 0.95, seed = NULL,
                                 return_draws = TRUE, max_output_mb = 256) {
      run_lrm_prediction_mcmc(self$chains, X_new, n_draws,
        self$run_settings$samples_gof, seed, level, return_draws,
        max_output_mb, include_noise = TRUE)
    }
  )
)
