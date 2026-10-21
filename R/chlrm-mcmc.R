# Exact integer counts accumulated in bounded draw blocks, then divided once.
# This avoids one groups-by-all-draws logical matrix per component and apply's
# full transpose, without subsampling or changing modal/tie conventions.
.chlrm_allocation_marginals <- function(z, K, block_size = 1024L) {
  S <- ncol(z); m <- nrow(z)
  counts <- matrix(0, m, K)
  occupied <- integer(S)
  names(occupied) <- colnames(z)
  for (first in seq.int(1L, S, by = block_size)) {
    index <- seq.int(first, min(S, first + block_size - 1L))
    block <- z[, index, drop = FALSE]
    for (k in seq_len(K)) counts[, k] <- counts[, k] + rowSums(block == k)
    occupied[index] <- vapply(seq_len(ncol(block)), function(s)
      length(unique(block[, s])), integer(1L))
  }
  list(probabilities = counts / S, occupied = occupied)
}

mcmc_chlrm <- R6::R6Class(
  "mcmc_chlrm",
  public = list(
    hyperparams     = NULL,
    run_settings    = NULL,
    chlrm_mcmc_object = NULL,
    chains          = NULL,
    relabeled_chains = NULL,
    label_reference = NULL,
    posterior_summary_result = NULL,
    rng_state      = NULL,
    fit_status     = "initialized",
    fit_metadata   = NULL,

    initialize = function(y, X, K, hyperparams, run_settings) {
      self$run_settings     <- run_settings
      dimensions <- validate_grouped_data(X, y, K)
      validate_mcmc_settings(run_settings)
      hyperparams <- validate_hlrm_hyperparams(
        hyperparams, dimensions$p, dimensions$m, "CHLRM", K)
      self$hyperparams <- hyperparams
      scope <- inference_rng_scope(seed = run_settings$seed)
      on.exit(scope$restore(), add = TRUE)
      X_use <- X
      y_use <- y

      # ---------------------------
      # Stats
      # ---------------------------
      self$chlrm_mcmc_object <- list(
        hyperparams = hyperparams,
        params = list(),
        stats = list(
          X   = X_use,
          y   = y_use,
          yty = unlist(lapply(y_use, crossprod)),
          XtX = lapply(X_use, crossprod),
          Xty = Map(crossprod, X_use, y_use),
          nj  = lengths(y_use),
          m   = length(y_use),
          p   = ncol(X_use[[1]]),
          K   = K
        )
      )
      self$chlrm_mcmc_object  <- initialize_chlrm_mcmc(self$chlrm_mcmc_object)
      self$chains           <- list()
      self$relabeled_chains <- list()
      self$label_reference <- NULL
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
      on.exit({
        if (self$fit_status == "running") self$fit_status <- "failed"
      }, add = TRUE)

      B <- n_samples_mcmc(self$run_settings)     # iter totales (incl. burn-in/thin)
      S <- self$run_settings$samples             # draws retained after thinning

      m <- self$chlrm_mcmc_object$stats$m
      p <- self$chlrm_mcmc_object$stats$p
      K <- self$chlrm_mcmc_object$stats$K

      # Cadenas a almacenar
      Betaj     <- array(NA_real_, dim = c(S, m, p))  # S × m × p
      sigma2j   <- matrix(NA_real_, nrow = S, ncol = m)   # S × m
      z_mat     <- matrix(NA_integer_, nrow = m, ncol = S) # m x S (groups x draws)
      omega_mat <- matrix(NA_real_, nrow = S, ncol = K)    # S × K
      Betak_arr <- array(NA_real_, dim = c(S, K, p))       # S × K × p
      sigma2k_m <- matrix(NA_real_, nrow = S, ncol = K)    # S × K
      Beta_g    <- matrix(NA_real_, nrow = S, ncol = p)
      Sigma_g   <- array(NA_real_, dim = c(S, p, p))
      xi2_g     <- numeric(S)

      ll        <- numeric(S)
      lj        <- numeric(S)
      iter_time <- numeric(S)

      t0 <- proc.time()[["elapsed"]]
      if (verbose) pb <- txtProgressBar(min = 0, max = B, width = 50, char = "=", style = 3)
      if (verbose) on.exit(close(pb), add = TRUE)

      for (b in seq_len(B)) {

        self$chlrm_mcmc_object <- sample_allocations_chlrm(self$chlrm_mcmc_object)
        self$chlrm_mcmc_object <- sample_omega_chlrm(  self$chlrm_mcmc_object)   # omega
        self$chlrm_mcmc_object <- sample_component_coefficients_chlrm(self$chlrm_mcmc_object)
        self$chlrm_mcmc_object <- sample_component_variances_chlrm(self$chlrm_mcmc_object)
        self$chlrm_mcmc_object <- sample_global_coefficients_chlrm(self$chlrm_mcmc_object)
        self$chlrm_mcmc_object <- sample_global_covariance_chlrm(self$chlrm_mcmc_object)
        self$chlrm_mcmc_object <- sample_xi2_chlrm(  self$chlrm_mcmc_object)     #xi2

        idx <- thin_chain_mcmc(b, self$run_settings)
        if (!is.null(idx)) {
          z_current <- self$chlrm_mcmc_object$params$z
          group_coefficients <- self$chlrm_mcmc_object$params$Betak[
            z_current, , drop = FALSE]
          group_variances <- self$chlrm_mcmc_object$params$sigma2k[z_current]

          Betaj[idx, , ]    <- group_coefficients                            # m × p
          sigma2j[idx, ]    <- group_variances                               # m
          z_mat[, idx]      <- z_current                                     # m (guardar como columna)
          omega_mat[idx, ]  <- self$chlrm_mcmc_object$params$omega           # K
          Betak_arr[idx, , ]<- self$chlrm_mcmc_object$params$Betak           # K × p
          sigma2k_m[idx, ]  <- self$chlrm_mcmc_object$params$sigma2k         # K
          Beta_g[idx, ]     <- self$chlrm_mcmc_object$params$Beta
          Sigma_g[idx, , ]  <- self$chlrm_mcmc_object$params$Sigma
          xi2_g[idx]        <- self$chlrm_mcmc_object$params$xi2

          iter_time[idx] <- proc.time()[["elapsed"]] - t0

          ll[idx] <- loglik_total_chlrm(
            y_list  = self$chlrm_mcmc_object$stats$y,
            X_list  = self$chlrm_mcmc_object$stats$X,
            Betaj   = group_coefficients,
            sigma2j = group_variances
          )
          lj[idx] <- logprior_chlrm(self$chlrm_mcmc_object) + ll[idx]
        }

        if (verbose) setTxtProgressBar(pb, b)
      }

      total_time <- proc.time()[["elapsed"]] - t0

      self$chains <- list(
        Betaj     = Betaj,       # S × m × p
        sigma2j   = sigma2j,     # S × m
        z         = z_mat,       # m x S (groups x draws)
        omega     = omega_mat,   # S × K
        Betak     = Betak_arr,   # S × K × p
        sigma2k   = sigma2k_m,   # S × K
        Beta      = Beta_g,
        Sigma     = Sigma_g,
        xi2       = xi2_g,
        loglik    = ll,
        logjoint  = lj
      )

      if (scope$isolated) self$rng_state <- rng_snapshot()
      self$fit_status <- "fitted"
      self$fit_metadata <- list(total_seconds = total_time,
        warmup_iterations = self$run_settings$burn_in,
        stored_samples = S, extensions = list())
      invisible(new_fit_result(
        "CHLRM", "MCMC", "posterior draws", self$chains,
        diagnostics = list(convergence = "not assessed for a single chain"),
        total_seconds = total_time, cumulative_seconds = iter_time,
        metadata = self$fit_metadata))
    },

    summarize_posterior = function(relabel = TRUE, label_reference = NULL,
                                   include_incidence = TRUE,
                                   max_relabel_iter = 100L) {
      if (self$fit_status != "fitted") stop("Call fit() successfully before summarize_posterior()")
      if (!is.logical(relabel) || length(relabel) != 1L || is.na(relabel))
        stop("relabel must be TRUE or FALSE")
      if (!is.logical(include_incidence) || length(include_incidence) != 1L ||
          is.na(include_incidence))
        stop("include_incidence must be TRUE or FALSE")
      validate_integer_scalar(max_relabel_iter, "max_relabel_iter", 1L)
      if (!isTRUE(relabel) && !is.null(label_reference))
        stop("label_reference requires relabel = TRUE")
      summary_chains <- if (isTRUE(relabel)) {
        if (is.null(label_reference)) label_reference <- fit_label_reference_chlrm(
          self$chains, max_iter = max_relabel_iter)
        align_to_reference_chlrm(self$chains, label_reference)
      } else self$chains
      self$relabeled_chains <- summary_chains
      self$label_reference <- if (isTRUE(relabel)) label_reference else NULL

      z_mat <- summary_chains$z
      m <- nrow(z_mat)
      S <- ncol(z_mat)
      K <- dim(summary_chains$Betak)[2]

      partition <- if (isTRUE(include_incidence))
        dahl_partition_chlrm(summary_chains$z) else NULL

      allocation <- .chlrm_allocation_marginals(z_mat, K)
      prob <- allocation$probabilities
      z_mode <- max.col(prob, ties.method = "first")
      assignment_entropy <- -rowSums(ifelse(prob > 0, prob * log(prob), 0))

      nc <- allocation$occupied
      tab_nc  <- table(nc)
      k_vals   <- as.integer(names(tab_nc))
      prob_nc <- as.numeric(tab_nc) / S
      names(prob_nc) <- names(tab_nc)

      parameters <- new_parameter_slots("CHLRM")
      parameters$population_coefficients <- summarize_draws_parameter(
        summary_chains$Beta, margin = 2)
      parameters$population_covariance <- summarize_draws_parameter(
        summary_chains$Sigma, margin = c(2, 3))
      parameters$variance_scale <- summarize_draws_parameter(summary_chains$xi2)
      parameters$component_coefficients <- summarize_draws_parameter(
        summary_chains$Betak, margin = c(2, 3))
      parameters$component_variances <- summarize_draws_parameter(
        summary_chains$sigma2k, margin = 2)
      parameters$mixture_weights <- summarize_draws_parameter(
        summary_chains$omega, margin = 2)
      allocations <- list(
        labels = z_mode, probabilities = prob, entropy = assignment_entropy,
        component_count = list(values = k_vals, probabilities = prob_nc))
      if (!is.null(partition)) {
        allocations$dahl_labels <- partition$z
        allocations$co_clustering <- partition$psm
      }
      relabeling_info <- list(applied = isTRUE(relabel))
      if (isTRUE(relabel)) relabeling_info <- c(relabeling_info, list(
        reference_draw = label_reference$initial_reference_draw,
        iterations = label_reference$iterations,
        algorithm = label_reference$algorithm))
      self$posterior_summary_result <- new_posterior_summary_result(
        "CHLRM", "MCMC", parameters, allocations, relabeling_info)
      invisible(self$posterior_summary_result)
    },

    evaluate_in_sample_fit = function(return_details = FALSE) {
      if (self$fit_status != "fitted")
        stop("Call fit() successfully before evaluate_in_sample_fit()")
      stats <- self$chlrm_mcmc_object$stats
      draws <- list(Betak = self$chains$Betak,
                    sigma2k = self$chains$sigma2k,
                    omega = self$chains$omega)
      invisible(evaluate_in_sample_fit_chlrm(
        stats$y, stats$X, draws, "MCMC", z = self$chains$z,
        return_details = return_details))
    },

    posterior_predictive_check = function(n_draws = NULL, seed = NULL,
                                           return_yrep = FALSE,
                                           max_output_mb = 256,
                                           discrepancies = ppc_discrepancies()) {
      if (self$fit_status != "fitted")
        stop("Call fit() successfully before posterior_predictive_check()")
      stats <- self$chlrm_mcmc_object$stats
      seeds <- evaluation_rng_streams(seed)
      prepared <- prepare_chlrm_mcmc_evaluation(
        self$chains, stats, n_draws, self$run_settings$samples_gof,
        seeds$draws, max_output_mb)
      result <- posterior_predictive_check_chlrm(
        stats$y, stats$X, prepared$draws, "MCMC", z = prepared$z,
        seed = seeds$replicates, return_yrep = return_yrep,
        max_output_mb = max_output_mb, discrepancies = discrepancies)
      result$metadata$draw_indices <- prepared$indices
      invisible(result)
    },

    predictive_criteria = function(n_draws = NULL, seed = NULL,
                                   predictive_target = "new_groups",
                                   return_log_likelihood = FALSE,
                                   max_output_mb = 256) {
      if (self$fit_status != "fitted")
        stop("Call fit() successfully before predictive_criteria()")
      if (!identical(predictive_target, "new_groups"))
        stop("CHLRM predictive_criteria() currently targets new_groups")
      stats <- self$chlrm_mcmc_object$stats
      prepared <- prepare_chlrm_mcmc_evaluation(
        self$chains, stats, n_draws, self$run_settings$samples_gof,
        seed, max_output_mb)
      result <- predictive_criteria_chlrm(
        stats$y, stats$X, prepared$draws, method = "MCMC",
        return_log_likelihood = return_log_likelihood,
        max_output_mb = max_output_mb)
      result$metadata$draw_indices <- prepared$indices
      invisible(result)
    },

    predict = function(X_new, group_ids = NULL, n_draws = NULL, level = 0.95,
                       seed = NULL, return_draws = FALSE, max_output_mb = 256) {
      run_chlrm_prediction_mcmc(self$chains, X_new, group_ids, n_draws,
        self$run_settings$samples_gof, seed, level, return_draws,
        max_output_mb, include_noise = FALSE)
    },

    posterior_predict = function(X_new, group_ids = NULL, n_draws = NULL,
                                 level = 0.95, seed = NULL,
                                 return_draws = TRUE, max_output_mb = 256) {
      run_chlrm_prediction_mcmc(self$chains, X_new, group_ids, n_draws,
        self$run_settings$samples_gof, seed, level, return_draws,
        max_output_mb, include_noise = TRUE)
    }
  )
)
