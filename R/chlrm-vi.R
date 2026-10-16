vi_chlrm <- R6::R6Class("vi_chlrm",
                       public = list(
                         hyperparams = NULL,
                         run_settings = NULL,
                         chlrm_vi_object = NULL,
                         results = NULL,
                         aligned_results = NULL,
                         label_reference = NULL,
                         posterior_summary_result = NULL,
                         rng_state = NULL,
                         fit_status = "initialized",
                         diagnostics = NULL,

                         initialize = function(y, X, K, hyperparams, run_settings){
                           self$run_settings <- run_settings
                           dimensions <- validate_grouped_data(X, y, K)
                           hyperparams <- validate_hlrm_hyperparams(
                             hyperparams, dimensions$p, dimensions$m, "CHLRM", K)
                           self$hyperparams <- hyperparams
                           vi_stopping_control(run_settings, length(y))
                           scope <- inference_rng_scope(seed = run_settings$seed)
                           on.exit(scope$restore(), add = TRUE)
                           X_use <- X
                           y_use <- y

                           XtX_use <- lapply(X_use, crossprod)
                           Xty_use <- Map(crossprod, X_use, y_use)
                           p_use <- ncol(X_use[[1]])

                           # ---------------------------
                           # Stats
                           # ---------------------------
                           self$chlrm_vi_object <- list(
                             hyperparams = hyperparams,
                             params = list(),
                             stats = list(
                               X   = X_use,
                               y   = y_use,
                               yty = unlist(lapply(y_use, crossprod)),
                               XtX = XtX_use,
                               Xty = Xty_use,
                               Xty_mat = t(matrix(vapply(Xty_use, as.numeric, numeric(p_use)),
                                                  nrow = p_use)),
                               XtX_mat = matrix(vapply(XtX_use, as.numeric, numeric(p_use * p_use)),
                                                nrow = p_use * p_use),
                               nj  = lengths(y_use),
                               m   = length(y_use),
                               p   = p_use,
                               K   = K,
                               Lambda0inv = hyperparams$Lambda0inv
                             )
                           )

                           self$chlrm_vi_object <- initialize_chlrm_vi(
                             self$chlrm_vi_object, run_settings$initialization)

                           self$results      <- list()
                           self$aligned_results <- list()
                           self$label_reference <- NULL
                           self$posterior_summary_result <- NULL
                           if (scope$isolated) self$rng_state <- rng_snapshot()
                         },


                         fit = function(verbose = TRUE, seed = NULL) {
                           # --------------------------------------------
                           # Run settings
                           # --------------------------------------------
                           if (self$fit_status != "initialized")
                             stop("This model has already been fitted or failed; construct a new model for a new run")
                           control <- vi_stopping_control(self$run_settings, self$chlrm_vi_object$stats$m)
                           scope <- inference_rng_scope(seed, self$rng_state)
                           on.exit(scope$restore(), add = TRUE)
                           self$fit_status <- "running"
                           on.exit({
                             if (self$fit_status == "running") {
                               self$fit_status <- "failed"
                               self$diagnostics <- list(converged = FALSE, stop_reason = "error")
                             }
                           }, add = TRUE)
                           max_iter <- control$iterations
                           stochastic <- control$stochastic
                           minibatch_m <- control$minibatch_m
                           tau <- control$tau
                           kappa <- control$kappa
                           elbo_every <- control$every

                           # --------------------------------------------
                           # Book-keeping
                           # --------------------------------------------
                           elbo <- rep(NA_real_, max_iter)
                           iter_time <- rep(NA_real_, max_iter)
                           elbo_time <- 0
                           t0 <- proc.time()[["elapsed"]]
                           if (verbose) pb <- txtProgressBar(min = 0, max = max_iter, width = 50, char = "=", style = 3)
                           if (verbose) on.exit(close(pb), add = TRUE)

                           m <- self$chlrm_vi_object$stats$m
                           monitor <- new_vi_monitor()

                           # --------------------------------------------
                           # Main loop (CAVI ⇔ stochastic=FALSE; SVI ⇔ TRUE)
                           # --------------------------------------------
                           for (iter in seq_len(max_iter)) {

                             if (stochastic) {
                               size  <- min(minibatch_m, m)
                               J     <- sample.int(m, size = size, replace = FALSE)
                               learning_rate <- step_size(iter, tau = tau,
                                                          kappa = kappa)
                             } else {
                               J     <- NULL
                               learning_rate <- 1
                             }

                             # ===== CHLRM updates =====
                             # 1) Responsabilidades locales q(z_j)
                             self$chlrm_vi_object <- update_q_allocations_chlrm(self$chlrm_vi_object, batch_idx = J)

                             # 2) Bloques por clúster con paso SVI y minibatch J
                             self$chlrm_vi_object <- update_q_component_coefficients_chlrm(self$chlrm_vi_object, batch_idx = J, learning_rate = learning_rate)
                             self$chlrm_vi_object <- update_q_component_variances_chlrm(self$chlrm_vi_object, batch_idx = J, learning_rate = learning_rate)
                             self$chlrm_vi_object <- update_q_mixture_weights_chlrm(self$chlrm_vi_object, batch_idx = J, learning_rate = learning_rate)

                             # 3) Globales (SVI en parámetros naturales, no dependen de J)
                             self$chlrm_vi_object <- update_q_global_coefficients_chlrm(self$chlrm_vi_object, learning_rate = learning_rate)
                             self$chlrm_vi_object <- update_q_global_covariance_chlrm(self$chlrm_vi_object, learning_rate = learning_rate)
                             self$chlrm_vi_object <- update_q_xi2_chlrm(self$chlrm_vi_object, learning_rate = learning_rate)

                             if (verbose) setTxtProgressBar(pb, iter)
                             scheduled <- iter %% elbo_every == 0L
                             should_eval <- iter == max_iter || scheduled
                             if (should_eval) {
                               te <- proc.time()[["elapsed"]]
                               elbo_new <- elbo_chlrm(self$chlrm_vi_object)
                               if (!is.finite(elbo_new)) stop("Non-finite ELBO at iteration ", iter)
                               elbo_time <- elbo_time + proc.time()[["elapsed"]] - te
                               elbo[iter] <- elbo_new
                               iter_time[iter] <- proc.time()[["elapsed"]] - t0

                               if (scheduled) monitor <- update_vi_monitor(
                                 monitor, elbo_new, iter, control,
                                 vi_convergence_parameters(
                                   self$chlrm_vi_object$params, "CHLRM"))
                               if (monitor$converged) break
                             }
                           }

                           total_time <- proc.time()[["elapsed"]] - t0
                           eval_idx <- which(!is.na(elbo))
                           elbo <- elbo[eval_idx]
                           iter_time <- iter_time[eval_idx]

                           # --------------------------------------------
                           # Save results (incluye parámetros clave CHLRM)
                           # --------------------------------------------
                           self$results <- list(
                             # mezcla
                             rho    = self$chlrm_vi_object$params$rho,
                             alpha  = self$chlrm_vi_object$params$alpha,
                             # clúster
                             MuBetak    = self$chlrm_vi_object$params$MuBetak,
                             SigmaBetak = self$chlrm_vi_object$params$SigmaBetak,
                             ak         = self$chlrm_vi_object$params$ak,
                             bk         = self$chlrm_vi_object$params$bk,
                             # globales
                             MuBeta     = self$chlrm_vi_object$params$MuBeta,
                             SigmaBeta  = self$chlrm_vi_object$params$SigmaBeta,
                             nuS        = self$chlrm_vi_object$params$nuS,
                             SS         = self$chlrm_vi_object$params$SS,
                             axi2       = self$chlrm_vi_object$params$axi2,
                             bxi2       = self$chlrm_vi_object$params$bxi2,
                             # tracking
                             elbo       = elbo,
                             elbo_iter  = eval_idx
                           )
                           if (scope$isolated) self$rng_state <- rng_snapshot()
                           self$fit_status <- "fitted"
                           self$diagnostics <- list(
                             converged = monitor$converged,
                             stop_reason = vi_convergence_reason(
                               control, monitor$converged),
                             iterations = iter, stable_checks = monitor$stable,
                             relative_change = monitor$relative_change,
                             parameter_change = monitor$parameter_change,
                             smoothed_elbo = monitor$smoothed_elbo,
                             history = vi_monitor_history(monitor),
                             control = control)

                           invisible(new_fit_result(
                             "CHLRM", inference_method(self$run_settings, "CAVI"),
                             "variational parameters", self$results,
                             self$diagnostics, total_time, iter_time,
                             evaluation_iterations = eval_idx,
                             update_seconds = total_time - elbo_time,
                             objective_seconds = elbo_time,
                             metadata = list(
                               sampling_unit = if (stochastic)
                                 "group" else "full_data",
                               minibatch_size = if (stochastic)
                                 minibatch_m else m)))
                         },

                         summarize_posterior = function(relabel = TRUE,
                                               label_reference = NULL,
                                               include_incidence = TRUE,
                                               max_relabel_iter = 100L) {
                           if (self$fit_status != "fitted") stop("Call fit() successfully before summarize_posterior()")
                           if (!is.logical(relabel) || length(relabel) != 1L || is.na(relabel))
                             stop("relabel must be TRUE or FALSE")
                           if (!is.logical(include_incidence) ||
                               length(include_incidence) != 1L || is.na(include_incidence))
                             stop("include_incidence must be TRUE or FALSE")
                           validate_integer_scalar(max_relabel_iter,
                                                   "max_relabel_iter", 1L)
                           if (!isTRUE(relabel) && !is.null(label_reference))
                             stop("label_reference requires relabel = TRUE")
                           # A multistart result may already carry the common
                           # MCMC/variational reference. Reuse it unless the caller
                           # explicitly requests raw labels with NULL.
                           if (isTRUE(relabel) && missing(label_reference))
                             label_reference <- self$label_reference
                           if (!isTRUE(relabel)) label_reference <- NULL
                           # Tamaños
                           K <- self$chlrm_vi_object$stats$K
                           p <- self$chlrm_vi_object$stats$p
                           m <- self$chlrm_vi_object$stats$m

                           summary_results <- if (is.null(label_reference)) self$results else
                             align_to_reference_chlrm(self$results, label_reference)
                           self$aligned_results <- summary_results
                           self$label_reference <- label_reference

                           rho <- summary_results$rho
                           z_pred <- max.col(rho, ties.method = "first")
                           assignment_entropy <- -rowSums(ifelse(rho > 0, rho * log(rho), 0))

                           A_soft <- NULL
                           if (isTRUE(include_incidence)) {
                             A_soft <- rho %*% t(rho)
                             diag(A_soft) <- 1
                           }

                           parameters <- new_parameter_slots("CHLRM")
                           parameters$population_coefficients <- summarize_normal_parameter(
                             summary_results$MuBeta, summary_results$SigmaBeta)
                           parameters$population_covariance <- summarize_inverse_wishart_mean(
                             summary_results$nuS, summary_results$SS)
                           parameters$variance_scale <- summarize_gamma_parameter(
                             summary_results$axi2, summary_results$bxi2)
                           parameters$component_coefficients <- summarize_independent_normals(
                             summary_results$MuBetak, summary_results$SigmaBetak)
                           parameters$component_variances <- summarize_inverse_gamma_parameter(
                             summary_results$ak, summary_results$bk)
                           parameters$mixture_weights <- summarize_dirichlet_parameter(
                             summary_results$alpha)
                           allocations <- list(
                             labels = z_pred, probabilities = rho,
                             entropy = assignment_entropy)
                           if (!is.null(A_soft))
                             allocations$co_clustering <- A_soft
                           relabeling_info <- list(
                             applied = !is.null(label_reference))
                           if (!is.null(label_reference)) {
                             relabeling_info <- c(relabeling_info, list(
                               reference_draw =
                                 label_reference$initial_reference_draw,
                               iterations = label_reference$iterations,
                               algorithm = label_reference$algorithm))
                             permutation <- attr(summary_results, "permutation")
                             if (!is.null(permutation))
                               relabeling_info$permutation <- permutation
                           }
                           self$posterior_summary_result <- new_posterior_summary_result(
                             "CHLRM", inference_method(self$run_settings, "CAVI"),
                             parameters, allocations, relabeling_info)
                           invisible(self$posterior_summary_result)
                         },

                         evaluate_in_sample_fit = function(return_details = FALSE) {
                           if (self$fit_status != "fitted")
                             stop("Call fit() successfully before evaluate_in_sample_fit()")
                           stats <- self$chlrm_vi_object$stats
                           invisible(evaluate_in_sample_fit_chlrm(
                             stats$y, stats$X, NULL,
                             inference_method(self$run_settings, "CAVI"),
                             rho = self$results$rho,
                             component_mean = self$results$MuBetak,
                             return_details = return_details))
                         },

                         posterior_predictive_check = function(n_draws = NULL, seed = NULL,
                                                                return_yrep = FALSE,
                                                                max_output_mb = 256,
                                                                discrepancies = ppc_discrepancies()) {
                           if (self$fit_status != "fitted")
                             stop("Call fit() successfully before posterior_predictive_check()")
                           stats <- self$chlrm_vi_object$stats
                           seeds <- evaluation_rng_streams(seed)
                           prepared <- prepare_chlrm_vi_evaluation(
                             self$chlrm_vi_object, self$results, n_draws,
                             self$run_settings$samples_gof, seeds$draws,
                             max_output_mb)
                           invisible(posterior_predictive_check_chlrm(
                             stats$y, stats$X, prepared$draws,
                             inference_method(self$run_settings, "CAVI"),
                             rho = prepared$rho, seed = seeds$replicates,
                             return_yrep = return_yrep,
                             max_output_mb = max_output_mb,
                             discrepancies = discrepancies))
                         },

                         predictive_criteria = function(n_draws = NULL, seed = NULL,
                                                        predictive_target = "new_groups",
                                                        return_log_likelihood = FALSE,
                                                        max_output_mb = 256) {
                           if (self$fit_status != "fitted")
                             stop("Call fit() successfully before predictive_criteria()")
                           if (!identical(predictive_target, "new_groups"))
                             stop("CHLRM predictive_criteria() currently targets new_groups")
                           stats <- self$chlrm_vi_object$stats
                           prepared <- prepare_chlrm_vi_evaluation(
                             self$chlrm_vi_object, self$results, n_draws,
                             self$run_settings$samples_gof, seed, max_output_mb)
                           invisible(predictive_criteria_chlrm(
                             stats$y, stats$X, prepared$draws,
                             method = inference_method(self$run_settings, "CAVI"),
                             return_log_likelihood = return_log_likelihood,
                             max_output_mb = max_output_mb))
                         },

                         predict = function(X_new, group_ids = NULL, n_draws = NULL,
                                            level = 0.95, seed = NULL,
                                            return_draws = FALSE, max_output_mb = 256) {
                           run_chlrm_prediction_vi(self$chlrm_vi_object, self$results,
                             X_new, group_ids, n_draws, self$run_settings$samples_gof,
                             seed, level, return_draws, max_output_mb,
                             include_noise = FALSE,
                             stochastic = self$run_settings$stochastic)
                         },

                         posterior_predict = function(X_new, group_ids = NULL, n_draws = NULL,
                                                      level = 0.95, seed = NULL,
                                                      return_draws = TRUE,
                                                      max_output_mb = 256) {
                           run_chlrm_prediction_vi(self$chlrm_vi_object, self$results,
                             X_new, group_ids, n_draws, self$run_settings$samples_gof,
                             seed, level, return_draws, max_output_mb,
                             include_noise = TRUE,
                             stochastic = self$run_settings$stochastic)
                         }
                       )
)
