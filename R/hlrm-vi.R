vi_hlrm <- R6::R6Class("vi_hlrm",
                        public = list(
                          hyperparams = NULL,
                          run_settings = NULL,
                          hlrm_vi_object = NULL,
                          results = NULL,
                          posterior_summary_result = NULL,
                          rng_state = NULL,
                          fit_status = "initialized",
                          diagnostics = NULL,

                          initialize = function(y, X, hyperparams, run_settings){
                            self$run_settings <- run_settings
                            dimensions <- validate_grouped_data(X, y)
                            hyperparams <- validate_hlrm_hyperparams(
                              hyperparams, dimensions$p, dimensions$m, "HLRM")
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
                            self$hlrm_vi_object <- list(
                              hyperparams = hyperparams,
                              params = list(),
                              stats = list(
                                X   = X_use,
                                y   = y_use,
                                yty = unlist(lapply(y_use, crossprod)),
                                XtX = XtX_use,
                                Xty = Xty_use,
                                Xty_mat = t(vapply(Xty_use, as.numeric, numeric(p_use))),
                                XtX_mat = vapply(XtX_use, as.numeric, numeric(p_use * p_use)),
                                nj  = lengths(y_use),
                                m   = length(y_use),
                                p   = p_use,
                                Lambda0inv = hyperparams$Lambda0inv
                              )
                            )

                            self$hlrm_vi_object <- initialize_hlrm_vi(self$hlrm_vi_object)

                            self$results <- list()
                            self$posterior_summary_result <- NULL
                            self$diagnostics <- NULL
                            if (scope$isolated) self$rng_state <- rng_snapshot()
                          },

                          fit = function(verbose = TRUE, seed = NULL) {
                            if (self$fit_status != "initialized")
                              stop("This model has already been fitted or failed; construct a new model for a new run")
                            control <- vi_stopping_control(self$run_settings, self$hlrm_vi_object$stats$m)
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

                            elbo <- rep(NA_real_, max_iter)
                            iter_time <- rep(NA_real_, max_iter)
                            elbo_time <- 0
                            t0 <- proc.time()[["elapsed"]]

                            if (verbose) pb <- txtProgressBar(min = 0, max = max_iter, width = 50, char = "=", style = 3)
                            if (verbose) on.exit(close(pb), add = TRUE)

                            m <- self$hlrm_vi_object$stats$m
                            monitor <- new_vi_monitor()

                            for (iter in seq_len(max_iter)) {

                              if (stochastic) {
                                size  <- min(minibatch_m, m)
                                J     <- sample.int(m, size = size, replace = FALSE)
                                learning_rate <- step_size(
                                  iter, offset = tau, exponent = kappa)
                              } else {
                                J     <- NULL
                                learning_rate <- 1
                              }

                              self$hlrm_vi_object <- update_q_group_coefficients_hlrm(self$hlrm_vi_object, batch_idx = J)
                              self$hlrm_vi_object <- update_q_group_variances_hlrm(self$hlrm_vi_object, batch_idx = J)
                              self$hlrm_vi_object <- update_q_global_coefficients_hlrm(self$hlrm_vi_object, learning_rate = learning_rate, batch_idx = J)
                              self$hlrm_vi_object <- update_q_global_covariance_hlrm(self$hlrm_vi_object, learning_rate = learning_rate, batch_idx = J)
                              self$hlrm_vi_object <- update_q_xi2_hlrm(self$hlrm_vi_object, learning_rate = learning_rate, batch_idx = J)

                              if (verbose) setTxtProgressBar(pb, iter)
                              scheduled <- iter %% elbo_every == 0L
                              should_eval <- iter == max_iter || scheduled
                              if (should_eval) {
                                te <- proc.time()[["elapsed"]]
                                elbo_new <- elbo_hlrm(self$hlrm_vi_object)
                                if (!is.finite(elbo_new)) stop("Non-finite ELBO at iteration ", iter)
                                elbo_time <- elbo_time + proc.time()[["elapsed"]] - te
                                elbo[iter] <- elbo_new
                                iter_time[iter] <- proc.time()[["elapsed"]] - t0
                                if (scheduled) monitor <- update_vi_monitor(
                                  monitor, elbo_new, iter, control,
                                  vi_convergence_parameters(
                                    self$hlrm_vi_object$params, "HLRM"))
                                if (monitor$converged) break
                              }
                            }

                            total_time <- proc.time()[["elapsed"]] - t0
                            eval_idx <- which(!is.na(elbo))
                            elbo <- elbo[eval_idx]
                            iter_time <- iter_time[eval_idx]

                            self$results <- list(
                              MuBetaj    = self$hlrm_vi_object$params$MuBetaj,
                              SigmaBetaj = self$hlrm_vi_object$params$SigmaBetaj,
                              aj         = self$hlrm_vi_object$params$aj,
                              bj         = self$hlrm_vi_object$params$bj,
                              MuBeta     = self$hlrm_vi_object$params$MuBeta,
                              SigmaBeta  = self$hlrm_vi_object$params$SigmaBeta,
                              nuS        = self$hlrm_vi_object$params$nuS,
                              SS         = self$hlrm_vi_object$params$SS,
                              axi2       = self$hlrm_vi_object$params$axi2,
                              bxi2       = self$hlrm_vi_object$params$bxi2,
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
                              "HLRM", inference_method(self$run_settings, "CAVI"),
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

                          summarize_posterior = function() {
                            if (self$fit_status != "fitted") stop("Call fit() successfully before summarize_posterior()")
                            p <- self$hlrm_vi_object$stats$p
                            parameters <- new_parameter_slots("HLRM")
                            parameters$population_coefficients <- summarize_normal_parameter(
                              self$results$MuBeta, self$results$SigmaBeta)
                            parameters$population_covariance <- summarize_inverse_wishart_mean(
                              self$results$nuS, self$results$SS)
                            parameters$variance_scale <- summarize_gamma_parameter(
                              self$results$axi2, self$results$bxi2)
                            parameters$group_coefficients <- summarize_independent_normals(
                              self$results$MuBetaj, self$results$SigmaBetaj)
                            parameters$group_variances <- summarize_inverse_gamma_parameter(
                              self$results$aj, self$results$bj)
                            self$posterior_summary_result <- new_posterior_summary_result(
                              "HLRM", inference_method(self$run_settings, "CAVI"), parameters)
                            invisible(self$posterior_summary_result)
                          },

    evaluate_in_sample_fit = function(return_details = FALSE) {
                            if (self$fit_status != "fitted")
                              stop("Call fit() successfully before evaluate_in_sample_fit()")
                            stats <- self$hlrm_vi_object$stats
                            invisible(evaluate_in_sample_fit_gaussian(
                              stats$y, stats$X, self$results$MuBetaj, "HLRM",
                              inference_method(self$run_settings, "CAVI"),
                              return_details))
                          },

    posterior_predictive_check = function(n_draws = NULL, seed = NULL,
                                           return_yrep = FALSE,
                                           max_output_mb = 256,
                                           discrepancies = ppc_discrepancies()) {
                            if (self$fit_status != "fitted")
                              stop("Call fit() successfully before posterior_predictive_check()")
                            stats <- self$hlrm_vi_object$stats
                            seeds <- evaluation_rng_streams(seed)
                            prepared <- prepare_gaussian_vi_evaluation(
                              "HLRM", self$hlrm_vi_object, self$results, n_draws,
                              self$run_settings$samples_gof, seeds$draws,
                              max_output_mb)
                            invisible(posterior_predictive_check_gaussian(
                              stats$y, stats$X, prepared$draws, "HLRM",
                              inference_method(self$run_settings, "CAVI"),
                              seeds$replicates,
                              return_yrep, max_output_mb, discrepancies))
                          },

    predictive_criteria = function(n_draws = NULL, seed = NULL,
                                   predictive_target = c("observed_groups", "new_groups"),
                                   return_log_likelihood = FALSE,
                                   max_output_mb = 256,
                                   marginal_control = list()) {
                            if (self$fit_status != "fitted")
                              stop("Call fit() successfully before predictive_criteria()")
                            predictive_target <- match.arg(predictive_target)
                            stats <- self$hlrm_vi_object$stats
                            method <- inference_method(self$run_settings, "CAVI")
                            if (predictive_target == "new_groups") {
                              draw_count <- prediction_draw_count(
                                n_draws, self$run_settings$samples_gof)
                              check_gaussian_evaluation_budget(
                                lengths(stats$y), stats$p, draw_count,
                                max_output_mb)
                              scope <- inference_rng_scope(seed)
                              on.exit(scope$restore(), add = TRUE)
                              global <- vi_global_draws_hlrm(
                                self$results, draw_count)
                              return(invisible(predictive_criteria_hlrm_marginal(
                                stats$y, stats$X, global,
                                self$hlrm_vi_object$hyperparams$nu0,
                                method = method, control = marginal_control,
                                return_log_likelihood = return_log_likelihood,
                                max_output_mb = max_output_mb)))
                            }
                            prepared <- prepare_gaussian_vi_evaluation(
                              "HLRM", self$hlrm_vi_object, self$results, n_draws,
                              self$run_settings$samples_gof, seed, max_output_mb)
                            invisible(predictive_criteria_gaussian(
                              stats$y, stats$X, prepared$draws, model = "HLRM",
                              method = method,
                              predictive_target = predictive_target,
                              return_log_likelihood = return_log_likelihood,
                              max_output_mb = max_output_mb))
                          },

                          predict = function(X_new, group_ids = NULL, n_draws = NULL,
                                             level = 0.95, seed = NULL,
                                             return_draws = FALSE, max_output_mb = 256) {
                            run_hlrm_prediction_vi(self$hlrm_vi_object, self$results,
                              X_new, group_ids, self$hlrm_vi_object$hyperparams$nu0,
                              n_draws, self$run_settings$samples_gof, seed, level,
                              return_draws, max_output_mb, include_noise = FALSE,
                              stochastic = self$run_settings$stochastic)
                          },

                          posterior_predict = function(X_new, group_ids = NULL, n_draws = NULL,
                                                       level = 0.95, seed = NULL,
                                                       return_draws = TRUE,
                                                       max_output_mb = 256) {
                            run_hlrm_prediction_vi(self$hlrm_vi_object, self$results,
                              X_new, group_ids, self$hlrm_vi_object$hyperparams$nu0,
                              n_draws, self$run_settings$samples_gof, seed, level,
                              return_draws, max_output_mb, include_noise = TRUE,
                              stochastic = self$run_settings$stochastic)
                          }
                        )
)
