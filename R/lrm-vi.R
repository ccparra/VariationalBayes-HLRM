vi_lrm <- R6::R6Class("vi_lrm",
                        public = list(
                          hyperparams = NULL,
                          run_settings = NULL,
                          lrm_vi_object = NULL,
                          results = NULL,
                          posterior_summary_result = NULL,
                          rng_state = NULL,
                          fit_status = "initialized",
                          diagnostics = NULL,

                          initialize = function(y, X, hyperparams, run_settings){
                            X <- as.matrix(X); y <- as.numeric(y)
                            if (!nrow(X) || !ncol(X) || nrow(X) != length(y) ||
                                any(!is.finite(X)) || any(!is.finite(y)))
                              stop("LRM requires finite y and X with compatible nonzero dimensions")
                            vi_stopping_control(run_settings, length(y))
                            hyperparams <- validate_lrm_hyperparams(hyperparams, ncol(X))
                            scope <- inference_rng_scope(seed = run_settings$seed)
                            on.exit(scope$restore(), add = TRUE)
                            self$run_settings <- run_settings
                            self$hyperparams <- hyperparams
                            self$lrm_vi_object <- list(hyperparams = hyperparams,
                                                       params = list(),
                                                       stats = list(
                                                         X = X,
                                                         y = y,
                                                         XtX = crossprod(X),
                                                         Xty = crossprod(X,y),
                                                         yty = drop(crossprod(y)),
                                                         n = length(y),
                                                         p =  ncol(X)
                                                       ))
                            self$lrm_vi_object <- initialize_lrm_vi(self$lrm_vi_object)
                            self$results <- list()
                            self$posterior_summary_result <- NULL
                            self$diagnostics <- NULL
                            if (scope$isolated) self$rng_state <- rng_snapshot()
                          },

                          fit = function(verbose = TRUE, seed = NULL) {
                            if (self$fit_status != "initialized")
                              stop("This model has already been fitted or failed; construct a new model for a new run")
                            control <- vi_stopping_control(
                              self$run_settings, self$lrm_vi_object$stats$n)
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
                            minibatch_n <- control$minibatch_m
                            tau <- control$tau
                            kappa <- control$kappa
                            elbo_every <- control$every

                            elbo     <- rep(NA_real_, max_iter)
                            monitor <- new_vi_monitor()

                            iter_time <- rep(NA_real_, max_iter)
                            elbo_time <- 0
                            t0 <- proc.time()[["elapsed"]]

                            if (verbose) pb <- txtProgressBar(min = 0, max = max_iter, width = 50, char = "=", style = 3)
                            if (verbose) on.exit(close(pb), add = TRUE)

                            n <- self$lrm_vi_object$stats$n

                            for (iter in seq_len(max_iter)) {

                              if (stochastic) {
                                B <- sample.int(
                                  n, size = min(minibatch_n, n),
                                  replace = FALSE)
                                learning_rate <- step_size(
                                  iter, tau = tau, kappa = kappa)
                              } else {
                                B <- NULL
                                learning_rate <- 1
                              }

                              sufficient_statistics <-
                                .lrm_scaled_sufficient_statistics(
                                  self$lrm_vi_object$stats, B)

                              self$lrm_vi_object <- update_q_coefficients_lrm(
                                self$lrm_vi_object, B, learning_rate,
                                sufficient_statistics)
                              self$lrm_vi_object <- update_q_residual_variance_lrm(
                                self$lrm_vi_object, B, learning_rate,
                                sufficient_statistics)

                              if (verbose) setTxtProgressBar(pb, iter)
                              scheduled <- iter %% elbo_every == 0L
                              should_eval <- iter == max_iter || scheduled
                              if (should_eval) {
                                te <- proc.time()[["elapsed"]]
                                elbo_new <- elbo_lrm(self$lrm_vi_object)
                                if (!is.finite(elbo_new))
                                  stop("Non-finite ELBO at iteration ", iter)
                                elbo_time <- elbo_time +
                                  proc.time()[["elapsed"]] - te
                                elbo[iter] <- elbo_new
                                iter_time[iter] <- proc.time()[["elapsed"]] - t0
                                if (scheduled) monitor <- update_vi_monitor(
                                  monitor, elbo_new, iter, control,
                                  vi_convergence_parameters(
                                    self$lrm_vi_object$params, "LRM"))
                                if (monitor$converged) break
                              }
                            }
                            total_time <- proc.time()[["elapsed"]] - t0
                            eval_idx <- which(!is.na(elbo))
                            elbo <- elbo[eval_idx]
                            iter_time <- iter_time[eval_idx]

                            self$results <- list(
                              MuBeta    = self$lrm_vi_object$params$MuBeta,
                              SigmaBeta = self$lrm_vi_object$params$SigmaBeta,
                              asigma2   = self$lrm_vi_object$params$asigma2,
                              bsigma2   = self$lrm_vi_object$params$bsigma2,
                              elbo      = elbo,
                              elbo_iter = eval_idx
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
                              "LRM", inference_method(self$run_settings, "CAVI"),
                              "variational parameters", self$results,
                              self$diagnostics, total_time, iter_time,
                              evaluation_iterations = eval_idx,
                              update_seconds = total_time - elbo_time,
                              objective_seconds = elbo_time,
                              metadata = list(
                                sampling_unit = if (stochastic)
                                  "observation" else "full_data",
                                minibatch_size = if (stochastic)
                                  minibatch_n else n)))
                          },

                          summarize_posterior = function() {
                            if (self$fit_status != "fitted") stop("Call fit() successfully before summarize_posterior()")
                            parameters <- new_parameter_slots("LRM")
                            parameters$coefficients <- summarize_normal_parameter(
                              self$results$MuBeta, self$results$SigmaBeta)
                            parameters$residual_variance <- summarize_inverse_gamma_parameter(
                              self$results$asigma2, self$results$bsigma2)
                            self$posterior_summary_result <- new_posterior_summary_result(
                              "LRM", inference_method(self$run_settings, "CAVI"),
                              parameters)
                            invisible(self$posterior_summary_result)
                          },

                          evaluate_in_sample_fit = function(return_details = FALSE) {
                            if (self$fit_status != "fitted")
                              stop("Call fit() successfully before evaluate_in_sample_fit()")
                            stats <- self$lrm_vi_object$stats
                            invisible(evaluate_in_sample_fit_gaussian(
                              list(stats$y), list(stats$X),
                              matrix(self$results$MuBeta, 1L), "LRM",
                              inference_method(self$run_settings, "CAVI"),
                              return_details))
                          },

                          posterior_predictive_check = function(n_draws = NULL, seed = NULL,
                                                                 return_yrep = FALSE,
                                                                 max_output_mb = 256,
                                                                 discrepancies = ppc_discrepancies()) {
                            if (self$fit_status != "fitted")
                              stop("Call fit() successfully before posterior_predictive_check()")
                            stats <- self$lrm_vi_object$stats
                            seeds <- evaluation_rng_streams(seed)
                            prepared <- prepare_gaussian_vi_evaluation(
                              "LRM", self$lrm_vi_object, self$results, n_draws,
                              self$run_settings$samples_gof, seeds$draws,
                              max_output_mb)
                            invisible(posterior_predictive_check_gaussian(
                              list(stats$y), list(stats$X), prepared$draws,
                              "LRM", inference_method(self$run_settings, "CAVI"),
                              seeds$replicates, return_yrep,
                              max_output_mb, discrepancies))
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
                            stats <- self$lrm_vi_object$stats
                            prepared <- prepare_gaussian_vi_evaluation(
                              "LRM", self$lrm_vi_object, self$results, n_draws,
                              self$run_settings$samples_gof, seed, max_output_mb)
                            invisible(predictive_criteria_gaussian(
                              list(stats$y), list(stats$X), prepared$draws,
                              model = "LRM",
                              method = inference_method(self$run_settings, "CAVI"),
                              predictive_target = predictive_target,
                              predictive_unit = predictive_unit,
                              return_log_likelihood = return_log_likelihood,
                              max_output_mb = max_output_mb, groups = groups))
                          },

                          predict = function(X_new, n_draws = NULL, level = 0.95, seed = NULL,
                                             return_draws = FALSE, max_output_mb = 256) {
                            run_lrm_prediction_vi(self$lrm_vi_object, self$results, X_new, n_draws,
                              self$run_settings$samples_gof, seed, level, return_draws,
                              max_output_mb, include_noise = FALSE,
                              stochastic = self$run_settings$stochastic)
                          },

                          posterior_predict = function(X_new, n_draws = NULL, level = 0.95,
                                                       seed = NULL, return_draws = TRUE,
                                                       max_output_mb = 256) {
                            run_lrm_prediction_vi(self$lrm_vi_object, self$results, X_new, n_draws,
                              self$run_settings$samples_gof, seed, level, return_draws,
                              max_output_mb, include_noise = TRUE,
                              stochastic = self$run_settings$stochastic)
                          }
                        )
)
