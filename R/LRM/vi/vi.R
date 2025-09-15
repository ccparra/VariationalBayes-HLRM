source('R/utils.R')
source('R/LRM/utils.R')
source('R/LRM/vi/utils.R')

library(R6)

vi_lrm <- R6::R6Class("vi_lrm",
                        public = list(
                          hyperparams = NULL,
                          run_settings = NULL,
                          lrm_vi_object = NULL,
                          results = NULL,
                          hat_values = NULL,
                          goodness_fit = NULL,
                          
                          initialize = function(y, X, hyperparams, run_settings){
                            self$run_settings <- run_settings
                            self$lrm_vi_object <- list(hyperparams = hyperparams,
                                                       params = list(),
                                                       stats = list(
                                                         X = as.matrix(X),
                                                         y = drop(y),
                                                         XtX = crossprod(X),
                                                         Xty = crossprod(X,y),
                                                         n = length(y),
                                                         p =  ncol(X)
                                                       ))
                            self$lrm_vi_object <- Initialize_lrm_vi(self$lrm_vi_object)
                            self$results <- list()
                            self$hat_values       <- list()
                            self$goodness_fit     <- list()
                          },

                          fit = function(verbose = TRUE) {
                            max_iter <- self$run_settings$iterations        
                            eps <- self$run_settings$epsilon
                            
                            elbo     <- numeric(max_iter)
                            elbo_old <- -Inf
                            
                            iter_time <- numeric(max_iter)
                            t0 <- proc.time()[["elapsed"]]
                            
                            if (verbose) pb <- txtProgressBar(min = 0, max = max_iter, width = 50, char = "=", style = 3)
                            
                            
                            for (iter in seq_len(max_iter)) {
                              
                              self$lrm_vi_object <- update_qBeta_lrm(self$lrm_vi_object)
                              self$lrm_vi_object <- update_qsigma2_lrm(self$lrm_vi_object)
                              
                              elbo_new <- ELBO_lrm(self$lrm_vi_object)
                              elbo[iter] <- elbo_new
                              iter_time[iter] <- proc.time()[["elapsed"]] - t0
                              
                              if (verbose) setTxtProgressBar(pb, iter)
                              
                              if (is.finite(elbo_old) && !is.null(eps) && abs(elbo_new - elbo_old) < eps) {
                                if (verbose) {
                                  cat(sprintf("\nConverged at iter %d | ELBO = %.6f | Δ = %.3e\n",
                                              iter, elbo_new, abs(elbo_new - elbo_old)))
                                }
                                elbo      <- elbo[seq_len(iter)]
                                iter_time <- iter_time[seq_len(iter)]
                                break
                              }
                              elbo_old <- elbo_new
                            }
                            if (verbose) close(pb)
                            total_time <- proc.time()[["elapsed"]] - t0
                            
                            self$results <- list(
                              MuBeta    = self$lrm_vi_object$params$MuBeta,
                              SigmaBeta = self$lrm_vi_object$params$SigmaBeta,
                              asigma2   = self$lrm_vi_object$params$asigma2,
                              bsigma2   = self$lrm_vi_object$params$bsigma2,
                              elbo      = elbo
                            )
                            invisible(list(params = self$results, 
                                           timing = list(up_to_iter = iter_time, 
                                                         total = total_time)))
                          },
                          
                          inferences = function() {
                            stopifnot(length(self$results) > 0)
                            
                            MuBeta    <- self$results$MuBeta   
                            SigmaBeta <- self$results$SigmaBeta 
                            asigma2   <- self$results$asigma2
                            bsigma2   <- self$results$bsigma2
                            
                            sd_beta <- sqrt(diag(SigmaBeta))
                            z975    <- stats::qnorm(0.975)
                            
                            CI_LI_b <- drop(MuBeta - z975 * sd_beta)
                            CI_LS_b <- drop(MuBeta + z975 * sd_beta)
                            
                            vals_beta <- list(
                              mean           = drop(MuBeta),
                              sd             = drop(sd_beta),
                              CI_LI          = CI_LI_b,
                              CI_LS          = CI_LS_b,
                              contains_zero  = (CI_LI_b <= 0 & 0 <= CI_LS_b),
                              contains_mean  = (CI_LI_b <= drop(MuBeta) & drop(MuBeta) <= CI_LS_b)
                            )
                            
                            mean_s2 <- if (asigma2 > 1) bsigma2 / (asigma2 - 1) else NA_real_
                            sd_s2   <- if (asigma2 > 2) bsigma2 / ((asigma2 - 1) * sqrt(asigma2 - 2)) else NA_real_
                            CI_LI_s2 <- 1 / stats::qgamma(0.975, shape = asigma2, rate = bsigma2) 
                            CI_LS_s2 <- 1 / stats::qgamma(0.025, shape = asigma2, rate = bsigma2) 
                            
                            vals_sigma2 <- list(
                              mean           = as.numeric(mean_s2),
                              sd             = as.numeric(sd_s2),
                              CI_LI          = as.numeric(CI_LI_s2),
                              CI_LS          = as.numeric(CI_LS_s2),
                              contains_zero  = (CI_LI_s2 <= 0 && 0 <= CI_LS_s2), 
                              contains_mean  = (!is.na(mean_s2) && CI_LI_s2 <= mean_s2 && mean_s2 <= CI_LS_s2)
                            )
                            
                            self$hat_values <- list(
                              Beta   = vals_beta,
                              sigma2 = vals_sigma2
                            )
                            
                            invisible(self$hat_values)
                          },
                          
                          goodness_of_fit = function(S = NULL, seed = NULL){
                            if (is.null(S)) S <- self$run_settings$samples_gof
                            y <- self$lrm_vi_object$stats$y
                            X <- self$lrm_vi_object$stats$X
                            
                            draws <- vi_sample_draws_lrm(self$lrm_vi_object, S = S, seed = seed)
                            chains <- list(Beta = draws$Beta, sigma2 = draws$sigma2)
                            
                            waic <- waic_lrm(chains = chains, y = y, X = X)
                            dic  <- dic_lrm( chains = chains, y = y, X = X)
                            
                            yhat <- yhat_vi_lrm(X, self$hat_values)
                            mse   <- mse(y, yhat)
                            r2  <- r2(y,yhat)
                
                            ppvalues <- ppval_lrm(y = y, 
                                                  X = X,
                                                  chains = chains,
                                                  return_yrep = FALSE,
                                                  seed = seed)
                            
                            self$goodness_fit = list(
                              waic = waic,
                              dic = dic,
                              mse = mse,
                              r2 = r2,
                              ppval = ppvalues
                            )
                            
                            invisible(self$goodness_fit)
                          },
                          
                          predict = function(X_new) {
                            stopifnot(length(self$hat_values) > 0)
                            X_new <- as.matrix(X_new)
                            invisible(yhat_vi_lrm(X = X_new, hat_values = self$hat_values))
                          }
                        )
)
