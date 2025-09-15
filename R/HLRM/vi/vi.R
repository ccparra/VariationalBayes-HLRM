source('R/utils.R')
source('R/HLRM/utils.R')
source('R/HLRM/vi/utils.R')

library(R6)

vi_hlrm <- R6::R6Class("vi_hlrm",
                        public = list(
                          hyperparams = NULL,
                          run_settings = NULL,
                          hlrm_vi_object = NULL,
                          results = NULL,
                          hat_values = NULL,
                          goodness_fit = NULL,
                          
                          initialize = function(y, X, standarize, hyperparams, run_settings){
                            self$run_settings <- run_settings
                            
                            std_flag <- isTRUE(standarize)
                            
                            if (std_flag) {
                              std <- standardize_Xy_list(X, y)
                              X_use <- std$X
                              y_use <- std$y
                              scale <- std$scale
                            } else {
                              X_use <- X
                              y_use <- y
                              scale <- NULL
                            }
                            
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
                                XtX = lapply(X_use, crossprod),
                                Xty = Map(crossprod, X_use, y_use),
                                nj  = lengths(y_use),
                                m   = length(y_use),
                                p   = ncol(X_use[[1]]),
                                standardized = std_flag,
                                scale = scale
                              )
                            )
                            
                            self$hlrm_vi_object <- Initialize_hlrm_vi(self$hlrm_vi_object)
                            
                            self$results <- list()
                            self$hat_values       <- list()
                            self$goodness_fit     <- list()
                          },

                          fit = function(verbose = TRUE) {
                            max_iter <- self$run_settings$iterations
                            eps      <- self$run_settings$epsilon
                            tau      <- self$run_settings$tau
                            kappa      <- self$run_settings$kappa
                            
                            
                            stochastic  <- isTRUE(self$run_settings$stochastic)
                            minibatch_m <- if (!is.null(self$run_settings$minibatch_m)) self$run_settings$minibatch_m else self$hlrm_vi_object$stats$m
                            tau         <- if (!is.null(self$run_settings$tau))   self$run_settings$tau   else 1
                            kappa       <- if (!is.null(self$run_settings$kappa)) self$run_settings$kappa else 0.6
                            
                            elbo <- numeric(max_iter)
                            iter_time <- numeric(max_iter)
                            t0 <- proc.time()[["elapsed"]]
                            
                            if (verbose) pb <- txtProgressBar(min = 0, max = max_iter, width = 50, char = "=", style = 3)
                            
                            m <- self$hlrm_vi_object$stats$m
                            elbo_old <- -Inf
                            
                            for (iter in seq_len(max_iter)) {
                              
                              if (stochastic) {
                                size  <- min(minibatch_m, m)
                                J     <- sample.int(m, size = size, replace = FALSE)
                                rho_t <- step_size(iter,tau=tau,kappa=kappa)
                              } else {
                                J     <- NULL   
                                rho_t <- 1
                              }
                              
                              self$hlrm_vi_object <- update_qbetaj_hlrm(self$hlrm_vi_object, batch_idx = J)
                              self$hlrm_vi_object <- update_qsigma2j_hlrm(self$hlrm_vi_object, batch_idx = J)
                              self$hlrm_vi_object <- update_qBeta_hlrm(self$hlrm_vi_object,  rho_t = rho_t, batch_idx = J)
                              self$hlrm_vi_object <- update_qSigma_hlrm(self$hlrm_vi_object, rho_t = rho_t, batch_idx = J)
                              self$hlrm_vi_object <- update_qxi2_hlrm(self$hlrm_vi_object,   rho_t = rho_t, batch_idx = J)
                              
                              elbo_new <- ELBO_hlrm(self$hlrm_vi_object)
                              elbo[iter] <- elbo_new
                              iter_time[iter] <- proc.time()[["elapsed"]] - t0
                              
                              if (verbose) setTxtProgressBar(pb, iter)
                              
                              elbo      <- elbo[seq_len(iter)]
                              iter_time <- iter_time[seq_len(iter)]
                              if (is.null(eps) && is.finite(elbo_old) && abs(elbo_new - elbo_old) < eps) {
                                if (verbose) {
                                  cat(sprintf("\nConverged at iter %d | ELBO = %.6f | Δ = %.3e\n",
                                              iter, elbo_new, abs(elbo_new - elbo_old)))
                                }
                                break
                              }
                              elbo_old <- elbo_new
                            }
                            
                            if (verbose) close(pb)
                            total_time <- proc.time()[["elapsed"]] - t0
                            
                            self$results <- list(
                              MuBetaj    = self$hlrm_vi_object$params$MuBetaj,
                              SigmaBetaj = self$hlrm_vi_object$params$SigmaBetaj,
                              aj         = self$hlrm_vi_object$params$aj,
                              bj         = self$hlrm_vi_object$params$bj,
                              elbo       = elbo
                            )
                            
                            invisible(list(
                              params = self$results,
                              timing = list(up_to_iter = iter_time, total = total_time)
                            ))
                          },
                          
                          inferences = function() {
                            MuBetaj    <- self$results$MuBetaj        
                            SigmaBetaj <- self$results$SigmaBetaj     
                            aj         <- self$results$aj             
                            bj         <- self$results$bj
                            
                            m <- self$hlrm_vi_object$stats$m
                            p <- self$hlrm_vi_object$stats$p
                            
                            z975 <- stats::qnorm(0.975)
                            
                            sd_bj   <- matrix(NA_real_, nrow = m, ncol = p)
                            CI_LI_bj<- matrix(NA_real_, nrow = m, ncol = p)
                            CI_LS_bj<- matrix(NA_real_, nrow = m, ncol = p)
                            
                            for (j in seq_len(m)) {
                              sd_j <- sqrt(diag(SigmaBetaj[[j]]))
                              sd_bj[j, ]    <- sd_j
                              CI_LI_bj[j, ] <- MuBetaj[j, ] - z975 * sd_j
                              CI_LS_bj[j, ] <- MuBetaj[j, ] + z975 * sd_j
                            }
                            
                            contains_zero_bj <- (CI_LI_bj <= 0 & 0 <= CI_LS_bj)
                            contains_mean_bj <- (CI_LI_bj <= MuBetaj & MuBetaj <= CI_LS_bj)
                            
                            vals_betaj <- list(
                              mean           = MuBetaj,
                              sd             = sd_bj,
                              CI_LI          = CI_LI_bj,
                              CI_LS          = CI_LS_bj,
                              contains_zero  = contains_zero_bj,
                              contains_mean  = contains_mean_bj
                            )
                            
                            mean_s2j <- ifelse(aj > 1, bj / (aj - 1), NA_real_)
                            sd_s2j   <- ifelse(aj > 2, bj / ((aj - 1) * sqrt(aj - 2)), NA_real_)  
                            
                         
                            CI_LI_s2j <- 1 / stats::qgamma(0.975, shape = aj, rate = bj)  
                            CI_LS_s2j <- 1 / stats::qgamma(0.025, shape = aj, rate = bj)  
                            
                            contains_zero_s2j <- (CI_LI_s2j <= 0 & 0 <= CI_LS_s2j)
                            contains_mean_s2j <- (CI_LI_s2j <= mean_s2j & mean_s2j <= CI_LS_s2j)
                            
                            vals_sigma2j <- list(
                              mean           = mean_s2j,
                              sd             = sd_s2j,
                              CI_LI          = CI_LI_s2j,
                              CI_LS          = CI_LS_s2j,
                              contains_zero  = contains_zero_s2j,
                              contains_mean  = contains_mean_s2j
                            )
                            
                            self$hat_values <- list(
                              Betaj   = vals_betaj,
                              sigma2j = vals_sigma2j
                            )
                            invisible(self$hat_values)
                          },
                          
                          goodness_of_fit = function(S = NULL, seed = NULL){
                            if (is.null(S)) S <- self$run_settings$samples_gof
                            y <- self$hlrm_vi_object$stats$y
                            X <- self$hlrm_vi_object$stats$X
                            
                            draws <- vi_sample_draws_hlrm(self$hlrm_vi_object, self$results,S = S, seed = seed)
                            chains <- list(Betaj = draws$Betaj, sigma2j = draws$sigma2j)
                            
                            waic <- waic_hlrm(chains = chains, y_list = y, X_list = X)
                            dic  <- dic_hlrm( chains = chains, y_list = y, X_list = X)
                            
                            yhat <- yhat_vi_hlrm(X, self$hat_values)
                            mse   <- mse(unlist(y), unlist(yhat))
                            r2   <- r2(unlist(y), unlist(yhat))
                            
                            ppvalues <- ppval_hlrm(y_list = y,
                                                   X_list= X,
                                                  chains = chains,
                                                  return_yrep = FALSE,
                                                  seed = seed)
                            
                            self$goodness_fit = list(waic = waic, 
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
                            invisible(yhat_vi_hlrm(X = X_new, hat_values = self$hat_values))
                          }
                        )
)
