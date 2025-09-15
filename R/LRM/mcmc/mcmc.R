source('R/utils.R')
source('R/LRM/utils.R')
source('R/LRM/mcmc/utils.R')

library(R6)

mcmc_lrm <- R6::R6Class(
  "mcmc_lrm",
  public = list(
    hyperparams     = NULL,
    run_settings    = NULL,
    lrm_mcmc_object = NULL,
    chains          = NULL,
    hat_values      = NULL,
    goodness_fit   = NULL,
    
    initialize = function(y, X, hyperparams, run_settings) {
      self$run_settings     <- run_settings
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
      self$lrm_mcmc_object  <- Initialize_lrm_mcmc(self$lrm_mcmc_object)
      self$chains           <- list()
      self$hat_values       <- list()
      self$goodness_fit     <- list() 
    },
    
    fit = function(verbose = TRUE, seed = NULL) {
      if (!is.null(seed)) set.seed(seed)
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
      
      for (b in seq_len(B)) {
        self$lrm_mcmc_object <- sample_Beta_lrm(self$lrm_mcmc_object)
        self$lrm_mcmc_object <- sample_sigma2_lrm(self$lrm_mcmc_object)
        
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
      if (verbose) close(pb)
      total_time <- proc.time()[["elapsed"]] - t0
      
      self$chains <- list(
        Beta     = Beta,
        sigma2   = sigma2,
        loglik   = ll,
        logjoint = lj
      )
      
      invisible(list(chains = self$chains, 
                     timing = list(up_to_iter = iter_time,
                                  total      = total_time)
                     )
      )
    },
    
    inferences = function() {
      stopifnot(length(self$chains) > 0)
      
      Beta   <- self$chains$Beta  
      s2_vec <- self$chains$sigma2      
      
      mu_b <- colMeans(Beta, na.rm = TRUE)    
      sd_b <- apply(Beta, 2, stats::sd, na.rm = TRUE) 
      qs_b <- apply(Beta, 2, stats::quantile, 
                    probs = c(0.025, 0.975), na.rm = TRUE)  
      CI_LI_b <- qs_b[1, ]
      CI_LS_b <- qs_b[2, ]
      
      contains_zero_b  <- (CI_LI_b <= 0 & 0 <= CI_LS_b) 
      contains_mean_b  <- (CI_LI_b <= mu_b & mu_b <= CI_LS_b)
      
      vals_beta <- list(
        mean           = drop(mu_b),
        sd             = drop(sd_b),
        CI_LI          = drop(CI_LI_b),
        CI_LS          = drop(CI_LS_b),
        contains_zero  = drop(contains_zero_b),
        contains_mean  = drop(contains_mean_b)
      )
      
      mean_s2 <- mean(s2_vec, na.rm = TRUE)
      sd_s2   <- stats::sd(s2_vec, na.rm = TRUE)
      qs_s2   <- stats::quantile(s2_vec, c(0.025, 0.975), na.rm = TRUE, names = FALSE)
      
      contains_zero_s2 <- (qs_s2[1] <= 0 & 0 <= qs_s2[2])
      contains_mean_s2 <- (qs_s2[1] <= mean_s2 & mean_s2 <= qs_s2[2])
      
      vals_sigma2 <- list(
        mean           = as.numeric(mean_s2),
        sd             = as.numeric(sd_s2),
        CI_LI          = qs_s2[1],
        CI_LS          = qs_s2[2],
        contains_zero  = contains_zero_s2,
        contains_mean  = contains_mean_s2
      )
      
      self$hat_values <- list(
        Beta   = vals_beta,
        sigma2 = vals_sigma2
      )
      invisible(self$hat_values)
    },
    
    goodness_of_fit = function(){
      waic <- waic_lrm(y = self$lrm_mcmc_object$stats$y, 
                            X = self$lrm_mcmc_object$stats$X,
                            chains = self$chains)
      dic <- dic_lrm(y = self$lrm_mcmc_object$stats$y,
                          X = self$lrm_mcmc_object$stats$X,
                          chains = self$chains)
      yhat <- yhat_mcmc_lrm(X = self$lrm_mcmc_object$stats$X,
                       hat_values = self$hat_values)
      mse <- mse(self$lrm_mcmc_object$stats$y,yhat)
      ppvalues <- ppval_lrm(y = self$lrm_mcmc_object$stats$y, 
                            X = self$lrm_mcmc_object$stats$X, 
                            chains = self$chains,
                            return_yrep = FALSE)
      r2 <- r2(self$lrm_mcmc_object$stats$y,yhat)
      
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
      invisible(yhat_mcmc_lrm(X = X_new, hat_values = self$hat_values))
    }
  )
)
