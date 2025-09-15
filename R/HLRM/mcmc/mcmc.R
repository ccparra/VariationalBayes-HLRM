source('R/utils.R')
source('R/HLRM/utils.R')
source('R/HLRM/mcmc/utils.R')

library(R6)

mcmc_hlrm <- R6::R6Class(
  "mcmc_hlrm",
  public = list(
    hyperparams     = NULL,
    run_settings    = NULL,
    hlrm_mcmc_object = NULL,
    chains          = NULL,
    hat_values      = NULL,
    goodness_fit   = NULL,
    
    initialize = function(y,X,hyperparams, run_settings) {
      self$run_settings     <- run_settings
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
      self$hlrm_mcmc_object  <- Initialize_hlrm_mcmc(self$hlrm_mcmc_object)
      self$chains           <- list()
      self$hat_values       <- list()
      self$goodness_fit     <- list() 
    },
    
    fit = function(verbose = TRUE, seed = NULL) {
      if (!is.null(seed)) set.seed(seed)
      B <- n_samples_mcmc(self$run_settings)
      S <- self$run_settings$samples
      m <- self$hlrm_mcmc_object$stats$m
      p <- self$hlrm_mcmc_object$stats$p
      
      Betaj <- array(NA_real_, dim = c(S, m, p))
      sigma2j  <- array(NA_real_, dim = c(S, m))
      ll     <- numeric(S)
      lj     <- numeric(S)
      iter_time <- numeric(S)
      
      t0 <- proc.time()[["elapsed"]]
      
      if (verbose) pb <- txtProgressBar(min = 0, max = B, width = 50, char = "=", style = 3)
      
      for (b in seq_len(B)) {
        self$hlrm_mcmc_object <- sample_Betaj_hlrm(self$hlrm_mcmc_object)
        self$hlrm_mcmc_object <- sample_sigma2j_hlrm(self$hlrm_mcmc_object)
        self$hlrm_mcmc_object <- sample_xi2_hlrm(self$hlrm_mcmc_object)
        self$hlrm_mcmc_object <- sample_Beta_hlrm(self$hlrm_mcmc_object)
        self$hlrm_mcmc_object <- sample_Sigma_hlrm(self$hlrm_mcmc_object)
        
        idx <- thin_chain_mcmc(b, self$run_settings)
        if (!is.null(idx)) {
          Betaj[idx, , ]      <- self$hlrm_mcmc_object$params$Betaj
          sigma2j[idx,]      <- self$hlrm_mcmc_object$params$sigma2j
          iter_time[idx]   <- proc.time()[["elapsed"]] - t0
          
          ll[idx] <- loglik_total_hlrm(
            y     = self$hlrm_mcmc_object$stats$y,
            X     = self$hlrm_mcmc_object$stats$X,
            betaj  = self$hlrm_mcmc_object$params$Betaj,
            sigma2j = self$hlrm_mcmc_object$params$sigma2j
          )
          lj[idx] <- logprior_hlrm(self$hlrm_mcmc_object) + ll[idx]
        }
        
        if (verbose) setTxtProgressBar(pb, b)
      }
      if (verbose) close(pb)
      total_time <- proc.time()[["elapsed"]] - t0
      
      self$chains <- list(
        Betaj     = Betaj,
        sigma2j   = sigma2j,
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
      Betaj   <- self$chains$Betaj      
      sigma2j <- self$chains$sigma2j    

      mu_bj <- apply(Betaj, c(2, 3), mean, na.rm = TRUE)  
      sd_bj <- apply(Betaj, c(2, 3), stats::sd, na.rm = TRUE)
      qs_bj <- apply(
        Betaj, c(2, 3),
        function(v) stats::quantile(v, probs = c(0.025, 0.975), na.rm = TRUE, names = FALSE)
      ) 
      CI_LI_bj <- qs_bj[1, , ]          
      CI_LS_bj <- qs_bj[2, , ]
      
      contains_zero_bj <- (CI_LI_bj <= 0 & 0 <= CI_LS_bj)    
      contains_mean_bj <- (CI_LI_bj <= mu_bj & mu_bj <= CI_LS_bj) 
      
      vals_betaj <- list(
        mean          = mu_bj,
        sd            = sd_bj,
        CI_LI         = CI_LI_bj,
        CI_LS         = CI_LS_bj,
        contains_zero = contains_zero_bj,
        contains_mean = contains_mean_bj
      )
      

      mu_s2j <- colMeans(sigma2j, na.rm = TRUE)           
      sd_s2j <- apply(sigma2j, 2, stats::sd, na.rm = TRUE)
      qs_s2j <- apply(
        sigma2j, 2,
        function(v) stats::quantile(v, probs = c(0.025, 0.975), na.rm = TRUE, names = FALSE)
      ) 
      CI_LI_s2j <- qs_s2j[1, ]
      CI_LS_s2j <- qs_s2j[2, ]
      
      contains_zero_s2j <- (CI_LI_s2j <= 0 & 0 <= CI_LS_s2j)
      contains_mean_s2j <- (CI_LI_s2j <= mu_s2j & mu_s2j <= CI_LS_s2j)
      
      vals_sigma2j <- list(
        mean          = mu_s2j,
        sd            = sd_s2j,
        CI_LI         = CI_LI_s2j,
        CI_LS         = CI_LS_s2j,
        contains_zero = contains_zero_s2j,
        contains_mean = contains_mean_s2j
      )
      
      self$hat_values <- list(
        Betaj   = vals_betaj,
        sigma2j = vals_sigma2j
      )
      invisible(self$hat_values)
    },
    
    goodness_of_fit = function(S = NULL,seed=NULL){
      if (is.null(S)) S <- self$run_settings$samples_gof
      if (!is.null(seed)) set.seed(seed)
      
      # datos
      y_list <- self$hlrm_mcmc_object$stats$y
      X_list <- self$hlrm_mcmc_object$stats$X
      
      # longitud común de las cadenas (n_iter): tomo la primera entrada
      first_el <- self$chains[[1]]
      n_iter <- if (is.null(dim(first_el))) length(first_el) else dim(first_el)[1]
      
      # índices aleatorios (sin reemplazo); si S > n_iter, tomo todo
      S_eff <- min(as.integer(S), n_iter)
      idx <- sort(sample.int(n_iter, size = S_eff, replace = FALSE))
      
      # subsetting genérico por 1ra dimensión
      chains_S <- lapply(self$chains, function(a) {
        if (is.null(dim(a))) {
          a[idx]
        } else {
          d <- dim(a)
          subs <- vector("list", length(d))
          subs[[1]] <- idx
          if (length(d) >= 2L) {
            for (k in 2:length(d)) subs[[k]] <- seq_len(d[k])
          }
          do.call(`[`, c(list(a), subs, list(drop = FALSE)))
        }
      })
      
      waic <- waic_hlrm(y_list = y_list, X_list = X_list, chains = chains_S)
      dic  <- dic_hlrm( y_list = y_list, X_list = X_list, chains = chains_S)
      
      yhat  <- yhat_mcmc_hlrm(X_list,self$hat_values)
      
      mse <- mse(unlist(y_list), unlist(yhat))
      r2 <- r2(unlist(y_list), unlist(yhat))
      
      ppvalues <- ppval_hlrm(
        y_list = y_list,
        X_list = X_list,
        chains = self$chains,
        return_yrep = FALSE
      )
      
      self$goodness_fit <- list(
        waic = waic,
        dic  = dic,
        mse  = mse,
        r2 = r2,
        ppval = ppvalues
      )
      invisible(self$goodness_fit)
    },

    predict = function(X_new, as_vector = FALSE) {
      stopifnot(length(self$hat_values) > 0)
      yhat_list <- yhat_mcmc_hlrm(X_list = X_new, hat_values = self$hat_values)
      if (isTRUE(as_vector)) return(unlist(yhat_list, use.names = FALSE))
     return( yhat_list)
    }
  )
)
