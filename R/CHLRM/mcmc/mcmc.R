source('R/utils.R')
source('R/CHLRM/utils.R')
source('R/CHLRM/mcmc/utils.R')

library(R6)

mcmc_chlrm <- R6::R6Class(
  "mcmc_chlrm",
  public = list(
    hyperparams     = NULL,
    run_settings    = NULL,
    chlrm_mcmc_object = NULL,
    chains          = NULL,
    hat_values      = NULL,
    goodness_fit   = NULL,
    
    initialize = function(y,X,K,standardize,hyperparams, run_settings) {
      self$run_settings     <- run_settings
      std_flag <- isTRUE(standardize)
      
      if (std_flag) {
        std <- standardize_Xy_chlrm(X, y)
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
          K   = K,
          standardized = std_flag,
          scale = scale
        )
      )
      self$chlrm_mcmc_object  <- Initialize_chlrm_mcmc(self$chlrm_mcmc_object)
      self$chains           <- list()
      self$hat_values       <- list()
      self$goodness_fit     <- list() 
    },
    
    fit = function(verbose = TRUE, seed = NULL) {
      if (!is.null(seed)) set.seed(seed)
      
      B <- n_samples_mcmc(self$run_settings)     # iter totales (incl. burn-in/thin)
      S <- self$run_settings$samples             # muestras guardadas tras thinning
      
      m <- self$chlrm_mcmc_object$stats$m
      p <- self$chlrm_mcmc_object$stats$p
      K <- self$chlrm_mcmc_object$stats$K
      
      # Cadenas a almacenar
      Betaj     <- array(NA_real_, dim = c(S, m, p))  # S × m × p
      sigma2j   <- matrix(NA_real_, nrow = S, ncol = m)   # S × m
      z_mat     <- matrix(NA_integer_, nrow = m, ncol = S) # m × S (grupos × draws)
      omega_mat <- matrix(NA_real_, nrow = S, ncol = K)    # S × K
      Betak_arr <- array(NA_real_, dim = c(S, K, p))       # S × K × p
      sigma2k_m <- matrix(NA_real_, nrow = S, ncol = K)    # S × K
      
      ll        <- numeric(S)
      lj        <- numeric(S)
      iter_time <- numeric(S)
      
      t0 <- proc.time()[["elapsed"]]
      if (verbose) pb <- txtProgressBar(min = 0, max = B, width = 50, char = "=", style = 3)
      
      for (b in seq_len(B)) {
        
        self$chlrm_mcmc_object <- sample_rhoj_chlrm(   self$chlrm_mcmc_object)   # z, n_k
        self$chlrm_mcmc_object <- sample_omega_chlrm(  self$chlrm_mcmc_object)   # omega 
        self$chlrm_mcmc_object <- sample_Betak_chlrm(  self$chlrm_mcmc_object)   # β_j
        self$chlrm_mcmc_object <- sample_sigma2k_chlrm(self$chlrm_mcmc_object)   # σ^2_k 
        self$chlrm_mcmc_object <- sample_Betak_chlrm(  self$chlrm_mcmc_object)   # β_k 
        self$chlrm_mcmc_object <- sample_Sigma_chlrm(  self$chlrm_mcmc_object)   # Σ
        self$chlrm_mcmc_object <- sample_xi2_chlrm(  self$chlrm_mcmc_object)     #xi2
        
        idx <- thin_chain_mcmc(b, self$run_settings)
        if (!is.null(idx)) {

          Betaj[idx, , ]    <- self$chlrm_mcmc_object$params$Betaj          # m × p
          sigma2j[idx, ]    <- self$chlrm_mcmc_object$params$sigma2j         # m
          z_mat[, idx]      <- self$chlrm_mcmc_object$params$z               # m (guardar como columna)
          omega_mat[idx, ]  <- self$chlrm_mcmc_object$params$omega           # K
          Betak_arr[idx, , ]<- self$chlrm_mcmc_object$params$Betak           # K × p
          sigma2k_m[idx, ]  <- self$chlrm_mcmc_object$params$sigma2k         # K
          
          iter_time[idx] <- proc.time()[["elapsed"]] - t0

          ll[idx] <- loglik_total_chlrm(
            y_list  = self$chlrm_mcmc_object$stats$y,
            X_list  = self$chlrm_mcmc_object$stats$X,
            Betaj   = self$chlrm_mcmc_object$params$Betaj,
            sigma2j = self$chlrm_mcmc_object$params$sigma2j
          )
          lj[idx] <- logprior_chlrm(self$chlrm_mcmc_object) + ll[idx]
        }
        
        if (verbose) setTxtProgressBar(pb, b)
      }
      
      if (verbose) close(pb)
      total_time <- proc.time()[["elapsed"]] - t0
      
      self$chains <- list(
        Betaj     = Betaj,       # S × m × p
        sigma2j   = sigma2j,     # S × m
        z         = z_mat,       # m × S (grupos × draws)
        omega     = omega_mat,   # S × K
        Betak     = Betak_arr,   # S × K × p
        sigma2k   = sigma2k_m,   # S × K
        loglik    = ll,
        logjoint  = lj
      )
      
      invisible(list(
        chains = self$chains,
        timing = list(up_to_iter = iter_time, total = total_time)
      ))
    },
    
    inferences = function() {
      Betak   <- self$chains$Betak    # S × K × p
      sigma2k <- self$chains$sigma2k  # S × K
      omega   <- self$chains$omega    # S × K
      z_mat   <- self$chains$z        # m × S
      
      ## ----- Betak (por cluster)
      mu_bk <- apply(Betak,  c(2, 3), mean, na.rm = TRUE)  # K × p
      sd_bk <- apply(Betak,  c(2, 3), stats::sd, na.rm = TRUE)
      qs_bk <- apply(Betak,  c(2, 3), function(v)
        stats::quantile(v, probs = c(0.025, 0.975), na.rm = TRUE, names = FALSE))
      CI_LI_bk <- qs_bk[1, , ]  # K × p
      CI_LS_bk <- qs_bk[2, , ]
      
      contains_zero_bk <- (CI_LI_bk <= 0 & 0 <= CI_LS_bk)
      contains_mean_bk <- (CI_LI_bk <= mu_bk & mu_bk <= CI_LS_bk)
      
      vals_Betak <- list(
        mean           = mu_bk,
        sd             = sd_bk,
        CI_LI          = CI_LI_bk,
        CI_LS          = CI_LS_bk,
        contains_zero  = contains_zero_bk,
        contains_mean  = contains_mean_bk
      )
      
      ## ----- sigma2k
      mu_s2k <- colMeans(sigma2k, na.rm = TRUE)           # K
      sd_s2k <- apply(sigma2k, 2, stats::sd, na.rm = TRUE)
      qs_s2k <- apply(sigma2k, 2, function(v)
        stats::quantile(v, probs = c(0.025, 0.975), na.rm = TRUE, names = FALSE))
      CI_LI_s2k <- qs_s2k[1, ]
      CI_LS_s2k <- qs_s2k[2, ]
      
      contains_zero_s2k <- (CI_LI_s2k <= 0 & 0 <= CI_LS_s2k)
      contains_mean_s2k <- (CI_LI_s2k <= mu_s2k & mu_s2k <= CI_LS_s2k)
      
      vals_sigma2k <- list(
        mean           = mu_s2k,
        sd             = sd_s2k,
        CI_LI          = CI_LI_s2k,
        CI_LS          = CI_LS_s2k,
        contains_zero  = contains_zero_s2k,
        contains_mean  = contains_mean_s2k
      )
      
      ## ----- omega
      mu_w <- colMeans(omega, na.rm = TRUE)               # K
      sd_w <- apply(omega, 2, stats::sd, na.rm = TRUE)
      qs_w <- apply(omega, 2, function(v)
        stats::quantile(v, probs = c(0.025, 0.975), na.rm = TRUE, names = FALSE))
      CI_LI_w <- qs_w[1, ]
      CI_LS_w <- qs_w[2, ]
      
      contains_zero_w <- (CI_LI_w <= 0 & 0 <= CI_LS_w)
      contains_mean_w <- (CI_LI_w <= mu_w & mu_w <= CI_LS_w)
      
      vals_omega <- list(
        mean           = mu_w,
        sd             = sd_w,
        CI_LI          = CI_LI_w,
        CI_LS          = CI_LS_w,
        contains_zero  = contains_zero_w,
        contains_mean  = contains_mean_w
      )
      
      #Cluster assignment
      m <- nrow(z_mat)
      S <- ncol(z_mat)
      K <- dim(Betak)[2]
      
      A <- incidence_matrix_mcmc_chlrm(t(z_mat))
      
      prob <- matrix(0, nrow = m, ncol = K)
      for (k in 1:K) prob[, k] <- rowMeans(z_mat == k)
      z_mode <- max.col(prob, ties.method = "first")
      
      nc <- apply(z_mat, 2, function(z) length(unique(z)))
      tab_nc  <- table(nc)
      k_vals   <- as.integer(names(tab_nc))
      prob_nc <- as.numeric(tab_nc) / S
      names(prob_nc) <- names(tab_nc)
      
      self$hat_values <- list(
        Betak   = vals_Betak,
        sigma2k = vals_sigma2k,
        omega   = vals_omega,
        z       = list(z_pred=z_mode,k_vals=k_vals,prob_nc=prob_nc,incidence_matrix=A)
      )
      invisible(self$hat_values)
    },
    
    goodness_of_fit = function(S = NULL,seed=NULL){
      if (is.null(S)) S <- self$run_settings$samples_gof
      
      # datos
      y_list <- self$chlrm_mcmc_object$stats$y
      X_list <- self$chlrm_mcmc_object$stats$X
      nj <- self$chlrm_mcmc_object$stats$nj
      
      # Dentro de goodness_of_fit(), justo antes de calcular waic/dic:
      chains <- self$chains
      
      # n_iter de referencia (usa loglik si está; si no, el primero)
      n_iter <- if (!is.null(chains$loglik)) length(chains$loglik) else {
        a1 <- chains[[1]]
        if (is.null(dim(a1))) length(a1) else dim(a1)[1L]
      }
      
      if (!is.null(seed)) set.seed(seed)
      S_eff <- min(as.integer(S), n_iter)
      idx <- sort(sample.int(n_iter, size = S_eff, replace = FALSE))
      
      subset_iter <- function(a, idx, n_iter){
        if (is.null(dim(a))) {
          if (length(a) == n_iter) return(a[idx])
          return(a)  # escalar u otra cosa: déjalo igual
        }
        d <- dim(a)
        if (d[1L] == n_iter) {
          # iteraciones ya en eje 1
          subs <- lapply(seq_along(d), function(k) if (k == 1L) idx else seq_len(d[k]))
          return(do.call(`[`, c(list(a), subs, list(drop = FALSE))))
        }
        if (any(d == n_iter)) {
          # mueve el eje de iteraciones a la 1ª posición, submuestrea y revierte
          ax <- which(d == n_iter)[1L]
          perm <- c(ax, setdiff(seq_along(d), ax))
          invp <- match(seq_along(perm), perm)
          a_perm <- aperm(a, perm)
          subs <- lapply(seq_along(dim(a_perm)), function(k) if (k == 1L) idx else seq_len(dim(a_perm)[k]))
          a_sub  <- do.call(`[`, c(list(a_perm), subs, list(drop = FALSE)))
          return(aperm(a_sub, invp))
        }
        # si no coincide ninguna dimensión con n_iter, déjalo como está
        a
      }
      
      chains_S <- lapply(chains, subset_iter, idx = idx, n_iter = n_iter)
      
      waic <- waic_chlrm(y_list = y_list, X_list = X_list,
                         chains = chains_S, nj = nj)
      dic  <- dic_chlrm(y_list = y_list, X_list = X_list,
                        chains = chains_S, nj = nj)
      
      yhat  <- yhat_mcmc_chlrm(X_list,self$hat_values)
      
      mse <- mse(unlist(y_list), unlist(yhat))
      r2 <- r2(unlist(y_list), unlist(yhat))
      
      ppvalues <- ppval_chlrm(
        y_list = y_list,
        X_list = X_list,
        chains = chains_S,
        nj = nj, 
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
      yhat_list <- yhat_mcmc_chlrm(X_list = X_new, hat_values = self$hat_values)
      if (isTRUE(as_vector)) return(unlist(yhat_list, use.names = FALSE))
      return(yhat_list)
    }
  )
)
