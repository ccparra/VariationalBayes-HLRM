source('R/utils.R')
source('R/CHLRM/utils.R')
source('R/CHLRM/vi/utils.R')

library(R6)

vi_chlrm <- R6::R6Class("vi_chlrm",
                       public = list(
                         hyperparams = NULL,
                         run_settings = NULL,
                         chlrm_vi_object = NULL,
                         results = NULL,
                         hat_values = NULL,
                         goodness_fit = NULL,
                         
                         initialize = function(y, X, K, standarize, hyperparams, run_settings){
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
                           self$chlrm_vi_object <- list(
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
                           
                           self$chlrm_vi_object <- Initialize_chlrm_vi(self$chlrm_vi_object)
                           
                           self$results      <- list()
                           self$hat_values   <- list()
                           self$goodness_fit <- list()
                         },
                         
                         
                         fit = function(verbose = TRUE) {
                           # --------------------------------------------
                           # Run settings
                           # --------------------------------------------
                           max_iter    <- self$run_settings$iterations
                           eps         <- self$run_settings$epsilon
                           stochastic  <- isTRUE(self$run_settings$stochastic)
                           minibatch_m <- if (!is.null(self$run_settings$minibatch_m)) self$run_settings$minibatch_m else self$chlrm_vi_object$stats$m
                           tau         <- if (!is.null(self$run_settings$tau))   self$run_settings$tau   else 100
                           kappa       <- if (!is.null(self$run_settings$kappa)) self$run_settings$kappa else 0.7
                           
                           # --------------------------------------------
                           # Book-keeping
                           # --------------------------------------------
                           elbo <- numeric(max_iter)
                           iter_time <- numeric(max_iter)
                           t0 <- proc.time()[["elapsed"]]
                           if (verbose) pb <- txtProgressBar(min = 0, max = max_iter, width = 50, char = "=", style = 3)
                           
                           m <- self$chlrm_vi_object$stats$m
                           elbo_old <- -Inf
                           
                           # --------------------------------------------
                           # Main loop (CAVI ⇔ stochastic=FALSE; SVI ⇔ TRUE)
                           # --------------------------------------------
                           for (iter in seq_len(max_iter)) {
                             
                             if (stochastic) {
                               size  <- min(minibatch_m, m)
                               J     <- sample.int(m, size = size, replace = FALSE)
                               rho_t <- step_size(iter, tau = tau, kappa = kappa)
                             } else {
                               J     <- NULL
                               rho_t <- 1
                             }
                             
                             # ===== CHLRM updates =====
                             # 1) Responsabilidades locales q(z_j)
                             self$chlrm_vi_object <- update_qgammaj_chlrm(self$chlrm_vi_object, batch_idx = J)
                             
                             # 2) Bloques por clúster (SVI con rho_t y minibatch J)
                             self$chlrm_vi_object <- update_qBetak_chlrm(self$chlrm_vi_object, batch_idx = J, rho_t = rho_t)
                             self$chlrm_vi_object <- update_qsigma2k_chlrm(self$chlrm_vi_object, batch_idx = J, rho_t = rho_t)
                             self$chlrm_vi_object <- update_qomega_chlrm(self$chlrm_vi_object,  batch_idx = J, rho_t = rho_t)
                             
                             # 3) Globales (SVI en parámetros naturales, no dependen de J)
                             self$chlrm_vi_object <- update_qBeta_chlrm(self$chlrm_vi_object,  rho_t = rho_t)
                             self$chlrm_vi_object <- update_qSigma_chlrm(self$chlrm_vi_object, rho_t = rho_t)
                             self$chlrm_vi_object <- update_qxi2_chlrm(self$chlrm_vi_object,   rho_t = rho_t)
                             
                             # ===== ELBO & timing =====
                             elbo_new <- ELBO_chlrm(self$chlrm_vi_object)
                             elbo[iter] <- elbo_new
                             iter_time[iter] <- proc.time()[["elapsed"]] - t0
                             
                             if (verbose) setTxtProgressBar(pb, iter)
                             
                             # Convergencia (solo CAVI determinístico)
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
                             elbo       = elbo
                           )
                           
                           invisible(list(
                             params = self$results,
                             timing = list(up_to_iter = iter_time, total = total_time)
                           ))
                         },
                         
                         inferences = function() {
                           # Tamaños
                           K <- self$chlrm_vi_object$stats$K
                           p <- self$chlrm_vi_object$stats$p
                           m <- self$chlrm_vi_object$stats$m
                           
                           # Resultados del fit() (VI)
                           MuBetak    <- self$results$MuBetak        # K × p
                           SigmaBetak <- self$results$SigmaBetak     # list length K de p×p
                           ak         <- self$results$ak             # K (IG shape)
                           bk         <- self$results$bk             # K (IG scale)
                           alpha      <- self$results$alpha          # K (Dirichlet)
                           rho        <- self$results$rho            # m × K
                           
                           # -------------------------------
                           # Betak: sd y CIs normales
                           # -------------------------------
                           z975 <- stats::qnorm(0.975)
                           sd_bk <- t(vapply(seq_len(K), function(k) sqrt(diag(SigmaBetak[[k]])), numeric(p)))
                           CI_LI_bk <- MuBetak - z975 * sd_bk
                           CI_LS_bk <- MuBetak + z975 * sd_bk
                           contains_zero_bk <- (CI_LI_bk <= 0 & 0 <= CI_LS_bk)
                           contains_mean_bk <- (CI_LI_bk <= MuBetak & MuBetak <= CI_LS_bk)
                           
                           vals_Betak <- list(
                             mean           = MuBetak,            # K × p
                             sd             = sd_bk,              # K × p
                             CI_LI          = CI_LI_bk,           # K × p
                             CI_LS          = CI_LS_bk,           # K × p
                             contains_zero  = contains_zero_bk,   # K × p (lógico)
                             contains_mean  = contains_mean_bk    # K × p (lógico)
                           )
                           
                           # -------------------------------
                           # sigma2k ~ IG(a_k, b_k) (shape–scale)
                           #  E[s2] = b/(a-1) (a>1)
                           #  SD[s2]= b/((a-1)*sqrt(a-2)) (a>2)
                           #  CI via 1/qgamma(·; shape=a, rate=b)
                           # -------------------------------
                           mean_s2k  <- ifelse(ak > 1, bk / (ak - 1), NA_real_)
                           sd_s2k    <- ifelse(ak > 2, bk / ((ak - 1) * sqrt(ak - 2)), NA_real_)
                           CI_LI_s2k <- 1 / stats::qgamma(0.975, shape = ak, rate = bk)
                           CI_LS_s2k <- 1 / stats::qgamma(0.025, shape = ak, rate = bk)
                           
                           vals_sigma2k <- list(
                             mean           = mean_s2k,           # K
                             sd             = sd_s2k,             # K
                             CI_LI          = CI_LI_s2k,          # K
                             CI_LS          = CI_LS_s2k,          # K
                             contains_zero  = (CI_LI_s2k <= 0 & 0 <= CI_LS_s2k),
                             contains_mean  = (CI_LI_s2k <= mean_s2k & mean_s2k <= CI_LS_s2k)
                           )
                           
                           # -------------------------------
                           # omega ~ Dirichlet(alpha)
                           #  E[ω_k] = α_k / α0
                           #  Var(ω_k) = α_k(α0-α_k) / [α0^2(α0+1)]
                           #  Marginal CI para ω_k via Beta(α_k, α0-α_k)
                           # -------------------------------
                           alpha0 <- sum(alpha)
                           omega_mean <- alpha / alpha0
                           sd_w <- sqrt( (alpha * (alpha0 - alpha)) / (alpha0^2 * (alpha0 + 1)) )
                           CI_LI_w <- mapply(function(a, b) stats::qbeta(0.025, shape1 = a, shape2 = b),
                                             alpha, alpha0 - alpha)
                           CI_LS_w <- mapply(function(a, b) stats::qbeta(0.975, shape1 = a, shape2 = b),
                                             alpha, alpha0 - alpha)
                           
                           vals_omega <- list(
                             mean           = omega_mean,         # K
                             sd             = sd_w,               # K
                             CI_LI          = CI_LI_w,            # K
                             CI_LS          = CI_LS_w,            # K
                             contains_zero  = (CI_LI_w <= 0 & 0 <= CI_LS_w),
                             contains_mean  = (CI_LI_w <= omega_mean & omega_mean <= CI_LS_w)
                           )
                           
                           # -------------------------------
                           # z: MAP por fila de rho + matriz de incidencia 0/1
                           # -------------------------------
                           z_pred <- max.col(rho, ties.method = "first")    # m-vector (1..K)
                           
                           A_soft <- rho %*% t(rho)
                           diag(A_soft) <- 1
                           
                           vals_z <- list(
                             prob = rho,                   # m × K (responsabilidades)
                             z_pred = z_pred,              # m
                             incidence_matrix = A_soft          # m × K (0/1)
                           )
                           
                           # -------------------------------
                           # Guardar hat-values
                           # -------------------------------
                           self$hat_values <- list(
                             Betak   = vals_Betak,
                             sigma2k = vals_sigma2k,
                             omega   = vals_omega,
                             z       = vals_z
                           )
                           invisible(self$hat_values)
                         },
                         
                         goodness_of_fit = function(S = NULL, seed = NULL){
                           if (is.null(S)) S <- self$run_settings$samples_gof
                           y <- self$chlrm_vi_object$stats$y
                           X <- self$chlrm_vi_object$stats$X
                           nj <- self$chlrm_vi_object$stats$nj
                           
                           draws <- vi_sample_draws_chlrm(self$chlrm_vi_object, self$results, S = S, seed = seed)
                           chains <- list(Betaj = draws$Betaj, sigma2j = draws$sigma2j)
                           
                           waic <- waic_chlrm(chains = chains, nj = nj,
                                              y_list = y, X_list = X)
                           dic  <- dic_chlrm( chains = chains, nj = nj,
                                              y_list = y, X_list = X)
                           
                           yhat <- yhat_vi_chlrm(X, self$hat_values)
                           mse   <- mse(unlist(y), unlist(yhat))
                           r2   <- r2(unlist(y), unlist(yhat))
                           
                           ppvalues <- ppval_chlrm(y_list = y,
                                                  X_list= X,
                                                  chains = chains,
                                                  nj = nj,
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
                         
                         predict = function(X_new,as_vector=FALSE) {
                           stopifnot(length(self$hat_values) > 0)
                           yhat_list <- yhat_vi_chlrm(X_list = X_new, hat_values = self$hat_values)
                           if (isTRUE(as_vector)) return(unlist(yhat_list, use.names = FALSE))
                           return(yhat_list)
                         }
                       )
)
