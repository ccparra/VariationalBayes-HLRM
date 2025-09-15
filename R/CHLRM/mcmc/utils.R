Initialize_chlrm_mcmc <- function(chlrm_mcmc_object) {
  Mu0     <- chlrm_mcmc_object$hyperparams$Mu0
  Lambda0 <- chlrm_mcmc_object$hyperparams$Lambda0
  S0      <- chlrm_mcmc_object$hyperparams$S0
  nu0     <- chlrm_mcmc_object$hyperparams$nu0
  n0      <- chlrm_mcmc_object$hyperparams$n0 
  a0      <- chlrm_mcmc_object$hyperparams$a0 
  b0      <- chlrm_mcmc_object$hyperparams$b0
  alpha0  <- chlrm_mcmc_object$hyperparams$alpha0
  
  p <- chlrm_mcmc_object$stats$p
  m <- chlrm_mcmc_object$stats$m
  K <- chlrm_mcmc_object$stats$K
  
  omega <- drop(gtools::rdirichlet(1, alpha0))
  z <- sample.int(K, size = m, replace = TRUE, prob = omega)
  nk <- tabulate(z, nbins = K) 
  Beta  <- drop(mvtnorm::rmvnorm(n = 1, mean = Mu0, sigma = Lambda0))
  Sigma <- LaplacesDemon::rinvwishart(nu = n0, S = S0)
  xi2   <- stats::rgamma(1, shape = a0, rate = b0)

  Betak  <- mvtnorm::rmvnorm(n = K, mean = Beta, sigma = Sigma)
  sigma2k <- 1 / stats::rgamma(n = K, shape = 0.5 * nu0, rate = 0.5 * nu0 * xi2)

  Betaj   <- matrix(NA_real_, nrow = m, ncol = p)
  sigma2j <- numeric(m)
  for (j in seq_len(m)) {
    k <- z[j]
    Betaj[j, ] <- drop(mvtnorm::rmvnorm(1, mean = Betak[k, ], sigma = Sigma))
    sigma2j[j] <- sigma2k[k]
  }  
  
  chlrm_mcmc_object$params$omega   <- omega
  chlrm_mcmc_object$params$z       <- z 
  chlrm_mcmc_object$params$nk       <- nk 
  chlrm_mcmc_object$params$Beta    <- Beta
  chlrm_mcmc_object$params$Sigma   <- Sigma
  chlrm_mcmc_object$params$xi2     <- xi2
  chlrm_mcmc_object$params$Betak   <- Betak
  chlrm_mcmc_object$params$sigma2k <- sigma2k
  chlrm_mcmc_object$params$Betaj   <- Betaj
  chlrm_mcmc_object$params$sigma2j <- sigma2j
  
  return(chlrm_mcmc_object)
}

sample_rhoj_chlrm <- function(chlrm_mcmc_object) {
  K    <- chlrm_mcmc_object$stats$K
  m    <- chlrm_mcmc_object$stats$m
  nj   <- chlrm_mcmc_object$stats$nj
  XtX  <- chlrm_mcmc_object$stats$XtX    # lista p×p
  Xty  <- chlrm_mcmc_object$stats$Xty    # lista p×1
  y    <- chlrm_mcmc_object$stats$y
  yty  <- chlrm_mcmc_object$stats$yty

  
  Betak    <- chlrm_mcmc_object$params$Betak   # K×p
  omega    <- chlrm_mcmc_object$params$omega   # K
  sigma2k  <- chlrm_mcmc_object$params$sigma2k # K
  
  log_2pi      <- log(2 * pi)
  log_sigma2k  <- log(sigma2k)
  inv_sigma2k  <- 1 / sigma2k
  
  z  <- integer(m)
  B  <- t(Betak)  # p×K (para multiplicaciones vectorizadas)
  
  for (j in seq_len(m)) {
    XtXj <- XtX[[j]]      # p×p
    Xtyj <- Xty[[j]]      # p×1
    ytyj <- yty[j]
    
    # l_vec_k = β_k^T X_j^T y_j  (K×1)
    l_vec <- as.numeric(Betak %*% Xtyj)
    
    # q_vec_k = β_k^T X_j^T X_j β_k para todos los k a la vez:
    # M = XtXj %*% B (p×K); q_k = colSums(B * M)
    M     <- XtXj %*% B
    q_vec <- colSums(B * M)
    
    SSR_vec <- ytyj - 2 * l_vec + q_vec
    
    llk <- -0.5 * nj[j] * (log_2pi + log_sigma2k) - 0.5 * SSR_vec * inv_sigma2k
    log_post <- log(omega) + llk
    
    # estabilidad numérica
    log_post <- log_post - max(log_post)
    pk <- exp(log_post); pk <- pk / sum(pk)
    
    z[j] <- sample.int(K, size = 1, prob = pk)
  }
  
  nk <- tabulate(z, nbins = K)
  chlrm_mcmc_object$params$z       <- z
  chlrm_mcmc_object$params$nk      <- nk
  chlrm_mcmc_object$params$sigma2j <- chlrm_mcmc_object$params$sigma2k[z]
  return(chlrm_mcmc_object)
}

sample_omega_chlrm <- function(chlrm_mcmc_object) {

  alpha0 <- chlrm_mcmc_object$hyperparams$alpha0
  nk     <- chlrm_mcmc_object$params$nk  # conteos por clúster (longitud K)
  
  post_alpha <- alpha0 + nk
  chlrm_mcmc_object$params$omega <- drop(gtools::rdirichlet(1, post_alpha))
  return(chlrm_mcmc_object)
}

sample_Betak_chlrm <- function(chlrm_mcmc_object) {
  K       <- chlrm_mcmc_object$stats$K
  p       <- chlrm_mcmc_object$stats$p
  XtX     <- chlrm_mcmc_object$stats$XtX   # lista de matrices p×p
  Xty     <- chlrm_mcmc_object$stats$Xty   # lista de vectores/matrices p×1
  z       <- chlrm_mcmc_object$params$z    # asignaciones de grupo
  Sigma   <- chlrm_mcmc_object$params$Sigma
  SigInv  <- solve(Sigma)
  Beta    <- chlrm_mcmc_object$params$Beta
  sigma2k <- chlrm_mcmc_object$params$sigma2k
  Betak   <- chlrm_mcmc_object$params$Betak
  
  bcol <- matrix(Beta, ncol = 1)           # p×1
  
  for (k in seq_len(K)) {
    idx <- which(z == k)
    if (length(idx) == 0L) {
      XtXk <- matrix(0, p, p)
      Xtyk <- matrix(0, p, 1)
    } else if (length(idx) == 1L) {
      XtXk <- chlrm_mcmc_object$stats$XtX[[idx]]
      Xtyk <- chlrm_mcmc_object$stats$Xty[[idx]]
    } else {
      XtXk <- Reduce(`+`, chlrm_mcmc_object$stats$XtX[idx])
      Xtyk <- Reduce(`+`, chlrm_mcmc_object$stats$Xty[idx])
    }
    
    Prec_k <- (1 / sigma2k[k]) * XtXk + SigInv
    Vk     <- solve(Prec_k); Vk <- .sym(Vk)
    muk    <- Vk %*% ((1 / sigma2k[k]) * Xtyk + SigInv %*% bcol)
    
    Betak[k, ] <- drop(mvtnorm::rmvnorm(1, mean = drop(muk), sigma = Vk))
  }
  
  Betaj <- Betak[z,]
  
  chlrm_mcmc_object$params$Betak <- Betak
  chlrm_mcmc_object$params$Betaj <- Betaj
  
  return(chlrm_mcmc_object)
}

sample_sigma2k_chlrm <- function(chlrm_mcmc_object) {
  # Stats
  X   <- chlrm_mcmc_object$stats$X      # lista X_j (n_j x p) — no la usamos aquí
  y   <- chlrm_mcmc_object$stats$y      # lista y_j
  nj  <- chlrm_mcmc_object$stats$nj     # longitudes por grupo
  K   <- chlrm_mcmc_object$stats$K
  XtX <- chlrm_mcmc_object$stats$XtX    # lista p x p
  Xty <- chlrm_mcmc_object$stats$Xty    # lista p x 1
  yty <- chlrm_mcmc_object$stats$yty    #1x p 
  
  # Parámetros
  z     <- chlrm_mcmc_object$params$z
  Betaj <- chlrm_mcmc_object$params$Betaj
  xi2   <- chlrm_mcmc_object$params$xi2
  nu0   <- chlrm_mcmc_object$hyperparams$nu0
  
  sigma2k <- numeric(K)
  
  for (k in seq_len(K)) {
    idx <- which(z == k)
    if (length(idx) > 0) {
      n_k_obs <- sum(nj[idx])
      ssr_k <- 0
      for (j in idx) {
        bj   <- Betaj[j, ]
        xtyj <- Xty[[j]]          # p x 1
        xtxj <- XtX[[j]]          # p x p
        
        l_j  <- sum(bj * xtyj)                 # β_j^T X_j^T y_j
        q_j  <- sum(bj * (xtxj %*% bj))        # β_j^T X_j^T X_j β_j
        
        ssr_k <- ssr_k + (yty[j] - 2 * l_j + q_j)
      }
      
      a_k <- 0.5 * (nu0 + n_k_obs)
      b_k <- 0.5 * (nu0 * xi2 + ssr_k)
    } else {
      # sin grupos en el cluster: prior
      a_k <- 0.5 * nu0
      b_k <- 0.5 * nu0 * xi2
    }
    
    sigma2k[k] <- 1 / stats::rgamma(1, shape = a_k, rate = b_k)
  }
  
  chlrm_mcmc_object$params$sigma2k <- sigma2k
  chlrm_mcmc_object$params$sigma2j <- sigma2k[z]
  
  return(chlrm_mcmc_object)
}

sample_xi2_chlrm <- function(chlrm_mcmc_object){
  a0   <- chlrm_mcmc_object$hyperparams$a0
  b0   <- chlrm_mcmc_object$hyperparams$b0
  nu0  <- chlrm_mcmc_object$hyperparams$nu0
  
  sigma2k <- chlrm_mcmc_object$params$sigma2k 
  nk <- chlrm_mcmc_object$params$nk

  active <- nk > 0
  kappa  <- sum(active)
  invsum <- sum(1 / sigma2k[active])
  
  shape <- a0 + 0.5 * nu0 * kappa
  rate  <- b0 + 0.5 * nu0 * invsum
  
  chlrm_mcmc_object$params$xi2 <- stats::rgamma(1, shape = shape, rate = rate)
  return(chlrm_mcmc_object)
}

sample_Beta_chlrm <- function(chlrm_mcmc_object){
  # Hiperparámetros
  Mu0     <- chlrm_mcmc_object$hyperparams$Mu0      # μ0 (p)
  Lambda0 <- chlrm_mcmc_object$hyperparams$Lambda0  # Λ0 (p×p)
  
  # Parámetros actuales
  Sigma <- chlrm_mcmc_object$params$Sigma           # Σ (p×p)
  Betak <- chlrm_mcmc_object$params$Betak           # matriz K×p
  nk <- chlrm_mcmc_object$params$nk
  
  active <- nk > 0                                  # nk(K)
  kappa  <- sum(active)
  

  Lam0Inv <- solve(Lambda0)
  SigInv  <- solve(Sigma)
  
  sumBetak <- colSums(Betak[active, , drop = FALSE])  # ∑_{k:nk>0} β_k (p)
  
  Prec <- Lam0Inv + kappa * SigInv
  V    <- solve(Prec); V <- 0.5 * (V + t(V))      
  bvec <- matrix(sumBetak, ncol = 1)
  mu   <- V %*% (Lam0Inv %*% Mu0 + SigInv %*% bvec)
  
  chlrm_mcmc_object$params$Beta <- drop(
    mvtnorm::rmvnorm(n = 1, mean = drop(mu), sigma = V)
  )
  return(chlrm_mcmc_object)
}

sample_Sigma_chlrm <- function(chlrm_mcmc_object){
  # Hiperparámetros
  S0 <- chlrm_mcmc_object$hyperparams$S0   # p×p
  n0 <- chlrm_mcmc_object$hyperparams$n0
  
  # Parámetros actuales
  Betak <- chlrm_mcmc_object$params$Betak  # K×p
  Beta  <- chlrm_mcmc_object$params$Beta   # p
  nk <- chlrm_mcmc_object$params$nk       #K

  active <- nk > 0
  kappa  <- sum(active)
  
  #  sum_k (βk-β)(βk-β)^T
  D  <- sweep(Betak[active, , drop = FALSE], 2, Beta, "-")  # kappa×p
  SS <- S0 + crossprod(D)                                   # p×p
  nuS <- n0 + kappa   

  # Muestreo Inv-Wishart
  Sigma_draw <- LaplacesDemon::rinvwishart(nu = nuS, S = SS)
  Sigma_draw <- .sym(Sigma_draw)  
  
  chlrm_mcmc_object$params$Sigma <- Sigma_draw
  return(chlrm_mcmc_object)
}

yhat_mcmc_chlrm <- function(X_list, hat_values) {
  
  Betak <- hat_values$Betak$mean
  z     <- hat_values$z$z_pred
  
  m <- length(X_list)
  yhat_list <- vector("list", m)
  
  for (j in seq_len(m)) {
    k <- z[j]
    yhat_list[[j]] <- as.vector(X_list[[j]] %*% Betak[k, ])
  }
  yhat_list
}

logprior_chlrm <- function(chlrm_mcmc_object) {
  # --- Hyper ---
  Mu0      <- chlrm_mcmc_object$hyperparams$Mu0        # p
  Lambda0  <- chlrm_mcmc_object$hyperparams$Lambda0    # p×p
  n0       <- chlrm_mcmc_object$hyperparams$n0
  S0       <- chlrm_mcmc_object$hyperparams$S0         # p×p
  nu0      <- chlrm_mcmc_object$hyperparams$nu0
  a0       <- chlrm_mcmc_object$hyperparams$a0
  b0       <- chlrm_mcmc_object$hyperparams$b0
  alpha0   <- chlrm_mcmc_object$hyperparams$alpha0     # K-vector (Dirichlet)
  
  # --- Params ---
  Betak    <- chlrm_mcmc_object$params$Betak           # K×p
  Beta     <- chlrm_mcmc_object$params$Beta            # p
  Sigma    <- chlrm_mcmc_object$params$Sigma           # p×p
  sigma2k  <- chlrm_mcmc_object$params$sigma2k         # K
  xi2      <- chlrm_mcmc_object$params$xi2             # scalar
  omega    <- chlrm_mcmc_object$params$omega           # K
  z        <- chlrm_mcmc_object$params$z               # length m (1..K)
  
  K <- chlrm_mcmc_object$stats$K
  
  # 1) β_k | β, Σ  (producto sobre k)
  lp_betak <- sum(apply(Betak, 1L, function(bk)
    mvtnorm::dmvnorm(x = bk, mean = Beta, sigma = Sigma, log = TRUE)
  ))
  
  # 2) β | μ0, Λ0
  lp_beta <- mvtnorm::dmvnorm(x = Beta, mean = Mu0, sigma = Lambda0, log = TRUE)
  
  # 3) Σ ~ IW(n0, S0)
  lp_Sigma <- LaplacesDemon::dinvwishart(Sigma = Sigma, nu = n0, S = S0, log = TRUE)
  
  # 4) σ_k^2 | ξ^2  ~ IG(ν0/2, ν0 ξ^2/2)
  #    usando 1/σ_k^2 ~ Gamma(ν0/2, ν0 ξ^2/2) + jacobiano log(1/σ^2) = -2 log σ
  lp_s2k <- sum(
    stats::dgamma(x = 1 / sigma2k,
                  shape = 0.5 * nu0,
                  rate  = 0.5 * nu0 * xi2,
                  log   = TRUE) +
      log(1 / (sigma2k^2))
  )
  
  # 5) ξ^2 ~ G(a0, b0)
  lp_xi2 <- stats::dgamma(x = xi2, shape = a0, rate = b0, log = TRUE)
  
  # 6) ω ~ Dir(α0)
  #    (usa forma cerrada para evitar dependencias)
  lp_omega <- lgamma(sum(alpha0)) - sum(lgamma(alpha0)) +
    sum((alpha0 - 1) * log(omega))
  
  # 7) z | ω  (producto de ω_{z_j})
  lp_z <- sum(log(omega[z]))
  
  as.numeric(lp_betak + lp_beta + lp_Sigma + lp_s2k + lp_xi2 + lp_omega + lp_z)
}
incidence_matrix_mcmc_chlrm <- function(Z) {
  # Z: S x m  (filas = draws, columnas = grupos)
  S <- nrow(Z); n <- ncol(Z)
  A <- matrix(0, n, n)
  for (s in seq_len(S)) {
    labs <- Z[s, ]
    A <- A + (outer(labs, labs, "==") * 1.0)  # co-clustering por draw
  }
  A <- A / S
  diag(A) <- 1
  A
}
