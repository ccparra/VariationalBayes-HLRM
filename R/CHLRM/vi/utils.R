source('R/utils.R')

Initialize_chlrm_vi <- function(chlrm_vi_object) {
  # tamaños
  p <- chlrm_vi_object$stats$p
  m <- chlrm_vi_object$stats$m
  K <- chlrm_vi_object$stats$K
  
  # hiperparámetros
  Mu0     <- chlrm_vi_object$hyperparams$Mu0
  Lambda0 <- chlrm_vi_object$hyperparams$Lambda0
  n0      <- chlrm_vi_object$hyperparams$n0
  S0      <- chlrm_vi_object$hyperparams$S0
  nu0     <- chlrm_vi_object$hyperparams$nu0
  a0      <- chlrm_vi_object$hyperparams$a0
  b0      <- chlrm_vi_object$hyperparams$b0
  alpha0  <- chlrm_vi_object$hyperparams$alpha0
  B_ols      <- chlrm_vi_object$hyperparams$Betaj_ols
  
  # ---------------------------
  # K-means sobre Betaj_ols
  # ---------------------------
  cl <- stats::kmeans(B_ols, centers = K)
  
  # responsabilidades (duras) y Dirichlet
  rho <- matrix(0, m, K); rho[cbind(seq_len(m), cl$cluster)] <- 1
  chlrm_vi_object$params$rho   <- rho
  chlrm_vi_object$params$alpha <- alpha0 + colSums(rho)
  
  # betas por clúster (medias = centros Kmeans)
  chlrm_vi_object$params$MuBetak <- cl$centers
  
  # ---------------------------
  # Inicializar SigmaBetak desde la prior IW de Sigma
  # Prior: Sigma ~ IW(n0, S0^{-1})  => rinvwishart(nu=n0, S=solve(S0))
  # ---------------------------
  Sigma0_draw <- LaplacesDemon::rinvwishart(nu = n0, S = solve(S0))
  Sigma0_draw <-.sym(Sigma0_draw) 
  
  chlrm_vi_object$params$SigmaBetak <- replicate(K, Sigma0_draw, simplify = FALSE)
  
  # Beta global ~ N(Mu0, Lambda0) 
  Sigma0_safe <- Lambda0 + 1e-8 * diag(p)
  mu_draw <- as.numeric(mvtnorm::rmvnorm(n = 1, mean = Mu0, sigma = .sym(Sigma0_safe)))
  chlrm_vi_object$params$MuBeta    <- mu_draw
  chlrm_vi_object$params$SigmaBeta <- Lambda0
  
  # Sigma global (q(Sigma) arranca en el prior)
  chlrm_vi_object$params$nuS <- n0 + K
  chlrm_vi_object$params$SS  <- S0
  
  # xi^2 ~ Gamma(a0, b0), y σ_k^2 | xi^2 ~ IG(nu0/2, nu0*xi^2/2)
  xi2_draw <- stats::rgamma(1, shape = a0, rate = b0)
  chlrm_vi_object$params$ak <- rep(0.5 * nu0, K)
  chlrm_vi_object$params$bk <- rep(0.5 * nu0 * xi2_draw, K)
  
  # q(xi^2) (forma cerrada inicial)
  chlrm_vi_object$params$axi2 <- a0 + 0.5 * nu0 * K
  Einvs2_sum <- sum(E_inv_sigma2(chlrm_vi_object$params$ak, chlrm_vi_object$params$bk))
  chlrm_vi_object$params$bxi2 <- b0 + 0.5 * nu0 * Einvs2_sum
  
  return(chlrm_vi_object)
}

# SSE_{k,j} con q(beta_k) = N(mu_k, Sigma_k)
.SSE_kj <- function(yj, Xj, XtXj, mu_k, Sig_k, yty_j) {
  term1 <- yty_j
  term2 <- -2 * drop(crossprod(mu_k, crossprod(Xj, yj)))
  term3 <- sum(XtXj * Sig_k)
  term4 <- drop(crossprod(mu_k, XtXj %*% mu_k))
  term1 + term2 + term3 + term4
}

# ================================================================
# 1) q(gamma_j) = Categorical(rho_j•)
#    log rho_{jk} ∝ E[log ω_k] - ½ n_j E[log σ_k^2] - ½ E[σ_k^{-2}] SSE_{k,j}
# ================================================================
update_qgammaj_chlrm <- function(chlrm_vi_object, batch_idx = NULL) {
  m  <- chlrm_vi_object$stats$m
  K  <- chlrm_vi_object$stats$K
  J  <- if (is.null(batch_idx)) seq_len(m) else batch_idx
  
  Elog_omega <- E_log_omega(chlrm_vi_object$params$alpha)
  Einv_s2_k  <- E_inv_sigma2(chlrm_vi_object$params$ak, chlrm_vi_object$params$bk)
  Elog_s2_k  <- E_log_sigma2(chlrm_vi_object$params$ak, chlrm_vi_object$params$bk)
  
  for (j in J) {
    Xj   <- chlrm_vi_object$stats$X[[j]]
    yj   <- chlrm_vi_object$stats$y[[j]]
    XtXj <- chlrm_vi_object$stats$XtX[[j]]
    ytyj <- chlrm_vi_object$stats$yty[j]
    nj   <- chlrm_vi_object$stats$nj[j]
    
    log_rho <- numeric(K)
    for (k in seq_len(K)) {
      mu_k  <- chlrm_vi_object$params$MuBetak[k, ]
      Sig_k <- chlrm_vi_object$params$SigmaBetak[[k]]
      sse_kj <- .SSE_kj(yj, Xj, XtXj, mu_k, Sig_k, ytyj)
      log_rho[k] <- Elog_omega[k] -
        0.5 * (nj * Elog_s2_k[k] + Einv_s2_k[k] * sse_kj)
    }
    lse <- logsumexp(log_rho)
    chlrm_vi_object$params$rho[j, ] <- exp(log_rho - lse)
  }
  return(chlrm_vi_object)
}

# ================================================================
# 2) q(beta_k) = N(mu_k, Sigma_k) con SVI (rho_t; CAVI ⇔ rho_t=1)
#    V_k = (E[Σ^{-1}] + E[σ_k^{-2}] * Σ_j ρ_{jk} Xj'Xj)^{-1}
#    m_k = V_k (E[Σ^{-1}] MuBeta + E[σ_k^{-2}] * Σ_j ρ_{jk} Xj'yj)
# ================================================================
update_qBetak_chlrm <- function(chlrm_vi_object, batch_idx = NULL, rho_t = 1) {
  m  <- chlrm_vi_object$stats$m
  p  <- chlrm_vi_object$stats$p
  K  <- chlrm_vi_object$stats$K
  J  <- if (is.null(batch_idx)) seq_len(m) else batch_idx
  w  <- m / length(J)
  

  E_SigInv <- .sym(chlrm_vi_object$params$nuS * solve(chlrm_vi_object$params$SS))
  EinvS    <- E_inv_sigma2(chlrm_vi_object$params$ak, chlrm_vi_object$params$bk)
  
  for (k in seq_len(K)) {
    Ak <- matrix(0, p, p)
    bvec <- rep(0, p)
    for (j in J) {
      rjk <- chlrm_vi_object$params$rho[j, k]
      if (rjk > 0) {
        Ak   <- Ak   + w * rjk * chlrm_vi_object$stats$XtX[[j]]
        bvec <- bvec + w * rjk * chlrm_vi_object$stats$Xty[[j]]
      }
    }
    V_hat <- solve(E_SigInv + EinvS[k] * Ak); V_hat <- .sym(V_hat)
    m_hat <- as.numeric(V_hat %*% (E_SigInv %*% chlrm_vi_object$params$MuBeta + EinvS[k] * bvec))
    
    # SVI
    Sig_old <- chlrm_vi_object$params$SigmaBetak[[k]]
    mu_old  <- chlrm_vi_object$params$MuBetak[k, ]
    chlrm_vi_object$params$SigmaBetak[[k]] <- .sym((1 - rho_t) * Sig_old + rho_t * V_hat)
    chlrm_vi_object$params$MuBetak[k, ]    <- (1 - rho_t) * mu_old + rho_t * m_hat
  }
  chlrm_vi_object
}


# ================================================================
# 3) q(σ_k^2) = IG(a_k,b_k) con SVI (rho_t; CAVI ⇔ rho_t=1)
#    a_hat = (ν0 + Σ_j ρ_{jk} n_j)/2
#    b_hat = 1/2 * [ ν0 * E[ξ^2] + Σ_j ρ_{jk} SSE_{k,j} ]
# ================================================================
update_qsigma2k_chlrm <- function(chlrm_vi_object, batch_idx = NULL, rho_t = 1) {
  m   <- chlrm_vi_object$stats$m
  K   <- chlrm_vi_object$stats$K
  J   <- if (is.null(batch_idx)) seq_len(m) else batch_idx
  w   <- m / length(J)
  
  nu0 <- chlrm_vi_object$hyperparams$nu0
  Exi2 <- chlrm_vi_object$params$axi2 / chlrm_vi_object$params$bxi2
  
  for (k in seq_len(K)) {
    Nobs_k <- 0
    SSE_k  <- 0
    mu_k   <- chlrm_vi_object$params$MuBetak[k, ]
    Sig_k  <- chlrm_vi_object$params$SigmaBetak[[k]]
    
    for (j in J) {
      rjk <- chlrm_vi_object$params$rho[j, k]
      if (rjk > 0) {
        Xj   <- chlrm_vi_object$stats$X[[j]]
        yj   <- chlrm_vi_object$stats$y[[j]]
        XtXj <- chlrm_vi_object$stats$XtX[[j]]
        ytyj <- chlrm_vi_object$stats$yty[j]
        
        sse_kj <- .SSE_kj(yj, Xj, XtXj, mu_k, Sig_k, ytyj)
        Nobs_k <- Nobs_k + w * rjk * chlrm_vi_object$stats$nj[j]
        SSE_k  <- SSE_k  + w * rjk * sse_kj
      }
    }
    a_hat <- 0.5 * (nu0 + Nobs_k)
    b_hat <- 0.5 * (nu0 * Exi2 + SSE_k)
    
    chlrm_vi_object$params$ak[k] <- (1 - rho_t) * chlrm_vi_object$params$ak[k] + rho_t * a_hat
    chlrm_vi_object$params$bk[k] <- (1 - rho_t) * chlrm_vi_object$params$bk[k] + rho_t * b_hat
  }
  return(chlrm_vi_object)
}

# ================================================================
# 4) q(ω) = Dir(alpha) con SVI (rho_t; CAVI ⇔ rho_t=1)
#    alpha_hat = alpha0 + Σ_j ρ_{jk}  (minibatch escalado)
# ================================================================
update_qomega_chlrm <- function(chlrm_vi_object, batch_idx = NULL, rho_t = 1) {
  m  <- chlrm_vi_object$stats$m
  K  <- chlrm_vi_object$stats$K
  J  <- if (is.null(batch_idx)) seq_len(m) else batch_idx
  w  <- m / length(J)
  
  alpha0 <- chlrm_vi_object$hyperparams$alpha0
  Nk_hat <- rep(0, K)
  for (j in J) Nk_hat <- Nk_hat + w * chlrm_vi_object$params$rho[j, ]
  alpha_hat <- alpha0 + Nk_hat
  
  chlrm_vi_object$params$alpha <- (1 - rho_t) * chlrm_vi_object$params$alpha + rho_t * alpha_hat
  return(chlrm_vi_object)
}

# ================================================================
# 5) q(β) global = N(MuBeta, SigmaBeta) con SVI
#    J_hat = Λ0^{-1} + K * E[Σ^{-1}]
#    h_hat = Λ0^{-1} μ0 + E[Σ^{-1}] * Σ_k μ_k
# ================================================================
update_qBeta_chlrm <- function(chlrm_vi_object, rho_t = 1) {
  p  <- chlrm_vi_object$stats$p
  K  <- chlrm_vi_object$stats$K
  
  mu0  <- chlrm_vi_object$hyperparams$Mu0
  Lam0 <- chlrm_vi_object$hyperparams$Lambda0
  J0   <- solve(Lam0)
  
  E_SigInv <- .sym(chlrm_vi_object$params$nuS * solve(chlrm_vi_object$params$SS))
  
  J_hat <- J0 + K * E_SigInv
  h_hat <- as.numeric(J0 %*% mu0 + E_SigInv %*% colSums(chlrm_vi_object$params$MuBetak))
  
  Sig_old <- chlrm_vi_object$params$SigmaBeta
  mu_old  <- chlrm_vi_object$params$MuBeta
  
  # Natural-param mixing
  J_old <- solve(Sig_old)
  h_old <- as.numeric(J_old %*% mu_old)
  
  J_new <- (1 - rho_t) * J_old + rho_t * J_hat
  h_new <- (1 - rho_t) * h_old + rho_t * h_hat
  
  Sig_new <- .sym(solve(J_new))
  mu_new  <- as.numeric(Sig_new %*% h_new)
  
  chlrm_vi_object$params$SigmaBeta <- Sig_new
  chlrm_vi_object$params$MuBeta    <- mu_new
  return(chlrm_vi_object)
}

# ================================================================
# 6) q(Σ) = IW(nuS, SS) con SVI (guardamos (nuS, SS) tal que E[Σ^{-1}] = nuS * solve(SS))
#    nu_hat = n0 + K
#    SS_hat = S0 + Σ_k [ Σ_beta_k + Σ_beta + (μ_k - μ_beta)(μ_k - μ_beta)' ]
# ================================================================
update_qSigma_chlrm <- function(chlrm_vi_object, rho_t = 1) {
  p  <- chlrm_vi_object$stats$p
  K  <- chlrm_vi_object$stats$K
  
  n0 <- chlrm_vi_object$hyperparams$n0
  S0 <- chlrm_vi_object$hyperparams$S0
  
  SigmaBeta <- chlrm_vi_object$params$SigmaBeta
  MuBeta    <- chlrm_vi_object$params$MuBeta
  
  SS_hat <- S0
  for (k in seq_len(K)) {
    Sig_k <- chlrm_vi_object$params$SigmaBetak[[k]]
    mu_k  <- chlrm_vi_object$params$MuBetak[k, ]
    diff  <- mu_k - MuBeta
    Mk    <- Sig_k + SigmaBeta + tcrossprod(diff, diff)
    SS_hat <- SS_hat + Mk
  }
  nu_hat <- n0 + K
  
  chlrm_vi_object$params$nuS <- (1 - rho_t) * chlrm_vi_object$params$nuS + rho_t * nu_hat
  chlrm_vi_object$params$SS  <- .sym((1 - rho_t) * chlrm_vi_object$params$SS + rho_t * SS_hat)
  
  return(chlrm_vi_object)
}

# ================================================================
# 7) q(ξ^2) = Gamma(a,b) (shape–rate) con SVI
#    a_hat = a0 + (K * ν0)/2
#    b_hat = b0 + (ν0/2) * Σ_k E[σ_k^{-2}]
# ================================================================
update_qxi2_chlrm <- function(chlrm_vi_object, rho_t = 1) {
  K   <- chlrm_vi_object$stats$K
  nu0 <- chlrm_vi_object$hyperparams$nu0
  a0  <- chlrm_vi_object$hyperparams$a0
  b0  <- chlrm_vi_object$hyperparams$b0
  
  EinvS_k <- E_inv_sigma2(chlrm_vi_object$params$ak, chlrm_vi_object$params$bk)
  
  a_hat <- a0 + 0.5 * K * nu0
  b_hat <- b0 + 0.5 * nu0 * sum(EinvS_k)
  
  chlrm_vi_object$params$axi2 <- (1 - rho_t) * chlrm_vi_object$params$axi2 + rho_t * a_hat
  chlrm_vi_object$params$bxi2 <- (1 - rho_t) * chlrm_vi_object$params$bxi2 + rho_t * b_hat
  return(chlrm_vi_object)
}


ELBO_chlrm <- function(chlrm_vi_object) {
  # -----------------------------
  # Stats y tamaños
  # -----------------------------
  X_list <- chlrm_vi_object$stats$X
  y_list <- chlrm_vi_object$stats$y
  XtX    <- chlrm_vi_object$stats$XtX
  Xty    <- chlrm_vi_object$stats$Xty
  yty    <- chlrm_vi_object$stats$yty
  nj     <- chlrm_vi_object$stats$nj
  m      <- chlrm_vi_object$stats$m
  p      <- chlrm_vi_object$stats$p
  K      <- chlrm_vi_object$stats$K
  
  # -----------------------------
  # Hiperparámetros
  # -----------------------------
  Mu0       <- chlrm_vi_object$hyperparams$Mu0
  Lambda0   <- chlrm_vi_object$hyperparams$Lambda0
  Lambda0inv <- solve(Lambda0)
  n0        <- chlrm_vi_object$hyperparams$n0
  S0        <- chlrm_vi_object$hyperparams$S0
  nu0       <- chlrm_vi_object$hyperparams$nu0
  a0        <- chlrm_vi_object$hyperparams$a0
  b0        <- chlrm_vi_object$hyperparams$b0
  alpha0    <- chlrm_vi_object$hyperparams$alpha0   # Dirichlet prior (K)
  
  # -----------------------------
  # Parámetros variacionales
  # -----------------------------
  # Mezcla
  rho    <- chlrm_vi_object$params$rho       # m x K
  alpha  <- chlrm_vi_object$params$alpha     # K (Dirichlet)
  
  # Clúster betas y varianzas
  MuBetak     <- chlrm_vi_object$params$MuBetak      # K x p
  SigmaBetak  <- chlrm_vi_object$params$SigmaBetak   # list K de p×p
  ak          <- chlrm_vi_object$params$ak           # K (IG shape)
  bk          <- chlrm_vi_object$params$bk           # K (IG scale)
  
  # Global beta
  MuBeta     <- chlrm_vi_object$params$MuBeta
  SigmaBeta  <- chlrm_vi_object$params$SigmaBeta
  
  # Global Sigma (IW) parametrizada por (nuS, SS) tal que E[Sigma^{-1}] = nuS * solve(SS)
  nuS        <- chlrm_vi_object$params$nuS
  SS         <- chlrm_vi_object$params$SS
  
  # Global xi^2 ~ Gamma(axi2, bxi2) (shape–rate)
  axi2       <- chlrm_vi_object$params$axi2
  bxi2       <- chlrm_vi_object$params$bxi2
  
  # -----------------------------
  # Expectativas útiles
  # -----------------------------
  Elog_sigma2_k <- E_log_sigma2(ak, bk)            # vector K
  Einv_sigma2_k <- E_inv_sigma2(ak, bk)            # vector K
  ESigmainv     <- .sym(nuS * solve(SS))                 # E[Σ^{-1}]
  ElogSigma     <- ElogdetIW(nuS, SS)              # E[log|Σ|]
  Exi2          <- axi2 / bxi2
  Elogxi2       <- digamma(axi2) - log(bxi2)
  Elog_omega    <- E_log_omega(alpha)
  logB_alpha0   <- lgamma(sum(alpha0)) - sum(lgamma(alpha0))  # const Dir prior
  logB_alpha    <- lgamma(sum(alpha))  - sum(lgamma(alpha))   # const Dir var

  # -----------------------------
  # 1) Verosimilitud: E_q[ log p(y | z, beta_k, sigma_k^2) ]
  #     = sum_j sum_k rho_jk [ -1/2 * n_j (log 2π + E log σ_k^2) - 1/2 * E σ_k^{-2} * SSE_{k,j} ]
  # -----------------------------
  Elogpy <- 0
  for (j in seq_len(m)) {
    XtXj <- XtX[[j]]
    Xtyj <- Xty[[j]]
    ytyj <- yty[j]
    nj_j <- nj[j]
    
    # para cada k, SSE_{k,j} con productos precomputados
    for (k in seq_len(K)) {
      mu_k  <- MuBetak[k, ]
      Sig_k <- SigmaBetak[[k]]
      
      # SSE_{k,j} = y'y - 2 mu' X'y + tr(X'X Sig_k) + mu' X'X mu
      term1 <- ytyj
      term2 <- -2 * sum(mu_k * Xtyj)
      term3 <- sum(XtXj * Sig_k)
      term4 <- drop(t(mu_k) %*% XtXj %*% mu_k)
      sse_kj <- term1 + term2 + term3 + term4
      
      Elogpy <- Elogpy +
        rho[j, k] * ( -0.5 * nj_j * (log2pi + Elog_sigma2_k[k]) - 0.5 * Einv_sigma2_k[k] * sse_kj )
    }
  }
  
  # -----------------------------
  # 2) Capa jerárquica: E_q[ log p(beta_k | Beta, Sigma) ]
  #     = sum_k [ -1/2 ( p log 2π + E log|Σ| ) - 1/2 tr( E[Σ^{-1}] * M_k ) ]
  #     donde M_k = Var_q(beta_k) + Var_q(Beta) + (E beta_k - E Beta)(...)
  # -----------------------------
  ElogpBetak <- 0
  for (k in seq_len(K)) {
    mu_k  <- MuBetak[k, ]
    Sig_k <- SigmaBetak[[k]]
    d_k   <- mu_k - MuBeta
    Mk    <- Sig_k + SigmaBeta + tcrossprod(d_k)
    ElogpBetak <- ElogpBetak + ( -0.5 * (p * log2pi + ElogSigma) - 0.5 * sum(diag(ESigmainv %*% Mk)) )
  }
  
  # -----------------------------
  # 3) Prior global Beta: E_q[ log p(Beta) ]
  #     = -1/2 ( p log 2π + log|Λ0| + tr(Λ0^{-1}(Σ_beta + (μ_beta-μ0)(...) )) )
  # -----------------------------
  dB <- MuBeta - Mu0
  ElogpBeta <- -0.5 * ( p * log2pi + logdet(Lambda0) +
                          sum(diag(Lambda0inv %*% (SigmaBeta + tcrossprod(dB)))) )
  
  # -----------------------------
  # 4) Prior global Sigma ~ IW(n0, S0^{-1})
  #     log C_IW + ...
  # -----------------------------
  logC_IW <- 0.5 * n0 * logdet(S0) - 0.5 * n0 * p * log(2) - log_multigamma(0.5 * n0, p)
  ElogpSigma <- logC_IW - 0.5 * (n0 + p + 1) * ElogSigma - 0.5 * sum(diag(S0 %*% ESigmainv))
  
  # -----------------------------
  # 5) Prior de σ_k^2 | ξ^2  (sin σ0^2):
  #     α = ν0/2, β(ξ) = ν0 ξ^2 / 2
  #     E_q[log p(σ_k^2 | ξ^2)] = α[log(ν0/2)+Elogξ^2] - lgamma(α)
  #                              - (α+1) Elog σ_k^2 - (ν0/2) Eξ^2 * E σ_k^{-2}
  # -----------------------------
  alpha_ig <- 0.5 * nu0
  Elogps2k <- 0
  for (k in seq_len(K)) {
    Elogps2k <- Elogps2k +
      ( alpha_ig * (log(0.5 * nu0) + Elogxi2) - lgamma(alpha_ig) -
          (alpha_ig + 1) * Elog_sigma2_k[k] - 0.5 * nu0 * Exi2 * Einv_sigma2_k[k] )
  }
  
  # -----------------------------
  # 6) Prior de ξ^2 ~ Gamma(a0, b0) (shape–rate):
  #     a0 log b0 - lgamma(a0) + (a0-1) Elogξ^2 - b0 Eξ^2
  # -----------------------------
  Elogpxi2 <- a0 * log(b0) - lgamma(a0) + (a0 - 1) * Elogxi2 - b0 * Exi2
  
  # -----------------------------
  # 7) Mezcla: E_q[ log p(z | ω) ] + E_q[ log p(ω) ]
  # -----------------------------
  # E[log p(z|ω)] = Σ_j Σ_k ρ_{jk} E[log ω_k]
  Elogpz <- sum(drop(rho %*% Elog_omega))
  
  # E[log p(ω)] = log B(α0) + Σ_k (α0_k - 1) E[log ω_k]
  Elogpomega <- logB_alpha0 + sum( (alpha0 - 1) * Elog_omega )
  
  # -----------------------------
  # 8) Entropías (E_q[log q(·)])
  # -----------------------------
  # q(beta_k) normales
  ElogqBetak <- 0
  for (k in seq_len(K)) {
    ElogqBetak <- ElogqBetak + ( -0.5 * ( p * (1 + log2pi) + logdet(SigmaBetak[[k]]) ) )
  }
  # q(Beta) normal
  ElogqBeta <- -0.5 * ( p * (1 + log2pi) + logdet(SigmaBeta) )
  
  # q(Sigma) inverse-Wishart: misma forma que prior con (nuS, SS)
  logC_qIW  <- 0.5 * nuS * logdet(SS) - 0.5 * nuS * p * log(2) - log_multigamma(0.5 * nuS, p)
  ElogqSigma <- logC_qIW - 0.5 * (nuS + p + 1) * ElogSigma - 0.5 * sum(diag(SS %*% ESigmainv))
  
  # q(σ_k^2) inverse-gamma (shape–scale): Σ_k [ a_k log b_k - lgamma(a_k) - (a_k+1) Elog σ_k^2 - b_k E σ_k^{-2} ]
  Elogqs2k <- 0
  for (k in seq_len(K)) {
    Elogqs2k <- Elogqs2k + ( ak[k] * log(bk[k]) - lgamma(ak[k]) -
                               (ak[k] + 1) * Elog_sigma2_k[k] - bk[k] * Einv_sigma2_k[k] )
  }
  
  # q(ξ^2) gamma (shape–rate): a log b - lgamma(a) + (a-1) Elogξ^2 - b Eξ^2
  Elogqxi2 <- axi2 * log(bxi2) - lgamma(axi2) + (axi2 - 1) * Elogxi2 - bxi2 * Exi2
  
  # q(z): Σ_j Σ_k ρ_{jk} log ρ_{jk}
  Elogqz <- sum(rho * log(pmax(rho, .Machine$double.eps)))
  
  # q(ω) Dirichlet: log B(α) + Σ_k (α_k - 1) E[log ω_k]
  Elogqomega <- logB_alpha + sum( (alpha - 1) * Elog_omega )
  
  # -----------------------------
  # ELBO = Σ E_q[log p(·)]  -  Σ E_q[log q(·)]
  # -----------------------------
  Elogp_total <- Elogpy + ElogpBetak + ElogpBeta + ElogpSigma + Elogps2k + Elogpxi2 + Elogpz + Elogpomega
  Elogq_total <- ElogqBetak + ElogqBeta + ElogqSigma + Elogqs2k + Elogqxi2 + Elogqz + Elogqomega
  
  elbo <- Elogp_total - Elogq_total
  return(as.numeric(elbo))
}

vi_sample_draws_chlrm <- function(chlrm_vi_object, results, S, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  # Tamaños
  m <- chlrm_vi_object$stats$m
  p <- chlrm_vi_object$stats$p
  K <- chlrm_vi_object$stats$K
  
  # Acepta list(params=...) o list(...) plano
  params <- if (!is.null(results$params)) results$params else results
  
  # Parámetros VI necesarios
  MuBetak     <- params$MuBetak                 # K x p
  SigmaBetak  <- lapply(params$SigmaBetak, as.matrix)  # list K de p x p
  ak          <- as.numeric(params$ak)          # K
  bk          <- as.numeric(params$bk)          # K
  rho         <- params$rho                     # m x K  (responsabilidades por grupo)
  
  # Helpers mínimos
  .rcat <- function(prob) sample.int(length(prob), size = 1L, prob = prob)
  
  # Prealocación
  Betak_arr <- array(NA_real_, dim = c(S, K, p))  # S×K×p
  s2k_mat   <- matrix(NA_real_, nrow = S, ncol = K)
  z_mat     <- matrix(NA_integer_, nrow = m, ncol = S)  # m×S
  Betaj     <- array(NA_real_, dim = c(S, m, p))        # S×m×p
  s2j_mat   <- matrix(NA_real_, nrow = S, ncol = m)     # S×m
  
  # Muestreo por clúster (independiente en s)
  for (k in seq_len(K)) {
    Betak_arr[, k, ] <- mvtnorm::rmvnorm(
      n      = S,
      mean   = MuBetak[k, ],
      sigma  = SigmaBetak[[k]],
      method = "chol"
    )
    # σ_k^2 ~ IG(a_k, b_k)  (shape–rate) ⇒ σ^2 = 1 / Gamma(a_k, rate=b_k)
    s2k_mat[, k] <- 1 / stats::rgamma(n = S, shape = ak[k], rate = bk[k])
  }
  
  # Asignaciones por grupo con rho_j y parámetros inducidos
  for (s in seq_len(S)) {
    for (j in seq_len(m)) {
      z_js <- .rcat(rho[j, ])       # z_j | y ~ Cat(rho_j)
      z_mat[j, s] <- z_js
      Betaj[s, j, ] <- Betak_arr[s, z_js, ]
      s2j_mat[s, j] <- s2k_mat[s, z_js]
    }
  }
  
  list(
    Betaj   = Betaj,      # S × m × p
    sigma2j = s2j_mat,    # S × m
    z       = z_mat,      # m × S
    Betak   = Betak_arr,  # S × K × p
    sigma2k = s2k_mat     # S × K
  )
}




yhat_vi_chlrm <- function(X_list, hat_values) {

  Betak <- hat_values$Betak$mean
  z     <- hat_values$z$z_pred
  
  m <- length(X_list)
  yhat_list <- vector("list", m)
  
  for (j in seq_len(m)) {
    k <- z[j]
    yhat_list[[j]] <- as.vector(as.matrix(X_list[[j]]) %*% Betak[k, ])
  }
  yhat_list
}
