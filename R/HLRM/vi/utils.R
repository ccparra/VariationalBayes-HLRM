source('R/utils.R')

Initialize_hlrm_vi <- function(hlrm_vi_object) {
  # tamaños
  p <- hlrm_vi_object$stats$p
  m <- hlrm_vi_object$stats$m
  
  # hiperparámetros
  Mu0     <- hlrm_vi_object$hyperparams$Mu0
  Lambda0 <- hlrm_vi_object$hyperparams$Lambda0
  n0      <- hlrm_vi_object$hyperparams$n0
  S0      <- hlrm_vi_object$hyperparams$S0
  nu0     <- hlrm_vi_object$hyperparams$nu0
  a0      <- hlrm_vi_object$hyperparams$a0
  b0      <- hlrm_vi_object$hyperparams$b0
  B_ols   <- hlrm_vi_object$hyperparams$Betaj_ols  # siempre viene
  
  # ---------------------------
  # Inicialización local por grupo (j=1..m)
  # ---------------------------
  # Medias q(β_j): usar directamente OLS de cada grupo
  hlrm_vi_object$params$MuBetaj    <- B_ols                           # m × p
  # Covarianzas q(β_j): arrancar con la misma Λ0 para todos
  Sigma0_draw <- LaplacesDemon::rinvwishart(nu = n0, S = solve(S0))
  Sigma0_draw <- .sym(Sigma0_draw)
  hlrm_vi_object$params$SigmaBetaj <- replicate(m, Sigma0_draw, simplify = FALSE)
  
  # ---------------------------
  # Globales: Beta y Sigma jerárquica
  # ---------------------------
  # q(Beta) ~ N(MuBeta, SigmaBeta)
  Sigma0_safe <- Lambda0 + 1e-8 * diag(p)
  MuBeta_draw <- as.numeric(mvtnorm::rmvnorm(
    n     = 1,
    mean  = Mu0,
    sigma = .sym(Sigma0_safe)
  ))
  hlrm_vi_object$params$MuBeta    <- MuBeta_draw                 # p
  hlrm_vi_object$params$SigmaBeta <- Lambda0                          # p×p
  
  # q(Sigma) ~ IW(nuS, SS) 
  hlrm_vi_object$params$nuS <- n0 + m
  hlrm_vi_object$params$SS  <- S0
  
  # ---------------------------
  # Varianzas por grupo y escala global ξ²
  # ---------------------------
  # σ_j^2 | ξ^2 ~ IG(ν0/2, (ν0 ξ^2)/2)
  xi2_draw <- stats::rgamma(1, shape = a0, rate = b0)
  
  hlrm_vi_object$params$aj <- rep(0.5 * nu0, m)
  hlrm_vi_object$params$bj <- rep(0.5 * nu0 * xi2_draw, m)
  
  # q(ξ^2) ~ Gamma(axi2, bxi2): inicial en forma cerrada con E[σ_j^{-2}]
  hlrm_vi_object$params$axi2 <- a0 + 0.5 * nu0 * m
  Einvs2_sum <- sum(E_inv_sigma2(hlrm_vi_object$params$aj, hlrm_vi_object$params$bj))
  hlrm_vi_object$params$bxi2 <- b0 + 0.5 * nu0 * Einvs2_sum
  
  return(hlrm_vi_object)
}

update_qbetaj_hlrm <- function(hlrm_vi_object, batch_idx = NULL) {
  J <- if (is.null(batch_idx)) seq_len(hlrm_vi_object$stats$m) else batch_idx
  E_SigInv <- hlrm_vi_object$params$nuS * solve(hlrm_vi_object$params$SS)
  
  for (j in J) {
    Einvs2j <- hlrm_vi_object$params$aj[j] / hlrm_vi_object$params$bj[j]
    Vj   <- solve(E_SigInv + Einvs2j * hlrm_vi_object$stats$XtX[[j]])
    Vj   <- .sym(Vj)
    mj   <- Vj %*% (E_SigInv %*% hlrm_vi_object$params$MuBeta +
                      Einvs2j * crossprod(hlrm_vi_object$stats$X[[j]], hlrm_vi_object$stats$y[[j]]))
    hlrm_vi_object$params$SigmaBetaj[[j]] <- Vj
    hlrm_vi_object$params$MuBetaj[j, ]    <- as.numeric(mj)
  }
  return(hlrm_vi_object)
}

update_qsigma2j_hlrm <- function(hlrm_vi_object, batch_idx = NULL) {
  J   <- if (is.null(batch_idx)) seq_len(hlrm_vi_object$stats$m) else batch_idx
  Exi2 <- hlrm_vi_object$params$axi2 / hlrm_vi_object$params$bxi2
  nu0   <- hlrm_vi_object$hyperparams$nu0
  
  for (j in J) {
    Sj   <- hlrm_vi_object$params$SigmaBetaj[[j]]
    mj   <- hlrm_vi_object$params$MuBetaj[j, ]
    yj   <- hlrm_vi_object$stats$y[[j]]
    Xj   <- hlrm_vi_object$stats$X[[j]]
    XtXj <- hlrm_vi_object$stats$XtX[[j]]
    
    ssr_mean <- drop(crossprod(yj - Xj %*% mj))
    tr_term  <- sum(diag(XtXj %*% (Sj + tcrossprod(mj))))
    Essr     <- ssr_mean + tr_term
    
    hlrm_vi_object$params$aj[j] <- 0.5 * (nu0 + hlrm_vi_object$stats$nj[j])
    hlrm_vi_object$params$bj[j] <- 0.5 * (nu0 * Exi2 + Essr)
  }
  return(hlrm_vi_object)
}

update_qBeta_hlrm <- function(hlrm_vi_object, rho_t = 1, batch_idx = NULL) {
  m        <- hlrm_vi_object$stats$m
  J        <- if (is.null(batch_idx)) seq_len(m) else batch_idx
  cscale   <- m / length(J)
  
  Lambda0  <- hlrm_vi_object$hyperparams$Lambda0
  Lambda0inv <- solve(Lambda0)
  Mu0      <- hlrm_vi_object$hyperparams$Mu0
  E_SigInv <- hlrm_vi_object$params$nuS * solve(hlrm_vi_object$params$SS)
  
  Prec_hat <- Lambda0inv + cscale * length(J) * E_SigInv
  h_hat    <- Lambda0inv %*% Mu0 + E_SigInv %*% (cscale * colSums(hlrm_vi_object$params$MuBetaj[J, , drop = FALSE]))
  
  Prec_old <- solve(hlrm_vi_object$params$SigmaBeta)
  h_old    <- Prec_old %*% hlrm_vi_object$params$MuBeta
  
  Prec_new <- (1 - rho_t) * Prec_old + rho_t * Prec_hat
  h_new    <- (1 - rho_t) * h_old    + rho_t * h_hat
  
  SigmaB <- solve(Prec_new); SigmaB <- .sym(SigmaB)
  MuB    <- SigmaB %*% h_new
  
  hlrm_vi_object$params$SigmaBeta <- SigmaB
  hlrm_vi_object$params$MuBeta    <- as.numeric(MuB)
  return(hlrm_vi_object)
}

update_qSigma_hlrm <- function(hlrm_vi_object, rho_t = 1, batch_idx = NULL) {
  m      <- hlrm_vi_object$stats$m
  J      <- if (is.null(batch_idx)) seq_len(m) else batch_idx
  cscale <- m / length(J)
  
  S0 <- hlrm_vi_object$hyperparams$S0
  MuB <- hlrm_vi_object$params$MuBeta
  SB  <- hlrm_vi_object$params$SigmaBeta
  
  SumT <- Reduce(`+`, lapply(J, function(j) {
    mj <- hlrm_vi_object$params$MuBetaj[j, ]
    Sj <- hlrm_vi_object$params$SigmaBetaj[[j]]
    (Sj + tcrossprod(mj)) - tcrossprod(mj, MuB) - tcrossprod(MuB, mj) + (SB + tcrossprod(MuB))
  }))
  
  S_hat  <- S0 + cscale * SumT
  S_old  <- hlrm_vi_object$params$SS
  
  hlrm_vi_object$params$SS  <- .sym((1 - rho_t) * S_old  + rho_t * S_hat)
  return(hlrm_vi_object)
}

update_qxi2_hlrm <- function(hlrm_vi_object, rho_t = 1, batch_idx = NULL) {
  m      <- hlrm_vi_object$stats$m
  J      <- if (is.null(batch_idx)) seq_len(m) else batch_idx
  cscale <- m / length(J)
  
  a0 <- hlrm_vi_object$hyperparams$a0
  b0 <- hlrm_vi_object$hyperparams$b0
  nu0 <- hlrm_vi_object$hyperparams$nu0
  
  Einvs2j_sum <- sum(hlrm_vi_object$params$aj[J] / hlrm_vi_object$params$bj[J])
  
  b_hat <- b0 + 0.5 * nu0 * cscale * Einvs2j_sum
  eta2_old <- - hlrm_vi_object$params$bxi
  eta2_hat <- - b_hat
  
  eta2_new <- (1 - rho_t) * eta2_old + rho_t * eta2_hat
  
  hlrm_vi_object$params$bxi2 <- - eta2_new
  return(hlrm_vi_object)
}

ELBO_hlrm <- function(hlrm_vi_object) {
  
  X_list <- hlrm_vi_object$stats$X
  y_list <- hlrm_vi_object$stats$y
  XtX    <- hlrm_vi_object$stats$XtX
  nj     <- hlrm_vi_object$stats$nj
  m      <- hlrm_vi_object$stats$m
  p      <- hlrm_vi_object$stats$p
  
  Mu0     <- hlrm_vi_object$hyperparams$Mu0
  Lambda0 <- hlrm_vi_object$hyperparams$Lambda0
  Lambda0inv <- solve(Lambda0)
  n0      <- hlrm_vi_object$hyperparams$n0
  S0      <- hlrm_vi_object$hyperparams$S0
  nu0     <- hlrm_vi_object$hyperparams$nu0
  a0      <- hlrm_vi_object$hyperparams$a0
  b0      <- hlrm_vi_object$hyperparams$b0
  
  MuBetaj    <- hlrm_vi_object$params$MuBetaj      
  SigmaBetaj <- hlrm_vi_object$params$SigmaBetaj     
  aj         <- hlrm_vi_object$params$aj             
  bj         <- hlrm_vi_object$params$bj             
  
  MuBeta     <- hlrm_vi_object$params$MuBeta         
  SigmaBeta  <- hlrm_vi_object$params$SigmaBeta     
  
  nuS        <- hlrm_vi_object$params$nuS           
  SS         <- hlrm_vi_object$params$SS             
  
  axi2        <- hlrm_vi_object$params$axi2 
  bxi2        <- hlrm_vi_object$params$bxi2
  
  
  Elogpy <- 0
  for (j in seq_len(m)) {
    mj   <- MuBetaj[j, ]
    Sj   <- SigmaBetaj[[j]]
    yj   <- y_list[[j]]
    Xj   <- X_list[[j]]
    XtXj <- XtX[[j]]
    
    Elogs2j <- E_log_sigma2(aj[j], bj[j])
    Einvs2j <- E_inv_sigma2(aj[j], bj[j])
    
    ssrj_mean <- drop(crossprod(yj - Xj %*% mj))
    tr_Sj     <- sum(diag(XtXj %*% Sj))
    Essrj   <- ssrj_mean + tr_Sj
    
    Elogpy <- Elogpy + ( -0.5 * nj[j] * (log2pi + Elogs2j) - 0.5 * Einvs2j * Essrj )
    
  }
  
  ESigmainv   <- nuS * solve(SS)
  ElogSigma <- ElogdetIW(nuS, SS) 
  ElogpBetaj <- 0
  for (j in seq_len(m)) {
    dj <- MuBetaj[j, ] - MuBeta
    Cj <- SigmaBetaj[[j]] + SigmaBeta + tcrossprod(dj)
    ElogpBetaj <- ElogpBetaj +
      ( -0.5 * (p * log2pi + ElogSigma) - 0.5 * sum(diag(ESigmainv %*% Cj)) )
  }
  
  dB <- MuBeta - Mu0
  ElogpBeta <- -0.5 * ( p * log2pi + logdet(Lambda0) +
                           sum(diag(Lambda0inv %*% (SigmaBeta + tcrossprod(dB)))))
  
  
  logC_IW <- 0.5 * n0 * logdet(S0) - 0.5 * n0 * p * log(2) - log_multigamma(0.5 * n0, p)
  ElogpSigma <- logC_IW - 0.5 * (n0 + p + 1) * ElogSigma - 0.5 * sum(diag(S0 %*% ESigmainv))
  
  
  alpha <- 0.5 * nu0
  Exi2      <- axi2 / bxi2
  Elogxi2   <- digamma(axi2) - log(bxi2)
  Elogps2j <- 0
  for (j in seq_len(m)) {
    Elogs2j <- E_log_sigma2(aj[j], bj[j])
    Einvs2j <- E_inv_sigma2(aj[j], bj[j])
    Elogps2j <- Elogps2j +
      ( alpha * (log(0.5 * nu0) + Elogxi2) - lgamma(alpha) -
          (alpha + 1) * Elogs2j - 0.5 * nu0 * Exi2 * Einvs2j )
  }
  

  Elogpxi2 <- a0 * log(b0) - lgamma(a0) + (a0 - 1) * Elogxi2 - b0 * Exi2
  
  ElogqBetaj <- 0
  for (j in seq_len(m)) {
    ElogqBetaj <- ElogqBetaj + ( -0.5 * ( p * (1 + log2pi) + logdet(SigmaBetaj[[j]])))
  }
  
  ElogqBeta <- -0.5 * ( p * (1 + log2pi) + logdet(SigmaBeta) )
  
  logC_qIW <- 0.5 * nuS * logdet(SS) - 0.5 * nuS * p * log(2) - log_multigamma(0.5 * nuS, p)
  ElogqSigma <- logC_qIW - 0.5 * (nuS + p + 1) * ElogSigma - 0.5 * sum(diag(SS %*% ESigmainv))
  
  
  Elogqs2j <- 0
  for (j in seq_len(m)) {
    Elogs2j <- E_log_sigma2(aj[j], bj[j])
    Einvs2j <- E_inv_sigma2(aj[j], bj[j])
    Elogqs2j <- Elogqs2j + ( aj[j] * log(bj[j]) - lgamma(aj[j]) -
                                 (aj[j] + 1) * Elogs2j - bj[j] * Einvs2j )
  }
  
  Elogqxi2 <- axi2 * log(bxi2) - lgamma(axi2) + (axi2 - 1) * Elogxi2 - bxi2 * Exi2
  
  elbo <- (Elogpy + ElogpBetaj + ElogpBeta + ElogpSigma + Elogps2j + Elogpxi2) -
    (ElogqBetaj + ElogqBeta + ElogqSigma + Elogqs2j + Elogqxi2)
  
  return(as.numeric(elbo))
}

# Minimal VI draws: only Betaj and sigma2j
vi_sample_draws_hlrm <- function(hlrm_vi_object, results, S, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  params   <- if (!is.null(results$params)) results$params else results
  MuBetaj  <- params$MuBetaj        # m x p
  SigList  <- params$SigmaBetaj     # list length m of p x p
  aj       <- params$aj             # length m
  bj       <- params$bj             # length m
  
  m <- nrow(MuBetaj)
  p <- ncol(MuBetaj)
  
  Betaj   <- array(NA_real_, dim = c(S, m, p))  # S x m x p
  sigma2j <- matrix(NA_real_, nrow = S, ncol = m)
  
  for (j in seq_len(m)) {
    Betaj[, j, ] <- mvtnorm::rmvnorm(n = S, mean = MuBetaj[j, ], sigma = SigList[[j]], method = "chol")
    sigma2j[, j] <- 1 / stats::rgamma(n = S, shape = aj[j], rate = bj[j])  # IG(a_j, b_j)
  }
  
  list(Betaj = Betaj, sigma2j = sigma2j)
}



yhat_vi_hlrm <-function(X_list, hat_values) {
  B_hat <- hat_values$Betaj$mean
  lapply(seq_along(X_list), function(j) {
    as.vector(as.matrix(X_list[[j]]) %*% drop(B_hat[j, ]))
  })
}
