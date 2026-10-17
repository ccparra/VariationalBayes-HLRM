elbo_chlrm <- function(chlrm_vi_object) {
  # -----------------------------
  # Sufficient statistics and dimensions.
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
  # Prior parameters.
  # -----------------------------
  Mu0       <- chlrm_vi_object$hyperparams$Mu0
  Lambda0   <- chlrm_vi_object$hyperparams$Lambda0
  Lambda0inv <- if (!is.null(chlrm_vi_object$stats$Lambda0inv)) {
    chlrm_vi_object$stats$Lambda0inv
  } else solve_spd(Lambda0, name = "CHLRM Lambda0")
  n0        <- chlrm_vi_object$hyperparams$n0
  S0        <- chlrm_vi_object$hyperparams$S0
  nu0       <- chlrm_vi_object$hyperparams$nu0
  a0        <- chlrm_vi_object$hyperparams$a0
  b0        <- chlrm_vi_object$hyperparams$b0
  alpha0    <- chlrm_vi_object$hyperparams$alpha0   # Dirichlet prior (K)

  # -----------------------------
  # Variational parameters.
  # -----------------------------
  # Mezcla
  rho    <- chlrm_vi_object$params$rho       # m x K
  alpha  <- chlrm_vi_object$params$alpha     # K (Dirichlet)

  # Component coefficients and variances.
  MuBetak     <- chlrm_vi_object$params$MuBetak      # K x p
  SigmaBetak  <- chlrm_vi_object$params$SigmaBetak   # list K de p×p
  ak          <- chlrm_vi_object$params$ak           # K (IG shape)
  bk          <- chlrm_vi_object$params$bk           # K (IG scale)

  # Global beta
  MuBeta     <- chlrm_vi_object$params$MuBeta
  SigmaBeta  <- chlrm_vi_object$params$SigmaBeta

  # Global Sigma uses the IW parameters (nuS, SS), with
  # E[Sigma^{-1}] = nuS * solve_spd(SS)
  nuS        <- chlrm_vi_object$params$nuS
  SS         <- chlrm_vi_object$params$SS

  # Global xi^2 ~ Gamma(axi2, bxi2) (shape–rate)
  axi2       <- chlrm_vi_object$params$axi2
  bxi2       <- chlrm_vi_object$params$bxi2

  # -----------------------------
  # Required expectations.
  # -----------------------------
  Elog_sigma2_k <- expected_log_inverse_gamma(ak, bk)
  Einv_sigma2_k <- expected_inverse_inverse_gamma(ak, bk)
  ESigmainv     <- nuS * solve_spd(SS, name = "CHLRM variational IW scale")
  ElogSigma     <- expected_logdet_inverse_wishart(nuS, SS)
  Exi2          <- axi2 / bxi2
  Elogxi2       <- digamma(axi2) - log(bxi2)
  Elog_omega    <- expected_log_dirichlet(alpha)
  logB_alpha0   <- lgamma(sum(alpha0)) - sum(lgamma(alpha0))  # const Dir prior
  logB_alpha    <- lgamma(sum(alpha))  - sum(lgamma(alpha))   # const Dir var

  # -----------------------------
  # 1) Verosimilitud: E_q[ log p(y | z, beta_k, sigma_k^2) ]
  #     = sum_j sum_k rho_jk [ -1/2 * n_j (log 2π + E log σ_k^2) - 1/2 * E σ_k^{-2} * SSE_{k,j} ]
  # -----------------------------
  Elogpy <- 0
  all_groups <- seq_len(m)
  for (k in seq_len(K)) {
    sse_k <- .expected_sse_batch_chlrm(chlrm_vi_object$stats, all_groups,
                              MuBetak[k, ], SigmaBetak[[k]])
    contribution <- -0.5 * nj * (.log_2pi + Elog_sigma2_k[k]) -
      0.5 * Einv_sigma2_k[k] * sse_k
    Elogpy <- Elogpy + sum(rho[, k] * contribution)
  }

  # -----------------------------
  # 2) Hierarchical layer: E_q[log p(beta_k | Beta, Sigma)].
  #     = sum_k [ -1/2 ( p log 2π + E log|Σ| ) - 1/2 tr( E[Σ^{-1}] * M_k ) ]
  #     donde M_k = Var_q(beta_k) + Var_q(Beta) + (E beta_k - E Beta)(...)
  # -----------------------------
  ElogpBetak <- 0
  for (k in seq_len(K)) {
    mu_k  <- MuBetak[k, ]
    Sig_k <- SigmaBetak[[k]]
    d_k   <- mu_k - MuBeta
    Mk    <- Sig_k + SigmaBeta + tcrossprod(d_k)
    ElogpBetak <- ElogpBetak + ( -0.5 * (p * .log_2pi + ElogSigma) - 0.5 * sum(diag(ESigmainv %*% Mk)) )
  }

  # -----------------------------
  # 3) Prior global Beta: E_q[ log p(Beta) ]
  #     = -1/2 ( p log 2π + log|Λ0| + tr(Λ0^{-1}(Σ_beta + (μ_beta-μ0)(...) )) )
  # -----------------------------
  dB <- MuBeta - Mu0
  ElogpBeta <- -0.5 * ( p * .log_2pi + logdet_spd(Lambda0) +
                          sum(diag(Lambda0inv %*% (SigmaBeta + tcrossprod(dB)))) )

  # -----------------------------
  # 4) Prior global Sigma ~ IW(n0, S0^{-1})
  #     log C_IW + ...
  # -----------------------------
  logC_IW <- 0.5 * n0 * logdet_spd(S0) - 0.5 * n0 * p * log(2) - log_multigamma(0.5 * n0, p)
  ElogpSigma <- logC_IW - 0.5 * (n0 + p + 1) * ElogSigma - 0.5 * sum(diag(S0 %*% ESigmainv))

  # -----------------------------
  # 5) Prior for sigma_k^2 | xi^2, without a separate sigma0^2.
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
  # 8) Entropy terms E_q[log q(.)].
  # -----------------------------
  # q(beta_k) normales
  ElogqBetak <- 0
  for (k in seq_len(K)) {
    ElogqBetak <- ElogqBetak + ( -0.5 * ( p * (1 + .log_2pi) + logdet_spd(SigmaBetak[[k]]) ) )
  }
  # q(Beta) normal
  ElogqBeta <- -0.5 * ( p * (1 + .log_2pi) + logdet_spd(SigmaBeta) )

  # q(Sigma) is inverse-Wishart with parameters (nuS, SS).
  logC_qIW  <- 0.5 * nuS * logdet_spd(SS) - 0.5 * nuS * p * log(2) - log_multigamma(0.5 * nuS, p)
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

