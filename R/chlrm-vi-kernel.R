initialize_chlrm_vi <- function(chlrm_vi_object, initialization = NULL) {
  # Dimensions.
  p <- chlrm_vi_object$stats$p
  m <- chlrm_vi_object$stats$m
  K <- chlrm_vi_object$stats$K

  # Prior parameters.
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
  # Initialization is an algorithm setting, not a statistical prior.
  # ---------------------------
  if (is.null(initialization)) {
    initialization <- if (is.null(B_ols)) "random_responsibilities" else "kmeans_ols"
  }
  if (identical(initialization, "random_responsibilities")) {
    rho <- gtools::rdirichlet(m, rep(1, K))
    centers <- mvtnorm::rmvnorm(K, mean = Mu0, sigma = Lambda0)
  } else if (identical(initialization, "kmeans_ols")) {
    if (is.null(B_ols)) stop("kmeans_ols initialization requires Betaj_ols")
    cl <- stats::kmeans(B_ols, centers = K)
    rho <- matrix(0, m, K); rho[cbind(seq_len(m), cl$cluster)] <- 1
    centers <- cl$centers
  } else {
    stop("initialization must be random_responsibilities or kmeans_ols")
  }

  # Initial responsibilities are separate from the data-independent prior.
  chlrm_vi_object$params$rho   <- rho
  chlrm_vi_object$params$alpha <- alpha0 + colSums(rho)

  # Initial component means, distinct from their prior hyperparameters.
  chlrm_vi_object$params$MuBetak <- centers

  # ---------------------------
  # Inicializar SigmaBetak desde la prior IW de Sigma
  # Prior: Sigma ~ IW(n0, S0), consistently with MCMC and the ELBO.
  # ---------------------------
  Sigma0_draw <- LaplacesDemon::rinvwishart(nu = n0, S = S0)
  Sigma0_draw <- symmetrize_numeric(Sigma0_draw, "CHLRM initial Sigma draw")

  chlrm_vi_object$params$SigmaBetak <- replicate(K, Sigma0_draw, simplify = FALSE)

  # Beta global ~ N(Mu0, Lambda0)
  Sigma0_safe <- Lambda0 + 1e-8 * diag(p)
  mu_draw <- as.numeric(mvtnorm::rmvnorm(
    n = 1, mean = Mu0,
    sigma = symmetrize_numeric(Sigma0_safe, "CHLRM initial beta covariance")))
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
  Einvs2_sum <- sum(expected_inverse_inverse_gamma(
    chlrm_vi_object$params$ak, chlrm_vi_object$params$bk))
  chlrm_vi_object$params$bxi2 <- b0 + 0.5 * nu0 * Einvs2_sum

  return(chlrm_vi_object)
}
.expected_sse_batch_chlrm <- function(stats, J, mu_k, Sig_k) {
  Xty <- stats$Xty_mat[J, , drop = FALSE]
  XtX <- stats$XtX_mat[, J, drop = FALSE]
  linear <- as.numeric(Xty %*% mu_k)
  second_moment <- as.vector(Sig_k + tcrossprod(mu_k))
  quadratic <- as.numeric(crossprod(second_moment, XtX))
  pmax(stats$yty[J] - 2 * linear + quadratic, 0)
}

# ================================================================
# 1) q(gamma_j) = Categorical(rho_j•)
#    log rho_{jk} ∝ E[log ω_k] - ½ n_j E[log σ_k^2] - ½ E[σ_k^{-2}] SSE_{k,j}
# ================================================================
update_q_allocations_chlrm <- function(chlrm_vi_object, batch_idx = NULL) {
  m  <- chlrm_vi_object$stats$m
  K  <- chlrm_vi_object$stats$K
  J  <- if (is.null(batch_idx)) seq_len(m) else batch_idx

  Elog_omega <- expected_log_dirichlet(chlrm_vi_object$params$alpha)
  Einv_s2_k  <- expected_inverse_inverse_gamma(
    chlrm_vi_object$params$ak, chlrm_vi_object$params$bk)
  Elog_s2_k  <- expected_log_inverse_gamma(
    chlrm_vi_object$params$ak, chlrm_vi_object$params$bk)

  sse <- matrix(vapply(seq_len(K), function(k) {
    .expected_sse_batch_chlrm(chlrm_vi_object$stats, J,
                     chlrm_vi_object$params$MuBetak[k, ],
                     chlrm_vi_object$params$SigmaBetak[[k]])
  }, numeric(length(J))), nrow = length(J), ncol = K)
  log_rho <- -0.5 * (outer(chlrm_vi_object$stats$nj[J], Elog_s2_k) +
                       sweep(sse, 2L, Einv_s2_k, `*`))
  log_rho <- sweep(log_rho, 2L, Elog_omega, `+`)
  log_rho <- log_rho - apply(log_rho, 1L, max)
  rho_J <- exp(log_rho)
  chlrm_vi_object$params$rho[J, ] <- rho_J / rowSums(rho_J)
  return(chlrm_vi_object)
}

# ================================================================
# 2) q(beta_k) = N(mu_k, Sigma_k); CAVI equivale a learning_rate = 1.
#    V_k = (E[Σ^{-1}] + E[σ_k^{-2}] * Σ_j ρ_{jk} Xj'Xj)^{-1}
#    m_k = V_k (E[Σ^{-1}] MuBeta + E[σ_k^{-2}] * Σ_j ρ_{jk} Xj'yj)
# ================================================================
update_q_component_coefficients_chlrm <- function(chlrm_vi_object, batch_idx = NULL,
                                                  learning_rate = 1) {
  m  <- chlrm_vi_object$stats$m
  p  <- chlrm_vi_object$stats$p
  K  <- chlrm_vi_object$stats$K
  J  <- if (is.null(batch_idx)) seq_len(m) else batch_idx
  w  <- m / length(J)


  E_SigInv <- chlrm_vi_object$params$nuS *
    solve_spd(chlrm_vi_object$params$SS, name = "CHLRM variational IW scale")
  EinvS    <- expected_inverse_inverse_gamma(
    chlrm_vi_object$params$ak, chlrm_vi_object$params$bk)

  for (k in seq_len(K)) {
    weights <- w * chlrm_vi_object$params$rho[J, k]
    Ak <- matrix(chlrm_vi_object$stats$XtX_mat[, J, drop = FALSE] %*% weights,
                 nrow = p, ncol = p)
    bvec <- as.numeric(crossprod(weights,
                                chlrm_vi_object$stats$Xty_mat[J, , drop = FALSE]))
    P_hat <- symmetrize_numeric(E_SigInv + EinvS[k] * Ak,
                                "CHLRM component beta precision")
    V_hat <- solve_spd(P_hat, name = "CHLRM component beta precision")
    m_hat <- as.numeric(V_hat %*% (E_SigInv %*% chlrm_vi_object$params$MuBeta + EinvS[k] * bvec))

    # SVI must average Gaussian natural parameters (precision and
    # precision-weighted mean), not means and covariance matrices directly.
    Sig_old <- chlrm_vi_object$params$SigmaBetak[[k]]
    mu_old  <- chlrm_vi_object$params$MuBetak[k, ]
    P_old <- solve_spd(Sig_old, name = "CHLRM component beta covariance")
    h_old <- as.numeric(P_old %*% mu_old)
    h_hat <- as.numeric(P_hat %*% m_hat)
    P_new <- symmetrize_numeric((1 - learning_rate) * P_old +
      learning_rate * P_hat,
                                "CHLRM mixed component precision")
    h_new <- (1 - learning_rate) * h_old + learning_rate * h_hat
    Sig_new <- solve_spd(P_new, name = "CHLRM mixed component precision")
    chlrm_vi_object$params$SigmaBetak[[k]] <- Sig_new
    chlrm_vi_object$params$MuBetak[k, ]    <- as.numeric(Sig_new %*% h_new)
  }
  chlrm_vi_object
}


# ================================================================
# 3) q(σ_k^2) = IG(a_k,b_k); CAVI equivale a learning_rate = 1.
#    a_hat = (ν0 + Σ_j ρ_{jk} n_j)/2
#    b_hat = 1/2 * [ ν0 * E[ξ^2] + Σ_j ρ_{jk} SSE_{k,j} ]
# ================================================================
update_q_component_variances_chlrm <- function(chlrm_vi_object, batch_idx = NULL,
                                               learning_rate = 1) {
  m   <- chlrm_vi_object$stats$m
  K   <- chlrm_vi_object$stats$K
  J   <- if (is.null(batch_idx)) seq_len(m) else batch_idx
  w   <- m / length(J)

  nu0 <- chlrm_vi_object$hyperparams$nu0
  Exi2 <- chlrm_vi_object$params$axi2 / chlrm_vi_object$params$bxi2

  for (k in seq_len(K)) {
    mu_k   <- chlrm_vi_object$params$MuBetak[k, ]
    Sig_k  <- chlrm_vi_object$params$SigmaBetak[[k]]
    weights <- w * chlrm_vi_object$params$rho[J, k]
    sse_kj <- .expected_sse_batch_chlrm(chlrm_vi_object$stats, J, mu_k, Sig_k)
    Nobs_k <- sum(weights * chlrm_vi_object$stats$nj[J])
    SSE_k  <- sum(weights * sse_kj)
    a_hat <- 0.5 * (nu0 + Nobs_k)
    b_hat <- 0.5 * (nu0 * Exi2 + SSE_k)

    chlrm_vi_object$params$ak[k] <- (1 - learning_rate) *
      chlrm_vi_object$params$ak[k] + learning_rate * a_hat
    chlrm_vi_object$params$bk[k] <- (1 - learning_rate) *
      chlrm_vi_object$params$bk[k] + learning_rate * b_hat
  }
  return(chlrm_vi_object)
}

# ================================================================
# 4) q(ω) = Dir(alpha); CAVI equivale a learning_rate = 1.
#    alpha_hat = alpha0 + Σ_j ρ_{jk}  (minibatch escalado)
# ================================================================
update_q_mixture_weights_chlrm <- function(chlrm_vi_object, batch_idx = NULL,
                                           learning_rate = 1) {
  m  <- chlrm_vi_object$stats$m
  K  <- chlrm_vi_object$stats$K
  J  <- if (is.null(batch_idx)) seq_len(m) else batch_idx
  w  <- m / length(J)

  alpha0 <- chlrm_vi_object$hyperparams$alpha0
  Nk_hat <- w * colSums(chlrm_vi_object$params$rho[J, , drop = FALSE])
  alpha_hat <- alpha0 + Nk_hat

  chlrm_vi_object$params$alpha <- (1 - learning_rate) *
    chlrm_vi_object$params$alpha + learning_rate * alpha_hat
  return(chlrm_vi_object)
}

# ================================================================
# 5) Global q(beta) = N(MuBeta, SigmaBeta) under SVI.
#    J_hat = Λ0^{-1} + K * E[Σ^{-1}]
#    h_hat = Λ0^{-1} μ0 + E[Σ^{-1}] * Σ_k μ_k
# ================================================================
update_q_global_coefficients_chlrm <- function(chlrm_vi_object,
                                               learning_rate = 1) {
  p  <- chlrm_vi_object$stats$p
  K  <- chlrm_vi_object$stats$K

  mu0  <- chlrm_vi_object$hyperparams$Mu0
  Lam0 <- chlrm_vi_object$hyperparams$Lambda0
  J0 <- if (!is.null(chlrm_vi_object$stats$Lambda0inv)) {
    chlrm_vi_object$stats$Lambda0inv
  } else solve_spd(Lam0, name = "CHLRM Lambda0")

  E_SigInv <- chlrm_vi_object$params$nuS *
    solve_spd(chlrm_vi_object$params$SS, name = "CHLRM variational IW scale")

  J_hat <- J0 + K * E_SigInv
  h_hat <- as.numeric(J0 %*% mu0 + E_SigInv %*% colSums(chlrm_vi_object$params$MuBetak))

  Sig_old <- chlrm_vi_object$params$SigmaBeta
  mu_old  <- chlrm_vi_object$params$MuBeta

  # Natural-param mixing
  J_old <- solve_spd(Sig_old, name = "CHLRM variational global beta covariance")
  h_old <- as.numeric(J_old %*% mu_old)

  J_new <- (1 - learning_rate) * J_old + learning_rate * J_hat
  h_new <- (1 - learning_rate) * h_old + learning_rate * h_hat

  Sig_new <- solve_spd(J_new, name = "CHLRM mixed global beta precision")
  mu_new  <- as.numeric(Sig_new %*% h_new)

  chlrm_vi_object$params$SigmaBeta <- Sig_new
  chlrm_vi_object$params$MuBeta    <- mu_new
  return(chlrm_vi_object)
}

# ================================================================
# 6) q(Sigma) = IW(nuS, SS) under SVI. Store (nuS, SS) so that
#    E[Σ^{-1}] = nuS * solve_spd(SS))
#    nu_hat = n0 + K
#    SS_hat = S0 + Σ_k [ Σ_beta_k + Σ_beta + (μ_k - μ_beta)(μ_k - μ_beta)' ]
# ================================================================
update_q_global_covariance_chlrm <- function(chlrm_vi_object,
                                             learning_rate = 1) {
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

  chlrm_vi_object$params$nuS <- (1 - learning_rate) *
    chlrm_vi_object$params$nuS + learning_rate * nu_hat
  chlrm_vi_object$params$SS  <- symmetrize_numeric(
    (1 - learning_rate) * chlrm_vi_object$params$SS +
      learning_rate * SS_hat,
    "CHLRM mixed IW scale")

  return(chlrm_vi_object)
}

# ================================================================
# 7) q(xi^2) = Gamma(a,b), using the shape--rate convention under SVI.
#    a_hat = a0 + (K * ν0)/2
#    b_hat = b0 + (ν0/2) * Σ_k E[σ_k^{-2}]
# ================================================================
update_q_xi2_chlrm <- function(chlrm_vi_object, learning_rate = 1) {
  K   <- chlrm_vi_object$stats$K
  nu0 <- chlrm_vi_object$hyperparams$nu0
  a0  <- chlrm_vi_object$hyperparams$a0
  b0  <- chlrm_vi_object$hyperparams$b0

  EinvS_k <- expected_inverse_inverse_gamma(
    chlrm_vi_object$params$ak, chlrm_vi_object$params$bk)

  a_hat <- a0 + 0.5 * K * nu0
  b_hat <- b0 + 0.5 * nu0 * sum(EinvS_k)

  chlrm_vi_object$params$axi2 <- (1 - learning_rate) *
    chlrm_vi_object$params$axi2 + learning_rate * a_hat
  chlrm_vi_object$params$bxi2 <- (1 - learning_rate) *
    chlrm_vi_object$params$bxi2 + learning_rate * b_hat
  return(chlrm_vi_object)
}
