initialize_chlrm_mcmc <- function(chlrm_mcmc_object) {
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
  Beta  <- drop(mvtnorm::rmvnorm(n = 1, mean = Mu0, sigma = Lambda0))
  Sigma <- LaplacesDemon::rinvwishart(nu = n0, S = S0)
  xi2   <- stats::rgamma(1, shape = a0, rate = b0)

  Betak  <- mvtnorm::rmvnorm(n = K, mean = Beta, sigma = Sigma)
  sigma2k <- 1 / stats::rgamma(n = K, shape = 0.5 * nu0, rate = 0.5 * nu0 * xi2)

  chlrm_mcmc_object$params$omega   <- omega
  chlrm_mcmc_object$params$z       <- z
  chlrm_mcmc_object$params$Beta    <- Beta
  chlrm_mcmc_object$params$Sigma   <- Sigma
  chlrm_mcmc_object$params$xi2     <- xi2
  chlrm_mcmc_object$params$Betak   <- Betak
  chlrm_mcmc_object$params$sigma2k <- sigma2k

  return(chlrm_mcmc_object)
}

sample_allocations_chlrm <- function(chlrm_mcmc_object) {
  K    <- chlrm_mcmc_object$stats$K
  m    <- chlrm_mcmc_object$stats$m
  group_sizes <- chlrm_mcmc_object$stats$nj
  XtX  <- chlrm_mcmc_object$stats$XtX    # lista p×p
  Xty  <- chlrm_mcmc_object$stats$Xty    # lista p×1
  yty  <- chlrm_mcmc_object$stats$yty


  Betak    <- chlrm_mcmc_object$params$Betak   # K×p
  omega    <- chlrm_mcmc_object$params$omega   # K
  sigma2k  <- chlrm_mcmc_object$params$sigma2k # K

  log_2pi      <- log(2 * pi)
  log_sigma2k  <- log(sigma2k)
  inv_sigma2k  <- 1 / sigma2k

  z  <- integer(m)
  B  <- t(Betak)  # p x K for vectorized multiplication

  for (j in seq_len(m)) {
    XtXj <- XtX[[j]]      # p×p
    Xtyj <- Xty[[j]]      # p×1
    ytyj <- yty[j]

    # l_vec_k = β_k^T X_j^T y_j  (K×1)
    l_vec <- as.numeric(Betak %*% Xtyj)

    # Compute beta_k' X_j' X_j beta_k for every component at once.
    # M = XtXj %*% B (p×K); q_k = colSums(B * M)
    M     <- XtXj %*% B
    q_vec <- colSums(B * M)

    SSR_vec <- ytyj - 2 * l_vec + q_vec

    llk <- -0.5 * group_sizes[j] * (log_2pi + log_sigma2k) -
      0.5 * SSR_vec * inv_sigma2k
    log_post <- log(omega) + llk

    # Numerical stabilization.
    log_post <- log_post - max(log_post)
    pk <- exp(log_post); pk <- pk / sum(pk)

    z[j] <- sample.int(K, size = 1, prob = pk)
  }

  chlrm_mcmc_object$params$z       <- z
  return(chlrm_mcmc_object)
}

sample_omega_chlrm <- function(chlrm_mcmc_object) {

  alpha0 <- chlrm_mcmc_object$hyperparams$alpha0
  nk <- tabulate(chlrm_mcmc_object$params$z,
                 nbins = chlrm_mcmc_object$stats$K)

  post_alpha <- alpha0 + nk
  chlrm_mcmc_object$params$omega <- drop(gtools::rdirichlet(1, post_alpha))
  return(chlrm_mcmc_object)
}

sample_component_coefficients_chlrm <- function(chlrm_mcmc_object) {
  K       <- chlrm_mcmc_object$stats$K
  p       <- chlrm_mcmc_object$stats$p
  XtX     <- chlrm_mcmc_object$stats$XtX   # lista de matrices p×p
  Xty     <- chlrm_mcmc_object$stats$Xty   # lista de vectores/matrices p×1
  z       <- chlrm_mcmc_object$params$z    # asignaciones de grupo
  Sigma   <- chlrm_mcmc_object$params$Sigma
  SigInv  <- solve_spd(Sigma, name = "CHLRM Sigma")
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
    Vk     <- solve_spd(Prec_k, name = "CHLRM component beta posterior precision")
    muk    <- Vk %*% ((1 / sigma2k[k]) * Xtyk + SigInv %*% bcol)

    Betak[k, ] <- drop(mvtnorm::rmvnorm(1, mean = drop(muk), sigma = Vk))
  }

  chlrm_mcmc_object$params$Betak <- Betak

  return(chlrm_mcmc_object)
}

sample_component_variances_chlrm <- function(chlrm_mcmc_object) {
  # Stats
  group_sizes <- chlrm_mcmc_object$stats$nj
  K   <- chlrm_mcmc_object$stats$K
  XtX <- chlrm_mcmc_object$stats$XtX    # lista p x p
  Xty <- chlrm_mcmc_object$stats$Xty    # lista p x 1
  yty <- chlrm_mcmc_object$stats$yty    #1x p

  # Current parameters.
  z     <- chlrm_mcmc_object$params$z
  Betak <- chlrm_mcmc_object$params$Betak
  xi2   <- chlrm_mcmc_object$params$xi2
  nu0   <- chlrm_mcmc_object$hyperparams$nu0

  sigma2k <- numeric(K)

  for (k in seq_len(K)) {
    idx <- which(z == k)
    if (length(idx) > 0) {
      n_k_obs <- sum(group_sizes[idx])
      ssr_k <- 0
      for (j in idx) {
        bj   <- Betak[k, ]
        xtyj <- Xty[[j]]          # p x 1
        xtxj <- XtX[[j]]          # p x p

        l_j  <- sum(bj * xtyj)                 # β_j^T X_j^T y_j
        q_j  <- sum(bj * (xtxj %*% bj))        # β_j^T X_j^T X_j β_j

        ssr_k <- ssr_k + (yty[j] - 2 * l_j + q_j)
      }

      a_k <- 0.5 * (nu0 + n_k_obs)
      b_k <- 0.5 * (nu0 * xi2 + ssr_k)
    } else {
      # Empty component: draw from the prior.
      a_k <- 0.5 * nu0
      b_k <- 0.5 * nu0 * xi2
    }

    sigma2k[k] <- 1 / stats::rgamma(1, shape = a_k, rate = b_k)
  }

  chlrm_mcmc_object$params$sigma2k <- sigma2k

  return(chlrm_mcmc_object)
}

sample_xi2_chlrm <- function(chlrm_mcmc_object){
  a0   <- chlrm_mcmc_object$hyperparams$a0
  b0   <- chlrm_mcmc_object$hyperparams$b0
  nu0  <- chlrm_mcmc_object$hyperparams$nu0

  sigma2k <- chlrm_mcmc_object$params$sigma2k
  K <- chlrm_mcmc_object$stats$K
  invsum <- sum(1 / sigma2k)

  shape <- a0 + 0.5 * nu0 * K
  rate  <- b0 + 0.5 * nu0 * invsum

  chlrm_mcmc_object$params$xi2 <- stats::rgamma(1, shape = shape, rate = rate)
  return(chlrm_mcmc_object)
}

sample_global_coefficients_chlrm <- function(chlrm_mcmc_object){
  # Prior parameters.
  Mu0     <- chlrm_mcmc_object$hyperparams$Mu0      # μ0 (p)
  Lambda0 <- chlrm_mcmc_object$hyperparams$Lambda0  # Λ0 (p×p)

  # Current parameters.
  Sigma <- chlrm_mcmc_object$params$Sigma           # Σ (p×p)
  Betak <- chlrm_mcmc_object$params$Betak           # matriz K×p
  K <- chlrm_mcmc_object$stats$K

  Lam0Inv <- chlrm_mcmc_object$hyperparams$Lambda0inv
  SigInv  <- solve_spd(Sigma, name = "CHLRM Sigma")

  sumBetak <- colSums(Betak)                         # fixed-K prior includes empty components

  Prec <- Lam0Inv + K * SigInv
  V    <- solve_spd(Prec, name = "CHLRM global beta posterior precision")
  bvec <- matrix(sumBetak, ncol = 1)
  mu   <- V %*% (Lam0Inv %*% Mu0 + SigInv %*% bvec)

  chlrm_mcmc_object$params$Beta <- drop(
    mvtnorm::rmvnorm(n = 1, mean = drop(mu), sigma = V)
  )
  return(chlrm_mcmc_object)
}

sample_global_covariance_chlrm <- function(chlrm_mcmc_object){
  # Prior parameters.
  S0 <- chlrm_mcmc_object$hyperparams$S0   # p×p
  n0 <- chlrm_mcmc_object$hyperparams$n0

  # Current parameters.
  Betak <- chlrm_mcmc_object$params$Betak  # K×p
  Beta  <- chlrm_mcmc_object$params$Beta   # p
  K <- chlrm_mcmc_object$stats$K

  #  sum_k (βk-β)(βk-β)^T
  D  <- sweep(Betak, 2, Beta, "-")                          # K×p
  SS <- S0 + crossprod(D)                                   # p×p
  nuS <- n0 + K

  # Muestreo Inv-Wishart
  Sigma_draw <- LaplacesDemon::rinvwishart(nu = nuS, S = SS)
  Sigma_draw <- symmetrize_numeric(Sigma_draw, "CHLRM Sigma draw")

  chlrm_mcmc_object$params$Sigma <- Sigma_draw
  return(chlrm_mcmc_object)
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
  # Use the closed form to avoid an additional dependency.
  lp_omega <- lgamma(sum(alpha0)) - sum(lgamma(alpha0)) +
    sum((alpha0 - 1) * log(omega))

  # 7) z | ω  (producto de ω_{z_j})
  lp_z <- sum(log(omega[z]))

  as.numeric(lp_betak + lp_beta + lp_Sigma + lp_s2k + lp_xi2 + lp_omega + lp_z)
}
incidence_matrix_mcmc_chlrm <- function(Z) {
  # Z is S x m: rows are draws and columns are groups.
  S <- nrow(Z); n <- ncol(Z)
  A <- matrix(0, n, n)
  for (s in seq_len(S)) {
    labs <- Z[s, ]
    A <- A + (outer(labs, labs, "==") * 1.0)  # co-clustering for this draw
  }
  A <- A / S
  diag(A) <- 1
  A
}

dahl_partition_chlrm <- function(z_mat) {
  # z_mat: m x S. Dahl's least-squares draw is label invariant.
  Z <- t(z_mat)
  psm <- incidence_matrix_mcmc_chlrm(Z)
  losses <- apply(Z, 1L, function(z) {
    sum((outer(z, z, "==") - psm)^2)
  })
  list(draw = which.min(losses), z = Z[which.min(losses), ], psm = psm)
}
