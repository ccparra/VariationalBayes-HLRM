initialize_hlrm_mcmc <- function(hlrm_mcmc_object) {
  Mu0 <- hlrm_mcmc_object$hyperparams$Mu0
  Lambda0 <- hlrm_mcmc_object$hyperparams$Lambda0
  S0 <- hlrm_mcmc_object$hyperparams$S0
  sigma20 <- hlrm_mcmc_object$hyperparams$sigma20
  nu0 <- hlrm_mcmc_object$hyperparams$nu0
  n0 <- hlrm_mcmc_object$hyperparams$n0
  a0 <- hlrm_mcmc_object$hyperparams$a0
  b0 <- hlrm_mcmc_object$hyperparams$b0

  p <- hlrm_mcmc_object$stats$p
  m <- hlrm_mcmc_object$stats$m

  hlrm_mcmc_object$params$Betaj <- mvtnorm::rmvnorm(
    n=m,
    mean=Mu0,
    sigma=Lambda0)

  hlrm_mcmc_object$params$xi2 <- stats::rgamma(1,shape=a0,rate = b0)

  hlrm_mcmc_object$params$sigma2j <- 1/stats::rgamma(
    n = m,
    shape = 0.5*nu0 ,
    rate = 0.5*nu0*hlrm_mcmc_object$params$xi2)

  hlrm_mcmc_object$params$Beta <- drop(mvtnorm::rmvnorm(
    n=1,
    mean= Mu0,
    sigma=Lambda0))

  hlrm_mcmc_object$params$Sigma <- as.matrix(LaplacesDemon::rinvwishart(
    nu = n0,
    S = S0))

  return(hlrm_mcmc_object)
}

logprior_hlrm <- function(hlrm_mcmc_object) {
  Betaj <- hlrm_mcmc_object$params$Betaj
  Beta  <- hlrm_mcmc_object$params$Beta
  Sigma <- hlrm_mcmc_object$params$Sigma
  xi2   <- hlrm_mcmc_object$params$xi2

  lp_bj <- sum(apply(Betaj, 1L, function(Bj) {
    mvtnorm::dmvnorm(x = Bj, mean = Beta, sigma = Sigma, log = TRUE)
  }))

  sigma2j<- hlrm_mcmc_object$params$sigma2j
  lp_s2j <- sum(
    stats::dgamma(
      x = 1 / sigma2j,
      shape = 0.5 * hlrm_mcmc_object$hyperparams$nu0,
      rate  = 0.5 * hlrm_mcmc_object$hyperparams$nu0 * xi2,
      log = TRUE
    ) + log(1 / (sigma2j^2))
  )

  lp_b <- drop(mvtnorm::dmvnorm(
    x     = Beta,
    mean  = hlrm_mcmc_object$hyperparams$Mu0,
    sigma = hlrm_mcmc_object$hyperparams$Lambda0,
    log   = TRUE
  ))

  lp_S <- LaplacesDemon::dinvwishart(
    Sigma = Sigma,
    nu = hlrm_mcmc_object$hyperparams$n0,
    S = hlrm_mcmc_object$hyperparams$S0,
    log = TRUE
    )

  lp_xi2 <- stats::dgamma(
    x = xi2,
    shape = hlrm_mcmc_object$hyperparams$a0,
    rate = hlrm_mcmc_object$hyperparams$b0,
    log = TRUE)

  return(as.numeric(lp_bj + lp_s2j + lp_b + lp_S + lp_xi2))
}

sample_group_coefficients_hlrm <- function(hlrm_mcmc_object){
  m <- hlrm_mcmc_object$stats$m
  p <- hlrm_mcmc_object$stats$p
  X <- hlrm_mcmc_object$stats$X
  XtX <- hlrm_mcmc_object$stats$XtX
  Xty <- hlrm_mcmc_object$stats$Xty
  sigma2j <- hlrm_mcmc_object$params$sigma2j
  Sigma <- hlrm_mcmc_object$params$Sigma
  Beta <- hlrm_mcmc_object$params$Beta
  Betaj <- hlrm_mcmc_object$params$Betaj

  SigInv  <- solve_spd(Sigma, name = "HLRM Sigma")

  for (j in seq_len(m)) {
    XtXj <- as.matrix(XtX[[j]])
    Xtyj <- as.matrix(Xty[[j]])

    Vj  <- solve_spd((1 / sigma2j[j]) * XtXj + SigInv,
                     name = "HLRM group beta posterior precision")
    muj <- Vj %*% (((1/sigma2j[j])*Xtyj) + (SigInv%*%Beta))

    Betaj[j, ] <- drop(mvtnorm::rmvnorm(n = 1, mean = drop(muj), sigma = Vj))
  }

  hlrm_mcmc_object$params$Betaj <- Betaj
  return(hlrm_mcmc_object)
}

sample_group_variances_hlrm <- function(hlrm_mcmc_object){
  X <- hlrm_mcmc_object$stats$X
  y <- hlrm_mcmc_object$stats$y
  m <- hlrm_mcmc_object$stats$m
  nj <- hlrm_mcmc_object$stats$nj
  nu0 <- hlrm_mcmc_object$hyperparams$nu0
  Betaj <- hlrm_mcmc_object$params$Betaj
  xi2 <- hlrm_mcmc_object$params$xi2
  sigma2j <- hlrm_mcmc_object$params$sigma2j

  for (j in seq_len(m)) {
      aj <- 0.5*(nu0+nj[j])
      ssrj <- drop(crossprod((y[[j]]-X[[j]]%*%Betaj[j,])))
      bj <- 0.5*(nu0*xi2+ssrj)
      sigma2j[j]<-1/stats::rgamma(n=1,shape=aj,rate=bj)
    }
  hlrm_mcmc_object$params$sigma2j <- drop(sigma2j)
  return(hlrm_mcmc_object)
}

sample_xi2_hlrm <- function(hlrm_mcmc_object){
  m <- hlrm_mcmc_object$stats$m
  a0 <- hlrm_mcmc_object$hyperparams$a0
  b0 <- hlrm_mcmc_object$hyperparams$b0
  nu0 <- hlrm_mcmc_object$hyperparams$nu0
  sigma2j <- hlrm_mcmc_object$params$sigma2j

  axi2 <- 0.5*m*nu0 + a0
  bxi2 <- b0 + 0.5*nu0*sum(1 / sigma2j)
  hlrm_mcmc_object$params$xi2 <- stats::rgamma(n=1,shape=axi2,rate=bxi2)
  return(hlrm_mcmc_object)
}

sample_global_coefficients_hlrm <- function(hlrm_mcmc_object){
  m <- hlrm_mcmc_object$stats$m
  Mu0 <- hlrm_mcmc_object$hyperparams$Mu0
  Lambda0 <- hlrm_mcmc_object$hyperparams$Lambda0
  Sigma <- hlrm_mcmc_object$params$Sigma
  Betaj <- hlrm_mcmc_object$params$Betaj

  Lambda0Inv <- hlrm_mcmc_object$hyperparams$Lambda0inv
  Sigmainv <- solve_spd(Sigma, name = "HLRM Sigma")
  VB <- solve_spd(Lambda0Inv + m * Sigmainv,
                  name = "HLRM global beta posterior precision")
  muB <- VB%*%(Lambda0Inv%*%Mu0+Sigmainv%*%colSums(Betaj))
  hlrm_mcmc_object$params$Beta <- drop(mvtnorm::rmvnorm(n=1,mean=muB,sigma = VB))

  return(hlrm_mcmc_object)
}

sample_global_covariance_hlrm <- function(hlrm_mcmc_object){
  m <- hlrm_mcmc_object$stats$m
  p <- hlrm_mcmc_object$stats$p
  S0 <- hlrm_mcmc_object$hyperparams$S0
  n0 <- hlrm_mcmc_object$hyperparams$n0
  Betaj <- hlrm_mcmc_object$params$Betaj
  Beta <- hlrm_mcmc_object$params$Beta

  diff <- Betaj-matrix(Beta, nrow = m, ncol = p, byrow = TRUE)
  sS <- S0+crossprod(diff)
  nuS <- n0+m
  Sigma <- LaplacesDemon::rinvwishart(nu = nuS, S=sS)
  hlrm_mcmc_object$params$Sigma <- symmetrize_numeric(Sigma, "HLRM Sigma draw")
  return(hlrm_mcmc_object)
}
