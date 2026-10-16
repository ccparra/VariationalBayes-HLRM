initialize_lrm_mcmc <- function(lrm_mcmc_object) {
  lrm_mcmc_object$params$Beta <- drop(mvtnorm::rmvnorm(
    n = 1,
    mean = lrm_mcmc_object$hyperparams$Beta0,
    sigma = lrm_mcmc_object$hyperparams$Sigma0))

  lrm_mcmc_object$params$sigma2 <- 1 / stats::rgamma(
    n = 1,
    shape = 0.5*lrm_mcmc_object$hyperparams$nu0,
    rate  = 0.5*lrm_mcmc_object$hyperparams$nu0*lrm_mcmc_object$hyperparams$sigma20
  )
  return(lrm_mcmc_object)
}

logprior_lrm <- function(lrm_mcmc_object) {
    lp_beta <- mvtnorm::dmvnorm(
    x     = lrm_mcmc_object$params$Beta,
    mean  = lrm_mcmc_object$hyperparams$Beta0,
    sigma = lrm_mcmc_object$hyperparams$Sigma0,
    log   = TRUE
  )

  lp_s2 <- stats::dgamma(
    1 / lrm_mcmc_object$params$sigma2,
    shape = 0.5 * lrm_mcmc_object$hyperparams$nu0,
    rate = 0.5 * lrm_mcmc_object$hyperparams$nu0 * lrm_mcmc_object$hyperparams$sigma20,
    log = TRUE)
  + log(1 / (lrm_mcmc_object$params$sigma2^2))

  return(lp_beta + lp_s2)
}

sample_coefficients_lrm <- function(lrm_mcmc_object){
  Xty <- lrm_mcmc_object$stats$Xty
  XtX <- lrm_mcmc_object$stats$XtX
  sigma2 <- lrm_mcmc_object$params$sigma2
  Sigma0inv <- lrm_mcmc_object$hyperparams$Sigma0inv
  Beta0 <- lrm_mcmc_object$hyperparams$Beta0
  SigmaBeta <- solve_spd((XtX / sigma2) + Sigma0inv,
                         name = "LRM beta posterior precision")
  MuBeta <- SigmaBeta%*%((Xty/sigma2)+(Sigma0inv%*%Beta0))

  lrm_mcmc_object$params$Beta <- drop(mvtnorm::rmvnorm(n=1,
                                                       mean= MuBeta,
                                                       sigma= SigmaBeta))
  return(lrm_mcmc_object)
}

sample_residual_variance_lrm <- function(lrm_mcmc_object){
  y <- lrm_mcmc_object$stats$y
  X <- lrm_mcmc_object$stats$X
  n <- lrm_mcmc_object$stats$n
  nu0 <- lrm_mcmc_object$hyperparams$nu0
  sigma20 <- lrm_mcmc_object$hyperparams$sigma20
  Beta <- lrm_mcmc_object$params$Beta
  ssr <- drop(crossprod((y-X%*%Beta)))
  asigma2   <- 0.5*(nu0+n)
  bsigma2   <- 0.5*(nu0*sigma20 + ssr)

  lrm_mcmc_object$params$sigma2 <- 1/rgamma(n = 1,
                                            shape = asigma2,
                                            rate = bsigma2)
  return(lrm_mcmc_object)
}
