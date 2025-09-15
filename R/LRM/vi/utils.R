source('R/utils.R')

Initialize_lrm_vi <- function(lrm_vi_object) {
  # Hiperparámetros
  Beta0   <- lrm_vi_object$hyperparams$Beta0
  Sigma0  <- lrm_vi_object$hyperparams$Sigma0
  nu0     <- lrm_vi_object$hyperparams$nu0
  sigma20 <- lrm_vi_object$hyperparams$sigma20
  
  # q(Beta): inicializar con un draw desde la prior normal
  beta_draw <- as.numeric(
    mvtnorm::rmvnorm(1, mean = Beta0, sigma = Sigma0 + 1e-8 * diag(length(Beta0)))
  )
  lrm_vi_object$params$MuBeta    <- beta_draw
  lrm_vi_object$params$SigmaBeta <- Sigma0
  
  # q(sigma^2): IG(a, b) con a = nu0/2, b = (nu0 * sigma20)/2 (shape–scale)
  lrm_vi_object$params$asigma2 <- 0.5 * nu0
  lrm_vi_object$params$bsigma2 <- 0.5 * nu0 * sigma20
  
  return(lrm_vi_object)
}


update_qsigma2_lrm <- function(lrm_vi_object){
  SigmaBeta <- lrm_vi_object$params$SigmaBeta 
  MuBeta <- lrm_vi_object$params$MuBeta
  nu0 <- lrm_vi_object$hyperparams$nu0
  sigma20 <- lrm_vi_object$hyperparams$sigma20
  XtX <- lrm_vi_object$stats$XtX
  Xty <- lrm_vi_object$stats$Xty
  X <- lrm_vi_object$stats$X
  y <- lrm_vi_object$stats$y
  ssr <- drop(crossprod(y-X%*%MuBeta))
  n <- lrm_vi_object$stats$n
  asigma2 <- 0.5*(n+nu0)
  bsigma2 <- 0.5*(nu0*sigma20+ssr+sum(diag(XtX%*%SigmaBeta)))
  lrm_vi_object$params$asigma2 = asigma2
  lrm_vi_object$params$bsigma2 = bsigma2
  return(lrm_vi_object)
}

update_qBeta_lrm <- function(lrm_vi_object){
  asigma2 <- lrm_vi_object$params$asigma2
  bsigma2 <- lrm_vi_object$params$bsigma2
  Beta0 <- lrm_vi_object$hyperparams$Beta0
  Sigma0inv <- lrm_vi_object$hyperparams$Sigma0inv
  XtX <- lrm_vi_object$stats$XtX
  Xty <- lrm_vi_object$stats$Xty
  X <- lrm_vi_object$stats$X
  y <- lrm_vi_object$stats$y
  SigmaBeta <- solve((((asigma2/bsigma2)*XtX)+(Sigma0inv)))
  SigmaBeta <- .sym(SigmaBeta)
  MuBeta <- SigmaBeta%*%(((asigma2/bsigma2)*Xty)+(Sigma0inv%*%Beta0))
  lrm_vi_object$params$SigmaBeta = SigmaBeta
  lrm_vi_object$params$MuBeta = MuBeta
  return(lrm_vi_object)
}

ELBO_lrm <- function(lrm_vi_object){
  
  X <- lrm_vi_object$stats$X
  y <- lrm_vi_object$stats$y
  XtX <- lrm_vi_object$stats$XtX
  Xty <- lrm_vi_object$stats$Xty
  n <- lrm_vi_object$stats$n
  p <- lrm_vi_object$stats$p
  
  Beta0 <- lrm_vi_object$hyperparams$Beta0
  Sigma0<- lrm_vi_object$hyperparams$Sigma0
  Sigma0inv <- lrm_vi_object$hyperparams$Sigma0inv
  nu0 <- lrm_vi_object$hyperparams$nu0
  sigma20 <- lrm_vi_object$hyperparams$sigma20
  a0   <- 0.5 * lrm_vi_object$hyperparams$nu0
  b0   <- 0.5 * lrm_vi_object$hyperparams$nu0 * lrm_vi_object$hyperparams$sigma20
  
  SigmaBeta <- lrm_vi_object$params$SigmaBeta 
  MuBeta <- lrm_vi_object$params$MuBeta
  asigma2 <- lrm_vi_object$params$asigma2
  bsigma2 <- lrm_vi_object$params$bsigma2
  
  ssr <- drop(crossprod(y-X%*%MuBeta))
  tr_XSigXt <- sum(diag(XtX %*% SigmaBeta))
  Essr <- ssr+tr_XSigXt
  Elogs2 <- E_log_sigma2(asigma2,bsigma2)
  Einvs2 <- E_inv_sigma2(asigma2,bsigma2)
  Elogpy <- -0.5*n*log(2*pi)-0.5*n*Elogs2-0.5*Einvs2*Essr
  
  diff <- MuBeta-Beta0
  ElogpBeta <- -0.5*p*log(2*pi)-0.5*logdet(Sigma0) -
    0.5*(sum(diag(Sigma0inv %*% SigmaBeta))+drop(crossprod(diff,Sigma0inv%*%diff)))
  
  Elogpsigma2 <- a0*log(b0)-lgamma(a0)-(a0+1)*Elogs2-b0*Einvs2
  
  ElogqBeta <- -0.5*(p*(1+log(2*pi))+logdet(SigmaBeta))
  
  Elogqsigma2 <- asigma2*log(bsigma2)-lgamma(asigma2)-(asigma2+1)*Elogs2-bsigma2*Einvs2
  
  elbo <- Elogpy + ElogpBeta + Elogpsigma2 - (ElogqBeta + Elogqsigma2)
  
  return(as.numeric(elbo))
}

yhat_vi_lrm <- function(X,hat_values) {
  Beta <- hat_values$Beta$mean
  yhat     <- as.vector(X %*% Beta)
  return(yhat)
}

vi_sample_draws_lrm <- function(lrm_vi_object, S, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  MuBeta  <- lrm_vi_object$params$MuBeta
  SigmaBeta <- lrm_vi_object$params$SigmaBeta
  asigma2   <- lrm_vi_object$params$asigma2
  bsigma2   <- lrm_vi_object$params$bsigma2
  
  Beta <- mvtnorm::rmvnorm(n = S, mean = MuBeta, sigma = SigmaBeta)  
  sigma2   <- 1 / stats::rgamma(n = S, shape = asigma2, rate = bsigma2)
  list(Beta = Beta, sigma2 = sigma2)
}