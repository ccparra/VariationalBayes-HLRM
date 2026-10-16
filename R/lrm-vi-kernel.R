initialize_lrm_vi <- function(lrm_vi_object) {
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

.lrm_scaled_sufficient_statistics <- function(stats, batch_idx = NULL) {
  n <- stats$n
  if (is.null(batch_idx)) return(list(
    XtX = stats$XtX, Xty = stats$Xty,
    yty = if (is.null(stats$yty)) drop(crossprod(stats$y)) else stats$yty,
    weight = 1, size = n))
  B <- as.integer(batch_idx)
  if (!length(B) || anyNA(B) || any(B < 1L | B > n) || anyDuplicated(B))
    stop("batch_idx must contain distinct valid observation indices")
  weight <- n / length(B)
  X_B <- stats$X[B, , drop = FALSE]
  y_B <- stats$y[B]
  list(XtX = weight * crossprod(X_B),
       Xty = weight * crossprod(X_B, y_B),
       yty = weight * drop(crossprod(y_B)),
       weight = weight, size = length(B))
}


update_q_residual_variance_lrm <- function(lrm_vi_object, batch_idx = NULL,
                                           learning_rate = 1,
                                           batch_statistics = NULL){
  SigmaBeta <- lrm_vi_object$params$SigmaBeta
  MuBeta <- lrm_vi_object$params$MuBeta
  nu0 <- lrm_vi_object$hyperparams$nu0
  sigma20 <- lrm_vi_object$hyperparams$sigma20
  n <- lrm_vi_object$stats$n
  suff <- if (is.null(batch_statistics))
    .lrm_scaled_sufficient_statistics(lrm_vi_object$stats, batch_idx) else
      batch_statistics
  expected_sse <- suff$yty - 2 * sum(MuBeta * suff$Xty) +
    drop(crossprod(MuBeta, suff$XtX %*% MuBeta)) +
    sum(diag(suff$XtX %*% SigmaBeta))
  expected_sse <- max(expected_sse, 0)
  a_hat <- 0.5 * (n + nu0)
  b_hat <- 0.5 * (nu0 * sigma20 + expected_sse)
  # Inverse-gamma natural parameters are (-(a+1), -b), so convex
  # combinations reduce to the same combination of a and b.
  lrm_vi_object$params$asigma2 <-
    (1 - learning_rate) * lrm_vi_object$params$asigma2 +
    learning_rate * a_hat
  lrm_vi_object$params$bsigma2 <-
    (1 - learning_rate) * lrm_vi_object$params$bsigma2 +
    learning_rate * b_hat
  return(lrm_vi_object)
}

update_q_coefficients_lrm <- function(lrm_vi_object, batch_idx = NULL,
                                      learning_rate = 1,
                                      batch_statistics = NULL){
  asigma2 <- lrm_vi_object$params$asigma2
  bsigma2 <- lrm_vi_object$params$bsigma2
  Beta0 <- lrm_vi_object$hyperparams$Beta0
  Sigma0inv <- lrm_vi_object$hyperparams$Sigma0inv
  suff <- if (is.null(batch_statistics))
    .lrm_scaled_sufficient_statistics(lrm_vi_object$stats, batch_idx) else
      batch_statistics
  expected_precision <- asigma2 / bsigma2
  P_hat <- symmetrize_numeric(
    expected_precision * suff$XtX + Sigma0inv,
    "LRM variational beta precision target")
  h_hat <- as.numeric(expected_precision * suff$Xty + Sigma0inv %*% Beta0)
  # Mix the Gaussian factor in natural coordinates (precision and
  # precision-weighted mean), never in mean/covariance coordinates.
  P_old <- solve_spd(lrm_vi_object$params$SigmaBeta,
                     name = "LRM variational beta covariance")
  h_old <- as.numeric(P_old %*% lrm_vi_object$params$MuBeta)
  P_new <- symmetrize_numeric(
    (1 - learning_rate) * P_old + learning_rate * P_hat,
    "LRM mixed variational beta precision")
  h_new <- (1 - learning_rate) * h_old + learning_rate * h_hat
  SigmaBeta <- solve_spd(P_new, name = "LRM mixed variational beta precision")
  MuBeta <- SigmaBeta %*% h_new
  lrm_vi_object$params$SigmaBeta <- SigmaBeta
  lrm_vi_object$params$MuBeta <- as.numeric(MuBeta)
  return(lrm_vi_object)
}

elbo_lrm <- function(lrm_vi_object){

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
  Elogs2 <- expected_log_inverse_gamma(asigma2, bsigma2)
  Einvs2 <- expected_inverse_inverse_gamma(asigma2, bsigma2)
  Elogpy <- -0.5*n*log(2*pi)-0.5*n*Elogs2-0.5*Einvs2*Essr

  diff <- MuBeta-Beta0
  ElogpBeta <- -0.5*p*log(2*pi)-0.5*logdet_spd(Sigma0) -
    0.5*(sum(diag(Sigma0inv %*% SigmaBeta))+drop(crossprod(diff,Sigma0inv%*%diff)))

  Elogpsigma2 <- a0*log(b0)-lgamma(a0)-(a0+1)*Elogs2-b0*Einvs2

  ElogqBeta <- -0.5*(p*(1+log(2*pi))+logdet_spd(SigmaBeta))

  Elogqsigma2 <- asigma2*log(bsigma2)-lgamma(asigma2)-(asigma2+1)*Elogs2-bsigma2*Einvs2

  elbo <- Elogpy + ElogpBeta + Elogpsigma2 - (ElogqBeta + Elogqsigma2)

  return(as.numeric(elbo))
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
