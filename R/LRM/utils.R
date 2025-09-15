source('R/utils.R')
source('R/LRM/mcmc/utils.R')

# =============== LRM: recta (d=1) o hiperplano (d>1) + jitter en X =================
sim_lrm <- function(n, beta, sigma2, x_range, seed = NULL,
                    jitter_frac = 0.01, jitter_sd = NULL, clip = TRUE) {
  if (!is.null(seed)) set.seed(seed)
  
  p  <- length(beta)
  d  <- p - 1
  XR <- matrix(x_range, nrow = d, ncol = 2, byrow = TRUE)  # d×2
  
  # Construcción base (recta si d=1; grid si d>1)
  if (d == 1L) {
    x1    <- seq(XR[1, 1], XR[1, 2], length.out = n)
    Xpred <- matrix(x1, ncol = 1)
  } else {
    base  <- floor(n^(1/d))
    sizes <- rep(base, d)
    while (prod(sizes) < n) {
      for (i in seq_len(d)) {
        if (prod(sizes) >= n) break
        sizes[i] <- sizes[i] + 1L
      }
    }
    axes <- lapply(seq_len(d), function(i) seq(XR[i,1], XR[i,2], length.out = sizes[i]))
    G    <- do.call(expand.grid, axes)
    #Xpred <- as.matrix(G)[sample(nrow(G), n), , drop=FALSE]
    Xpred <- as.matrix(G)[seq_len(n), , drop = FALSE]
  }
  
  # ---- JITTER en X (por dimensión) ----
  rng <- XR[,2] - XR[,1]
  sdv <- if (is.null(jitter_sd)) jitter_frac * rng else as.numeric(jitter_sd)
  if (d == 1L) sdv <- sdv[1]
  for (i in seq_len(d)) {
    Xpred[, i] <- Xpred[, i] + rnorm(nrow(Xpred), sd = sdv[i])
    if (clip) Xpred[, i] <- pmin(pmax(Xpred[, i], XR[i,1]), XR[i,2])
  }
  
  X      <- cbind(1, Xpred)
  y_mean <- as.numeric(X %*% beta)
  y      <- y_mean + rnorm(n, sd = sqrt(sigma2))
  
  list(
    X = X, y = y,
    signal = y_mean, eps = y - y_mean,
    real_parameters = list(beta = beta, sigma2 = sigma2, x_range = XR,
                           jitter_sd = sdv, clip = clip)
  )
}




hyperparams_lrm <- function(data_list,
                            type   = c("OLS", "NON_INFORMATIVE"),
                            jitter = 1e-6) {
  type <- match.arg(type)
  
  X <- as.matrix(data_list$X)
  y <- as.numeric(data_list$y)
  
  n <- nrow(X)
  p <- ncol(X)
  
  # OLS con regularización mínima para estabilidad
  XtX <- crossprod(X)
  Xty <- crossprod(X, y)
  Beta_OLS  <- solve(XtX + jitter * diag(p), Xty)
  r         <- y - as.vector(X %*% Beta_OLS)
  sigma2_OLS <- sum(r^2) / max(1, (n - p))
  
  if (type == "OLS") {
    # Prior centrado en OLS; var ~ (n * sigma2_OLS) (X'X)^{-1}
    Sigma0 <- n * sigma2_OLS * solve(XtX + jitter * diag(p))
    Beta0  <- as.numeric(Beta_OLS)
    nu0    <- 2
    sigma20 <- sigma2_OLS
  } else { # NON_INFORMATIVE
    Sigma0  <- n * sigma2_OLS * diag(p)
    Beta0   <- rep(0, p)
    nu0     <- 2
    sigma20 <- sigma2_OLS
  }
  
  list(
    Beta0     = Beta0,
    Sigma0    = Sigma0,
    Sigma0inv = solve(Sigma0),
    nu0       = nu0,
    sigma20   = sigma20,
    type      = type
  )
}


loglik_pointwise_lrm <- function(y, X, beta, sigma2) {
  mu <- as.vector(X %*% beta)
  stats::dnorm(y, mean = mu, sd = sqrt(sigma2), log = TRUE)
}

loglik_total_lrm <- function(y, X, beta, sigma2) {
  sum(loglik_pointwise_lrm(y, X, beta, sigma2))
}

waic_lrm <- function(y, X, chains) {
  Beta_mat <- chains$Beta
  s2_vec   <- chains$sigma2
  
  S <- nrow(Beta_mat)
  n <- nrow(X)
  
  ll_point_mat <- matrix(NA_real_, n, S)
  for (s in seq_len(S)) {
    ll_point_mat[, s] <- loglik_pointwise_lrm(y, X, beta = Beta_mat[s, ], sigma2 = s2_vec[s])
  }
  lppd <- sum(apply(ll_point_mat, 1, logmeanexp))
  p_waic <- sum(apply(ll_point_mat, 1, stats::var, na.rm = TRUE))
  waic <- -2 * (lppd - p_waic)
  return(list(WAIC = waic, lppd = lppd, p_waic = p_waic))
}

dic_lrm <- function(y, X, chains) {
  Beta_mat <- chains$Beta
  s2_vec   <- chains$sigma2
  
  S <- nrow(Beta_mat)
  
  D_draws <- numeric(S)
  for (s in seq_len(S)) {
    D_draws[s] <- -2 * loglik_total_lrm(y, X, beta = Beta_mat[s, ], sigma2 = s2_vec[s])
  }
  Dbar <- mean(D_draws, na.rm = TRUE)
  
  beta_bar <- colMeans(Beta_mat, na.rm = TRUE)
  s2_bar   <- mean(s2_vec, na.rm = TRUE)
  D_hat    <- -2 * loglik_total_lrm(y, X, beta = beta_bar, sigma2 = s2_bar)
  
  p_D <- Dbar - D_hat
  DIC <- D_hat + 2 * p_D
  
  list(DIC = DIC, Dbar = Dbar, D_hat = D_hat, p_D = p_D)
}


yrep_lrm <- function(X, chains, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  Beta   <- chains$Beta
  sigma2 <- chains$sigma2
  
  n <- nrow(X); S <- nrow(Beta)
  Mu <- X %*% t(Beta) 
  
  yrep_vec <- rnorm(
    n = n * S,
    mean = as.vector(Mu),
    sd   = rep(sqrt(sigma2), each = n)
  )
  matrix(yrep_vec, nrow = n, ncol = S)
}

ppval_lrm <- function(y, X, chains, return_yrep = FALSE, seed = NULL) {
  y <- drop(y)
  Yrep <- yrep_lrm(X, chains,seed)
  
  rep_min   <- apply(Yrep, 2, min)
  rep_max   <- apply(Yrep, 2, max)
  rep_IQR   <- apply(Yrep, 2, IQR)
  rep_mean  <- colMeans(Yrep)
  rep_median<- apply(Yrep, 2, median)
  rep_sd    <- apply(Yrep, 2, sd)
  
  repdata <- cbind(
    min   = rep_min,
    max   = rep_max,
    IQR   = rep_IQR,
    mean  = rep_mean,
    med   = rep_median,
    sd    = rep_sd
  )
  
  t0 <- c(
    min(y), max(y), IQR(y),
    mean(y), median(y), sd(y)
  )
  names(t0) <- colnames(repdata)
  
  ppp <- colMeans(sweep(repdata, 2, t0, `>=`))
  
  out <- list(
    repdata = as.data.frame(repdata),
    t0      = t0,
    ppp     = ppp
  )
  if (return_yrep) out$y_rep <- Yrep
  out
}








