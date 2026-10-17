hyperparams_hlrm <- function(data,
                             type   = c("OLS", "NON_INFORMATIVE"),
                             jitter = 1e-6) {
  type <- match.arg(type)

  X <- data$X
  y <- data$y
  validate_grouped_data(X, y)
  n <- sum(lengths(y))
  m <- length(X)
  p <- ncol(X[[1]])
  if (!is.numeric(jitter) || length(jitter) != 1L || !is.finite(jitter) || jitter <= 0)
    stop("jitter must be a positive finite scalar")

  # Group-specific OLS coefficients and residual variances.
  B_ols  <- matrix(NA_real_, m, p)
  s2_ols <- numeric(m)

  for (j in seq_len(m)) {
    XtX <- crossprod(X[[j]])
    Xty <- crossprod(X[[j]], y[[j]])
    bj  <- solve_spd(XtX + jitter * diag(p), Xty,
                     "regularized group XtX")
    B_ols[j, ] <- bj

    rj <- y[[j]] - X[[j]] %*% bj
    s2_ols[j] <- mean(drop(rj)^2)
  }

  # OLS statistics used by the explicit OLS initialization strategy.
  Beta_ols  <- colMeans(B_ols)
  Sigma_ols <- if (m > 1L) stats::cov(B_ols) + jitter * diag(p) else
    max(s2_ols, jitter) * diag(p)
  s2_prior  <- max(stats::median(s2_ols), jitter)

  if (type == "OLS") {
    Mu0     <- as.numeric(Beta_ols)
    Lambda0 <- Sigma_ols
    S0      <- Sigma_ols
    n0      <- p + 2
    nu0     <- 2
    a0      <- 2
    b0      <- a0 / s2_prior
  } else {  # NON_INFORMATIVE
    Mu0     <- rep(0, p)
    Lambda0 <- s2_prior*diag(p)
    S0      <- s2_prior*diag(p)
    n0      <- p + 2
    nu0     <- 2
    a0      <- 2
    b0      <- a0/s2_prior
  }

  validate_hlrm_hyperparams(list(
    Mu0      = Mu0,
    Lambda0  = Lambda0,
    n0       = n0,
    S0       = S0,
    nu0      = nu0,
    a0       = a0,
    b0       = b0,
    sigma20  = s2_prior,
    type     = type,
    Betaj_ols = B_ols
  ), p, m, "HLRM")
}


loglik_pointwise_hlrm <- function(y_list, X_list, Betaj, sigma2j) {
  lj <- lapply(seq_along(X_list), function(j) {
    muj <- as.vector(X_list[[j]] %*% Betaj[j, ])
    sj <- sqrt(sigma2j[j])
    stats::dnorm(y_list[[j]], mean = muj, sd = sj , log = TRUE)
  })
  return(unlist(lj,use.names = FALSE))
}

loglik_total_hlrm <- function(y_list, X_list, Betaj, sigma2j) {
  return(sum(loglik_pointwise_hlrm(y_list, X_list, Betaj, sigma2j)))
}
