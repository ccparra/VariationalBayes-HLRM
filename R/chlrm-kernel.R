hyperparams_chlrm <- function(data, K,
                              type = c("OLS","NON_INFORMATIVE"),
                              jitter  = 1e-6) {
  type <- match.arg(type)
  X <- data$X; y <- data$y
  validate_grouped_data(X, y, K)
  if (!is.numeric(jitter) || length(jitter) != 1L || !is.finite(jitter) || jitter <= 0)
    stop("jitter must be a positive finite scalar")
  m <- length(X); p <- ncol(X[[1]])
  n <- sum(lengths(y))

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
  Beta_ols <- colMeans(B_ols)
  Sigma_ols <- if (m > 1L) stats::cov(B_ols) + jitter * diag(p) else
    max(s2_ols, jitter) * diag(p)
  s2_prior <- max(stats::median(s2_ols), jitter)

  if (type == "OLS") {
    Mu0       <- Beta_ols
    Lambda0   <- Sigma_ols
    S0        <- Sigma_ols
    n0        <- p + 2
    a0        <- 2
    nu0       <- 2
    b0        <- a0/s2_prior
    alpha0 <- rep(1, K)
  } else { # NON_INFORMATIVE
    Mu0     <- rep(0, p)
    Lambda0 <- s2_prior * diag(p)
    S0      <- s2_prior * diag(p)
    n0      <- p + 2
    a0      <- 2
    nu0     <- 2
    b0      <- a0/s2_prior
    alpha0 <- rep(1, K)
  }

  validate_hlrm_hyperparams(list(
    Mu0 = as.numeric(Mu0),
    Lambda0 = Lambda0,
    n0 = n0,
    S0 = S0,
    nu0 = nu0,
    a0 = a0, b0 = b0,
    alpha0 = alpha0,
    sigma20 = s2_prior,
    type = type,
    Betaj_ols= B_ols
  ), p, m, "CHLRM", K)
}

# Conditional likelihood used by the Gibbs log joint, not by marginal scoring.
loglik_pointwise_chlrm <- function(y_list, X_list, Betaj, sigma2j) {
  lj <- lapply(seq_along(X_list), function(j) {
    mu_j <- as.vector(X_list[[j]] %*% Betaj[j, ])
    stats::dnorm(y_list[[j]], mean = mu_j, sd = sqrt(max(sigma2j[j], 1e-12)), log = TRUE)
  })
  unlist(lj, use.names = FALSE)
}

loglik_total_chlrm <- function(y_list, X_list, Betaj, sigma2j) {
  sum(loglik_pointwise_chlrm(y_list, X_list, Betaj, sigma2j))
}
