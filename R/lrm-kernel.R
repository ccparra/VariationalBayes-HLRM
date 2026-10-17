hyperparams_lrm <- function(data,
                            type   = c("OLS", "NON_INFORMATIVE"),
                            jitter = 1e-6) {
  type <- match.arg(type)

  X <- as.matrix(data$X)
  y <- as.numeric(data$y)

  n <- nrow(X)
  p <- ncol(X)
  if (!is.numeric(jitter) || length(jitter) != 1L || !is.finite(jitter) || jitter <= 0)
    stop("jitter must be a positive finite scalar")

  # OLS with minimal regularization for numerical stability.
  XtX <- crossprod(X)
  Xty <- crossprod(X, y)
  ridge <- XtX + jitter * diag(p)
  Beta_OLS  <- solve_spd(ridge, Xty, "regularized XtX")
  r         <- y - as.vector(X %*% Beta_OLS)
  sigma2_OLS <- max(sum(r^2) / max(1, (n - p)), jitter)

  if (type == "OLS") {
    # Prior centrado en OLS; var ~ (n * sigma2_OLS) (X'X)^{-1}
    Sigma0 <- n * sigma2_OLS * solve_spd(ridge, name = "regularized XtX")
    Beta0  <- as.numeric(Beta_OLS)
    nu0    <- 2
    sigma20 <- sigma2_OLS
  } else { # NON_INFORMATIVE
    Sigma0  <- n * sigma2_OLS * diag(p)
    Beta0   <- rep(0, p)
    nu0     <- 2
    sigma20 <- sigma2_OLS
  }

  validate_lrm_hyperparams(list(
    Beta0     = Beta0,
    Sigma0    = Sigma0,
    nu0       = nu0,
    sigma20   = sigma20,
    type      = type
  ), p)
}


loglik_pointwise_lrm <- function(y, X, beta, sigma2) {
  mu <- as.vector(X %*% beta)
  stats::dnorm(y, mean = mu, sd = sqrt(sigma2), log = TRUE)
}

loglik_total_lrm <- function(y, X, beta, sigma2) {
  sum(loglik_pointwise_lrm(y, X, beta, sigma2))
}
