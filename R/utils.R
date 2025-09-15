n_samples_mcmc <- function(run_settings){
  if (!is.null(run_settings$thin)) {
    B <- run_settings$samples*run_settings$thin + run_settings$burn_in
  } else {
    B <- run_settings$samples + run_settings$burn_in  
  }
  return(B) 
}

thin_chain_mcmc <- function(b,run_settings){
  if (b > run_settings$burn_in) {
    if (!is.null(run_settings$thin)) {
      if (b %% run_settings$thin == 0) {
        i <- (b - run_settings$burn_in) / run_settings$thin
      } else {
        i <- NULL
      }
    } else {
      i <- b
    }
  } else {
    i <- NULL
  }
  return(i)
}

logmeanexp <- function(x) {
  m <- max(x)
  m + log(mean(exp(x - m)))
}

logsumexp <- function(x) {
  x_max <- max(x)
  x_max + log(sum(exp(x - x_max)))
}

mse <- function(y,yhat){
  return(mean((y-yhat)^2))
}

logdet <- function(M) {
  M <- as.matrix(M)
  as.numeric(determinant(M, logarithm = TRUE)$modulus)
}

.sym <- function(M) 0.5 * (M + t(M))

step_size <- function(t, tau = 10, kappa = 0.7) (t + tau)^(-kappa)

log_multigamma <- function(a, p) (p * (p - 1) / 4) * log(pi) + sum(lgamma(a + (1 - seq_len(p)) / 2))

ElogdetIW <- function(nu, S) { 
  p <- nrow(S)
  return(-sum(digamma((nu + 1 - seq_len(p)) / 2)) - p * log(2) + logdet(S))
}

log2pi <- log(2 * pi)

E_log_sigma2 <- function(a, b) log(b) - digamma(a) 
E_inv_sigma2 <- function(a, b) a / b

E_log_omega   <- function(alpha) digamma(alpha) - digamma(sum(alpha))


r2 <- function(y, yhat) {
  y    <- as.numeric(y)
  yhat <- as.numeric(yhat)
  sse <- sum((y - yhat)^2)
  sst <- sum((y - mean(y))^2)
 return(1 - sse/sst)
}


standardize_Xy_list<- function(X, y) {
  m <- length(X); p <- ncol(X[[1]])
  X_stack <- do.call(rbind, X)
  y_stack <- unlist(y)
  
  mu_x <- colMeans(X_stack[, -1, drop = FALSE])
  sd_x <- apply(X_stack[, -1, drop = FALSE], 2, sd)
  sd_x[sd_x == 0] <- 1
  
  mu_y <- mean(y_stack)
  sd_y <- stats::sd(y_stack); if (sd_y == 0) sd_y <- 1
  
  X_std <- vector("list", m)
  y_std <- vector("list", m)
  for (j in seq_len(m)) {
    Xj <- X[[j]]
    Xj[, -1] <- sweep(Xj[, -1, drop = FALSE], 2, mu_x, "-")
    Xj[, -1] <- sweep(Xj[, -1, drop = FALSE], 2, sd_x, "/")
    X_std[[j]] <- Xj
    y_std[[j]] <- (y[[j]] - mu_y) / sd_y
  }
  
  list(
    X = X_std, y = y_std,
    scale = list(mu_x = mu_x, sd_x = sd_x, mu_y = mu_y, sd_y = sd_y)
  )
}