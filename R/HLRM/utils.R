source('R/utils.R')
source('R/HLRM/mcmc/utils.R')

# sim_hlrm: p = # (intercepto + covariables)
# - Construye X con intercepto en la primera columna.
# - 'Beta_j' es matriz (m x p), 'sigma2_j' vector de longitud m.
sim_hlrm <- function(m, nj, Beta_j, sigma2_j, p,
                     dummies = 0,
                     seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  if (length(nj) == 1L)        nj       <- rep(nj, m)
  if (length(sigma2_j) == 1L)  sigma2_j <- rep(sigma2_j, m)
  
  if (!is.matrix(Beta_j)) stop("Beta_j debe ser matriz m x p.")
  if (ncol(Beta_j) != p)  stop("ncol(Beta_j) debe ser igual a p (incluye intercepto).")
  if (nrow(Beta_j) != m)  stop("nrow(Beta_j) debe ser igual a m.")
  if (length(sigma2_j) != m) stop("sigma2_j debe tener longitud m.")
  if (dummies < 0 || dummies != as.integer(dummies)) stop("'dummies' debe ser entero >= 0.")
  
  num_covars <- p - 1L - dummies
  if (num_covars < 0L) stop("p - 1 - dummies < 0: reduce dummies o incrementa p.")
  
  dummy_prob <- function(d) (0.2 * d) / (0.8 + 0.2 * d)
  
  X_list <- vector("list", m)
  y_list <- vector("list", m)
  g_list <- vector("list", m)
  
  for (j in seq_len(m)) {
    n  <- nj[j]
    bj <- Beta_j[j, ]
    s2 <- sigma2_j[j]
    
    # Intercepto
    Xj <- matrix(1, nrow = n, ncol = 1L)
    
    # Dummies
    if (dummies > 0L) {
      Xd <- matrix(NA_real_, nrow = n, ncol = dummies)
      for (d in seq_len(dummies)) {
        Xd[, d] <- rbinom(n, size = 1L, prob = dummy_prob(d))
      }
      Xj <- cbind(Xj, Xd)
    }
    
    # Numéricas
    if (num_covars > 0L) {
      mu_x <- 30 + 5 * j
      Xnum <- matrix(rnorm(n * num_covars, mean = mu_x, sd = 10), nrow = n, ncol = num_covars)
      Xnum <- round(Xnum, 2)
      Xj   <- cbind(Xj, Xnum)
    }
    
    mu_j <- as.vector(Xj %*% bj)
    y_j  <- mu_j + rnorm(n, mean = 0, sd = sqrt(s2))
    
    X_list[[j]] <- Xj
    y_list[[j]] <- y_j
    g_list[[j]] <- rep(j, n)
  }
  
  list(
    X = X_list,
    y = y_list,
    group = g_list,
    cluster_info = NULL,
    real_parameters = list(
      Beta_j   = Beta_j,
      sigma2_j = as.numeric(sigma2_j),
      p = p, m = m, nj = as.integer(nj),
      dummies = as.integer(dummies),
      num_covars = as.integer(num_covars)
    )
  )
}


df_to_hlrm_lists <- function(data, y, x, group) {
  g_raw <- data[[group]]
  g_lab <- if (is.numeric(g_raw)) {
    format(g_raw, trim = TRUE, scientific = FALSE)
  } else {
    as.character(g_raw)
  }
  lev <- unique(g_lab)
  m   <- length(lev)
  idx_by_g <- lapply(lev, function(lbl) which(g_lab == lbl))
  
  Xcov  <- as.matrix(data[, x, drop = FALSE]); storage.mode(Xcov) <- "double"
  Xfull <- cbind(1, Xcov); dimnames(Xfull) <- NULL
  yvec  <- as.numeric(data[[y]])
  
  X_list <- lapply(idx_by_g, function(idx){ Xi <- Xfull[idx, , drop = FALSE]; dimnames(Xi) <- NULL; Xi })
  y_list <- lapply(idx_by_g, function(idx){ yi <- yvec[idx]; names(yi) <- NULL; yi })
  g_list <- lapply(seq_len(m), function(j) rep.int(j, length(idx_by_g[[j]])))
  
  names(X_list) <- names(y_list) <- names(g_list) <- NULL
  
  mapping <- data.frame(group_id = seq_len(m), group_label = lev, stringsAsFactors = FALSE)
  
  list(X = X_list, y = y_list, group = g_list, mapping = mapping)
}

hyperparams_hlrm <- function(data,
                             type   = c("OLS", "NON_INFORMATIVE"),
                             jitter = 1e-6) {
  type <- match.arg(type)
  
  X <- data$X
  y <- data$y
  n <- sum(lengths(y))
  m <- length(X)
  p <- ncol(X[[1]])
  
  # ----- OLS por grupo: Betas y varianzas residuales -----
  B_ols  <- matrix(NA_real_, m, p)
  s2_ols <- numeric(m)
  
  for (j in seq_len(m)) {
    XtX <- crossprod(X[[j]])
    Xty <- crossprod(X[[j]], y[[j]])
    bj  <- solve(XtX + jitter * diag(p), Xty)     # regularización mínima
    B_ols[j, ] <- bj
    
    rj <- y[[j]] - X[[j]] %*% bj
    s2_ols[j] <- mean(drop(rj)^2)
  }
  
  # Estadísticos OLS para fijar priors tipo "OLS"
  Beta_ols  <- colMeans(B_ols)
  Sigma_ols <- stats::cov(B_ols) + jitter * diag(p)
  s2_prior  <- stats::median(s2_ols)
  
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
  
  list(
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
  )
}


loglik_pointwise_hlrm <- function(y_list, X_list, betaj, sigma2j) {
  lj <- lapply(seq_along(X_list), function(j) {
    muj <- as.vector(X_list[[j]] %*% betaj[j, ])
    sj <- sqrt(sigma2j[j])
    stats::dnorm(y_list[[j]], mean = muj, sd = sj , log = TRUE)
  })
  return(unlist(lj,use.names = FALSE))
}

loglik_total_hlrm <- function(y_list, X_list, betaj, sigma2j) {
  return(sum(loglik_pointwise_hlrm(y_list, X_list, betaj, sigma2j)))
}

waic_hlrm <- function(y_list, X_list, chains) {
  Betaj   <- chains$Betaj     
  sigma2j <- chains$sigma2j  
  S <- dim(Betaj)[1]
  n <-sum(lengths(y_list))
  
  ll_point_mat <- matrix(NA_real_, n, S)
  for (s in seq_len(S)) {
    ll_point_mat[, s] <- loglik_pointwise_hlrm(
      y_list  = y_list,
      X_list  = X_list,
      betaj   = Betaj[s, , ],
      sigma2j = sigma2j[s, ]
    )
  }
  lppd   <- sum(apply(ll_point_mat, 1, logmeanexp))
  p_waic <- sum(apply(ll_point_mat, 1, stats::var, na.rm = TRUE))
  waic   <- -2 * (lppd - p_waic)
  return(list(WAIC = waic, lppd = lppd, p_waic = p_waic))
}

dic_hlrm <- function(y_list, X_list, chains) {
  Betaj   <- chains$Betaj    
  sigma2j <- chains$sigma2j  
  S <- dim(Betaj)[1]
  
  D_draws <- numeric(S)
  for (s in seq_len(S)) {
    D_draws[s] <- -2 * loglik_total_hlrm(
      y_list  = y_list,
      X_list  = X_list,
      betaj   = Betaj[s, , ],
      sigma2j = sigma2j[s, ]
    )
  }
  Dbar <- mean(D_draws, na.rm = TRUE)
  
  betaj_bar <- apply(Betaj, c(2, 3), mean, na.rm = TRUE) 
  s2j_bar   <- colMeans(sigma2j, na.rm = TRUE)
  D_hat <- -2 * loglik_total_hlrm(
    y_list  = y_list,
    X_list  = X_list,
    betaj   = betaj_bar,
    sigma2j = s2j_bar
  )
  
  p_D <- Dbar - D_hat
  DIC <- D_hat + 2 * p_D
  return(list(DIC = DIC, Dbar = Dbar, D_hat = D_hat, p_D = p_D))
}

yrep_hlrm <- function(X_list, chains, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  Betaj   <- chains$Betaj 
  sigma2j <- chains$sigma2j
  S <- dim(Betaj)[1]
  n_j <- vapply(X_list, nrow, integer(1L))
  n   <- sum(n_j)
  
  Yrep <- matrix(NA_real_, nrow = n, ncol = S)
  for (s in seq_len(S)) {
    mu_vec <- unlist(lapply(seq_along(X_list), function(j) {
      as.vector(X_list[[j]] %*% as.numeric(Betaj[s, j, ]))
    }), use.names = FALSE)
    sd_vec <- unlist(lapply(seq_along(X_list), function(j) {
      rep(sqrt(sigma2j[s, j]), nrow(X_list[[j]]))
    }), use.names = FALSE)
    Yrep[, s] <- stats::rnorm(n = n, mean = mu_vec, sd = sd_vec)
  }
  return(Yrep)
}

ppval_hlrm <- function(y_list, X_list, chains, return_yrep = FALSE, seed = NULL) {
  y <- unlist(y_list, use.names = FALSE)
  Yrep <- yrep_hlrm(X_list, chains, seed)
  
  rep_min    <- apply(Yrep, 2, min)
  rep_max    <- apply(Yrep, 2, max)
  rep_IQR    <- apply(Yrep, 2, IQR)
  rep_mean   <- colMeans(Yrep)
  rep_median <- apply(Yrep, 2, median)
  rep_sd     <- apply(Yrep, 2, sd)
  
  repdata <- cbind(
    min   = rep_min,
    max   = rep_max,
    IQR   = rep_IQR,
    mean  = rep_mean,
    med   = rep_median,
    sd    = rep_sd
  )
  
  t0 <- c(min(y), max(y), IQR(y), mean(y), median(y), sd(y))
  names(t0) <- colnames(repdata)
  
  ppp <- colMeans(sweep(repdata, 2, t0, `>=`))
  
  out <- list(
    repdata = as.data.frame(repdata),
    t0      = t0,
    ppp     = ppp
  )
  if (return_yrep) out$y_rep <- Yrep
  return(out)
}






