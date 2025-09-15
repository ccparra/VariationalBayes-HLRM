source('R/utils.R')
source('R/CHLRM/mcmc/utils.R')

sim_chlrm <- function(m, nj, Betak, sigma2k, p,
                      z = NULL, omega = NULL, K,
                      dummies = 0,
                      seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  # --- tamaños por grupo y por clúster
  if (length(nj) == 1L)      nj      <- rep(nj, m)
  if (length(sigma2k) == 1L) sigma2k <- rep(sigma2k, K)
  
  # --- asignaciones de clúster por grupo
  if (is.null(z)) {
    if (is.null(omega)) omega <- rep(1 / K, K) else omega <- omega / sum(omega)
    z <- sample.int(K, size = m, replace = TRUE, prob = omega)
  } else {
    if (is.null(omega)) {
      tab <- tabulate(z, nbins = K); omega <- tab / sum(tab)
      if (all(tab == 0)) omega <- rep(1 / K, K)
    } else {
      omega <- omega / sum(omega)
    }
  }
  
  # --- chequeos básicos
  if (!is.matrix(Betak)) stop("Betak debe ser matriz K x p.")
  if (ncol(Betak) != p)  stop("ncol(Betak) debe ser igual a p.")
  if (length(sigma2k) != K) stop("sigma2k debe tener longitud K.")
  if (dummies < 0 || dummies != as.integer(dummies)) stop("dummies debe ser entero >= 0.")
  
  num_covars <- p - 1L - dummies
  if (num_covars < 0L) stop("p - 1 - dummies < 0: reduce dummies o incrementa p.")
  
  # --- helpers: prob. para cada dummy d
  #     p_d = (0.2*d) / (0.8 + 0.2*d)
  dummy_prob <- function(d) (0.2 * d) / (0.8 + 0.2 * d)
  
  X_list <- vector("list", m)
  y_list <- vector("list", m)
  g_list <- vector("list", m)
  
  for (j in seq_len(m)) {
    n  <- nj[j]
    kj <- z[j]                    # clúster del grupo j
    beta_k  <- Betak[kj, ]
    s2_k    <- sigma2k[kj]
    
    # --- X: [intercepto | dummies | numéricas]
    # Intercepto
    Xj <- matrix(1, nrow = n, ncol = 1L)
    
    # Dummies (si aplica)
    if (dummies > 0L) {
      Xd <- matrix(NA_real_, nrow = n, ncol = dummies)
      for (d in seq_len(dummies)) {
        p_d <- dummy_prob(d)
        Xd[, d] <- rbinom(n, size = 1L, prob = p_d)
      }
      Xj <- cbind(Xj, Xd)
    }
    
    # Numéricas (si aplica) — media depende del grupo j
    if (num_covars > 0L) {
      mu_x  <- 30 + 5 * j
      Xnum  <- matrix(rnorm(n * num_covars, mean = mu_x, sd = 10), nrow = n, ncol = num_covars)
      # (opcional) redondeo como en tu ejemplo:
      Xnum  <- round(Xnum, 2)
      Xj    <- cbind(Xj, Xnum)
    }
    
    # --- y = X beta_k + eps, eps ~ N(0, sigma2_k)
    mu_j <- as.vector(Xj %*% beta_k)
    y_j  <- mu_j + rnorm(n, mean = 0, sd = sqrt(s2_k))
    
    X_list[[j]] <- Xj
    y_list[[j]] <- y_j
    g_list[[j]] <- rep(j, n)
  }
  
  # salida compatible con tu pipeline original
  list(
    X = X_list,
    y = y_list,
    group = g_list,
    cluster_info = list(
      K = K,
      z = z,
      omega = omega
    ),
    real_parameters = list(
      Betak   = Betak,
      sigma2k = as.numeric(sigma2k),
      p = p,
      m = m,
      nj = as.integer(nj),
      dummies = as.integer(dummies),
      num_covars = as.integer(num_covars)
    )
  )
}



df_to_chlrm_lists <- function(data, y, x, group) {
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

hyperparams_chlrm <- function(data, K,
                              type = c("OLS","NON_INFORMATIVE"),
                              jitter  = 1e-6) {
  type <- match.arg(type)
  X <- data$X; y <- data$y
  m <- length(X); p <- ncol(X[[1]])
  n <- sum(lengths(y))
  
  B_ols  <- matrix(NA_real_, m, p)
  s2_ols <- numeric(m)
  for (j in seq_len(m)) {
    XtX <- crossprod(X[[j]])
    Xty <- crossprod(X[[j]], y[[j]])
    bj  <- solve(XtX + jitter * diag(p), Xty)
    B_ols[j, ] <- bj
    rj <- y[[j]] - X[[j]] %*% bj
    s2_ols[j] <- mean(drop(rj)^2)
  }
  Beta_ols <- colMeans(B_ols)
  Sigma_ols <- stats::cov(B_ols) + jitter * diag(p)
  s2_prior <- stats::median(s2_ols)
  
  if (type == "OLS") {
    Mu0       <- Beta_ols
    Lambda0   <- n*Sigma_ols
    S0        <- n*Sigma_ols
    n0        <- p + 2
    a0        <- 2
    nu0       <- 2 
    b0        <- a0/s2_prior
    alpha0 <- rep(1, K)
  } else { # NON_INFORMATIVE
    Mu0     <- rep(0, p)
    Lambda0 <-n* s2_prior*diag(p)
    S0      <- diag(p)
    n0      <- p + 2
    a0      <- 2
    nu0     <- 2
    b0      <- a0/s2_prior
    alpha0 <- rep(1, K)
  }
  
  list(
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
  )
}

# =========================
# 2) WAIC punto-a-punto (robusto)
#    Funciona para VI o MCMC (solo requiere chains con Betaj y sigma2j)
# =========================
logmeanexp_safe <- function(x) {
  m <- max(x)
  if (!is.finite(m)) return(-Inf)
  m + log(mean(exp(x - m)))
}

waic_chlrm <- function(y_list, X_list, chains, nj) {
  Betaj   <- chains$Betaj       # S × m × p
  sigma2j <- chains$sigma2j     # S × m
  S <- dim(Betaj)[1]
  
  n <- sum(nj)
  idx_end <- cumsum(nj)
  idx_ini <- c(1L, head(idx_end, -1L) + 1L)
  
  ll_point_mat <- matrix(NA_real_, n, S)
  
  for (s in seq_len(S)) {
    for (j in seq_along(y_list)) {
      mu_j <- as.vector(X_list[[j]] %*% Betaj[s, j, ])
      s2   <- max(sigma2j[s, j], 1e-12)  # clamp numérico
      ll_j <- stats::dnorm(y_list[[j]], mean = mu_j, sd = sqrt(s2), log = TRUE)
      ll_point_mat[idx_ini[j]:idx_end[j], s] <- ll_j
    }
  }
  
  # lppd_i y dos formas de p_waic por punto (más estable en mezclas)
  lppd_i     <- apply(ll_point_mat, 1, logmeanexp_safe)
  var_i      <- apply(ll_point_mat, 1, stats::var)
  mean_ll_i  <- rowMeans(ll_point_mat)
  pwaic2_i   <- 2 * (lppd_i - mean_ll_i)
  
  p_i <- pmin(var_i, pwaic2_i)
  p_i[p_i < 0] <- 0
  
  lppd   <- sum(lppd_i)
  p_waic <- sum(p_i)
  waic   <- -2 * (lppd - p_waic)
  
  list(WAIC = waic, lppd = lppd, p_waic = p_waic)
}

# =========================
# 3) DIC clásico (devianza); también aplica a VI/MCMC
# =========================
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

dic_chlrm <- function(y_list, X_list, chains, nj) {
  Betaj   <- chains$Betaj       # S × m × p
  sigma2j <- chains$sigma2j     # S × m
  S <- dim(Betaj)[1]
  
  # D(θ^s) por draw
  D_draws <- numeric(S)
  for (s in seq_len(S)) {
    s2j_s <- pmax(sigma2j[s, ], 1e-12)
    D_draws[s] <- -2 * loglik_total_chlrm(
      y_list  = y_list,
      X_list  = X_list,
      Betaj   = Betaj[s, , ],
      sigma2j = s2j_s
    )
  }
  Dbar <- mean(D_draws)
  
  # D(θ̂) con promedios
  betaj_bar <- apply(Betaj, c(2, 3), mean)
  s2j_bar   <- pmax(colMeans(sigma2j), 1e-12)
  
  D_hat <- -2 * loglik_total_chlrm(
    y_list  = y_list,
    X_list  = X_list,
    Betaj   = betaj_bar,
    sigma2j = s2j_bar
  )
  
  p_D <- Dbar - D_hat
  DIC <- D_hat + 2 * p_D
  
  list(DIC = DIC, Dbar = Dbar, D_hat = D_hat, p_D = p_D)
}



yrep_chlrm <- function(X_list, chains, nj = NULL, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  Betaj   <- chains$Betaj   # S × m × p
  sigma2j <- chains$sigma2j # S × m
  S <- dim(Betaj)[1]
  n   <- sum(nj)
  
  # índices de pegado
  idx_end <- cumsum(nj)
  idx_ini <- c(1L, head(idx_end, -1L) + 1L)
  Yrep <- matrix(NA_real_, nrow = n, ncol = S)
  
  for (s in seq_len(S)) {
    for (j in seq_along(X_list)) {
      mu_j <- as.vector(X_list[[j]] %*% Betaj[s, j, ])
      sd_j <- rep(sqrt(sigma2j[s, j]), nj[j])
      Yrep[idx_ini[j]:idx_end[j], s] <- stats::rnorm(nj[j], mean = mu_j, sd = sd_j)
    }
  }
  Yrep
}

ppval_chlrm <- function(y_list, X_list, chains, nj, return_yrep = FALSE, seed = NULL) {
  y <- unlist(y_list, use.names = FALSE)
  Yrep <- yrep_chlrm(X_list, chains, nj = nj, seed = seed)
  rep_IQR <- apply(Yrep, 2, IQR, na.rm = TRUE)
  rep_min    <- apply(Yrep, 2, min)
  rep_max    <- apply(Yrep, 2, max)
  rep_mean   <- colMeans(Yrep)
  rep_median <- apply(Yrep, 2, median)
  rep_sd     <- apply(Yrep, 2, sd)
  
  repdata <- cbind(
    min  = rep_min,
    max  = rep_max,
    IQR  = rep_IQR,
    mean = rep_mean,
    med  = rep_median,
    sd   = rep_sd
  )
  
  t0 <- c(min(y), max(y), IQR(y), mean(y), median(y), sd(y))
  names(t0) <- colnames(repdata)
  
  ppp <- colMeans(sweep(repdata, 2, t0, `>=`))
  
  out <- list(repdata = as.data.frame(repdata), t0 = t0, ppp = ppp)
  if (return_yrep) out$y_rep <- Yrep
  out
}







