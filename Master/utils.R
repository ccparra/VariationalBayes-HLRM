# ============================================================
# Runner de simulaciones LRM: CAVI (VI) vs MCMC
# Requiere que existan: sim_lrm(), hyperparams_lrm(), vi_lrm, mcmc_lrm
# ============================================================

# Matriz S×p robusta a vector, p×S, data.frame/list, etc.
.coerce_B_mat <- function(B, p) {
  if (is.list(B) || is.data.frame(B)) B <- as.matrix(do.call(rbind, B))
  B <- as.matrix(B)
  if (is.null(dim(B))) {
    if (length(B) != p) stop(sprintf("Longitud de Beta=%d no coincide con p=%d", length(B), p))
    B <- matrix(B, nrow = 1)  # 1×p
  }
  if (ncol(B) == p) return(B)      # S×p
  if (nrow(B) == p) return(t(B))   # p×S -> S×p
  stop(sprintf("Dimensiones de Beta incompatibles con p=%d: %dx%d", p, nrow(B), ncol(B)))
}



# -- util: validar x_range contra beta (p = 1 + #covariables)
.validate_x_range <- function(beta, x_range) {
  p <- length(beta)
  if (p == 1L) {
    if (!is.null(x_range)) warning("x_range ignorado porque p=1 (solo intercepto).")
    return(NULL)
  }
  if (is.null(x_range)) {
    # default: [-3, 3] para cada covariable
    return(cbind(rep(-3, p - 1L), rep(3, p - 1L)))
  }
  if (!is.matrix(x_range) || ncol(x_range) != 2L || nrow(x_range) != (p - 1L)) {
    stop(sprintf("x_range debe ser matriz de tamaño (p-1)×2; recibido %s×%s para p=%s.",
                 nrow(x_range), ncol(x_range), p))
  }
  x_range
}

# -- corre UNA simulación
run_one_sim_lrm <- function(sim, vi_settings, mcmc_settings, hyper_type = "OLS", verbose = TRUE) {
  # sim: list(beta, sigma2, n=, x_range=, seed_data=, seed_vi=, seed_mcmc=, jitter_frac=)
  beta        <- as.numeric(sim$beta)
  sigma2      <- as.numeric(sim$sigma2)
  n           <- if (!is.null(sim$n)) as.integer(sim$n) else 200L
  jitter_frac <- if (!is.null(sim$jitter_frac)) as.numeric(sim$jitter_frac) else 0
  seed_data   <- sim$seed_data
  seed_vi     <- sim$seed_vi
  seed_mcmc   <- sim$seed_mcmc
  
  x_range <- .validate_x_range(beta, sim$x_range)
  
  # ---- Simular datos ----
  data <- sim_lrm(n = n, beta = beta, sigma2 = sigma2,
                  x_range = x_range, seed = seed_data, jitter_frac = jitter_frac)
  X <- data$X; y <- data$y
  
  # ---- Hiperparámetros (fijado a OLS) ----
  hypers <- hyperparams_lrm(list(X = X, y = y), type = hyper_type)
  
  # ================== VI (CAVI) ==================
  if (!is.null(seed_vi)) set.seed(seed_vi)
  vi <- vi_lrm$new(y, X, hypers, vi_settings)
  results_vi    <- vi$fit(verbose = verbose)
  hat_values_vi <- vi$inferences()
  
  # y_hat VI (media posterior)
  p <- ncol(X)
  y_hat_vi <- tryCatch({
    if (!is.null(results_vi$params$MuBeta)) {
      mu <- as.numeric(results_vi$params$MuBeta)     # longitud p
    } else if (!is.null(hat_values_vi$Beta$mean)) {
      mu <- as.numeric(hat_values_vi$Beta$mean)      # longitud p
    } else stop("VI: no encuentro la media de Beta.")
    as.vector(X %*% mu)
  }, error = function(e) { message("[VI predict] ", e$message); NULL })
  
  goodness_vi <- tryCatch(
    vi$goodness_of_fit(seed = 1),
    error = function(e) { message("[VI GOF] ", e$message); NULL }
  )
  
  # ================== MCMC ==================
  if (!is.null(seed_mcmc)) set.seed(seed_mcmc)
  mcmc <- mcmc_lrm$new(y, X, hypers, mcmc_settings)
  results_mcmc    <- mcmc$fit(verbose = verbose, seed = seed_mcmc)
  hat_values_mcmc <- mcmc$inferences()
  
  # y_hat MCMC (promedio sobre los draws)
  y_hat_mcmc <- tryCatch({
    B <- results_mcmc$chains$Beta         # S×p o p×S
    B <- .coerce_B_mat(B, p)              # S×p
    rowMeans(X %*% t(B))                  # n×1
  }, error = function(e) { message("[MCMC predict] ", e$message); NULL })
  
  goodness_mcmc <- tryCatch(
    mcmc$goodness_of_fit(),
    error = function(e) {
      message("[MCMC GOF] ", e$message, " | usando fallback para MSE y R2.")
      out <- list(waic = NULL, dic = NULL, mse = NA_real_, r2 = NA_real_, ppval = NULL)
      B <- try(.coerce_B_mat(results_mcmc$chains$Beta, p), silent = TRUE)
      if (!inherits(B, "try-error")) {
        mu   <- colMeans(B)
        yhat <- as.vector(X %*% mu)
        res  <- y - yhat
        out$mse <- mean(res^2)
        out$r2  <- 1 - sum(res^2) / sum((y - mean(y))^2)
      }
      out
    }
  )
  
  # ---- Empaquetar TODO y devolver ----
  out <- list(
    model = "LRM",
    data = list(
      X = X, y = y,
      beta_true = beta, sigma2_true = sigma2,
      x_range = x_range, n = n,
      seeds = list(data = seed_data, vi = seed_vi, mcmc = seed_mcmc)
    ),
    hypers = hypers,
    vi = list(
      results    = results_vi,
      hat_values = hat_values_vi,
      gof        = goodness_vi,
      y_hat      = y_hat_vi,
      timing     = results_vi$timing
    ),
    mcmc = list(
      results    = results_mcmc,
      hat_values = hat_values_mcmc,
      gof        = goodness_mcmc,
      y_hat      = y_hat_mcmc,
      timing     = results_mcmc$timing
    )
  )
  return(out)
}



run_simulations_lrm <- function(sim_specs,
                                vi_settings   = list(iterations = 200000, epsilon = 1e-5, samples_gof = 1000),
                                mcmc_settings = list(samples = 20000, burn_in = 10000, thin = 10),
                                hyper_type    = "OLS",
                                verbose = TRUE,
                                save_rds_path = NULL) {
  
  stopifnot(is.list(sim_specs), length(sim_specs) >= 1L)
  
  res_list <- vector("list", length(sim_specs))
  names(res_list) <- if (!is.null(names(sim_specs))) names(sim_specs) else paste0("sim_", seq_along(sim_specs))
  
  msg <- function(...) { cat(sprintf(...), "\n"); flush.console() }
  
  for (i in seq_along(sim_specs)) {
    if (verbose) {
      msg("\n--- Running simulation %d/%d ---", i, length(sim_specs))
      utils::str(sim_specs[[i]], max.level = 1)
    }
    
    # envoltorio que loguea hitos dentro de run_one_sim_lrm
    res_list[[i]] <- tryCatch({
      # marca previa
      msg("[runner] calling run_one_sim_lrm(%d)", i)
      ans <- run_one_sim_lrm(sim_specs[[i]], vi_settings, mcmc_settings, hyper_type, verbose)
      # marca posterior
      msg("[runner] returned ok. Out slots: %s", paste(names(ans), collapse = ", "))
      if (!is.null(ans$data$X) && !is.null(ans$data$y)) {
        msg("[runner] dim(X)=%dx%d, length(y)=%d",
            nrow(ans$data$X), ncol(ans$data$X), length(ans$data$y))
      }
      ans
    },
    error = function(e) {
      # NO silenciar: informamos todo lo que podamos
      call_txt <- tryCatch(paste(deparse(conditionCall(e)), collapse=" "), error = function(...) "NA")
      msg("[ERROR] simulation %d -> %s", i, conditionMessage(e))
      msg("[ERROR] call: %s", call_txt)
      
      # devolvemos un objeto “error” con más contexto
      structure(list(
        error = conditionMessage(e),
        call  = call_txt,
        spec  = sim_specs[[i]]
      ), class = "sim_error")
    },
    warning = function(w) {
      msg("[WARN] simulation %d -> %s", i, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  }
  
  out <- list(
    simulations = res_list,
    settings = list(vi = vi_settings, mcmc = mcmc_settings, hyper_type = hyper_type)
  )
  
  if (!is.null(save_rds_path)) {
    saveRDS(out, file = save_rds_path)
    message(sprintf("Saved simulations to: %s", save_rds_path))
  }
  out
}

# ============================================================
# Runner de simulaciones HLRM: MCMC + CAVI (VI) + SVI
# Requiere: sim_hlrm(), hyperparams_hlrm(), vi_hlrm, mcmc_hlrm
# ============================================================

# -- corre UNA simulación HLRM
run_one_sim_hlrm <- function(sim,
                             vi_settings,
                             svi_settings,
                             mcmc_settings,
                             hyper_type = "OLS",
                             verbose = TRUE) {
  # sim: list con (mínimo) m, nj, Beta_j (m×p), sigma2_j (m)
  # opcionales: p, dummies, seed_data, seed_vi, seed_svi, seed_mcmc, standarize
  m            <- as.integer(sim$m)
  nj           <- sim$nj                    # escalar o vector de longitud m
  Beta_j       <- as.matrix(sim$Beta_j)     # m×p
  sigma2_j     <- as.numeric(sim$sigma2_j)  # m
  p            <- if (!is.null(sim$p)) as.integer(sim$p) else ncol(Beta_j)
  dummies      <- if (!is.null(sim$dummies)) as.integer(sim$dummies) else 0L
  seed_data    <- sim$seed_data
  seed_vi      <- sim$seed_vi
  seed_svi     <- sim$seed_svi
  seed_mcmc    <- sim$seed_mcmc
  standarize   <- if (!is.null(sim$standarize)) isTRUE(sim$standarize) else FALSE
  
  # ---- Simular datos ----
  if (!is.null(seed_data)) set.seed(seed_data)
  data <- sim_hlrm(
    m        = m,
    nj       = nj,
    Beta_j   = Beta_j,
    p        = p,
    dummies  = dummies,
    sigma2_j = sigma2_j,
    seed     = seed_data
  )
  X <- data$X; y <- data$y; g <- data$group
  
  # ---- Hiperparámetros ----
  hypers <- hyperparams_hlrm(list(X = X, y = y), type = hyper_type)
  
  # ================== VI (CAVI) ==================
  if (!is.null(seed_vi)) set.seed(seed_vi)
  vi <- vi_hlrm$new(y, X, standarize = standarize, hypers, vi_settings)
  results_vi     <- vi$fit(verbose = verbose)
  hat_values_vi  <- vi$inferences()
  y_hat_vi       <- vi$predict(X)
  goodness_vi    <- vi$goodness_of_fit(seed = 1)
  
  # ================== SVI (stochastic VI) ==================
  if (!is.null(seed_svi)) set.seed(seed_svi)
  svi <- vi_hlrm$new(y, X, standarize = standarize, hypers, svi_settings)
  results_svi     <- svi$fit(verbose = verbose)
  hat_values_svi  <- svi$inferences()
  y_hat_svi       <- svi$predict(X)
  goodness_svi    <- svi$goodness_of_fit(seed = 1)
  
  # ================== MCMC ==================
  if (!is.null(seed_mcmc)) set.seed(seed_mcmc)
  # asegurar que mcmc_settings tenga samples_gof
  if (is.null(mcmc_settings$samples_gof)) mcmc_settings$samples_gof <- 1000
  mcmc <- mcmc_hlrm$new(y, X, hypers, mcmc_settings)
  results_mcmc    <- mcmc$fit(verbose = verbose, seed = seed_mcmc)
  hat_values_mcmc <- mcmc$inferences()
  y_hat_mcmc      <- mcmc$predict(X)
  goodness_mcmc   <- mcmc$goodness_of_fit()
  
  # ---- Empaquetar ----
  list(
    data = list(
      X = X, y = y, group = g,
      Beta_j_true = Beta_j, sigma2_j_true = sigma2_j,
      m = m, nj = nj, p = p, dummies = dummies,
      seeds = list(data = seed_data, vi = seed_vi, svi = seed_svi, mcmc = seed_mcmc)
    ),
    vi = list(
      results = results_vi, hat_values = hat_values_vi,
      gof = goodness_vi, y_hat = y_hat_vi, timing = results_vi$timing
    ),
    svi = list(
      results = results_svi, hat_values = hat_values_svi,
      gof = goodness_svi, y_hat = y_hat_svi, timing = results_svi$timing
    ),
    mcmc = list(
      results = results_mcmc, hat_values = hat_values_mcmc,
      gof = goodness_mcmc, y_hat = y_hat_mcmc, timing = results_mcmc$timing
    )
  )
}

# -- corre VARIAS simulaciones HLRM (con overrides por escenario)
run_simulations_hlrm <- function(sim_specs,
                                 vi_settings   = list(iterations = 20000, samples_gof = 1000),
                                 svi_settings  = list(iterations = 20000, epsilon = 1e-5,
                                                      samples_gof = 1000, stochastic = TRUE,
                                                      minibatch_m = 6, tau = 10, kappa = 1),
                                 mcmc_settings = list(samples = 20000, burn_in = 5000, thin = 10,
                                                      samples_gof = 1000),
                                 hyper_type    = "OLS",
                                 verbose = TRUE,
                                 save_rds_path = NULL) {
  stopifnot(is.list(sim_specs), length(sim_specs) >= 1L)
  res_list <- vector("list", length(sim_specs))
  names(res_list) <- if (!is.null(names(sim_specs))) names(sim_specs) else paste0("sim_", seq_along(sim_specs))
  
  for (i in seq_along(sim_specs)) {
    if (verbose) cat(sprintf("\n--- Running HLRM simulation %d/%d ---\n", i, length(sim_specs)))
    sim <- sim_specs[[i]]
    
    vi_set_i   <- if (!is.null(sim$vi_settings))   sim$vi_settings   else vi_settings
    svi_set_i  <- if (!is.null(sim$svi_settings))  sim$svi_settings  else svi_settings
    mcmc_set_i <- if (!is.null(sim$mcmc_settings)) sim$mcmc_settings else mcmc_settings
    
    res_list[[i]] <- tryCatch(
      run_one_sim_hlrm(sim, vi_set_i, svi_set_i, mcmc_set_i, hyper_type, verbose),
      error = function(e) {
        warning(sprintf("HLRM simulation %d failed: %s", i, e$message))
        list(error = e$message)
      }
    )
  }
  
  out <- list(
    model = "HLRM",
    simulations = res_list,
    settings = list(vi = vi_settings, svi = svi_settings, mcmc = mcmc_settings, hyper_type = hyper_type)
  )
  if (!is.null(save_rds_path)) {
    saveRDS(out, file = save_rds_path)
    message(sprintf("Saved HLRM simulations to: %s", save_rds_path))
  }
  out
}


# ============================================================
# Runner de simulaciones CHLRM: MCMC + CAVI (VI) + SVI
# Requiere: sim_chlrm(), hyperparams_chlrm(), vi_chlrm, mcmc_chlrm
# ============================================================

# -- corre UNA simulación CHLRM
run_one_sim_chlrm <- function(sim,
                              vi_settings,
                              svi_settings,
                              mcmc_settings,
                              hyper_type = "OLS",
                              verbose = TRUE) {
  # sim: list con (mínimo) m, nj, K, Betak (K×p), sigma2k (K)
  # opcionales: z, omega, x_lim, p, seed_data, seed_vi, seed_svi, seed_mcmc
  m         <- as.integer(sim$m)
  nj        <- sim$nj
  K         <- as.integer(sim$K)
  Betak     <- as.matrix(sim$Betak)        # K×p
  sigma2k   <- as.numeric(sim$sigma2k)     # K
  p         <- if (!is.null(sim$p)) as.integer(sim$p) else ncol(Betak)
  z         <- sim$z
  omega     <- sim$omega
  x_lim     <- sim$x_lim
  seed_data <- sim$seed_data
  seed_vi   <- sim$seed_vi
  seed_svi  <- sim$seed_svi
  seed_mcmc <- sim$seed_mcmc
  
  # ---- Simular datos ----
  if (!is.null(seed_data)) set.seed(seed_data)
  data <- sim_chlrm(
    m = m, nj = nj,
    Betak = Betak, sigma2k = sigma2k, p = p,
    z = z, omega = omega, K = K,
    x_lim = if (!is.null(x_lim)) x_lim else c(-3, 3),
    seed = seed_data
  )
  X <- data$X; y <- data$y; g <- data$group
  
  # ---- Hiperparámetros ----
  hypers <- hyperparams_chlrm(list(X = X, y = y), type = hyper_type)
  
  # ================== VI (CAVI) ==================
  if (!is.null(seed_vi)) set.seed(seed_vi)
  vi <- vi_chlrm$new(y, X, hypers, vi_settings)
  results_vi     <- vi$fit(verbose = verbose)
  hat_values_vi  <- vi$inferences()
  y_hat_vi       <- vi$predict(X)
  goodness_vi    <- vi$goodness_of_fit(seed = 1)
  
  # ================== SVI (stochastic VI) ==================
  if (!is.null(seed_svi)) set.seed(seed_svi)
  svi <- vi_chlrm$new(y, X, hypers, svi_settings)  # requiere stochastic=TRUE en settings
  results_svi     <- svi$fit(verbose = verbose)
  hat_values_svi  <- svi$inferences()
  y_hat_svi       <- svi$predict(X)
  goodness_svi    <- svi$goodness_of_fit(seed = 1)
  
  # ================== MCMC ==================
  if (!is.null(seed_mcmc)) set.seed(seed_mcmc)
  if (is.null(mcmc_settings$samples_gof)) mcmc_settings$samples_gof <- 1000
  mcmc <- mcmc_chlrm$new(y, X, hypers, mcmc_settings)
  results_mcmc    <- mcmc$fit(verbose = verbose, seed = seed_mcmc)
  hat_values_mcmc <- mcmc$inferences()
  y_hat_mcmc      <- mcmc$predict(X)
  goodness_mcmc   <- mcmc$goodness_of_fit()
  
  # ---- Empaquetar ----
  list(
    data = list(
      X = X, y = y, group = g,
      Betak_true = Betak, sigma2k_true = sigma2k, K = K,
      m = m, nj = nj, p = p,
      z_true = data$cluster_info$z, omega_true = data$cluster_info$omega,
      seeds = list(data = seed_data, vi = seed_vi, svi = seed_svi, mcmc = seed_mcmc)
    ),
    vi = list(
      results = results_vi, hat_values = hat_values_vi,
      gof = goodness_vi, y_hat = y_hat_vi, timing = results_vi$timing
    ),
    svi = list(
      results = results_svi, hat_values = hat_values_svi,
      gof = goodness_svi, y_hat = y_hat_svi, timing = results_svi$timing
    ),
    mcmc = list(
      results = results_mcmc, hat_values = hat_values_mcmc,
      gof = goodness_mcmc, y_hat = y_hat_mcmc, timing = results_mcmc$timing
    )
  )
}

# -- corre VARIAS simulaciones CHLRM
run_simulations_chlrm <- function(sim_specs,
                                  vi_settings   = list(iterations = 20000, samples_gof = 1000),
                                  svi_settings  = list(iterations = 20000, epsilon = 1e-5,
                                                       samples_gof = 1000, stochastic = TRUE,
                                                       minibatch_m = 6, tau = 10, kappa = 1),
                                  mcmc_settings = list(samples = 20000, burn_in = 5000, thin = 10,
                                                       samples_gof = 1000),
                                  hyper_type    = "OLS",
                                  verbose = TRUE,
                                  save_rds_path = NULL) {
  stopifnot(is.list(sim_specs), length(sim_specs) >= 1L)
  res_list <- vector("list", length(sim_specs))
  names(res_list) <- if (!is.null(names(sim_specs))) names(sim_specs) else paste0("sim_", seq_along(sim_specs))
  
  for (i in seq_along(sim_specs)) {
    if (verbose) cat(sprintf("\n--- Running CHLRM simulation %d/%d ---\n", i, length(sim_specs)))
    sim <- sim_specs[[i]]
    
    vi_set_i   <- if (!is.null(sim$vi_settings))   sim$vi_settings   else vi_settings
    svi_set_i  <- if (!is.null(sim$svi_settings))  sim$svi_settings  else svi_settings
    mcmc_set_i <- if (!is.null(sim$mcmc_settings)) sim$mcmc_settings else mcmc_settings
    
    res_list[[i]] <- tryCatch(
      run_one_sim_chlrm(sim, vi_set_i, svi_set_i, mcmc_set_i, hyper_type, verbose),
      error = function(e) {
        warning(sprintf("CHLRM simulation %d failed: %s", i, e$message))
        list(error = e$message)
      }
    )
  }
  
  out <- list(
    model = "CHLRM",
    simulations = res_list,
    settings = list(vi = vi_settings, svi = svi_settings, mcmc = mcmc_settings, hyper_type = hyper_type)
  )
  if (!is.null(save_rds_path)) {
    saveRDS(out, file = save_rds_path)
    message(sprintf("Saved CHLRM simulations to: %s", save_rds_path))
  }
  out
}





# ---------- helpers ----------
.safe_chol <- function(S) {
  out <- try(chol(S), silent = TRUE)
  if (inherits(out, "try-error")) {
    E <- eigen((S + t(S))/2, symmetric = TRUE)
    vals <- pmax(E$values, 1e-12)
    S2 <- E$vectors %*% diag(vals) %*% t(E$vectors)
    chol(S2)
  } else out
}
ellipse_points <- function(mu, Sigma, prob = 0.9, n = 400) {
  stopifnot(length(mu) == 2, all(dim(Sigma) == c(2, 2)))
  q <- qchisq(prob, df = 2)
  t <- seq(0, 2*pi, length.out = n)
  L <- t(.safe_chol(Sigma))         # Σ = L L^T
  U <- rbind(cos(t), sin(t))
  E <- mu + sqrt(q) * (L %*% U)
  cbind(x = E[1, ], y = E[2, ])
}
.subcov <- function(S, j1, j2) {
  if (is.null(S)) return(NULL)
  if (!is.matrix(S)) stop("Covariance must be a matrix.")
  S[c(j1, j2), c(j1, j2), drop = FALSE]
}
.pct_str <- function(pv) {
  # formatea vector de coberturas en "%", p.ej. c(0.5,0.95) -> "50%, 95%"
  paste0(format(round(100*pv, ifelse(any((pv*100) %% 1 != 0), 1, 0)), trim = TRUE),
         "%", collapse = ", ")
}

# ---------- función principal ----------
plot_betas_joint <- function(results_mcmc, hat_values_mcmc, hat_values_vi,
                             j1 = 1, j2 = 2,
                             probs_ellipses = c(0.5, 0.95),   # Joint CI levels
                             ci_label_level = 95,              # texto de CI univariados en leyenda
                             add_kde = FALSE, kde_n = 200,
                             legend_pos = "topleft",
                             # knobs de leyenda
                             legend_cex = 0.75, legend_pt_cex = 0.8,
                             legend_seg_len = 1.4, legend_inset = 0.02,
                             legend_ncol = 1, legend_box_lwd = 0,
                             legend_box_col = NA, legend_bg = NA) {
  
  # scatter para contexto
  x <- results_mcmc$chains$Beta[, j1]
  y <- results_mcmc$chains$Beta[, j2]
  
  # medias e intervalos (usamos los que traes en hat_values)
  m_mu <- c(hat_values_mcmc$Beta$mean[j1], hat_values_mcmc$Beta$mean[j2])
  m_li <- c(hat_values_mcmc$Beta$CI_LI[j1], hat_values_mcmc$Beta$CI_LI[j2])
  m_ls <- c(hat_values_mcmc$Beta$CI_LS[j1], hat_values_mcmc$Beta$CI_LS[j2])
  v_mu <- c(hat_values_vi$Beta$mean[j1],   hat_values_vi$Beta$mean[j2])
  v_li <- c(hat_values_vi$Beta$CI_LI[j1],  hat_values_vi$Beta$CI_LI[j2])
  v_ls <- c(hat_values_vi$Beta$CI_LS[j1],  hat_values_vi$Beta$CI_LS[j2])
  
  # sd por si hay que construir diag
  m_sd <- c(hat_values_mcmc$Beta$sd[j1], hat_values_mcmc$Beta$sd[j2])
  v_sd <- c(hat_values_vi$Beta$sd[j1],   hat_values_vi$Beta$sd[j2])
  
  # covarianzas: tomamos SIEMPRE de hat_values (si no hay, asumimos independencia)
  S_mcmc <- .subcov(hat_values_mcmc$Beta$cov, j1, j2)
  if (is.null(S_mcmc)) S_mcmc <- diag(pmax(m_sd, .Machine$double.eps)^2)
  
  S_vi <- .subcov(hat_values_vi$Beta$cov, j1, j2)
  if (is.null(S_vi)) S_vi <- diag(pmax(v_sd, .Machine$double.eps)^2)
  
  # límites: scatter + mu ± 4*sd (ambos métodos)
  xr <- range(c(x, m_mu[1] + c(-4, 4)*m_sd[1], v_mu[1] + c(-4, 4)*v_sd[1]), na.rm = TRUE)
  yr <- range(c(y, m_mu[2] + c(-4, 4)*m_sd[2], v_mu[2] + c(-4, 4)*v_sd[2]), na.rm = TRUE)
  xr <- xr + c(-0.08, 0.08)*diff(xr); yr <- yr + c(-0.08, 0.08)*diff(yr)
  
  # base
  plot(x, y, pch = 4, col = "darkgray",
       xlab = bquote(beta[.(j1)]), ylab = bquote(beta[.(j2)]),
       xlim = xr, ylim = yr)
  
  # opcional: contornos kernel de MCMC
  if (add_kde && requireNamespace("MASS", quietly = TRUE)) {
    kd <- MASS::kde2d(x, y, n = kde_n, lims = c(xr[1], xr[2], yr[1], yr[2]))
    contour(kd, add = TRUE, drawlabels = FALSE,
            col = adjustcolor("red", 0.45), lwd = 1.2, lty = 1, nlevels = 6)
  }
  
  # elipses conjuntas (mismos niveles para ambos)
  for (p in probs_ellipses) {
    em <- ellipse_points(m_mu, S_mcmc, prob = p); lines(em[,1], em[,2], col = "red",   lwd = 2, lty = 1)
    ev <- ellipse_points(v_mu, S_vi,   prob = p); lines(ev[,1], ev[,2], col = "black", lwd = 2, lty = 2)
  }
  
  # medias e CI univariados (95% por etiqueta)
  abline(v = m_mu[1], col = "red",   lwd = 2,   lty = 1)
  abline(h = m_mu[2], col = "red",   lwd = 2,   lty = 1)
  abline(v = c(m_li[1], m_ls[1]), col = "red",   lwd = 1.2, lty = 3)
  abline(h = c(m_li[2], m_ls[2]), col = "red",   lwd = 1.2, lty = 3)
  
  abline(v = v_mu[1], col = "black", lwd = 2,   lty = 2)
  abline(h = v_mu[2], col = "black", lwd = 2,   lty = 2)
  abline(v = c(v_li[1], v_ls[1]), col = "black", lwd = 1.2, lty = 3)
  abline(h = c(v_li[2], v_ls[2]), col = "black", lwd = 1.2, lty = 3)
  
  # ------ legend (English) ------
  joint_txt <- .pct_str(probs_ellipses)  # e.g., "50%, 95%"
  legend(legend_pos,
         legend = c(
           sprintf("MCMC joint CR: %s", joint_txt),
           sprintf("VI joint CR: %s",   joint_txt),
           "MCMC mean", sprintf("MCMC margin CI %d%%", ci_label_level),
           "VI mean",   sprintf("VI margin CI %d%%",   ci_label_level)
         ),
         col    = c('red', 'black', "red", "red", "black", "black"),
         lty    = c(1, 1, 1, 3, 2, 3),
         lwd    = c(2, 2, 2, 1.2, 2, 1.2),
         cex    = legend_cex, pt.cex = legend_pt_cex,
         seg.len= legend_seg_len, inset = legend_inset, ncol = legend_ncol,
         box.lwd= legend_box_lwd, box.col = legend_box_col, bg = legend_bg,
         bty    = if (is.na(legend_box_col)) "n" else "o")
}

# ---------- helpers ----------
dinv_gamma <- function(x, alpha, beta) {
  # IG(alpha, beta[scale])  f(x) = beta^alpha / Gamma(alpha) * x^(-alpha-1) * exp(-beta/x), x>0
  out <- rep(NA_real_, length(x))
  pos <- x > 0 & is.finite(x)
  if (any(pos)) {
    logf <- alpha*log(beta) - lgamma(alpha) + (-alpha-1)*log(x[pos]) - beta/x[pos]
    out[pos] <- exp(logf)
  }
  out[!is.finite(out)] <- NA_real_
  out
}
ig_from_mean_sd <- function(mu, sd) {
  # Para IG con parametrización de "scale" beta:
  # E[X] = beta/(alpha-1), Var[X] = beta^2 / ((alpha-1)^2 (alpha-2))
  v <- sd^2
  if (!is.finite(mu) || !is.finite(sd) || sd <= 0) return(NULL)
  alpha <- 2 + (mu^2 / v)
  beta  <- mu * (alpha - 1)
  list(alpha = alpha, beta = beta)
}

# ---------- main ----------
plot_sigma2_density_ig <- function(results_mcmc, hat_values_mcmc, hat_values_vi,
                                   adjust = 2, L_vi = 4,        # L_vi para fijar rango x con mu±L*sd
                                   ci_label_level = 95,
                                   main_title = expression(paste("Posterior for ", sigma^2)),
                                   legend_pos = "topright",
                                   # legend knobs
                                   legend_cex = 0.75, legend_pt_cex = 0.8,
                                   legend_seg_len = 1.2, legend_inset = 0.02,
                                   legend_ncol = 1, legend_box_lwd = 0,
                                   legend_box_col = NA, legend_bg = NA) {
  
  # --- MCMC samples ---
  s <- results_mcmc$chains$sigma2
  stopifnot(is.numeric(s))
  d_mcmc <- density(s, adjust = adjust, na.rm = TRUE)
  
  # --- hat-values ---
  mu_m <- hat_values_mcmc$sigma2$mean
  li_m <- hat_values_mcmc$sigma2$CI_LI
  ls_m <- hat_values_mcmc$sigma2$CI_LS
  
  mu_v <- hat_values_vi$sigma2$mean
  sd_v <- max(hat_values_vi$sigma2$sd, .Machine$double.eps)
  li_v <- hat_values_vi$sigma2$CI_LI
  ls_v <- hat_values_vi$sigma2$CI_LS
  
  # Parámetros IG para VI
  prm <- ig_from_mean_sd(mu_v, sd_v)
  
  # --- rango x positivo: unión densidad MCMC y mu_v ± L_vi*sd_v ---
  xr <- range(d_mcmc$x[d_mcmc$x > 0], mu_v + c(-L_vi, L_vi) * sd_v, na.rm = TRUE)
  xr[1] <- max(xr[1], 1e-8)
  xg <- seq(xr[1], xr[2], length.out = 512)
  
  # densidad VI (IG)
  y_vi <- if (!is.null(prm)) dinv_gamma(xg, prm$alpha, prm$beta) else rep(NA_real_, length(xg))
  
  # ylim desde 0
  ymax <- max(d_mcmc$y, y_vi, na.rm = TRUE) * 1.1
  if (!is.finite(ymax) || ymax <= 0) ymax <- 1
  
  # --- plot ---
  plot(NA, xlim = xr, ylim = c(0, ymax), yaxs = "i",
       xlab = expression(sigma^2), ylab = "Density", main = main_title)
  
  # curvas
  lines(d_mcmc,          col = "red",   lwd = 2, lty = 1)  # MCMC
  if (any(is.finite(y_vi))) lines(xg, y_vi, col = "black", lwd = 2, lty = 2)  # VI (IG)
  
  # líneas: medias y CI univariados
  abline(v = mu_m,          col = "red",   lwd = 2,   lty = 1)
  abline(v = c(li_m, ls_m), col = "red",   lwd = 1.2, lty = 3)
  
  abline(v = mu_v,          col = "black", lwd = 2,   lty = 2)
  abline(v = c(li_v, ls_v), col = "black", lwd = 1.2, lty = 3)
  
  # leyenda
  legend(legend_pos,
         legend = c("MCMC",
                    "VI",
                    "MCMC mean", sprintf("MCMC CI %d%%", ci_label_level),
                    "VI mean",   sprintf("VI CI %d%%",   ci_label_level)),
         col    = c("red", "black", "red", "red", "black", "black"),
         lty    = c(1, 2, 1, 3, 2, 3),
         lwd    = c(2, 2, 2, 1.2, 2, 1.2),
         cex    = legend_cex, pt.cex = legend_pt_cex,
         seg.len= legend_seg_len, inset = legend_inset, ncol = legend_ncol,
         box.lwd= legend_box_lwd, box.col = legend_box_col, bg = legend_bg,
         bty    = if (is.na(legend_box_col)) "n" else "o")
}

plot_vi_vs_mcmc_time <- function(
    results_vi, results_mcmc,
    main_title = "CAVI ELBO vs MCMC mean log-joint",
    xlab = "t (s)",
    ylab = expression(E[q[Theta]]*"["*log*phantom(.)*p(y,X,Theta)*"]"),
    # estilos VI (CAVI)
    col_vi = "red3", pch_vi = 16, lwd_vi = 2, lty_vi = 1,
    # línea horizontal MCMC
    col_mcmc = "#1f78b4", lty_mcmc = 2, lwd_mcmc = 2,
    # suavizado opcional del ELBO (media móvil centrada)
    smooth_ma = NULL, col_smooth = "red3", lwd_smooth = 2.5,
    # grid y leyenda
    grid_col = "gray90",
    legend_pos = "bottomright",
    legend_cex = 0.9, legend_pt_cex = 1.0,
    legend_seg_len = 1.4, legend_inset = 0.02,
    legend_ncol = 1, legend_box_lwd = 0,
    legend_box_col = NA, legend_bg = NA
){
  # --- datos: tiempos y ELBO de VI ---
  t_vi   <- as.numeric(results_vi$timing$up_to_iter)
  t_vi <- c(0,t_vi)
  elbo   <- as.numeric(results_vi$params$elbo)
  elbo <- c(-300,elbo)
  
  stopifnot(length(t_vi) == length(elbo))
  
  # --- referencia: media de log-joint de MCMC ---
  lj_mean <- mean(results_mcmc$chains$logjoint, na.rm = TRUE)
  
  # --- suavizado (opcional) ---
  elbo_smooth <- NULL
  if (!is.null(smooth_ma) && is.finite(smooth_ma) && smooth_ma > 1) {
    w <- rep(1/smooth_ma, smooth_ma)
    elbo_smooth <- as.numeric(stats::filter(elbo, filter = w, sides = 2))
  }
  
  # --- límites del plot ---
  t_all <- t_vi
  y_all <- c(elbo, lj_mean, elbo_smooth)
  xlim  <- range(t_all, na.rm = TRUE)
  ylim  <- range(y_all, na.rm = TRUE)
  pad   <- 0.04 * diff(ylim); if (!is.finite(pad)) pad <- 1
  ylim  <- c(ylim[1] - pad, ylim[2] + pad)
  
  # --- lienzo ---
  plot(t_vi, elbo, type = "n",
       xlab = xlab, ylab = ylab, main = main_title,
       xlim = xlim, ylim = ylim,
       panel.first = grid(col = grid_col))
  
  # CAVI 
  lines(t_vi, elbo,  col = col_vi, lwd = lwd_vi, lty = lty_vi) 
  
  # Suavizado (si se pidió)
  if (!is.null(elbo_smooth)) {
    lines(t_vi, elbo_smooth, col = col_smooth, lwd = lwd_smooth, lty = 1)
  }
  
  # Media log-joint MCMC (línea horizontal)
  abline(h = lj_mean, col = col_mcmc, lty = lty_mcmc, lwd = lwd_mcmc)
  
  # --- leyenda ---
  leg <- c("VI ELBO", "MCMC posterior mean log-joint")
  col <- c(col_vi, col_mcmc)
  lty <- c(lty_vi, lty_mcmc)
  lwd <- c(lwd_vi, lwd_mcmc)
  pch <- c(pch_vi, NA)
  if (!is.null(elbo_smooth)) {
    leg <- c(leg, sprintf("CAVI %d-pt moving average", smooth_ma))
    col <- c(col, col_smooth)
    lty <- c(lty, 1)
    lwd <- c(lwd, lwd_smooth)
    pch <- c(pch, NA)
  }
  
  legend(legend_pos, legend = leg,
         col = col, lty = lty, lwd = lwd, pch = pch,
         cex = legend_cex, pt.cex = legend_pt_cex,
         seg.len = legend_seg_len, inset = legend_inset, ncol = legend_ncol,
         box.lwd = legend_box_lwd, box.col = legend_box_col, bg = legend_bg,
         bty = if (is.na(legend_box_col)) "n" else "o")
}



# ================== HELPERS (vectorizados) ==================
.format_num_one <- function(x, digits) {
  if (is.null(x) || length(x) == 0) return("--")
  x <- as.numeric(x)[1L]
  if (!is.finite(x)) return("--")
  formatC(x, digits = digits, format = "f")
}
.format_num <- function(x, digits) {
  # vectoriza .format_num_one
  vapply(x, .format_num_one, character(1), digits = digits)
}
.ci_str_one <- function(li, ls, digits) {
  paste0("$[", .format_num_one(li, digits), ", ", .format_num_one(ls, digits), "]$")
}
.ci_str <- function(li, ls, digits) {
  mapply(.ci_str_one, li, ls, MoreArgs = list(digits = digits), USE.NAMES = FALSE)
}
.bool_str <- function(b) ifelse(isTRUE(b), "Yes", "No")
.bool_str_vec <- function(b) if (length(b) == 1L) .bool_str(b) else ifelse(b, "Yes", "No")

# Construcción de filas para Beta y sigma^2
.build_rows_from_hat <- function(hat_values, digits = 3) {
  # --- Beta ---
  mu  <- as.numeric(hat_values$Beta$mean)
  sd  <- as.numeric(hat_values$Beta$sd)
  li  <- as.numeric(hat_values$Beta$CI_LI)
  ls  <- as.numeric(hat_values$Beta$CI_LS)
  cz  <- as.logical(hat_values$Beta$contains_zero)
  
  p <- length(mu)
  par_names <- paste0("$\\beta_{", seq_len(p), "}$")
  
  df_beta <- data.frame(
    Parameter   = par_names,
    Mean        = .format_num(mu, digits),
    SD          = .format_num(sd, digits),
    CI          = .ci_str(li, ls, digits),
    `0\\in CI`  = .bool_str_vec(cz),
    check.names = FALSE
  )
  
  # --- sigma^2 ---
  mu_s <- hat_values$sigma2$mean
  sd_s <- hat_values$sigma2$sd
  li_s <- hat_values$sigma2$CI_LI
  ls_s <- hat_values$sigma2$CI_LS
  cz_s <- hat_values$sigma2$contains_zero
  
  df_sig2 <- data.frame(
    Parameter   = "$\\sigma^2$",
    Mean        = .format_num_one(mu_s, digits),
    SD          = .format_num_one(sd_s, digits),
    CI          = .ci_str_one(li_s, ls_s, digits),
    `0\\in CI`  = .bool_str(cz_s),
    check.names = FALSE
  )
  
  rbind(df_beta, df_sig2)
}

# ================== TABLA (un método) ==================
latex_param_table_single <- function(hat_values,
                                     method_name = "MCMC",
                                     digits = 3,
                                     caption = NULL,
                                     label   = NULL,
                                     use_booktabs = FALSE) {
  df <- .build_rows_from_hat(hat_values, digits)
  
  top <- if (use_booktabs) "\\toprule" else "\\hline"
  mid <- if (use_booktabs) "\\midrule" else "\\hline"
  bot <- if (use_booktabs) "\\bottomrule" else "\\hline"
  if (is.null(caption)) caption <- sprintf("Posterior summaries (%s)", method_name)
  if (is.null(label))   label   <- sprintf("tab:post_%s", tolower(method_name))
  
  lines <- c(
    "\\begin{table}[ht]",
    "\\centering",
    paste0("\\caption{", caption, "}"),
    paste0("\\label{", label, "}"),
    "\\begin{tabular}{lrrrr}",
    top,
    "Parameter & Mean & SD & CI & $0\\in$CI \\\\",
    mid
  )
  
  for (i in seq_len(nrow(df))) {
    lines <- c(lines, sprintf("%s & %s & %s & %s & %s \\\\",
                              df$Parameter[i], df$Mean[i], df$SD[i], df$CI[i], df$`0\\in CI`[i]))
  }
  
  lines <- c(lines, bot, "\\end{tabular}", "\\end{table}")
  paste(lines, collapse = "\n")
}

# ================== TABLA (MCMC vs VI) ==================
latex_param_table_both <- function(hat_values_mcmc, hat_values_vi,
                                   digits = 3,
                                   caption = "Posterior summaries (MCMC vs VI)",
                                   label   = "tab:post_mcmc_vi",
                                   use_booktabs = FALSE) {
  df_m <- .build_rows_from_hat(hat_values_mcmc, digits)
  df_v <- .build_rows_from_hat(hat_values_vi,   digits)
  
  stopifnot(identical(df_m$Parameter, df_v$Parameter))
  df <- data.frame(
    Parameter         = df_m$Parameter,
    `Mean (MCMC)`     = df_m$Mean, `SD (MCMC)`     = df_m$SD,
    `CI (MCMC)`       = df_m$CI,   `0\\in CI (MCMC)` = df_m$`0\\in CI`,
    `Mean (VI)`       = df_v$Mean, `SD (VI)`       = df_v$SD,
    `CI (VI)`         = df_v$CI,   `0\\in CI (VI)`   = df_v$`0\\in CI`,
    check.names = FALSE
  )
  
  top <- if (use_booktabs) "\\toprule" else "\\hline"
  mid <- if (use_booktabs) "\\midrule" else "\\hline"
  bot <- if (use_booktabs) "\\bottomrule" else "\\hline"
  
  lines <- c(
    "\\begin{table}[ht]",
    "\\centering",
    paste0("\\caption{", caption, "}"),
    paste0("\\label{", label, "}"),
    "\\begin{tabular}{lrrrrrrrr}",
    top,
    "Parameter & Mean (MCMC) & SD (MCMC) & CI (MCMC) & $0\\in$CI (MCMC) & Mean (VI) & SD (VI) & CI (VI) & $0\\in$CI (VI) \\\\",
    mid
  )
  
  for (i in seq_len(nrow(df))) {
    lines <- c(lines, sprintf("%s & %s & %s & %s & %s & %s & %s & %s & %s \\\\",
                              df$Parameter[i],
                              df$`Mean (MCMC)`[i], df$`SD (MCMC)`[i],
                              df$`CI (MCMC)`[i],   df$`0\\in CI (MCMC)`[i],
                              df$`Mean (VI)`[i],   df$`SD (VI)`[i],
                              df$`CI (VI)`[i],     df$`0\\in CI (VI)`[i]))
  }
  
  lines <- c(lines, bot, "\\end{tabular}", "\\end{table}")
  paste(lines, collapse = "\n")
}

# Devuelve un string LaTeX con WAIC, DIC, R^2 y MSE para MCMC vs VI
latex_fit_table <- function(goodness_mcmc, goodness_vi,results_mcmc,results_vi,
                            digits = 3,
                            caption = "Model fit: MCMC vs CAVI",
                            label   = "tab:fit_mcmc_vi",
                            use_booktabs = FALSE) {
  # formateo seguro
  fmt <- function(x) {
    if (is.null(x) || length(x) == 0 || !is.finite(x)) return("--")
    formatC(as.numeric(x), digits = digits, format = "f")
  }
  
  # extraer métricas
  mcmc_WAIC <- fmt(goodness_mcmc$waic$WAIC)
  mcmc_time <- fmt(results_mcmc$timing$total)
  mcmc_DIC  <- fmt(goodness_mcmc$dic$DIC)
  mcmc_R2   <- fmt(goodness_mcmc$r2)
  mcmc_MSE  <- fmt(goodness_mcmc$mse)
  
  vi_WAIC <- fmt(goodness_vi$waic$WAIC)
  vi_time <- fmt(results_vi$timing$total)
  vi_DIC  <- fmt(goodness_vi$dic$DIC)
  vi_R2   <- fmt(goodness_vi$r2)
  vi_MSE  <- fmt(goodness_vi$mse)
  
  # separadores según estilo
  top <- if (use_booktabs) "\\toprule" else "\\hline"
  mid <- if (use_booktabs) "\\midrule" else "\\hline"
  bot <- if (use_booktabs) "\\bottomrule" else "\\hline"
  
  # construir LaTeX
  lines <- c(
    "\\begin{table}[ht]",
    "\\centering",
    paste0("\\caption{", caption, "}"),
    paste0("\\label{", label, "}"),
    "\\begin{tabular}{lrr}",
    top,
    "Metric & MCMC & CAVI \\\\",
    mid,
    paste0("WAIC & ", mcmc_WAIC, " & ", vi_WAIC, " \\\\"),
    paste0("DIC & ",  mcmc_DIC,  " & ", vi_DIC,  " \\\\"),
    paste0("$R^2$ & ", mcmc_R2,  " & ", vi_R2,   " \\\\"),
    paste0("MSE & ",  mcmc_MSE,  " & ", vi_MSE,  " \\\\"),
    paste0("Time (sec) & ",  mcmc_time,  " & ", vi_time,  " \\\\"),
    bot,
    "\\end{tabular}",
    "\\end{table}"
  )
  
  out <- paste(lines, collapse = "\n")
  return(out)
}


plot_y_predictive_densities <- function(
    y,
    y_hat_vi   = NULL,
    y_hat_svi  = NULL,
    y_hat_mcmc = NULL,
    S=NULL,
    adjust = 2,
    breaks = "FD",
    main = "Observed vs Predictive Densities",
    xlab = "y",
    ylab = "Density",
    legend_pos = "topright",
    legend_cex = 0.9,
    # estilos
    hist_col = "gray90",
    hist_border = "white",
    col_obs = "blue",   lty_obs = 4, lwd_obs = 2,
    col_vi  = "green4", lty_vi  = 1, lwd_vi  = 2,
    col_svi = "red3",   lty_svi = 2, lwd_svi = 2,
    col_mcmc= "orange", lty_mcmc= 3, lwd_mcmc=2,
    add_rug = FALSE,
    ylim = NULL,    # si NULL, se calcula automáticamente
    xlim = NULL,    # si NULL, usa el rango combinado
    na_rm = TRUE
){
  as_vec <- function(x) {
    if (is.null(x)) return(NULL)
    if (is.list(x)) return(unlist(x, use.names = FALSE))
    as.numeric(x)
  }
  
  yv   <- as_vec(y)
  vvi  <- as_vec(y_hat_vi)
  vsv  <- as_vec(y_hat_svi)
  vmcm <- as_vec(y_hat_mcmc)
  
  # densidades disponibles
  dens_list <- list()
  if (!is.null(yv))   dens_list$obs   <- density(yv,   adjust = adjust, na.rm = na_rm)
  if (!is.null(vvi))  dens_list$vi    <- density(vvi,  adjust = adjust, na.rm = na_rm)
  if (!is.null(vsv))  dens_list$svi   <- density(vsv,  adjust = adjust, na.rm = na_rm)
  if (!is.null(vmcm)) dens_list$mcmc  <- density(vmcm, adjust = adjust, na.rm = na_rm)
  
  if (length(dens_list) == 0) stop("No hay datos para graficar.")
  
  # límites X/Y
  if (is.null(xlim)) {
    x_all <- c(
      if (!is.null(yv)) yv,
      if (!is.null(vvi)) vvi,
      if (!is.null(vsv)) vsv,
      if (!is.null(vmcm)) vmcm
    )
    xlim <- range(x_all, na.rm = TRUE)
  }
  
  ymax <- max(vapply(dens_list, function(d) max(d$y, na.rm = TRUE), numeric(1)))
  if (is.null(ylim)) ylim <- c(0, ymax * 1.05) else ylim[1] <- 0
  
  # Histograma de y como base (escala de densidad)
  hist(yv, freq = FALSE, breaks = breaks, col = hist_col, border = hist_border,
       main = main, xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim)
  
  # Curvas
  if (!is.null(dens_list$obs))  lines(dens_list$obs,  col = col_obs,  lty = lty_obs,  lwd = lwd_obs)
  if (!is.null(dens_list$vi))   lines(dens_list$vi,   col = col_vi,   lty = lty_vi,   lwd = lwd_vi)
  if (!is.null(dens_list$svi))  lines(dens_list$svi,  col = col_svi,  lty = lty_svi,  lwd = lwd_svi)
  if (!is.null(dens_list$mcmc)) lines(dens_list$mcmc, col = col_mcmc, lty = lty_mcmc, lwd = lwd_mcmc)
  
  # Rug opcional
  if (isTRUE(add_rug)) rug(yv, col = adjustcolor("black", 0.4))
  
  # Leyenda dinámica según lo disponible
  leg  <- c("Observed")
  cols <- c(col_obs)
  ltys <- c(lty_obs)
  lwds <- c(lwd_obs)
  
  if (!is.null(dens_list$vi))   { leg  <- c(leg,  "VI"); cols <- c(cols,  col_vi);  ltys <- c(ltys,  lty_vi);  lwds <- c(lwds,  lwd_vi) }
  if (!is.null(dens_list$svi))  { leg  <- c(leg,  sprintf("SVI  |S|=%d",S));  cols <- c(cols,  col_svi); ltys <- c(ltys,  lty_svi); lwds <- c(lwds,  lwd_svi) }
  if (!is.null(dens_list$mcmc)) { leg  <- c(leg,  "MCMC"); cols <- c(cols,  col_mcmc);ltys <- c(ltys,  lty_mcmc);lwds <- c(lwds,  lwd_mcmc) }
  
  legend(legend_pos, legend = leg, col = cols, lty = ltys, lwd = lwds,
         bty = "n", cex = legend_cex)
  
  invisible(list(densities = dens_list, xlim = xlim, ylim = ylim))
}

latex_methods_overview <- function(all_runs,
                                   caption = "Method comparison across simulations",
                                   label   = "tab:methods_overview",
                                   use_booktabs = TRUE,
                                   digits_time = 2) {
  stopifnot(is.list(all_runs), !is.null(all_runs$simulations))
  sims <- all_runs$simulations
  sim_names <- names(sims)
  if (is.null(sim_names)) sim_names <- paste0("sim_", seq_along(sims))
  
  # util: escapar underscores en LaTeX
  .tex_escape <- function(s) gsub("_", "\\\\_", s, fixed = TRUE)
  
  # formateador "c/total (xx.x%)"
  .fmt_cnt <- function(c, tot) sprintf("%d/%d (%.1f\\%%)", c, tot, 100 * c / max(1, tot))
  
  # extrae conteos por método
  .one_method_counts <- function(sim, method = c("mcmc", "vi")) {
    method <- match.arg(method)
    out <- list(zero_in_ci = NA_integer_, truth_in_ci = NA_integer_, time = NA_real_, ok = FALSE)
    
    # tiempo (si existe)
    tt <- tryCatch(sim[[method]]$timing$total, error = function(e) NA_real_)
    if (!is.null(tt) && is.numeric(tt) && length(tt) == 1) out$time <- tt
    
    # hat values
    hv <- tryCatch(sim[[method]]$hat_values, error = function(e) NULL)
    if (is.null(hv) || is.null(hv$Beta)) return(out)
    
    # p (número de betas) y total de parámetros
    p <- tryCatch(ncol(sim$data$X), error = function(e) NA_integer_)
    if (is.na(p)) return(out)
    total_params <- p + 1L
    
    # 0 in CI: de betas + (opcionalmente sigma2, casi siempre FALSE)
    cz_beta <- hv$Beta$contains_zero
    if (is.null(cz_beta)) {
      # si no vino precomputado, lo calculo
      li <- hv$Beta$CI_LI; ls <- hv$Beta$CI_LS
      if (length(li) == p && length(ls) == p) {
        cz_beta <- (li <= 0 & 0 <= ls)
      } else {
        cz_beta <- rep(NA, p)
      }
    }
    cz_s2 <- tryCatch(isTRUE(hv$sigma2$contains_zero), error = function(e) FALSE)
    zero_cnt <- sum(cz_beta, na.rm = TRUE) + as.integer(cz_s2)
    
    # Truth in CI: compara contra valores verdaderos de la simulación
    btrue <- tryCatch(sim$data$beta_true, error = function(e) NULL)
    s2true <- tryCatch(sim$data$sigma2_true, error = function(e) NULL)
    ti_beta <- NA_integer_
    ti_sigma2 <- NA_integer_
    
    if (!is.null(btrue)) {
      li <- hv$Beta$CI_LI; ls <- hv$Beta$CI_LS
      if (length(li) == p && length(ls) == p && length(btrue) == p) {
        ti_beta <- sum((li <= btrue) & (btrue <= ls), na.rm = TRUE)
      } else ti_beta <- NA_integer_
    }
    
    if (!is.null(s2true) && !is.null(hv$sigma2)) {
      li2 <- hv$sigma2$CI_LI; ls2 <- hv$sigma2$CI_LS
      if (is.finite(li2) && is.finite(ls2) && is.finite(s2true)) {
        ti_sigma2 <- as.integer((li2 <= s2true) & (s2true <= ls2))
      } else ti_sigma2 <- NA_integer_
    }
    
    truth_cnt <- suppressWarnings(as.integer(ti_beta) + as.integer(ti_sigma2))
    out$zero_in_ci <- zero_cnt
    out$truth_in_ci <- truth_cnt
    out$ok <- TRUE
    out
  }
  
  # armar filas
  rows <- lapply(seq_along(sims), function(i) {
    sim <- sims[[i]]
    nm  <- .tex_escape(sim_names[i])
    n   <- tryCatch(sim$data$n, error = function(e) NA_integer_)
    p   <- tryCatch(ncol(sim$data$X), error = function(e) NA_integer_)
    total_params <- if (is.finite(p)) p + 1L else NA_integer_
    
    m_mcmc <- .one_method_counts(sim, "mcmc")
    m_vi   <- .one_method_counts(sim, "vi")
    
    # strings formateados
    z_mcmc <- if (is.finite(total_params) && !is.na(m_mcmc$zero_in_ci))
      .fmt_cnt(m_mcmc$zero_in_ci, total_params) else "---"
    t_mcmc <- if (is.finite(total_params) && !is.na(m_mcmc$truth_in_ci))
      .fmt_cnt(m_mcmc$truth_in_ci, total_params) else "---"
    s_mcmc <- if (is.finite(m_mcmc$time)) sprintf(paste0("%.", digits_time, "f"), m_mcmc$time) else "---"
    
    z_vi <- if (is.finite(total_params) && !is.na(m_vi$zero_in_ci))
      .fmt_cnt(m_vi$zero_in_ci, total_params) else "---"
    t_vi <- if (is.finite(total_params) && !is.na(m_vi$truth_in_ci))
      .fmt_cnt(m_vi$truth_in_ci, total_params) else "---"
    s_vi <- if (is.finite(m_vi$time)) sprintf(paste0("%.", digits_time, "f"), m_vi$time) else "---"
    
    c(nm, n, p, total_params, z_mcmc, t_mcmc, s_mcmc, z_vi, t_vi, s_vi)
  })
  
  M <- do.call(rbind, rows)
  colnames(M) <- c("Scenario","n","p","#params",
                   "MCMC_0inCI","MCMC_TruthInCI","MCMC_Time",
                   "VI_0inCI","VI_TruthInCI","VI_Time")
  
  # construir LaTeX manualmente (más flexible para doble encabezado)
  align <- "l r r r c c c c c c"
  top   <- if (use_booktabs) "\\toprule" else "\\hline"
  mid   <- if (use_booktabs) "\\midrule" else "\\hline"
  bot   <- if (use_booktabs) "\\bottomrule" else "\\hline"
  cmid  <- function(from,to) {
    if (use_booktabs) sprintf("\\cmidrule(lr){%d-%d}", from, to) else ""
  }
  
  # encabezados
  hdr1 <- paste0("& \\multicolumn{3}{c}{\\textbf{MCMC}} & \\multicolumn{3}{c}{\\textbf{VI}} \\\\",
                 "\n", cmid(5,7), cmid(8,10))
  hdr0 <- paste("Scenario & n & p & \\#params & $0\\in CI$ & Truth$\\in$CI & Time (s) & $0\\in CI$ & Truth$\\in$CI & Time (s) \\\\")
  
  # filas
  body <- apply(M, 1, function(r)
    paste(r, collapse = " & ")
  )
  body <- paste0(body, " \\\n", collapse = "")
  
  # tabla completa
  out <- paste0(
    "\\begin{table}[ht]\n\\centering\n",
    "\\caption{", caption, "}\n",
    "\\label{", label, "}\n",
    "\\begin{tabular}{", align, "}\n",
    top, "\n",
    " & & & & \\multicolumn{3}{c}{\\textbf{MCMC}} & \\multicolumn{3}{c}{\\textbf{VI}} \\\\\n",
    cmid(5,7), cmid(8,10), "\n",
    hdr0, "\n",
    mid, "\n",
    body,
    bot, "\n",
    "\\end{tabular}\n\\end{table}\n"
  )
  
  # Nota: si usas booktabs, recuerda \\usepackage{booktabs} en el preámbulo.
  out
}

latex_ic_counts_table <- function(hat_mcmc,
                                  hat_vi  = NULL,
                                  hat_svi = NULL,
                                  caption = "Cobertura de IC por método",
                                  label   = "tab:ic_counts",
                                  row_labels = c("$\\beta_{1j}$","$\\beta_{2j}$","$\\sigma^2_{j}$"),
                                  show_totals = TRUE,   # muestra "x / J" en lugar de solo "x"
                                  methods_order = c("MCMC","VI","SVI")) {
  
  # --- helper: computa conteos por método ---
  counts_for <- function(hv) {
    stopifnot(is.list(hv), !is.null(hv$Betaj), !is.null(hv$sigma2j))
    J <- nrow(hv$Betaj$contains_zero)
    z1 <- sum(hv$Betaj$contains_zero[,1], na.rm = TRUE)
    m1 <- sum(hv$Betaj$contains_mean[,1] , na.rm = TRUE)
    z2 <- sum(hv$Betaj$contains_zero[,2], na.rm = TRUE)
    m2 <- sum(hv$Betaj$contains_mean[,2] , na.rm = TRUE)
    zs <- sum(hv$sigma2j$contains_zero   , na.rm = TRUE)
    ms <- sum(hv$sigma2j$contains_mean   , na.rm = TRUE)
    list(J = J,
         beta1 = c(z0 = z1, mean = m1),
         beta2 = c(z0 = z2, mean = m2),
         s2    = c(z0 = zs, mean = ms))
  }
  
  have <- list(MCMC = !is.null(hat_mcmc),
               VI   = !is.null(hat_vi),
               SVI  = !is.null(hat_svi))
  meths <- methods_order[methods_order %in% names(have)[unlist(have)]]
  K <- length(meths)
  if (K == 0) stop("Debes pasar al menos hat_mcmc, hat_vi o hat_svi.")
  
  C <- list()
  if ("MCMC" %in% meths) C$MCMC <- counts_for(hat_mcmc)
  if ("VI"   %in% meths) C$VI   <- counts_for(hat_vi)
  if ("SVI"  %in% meths) C$SVI  <- counts_for(hat_svi)
  
  # asumimos mismo J entre métodos presentes
  J <- unique(sapply(C, `[[`, "J"))[1]
  
  # formato de celda
  fmt <- function(x) if (show_totals) sprintf("%d / %d", x, J) else sprintf("%d", x)
  
  # --- header LaTeX dinámico ---
  cols_spec <- paste(rep("cc", K), collapse = "")
  header1 <- paste0(" & ",
                    paste(sprintf("\\multicolumn{2}{c}{\\textbf{%s}}", meths),
                          collapse = " & "),
                    " \\\\")
  # cmidrules
  cmids <- character(K)
  start <- 2
  for (i in seq_len(K)) {
    cmids[i] <- sprintf("\\cmidrule(lr){%d-%d}", start, start+1)
    start <- start + 2
  }
  header2 <- paste0("Parámetro & ",
                    paste(rep(c("$\\#\\,CI\\ni 0$","$\\#\\,CI\\,\\ni\\,\\bar{\\theta}$"), K),
                          collapse = " & "),
                    " \\\\")
  
  # --- filas ---
  row_line <- function(parm) {
    cells <- unlist(lapply(meths, function(m) {
      c(fmt(C[[m]][[parm]]["z0"]), fmt(C[[m]][[parm]]["mean"]))
    }))
    paste(c(row_labels[ c(beta1=1,beta2=2,s2=3)[[parm]] ], cells), collapse = " & ")
  }
  
  rows <- paste0(
    row_line("beta1"), " \\\n",
    row_line("beta2"), " \\\n",
    row_line("s2"),    " \\\n"
  )
  
  # --- ensamblar tabla completa ---
  latex <- sprintf(
    "\\begin{table}[ht]
\\centering
\\caption{%s}
\\label{%s}
\\begin{tabular}{l%s}
\\toprule
%s
%s
%s
\\midrule
%s
\\bottomrule
\\end{tabular}
\\end{table}",
    caption, label, paste0(rep(" ", 0), cols_spec),
    header1, paste(cmids, collapse=""), header2, rows
  )
  
  return(latex)
}




plot_hlrm_betaj_panels <- function(hat_values, which = c(1, 2),
                                   ref_values = NULL, ...) {
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  nr <- length(which)
  par(mfrow = c(nr, 1))
  if (is.null(ref_values)) ref_values <- rep(NA_real_, nr)
  for (i in seq_along(which)) {
    rv <- if (is.na(ref_values[i])) NULL else ref_values[i]
    plot_hlrm_betaj_by_group(hat_values, jcoef = which[i], ref_value = rv, ...)
  }
}



plot_hlrm_lines_by_group <- function(
    data, x, y, group,
    hat_values_mcmc,
    hat_values_vi  = NULL,
    hat_values_svi = NULL,
    
    # eje Y
    sync_ylim = TRUE,
    ylim      = NULL,
    y_pad_frac = 0.06,
    
    # layout y títulos
    layout = c(4, 6),
    main   = "Rectas por grupo: MCMC vs VI vs SVI",
    xlab   = NULL,
    ylab   = NULL,
    group_label = "Group",
    
    # >>> NUEVO: control de orden de grupos
    group_order = NULL,         # vector opcional con el orden deseado de etiquetas
    
    # estética de puntos
    point_cols = NULL,
    point_pch  = NULL,
    point_cex  = 0.8,
    
    # estética de líneas
    col_mcmc = "red",       lty_mcmc = "13", lwd_mcmc = 1.5,
    col_vi   = "black",     lty_vi   = "42", lwd_vi   = 1.5,
    col_svi  = "darkgreen", lty_svi  = "44", lwd_svi  = 1.5,
    
    # márgenes y leyenda
    mar = c(3, 3, 1.2, 0.5),
    oma = c(1.5, 1.5, 3.6, 1),
    cex_axis = 0.9,
    cex_lab  = 0.95,
    legend_outside = TRUE,
    legend_side    = "top",   # "top","bottom","left","right"
    legend_cex     = 0.9,
    legend_inset   = -0.10,   # escalar o c(dx,dy)
    ...
){
  .col_from_len <- function(J) grDevices::hcl.colors(J, palette = "Dark 3")
  .pch_from_len <- function(J){ base <- 0:25; base[(1:J - 1) %% length(base) + 1] }
  .get_B <- function(obj, tag){
    if (is.null(obj) || !is.list(obj) || is.null(obj$Betaj) || is.null(obj$Betaj$mean))
      stop(sprintf("‘%s’ debe tener $Betaj$mean (J×2).", tag))
    obj$Betaj$mean
  }
  .fix_lty <- function(v) as.character(v)
  
  x_name <- deparse(substitute(x))
  y_name <- deparse(substitute(y))
  g_name <- deparse(substitute(group))
  
  xv <- data[[x_name]]; yv <- data[[y_name]]; gv <- data[[g_name]]
  if (is.null(xv) || is.null(yv) || is.null(gv)) stop("Revisa x, y y group.")
  
  # ====== ORDEN DE GRUPOS (ARREGLADO) ======
  g_char <- as.character(gv)
  
  if (!is.null(group_order)) {
    # usar el orden que te pasen
    glab <- as.character(group_order)
    gid  <- match(g_char, glab)
    if (any(is.na(gid)))
      stop("Algunos grupos de 'data' no están en 'group_order'.")
  } else {
    # detectar si las etiquetas son numéricas (p.ej. "1","2","10"...)
    g_try_num <- suppressWarnings(as.numeric(g_char))
    if (all(!is.na(g_try_num))) {
      # ordenar numéricamente
      u  <- unique(g_try_num)
      ord <- order(u)
      u  <- u[ord]
      glab <- as.character(u)
      gid  <- match(g_try_num, u)
    } else if (is.factor(gv)) {
      # si es factor no numérico: usar niveles en el orden del factor
      glab <- levels(gv)
      gid  <- as.integer(gv)
    } else {
      # no numérico: mantener orden de aparición (no alfabético)
      u <- unique(g_char)
      glab <- u
      gid  <- match(g_char, u)
    }
  }
  J <- length(glab)
  
  B_mcmc <- .get_B(hat_values_mcmc, "hat_values_mcmc")
  if (ncol(B_mcmc) < 2) stop("MCMC: Betaj$mean debe ser J×2 (intercepto, pendiente).")
  have_vi  <- !is.null(hat_values_vi)
  have_svi <- !is.null(hat_values_svi)
  B_vi  <- if (have_vi)  .get_B(hat_values_vi,  "hat_values_vi")  else NULL
  B_svi <- if (have_svi) .get_B(hat_values_svi, "hat_values_svi") else NULL
  if (have_vi  && ncol(B_vi)  < 2) stop("VI: Betaj$mean debe ser J×2.")
  if (have_svi && ncol(B_svi) < 2) stop("SVI: Betaj$mean debe ser J×2.")
  
  if (is.null(point_cols)) point_cols <- .col_from_len(J)
  if (length(point_cols) < J) point_cols <- rep(point_cols, length.out = J)
  if (is.null(point_pch))  point_pch  <- .pch_from_len(J)
  if (length(point_pch) < J) point_pch <- rep(point_pch, length.out = J)
  
  # límites Y
  if (!is.null(ylim)) {
    ylim_all <- ylim
  } else if (sync_ylim) {
    yr <- range(yv, na.rm = TRUE)
    if (!is.finite(diff(yr)) || diff(yr) == 0) yr <- yr + c(-1, 1)
    pad <- diff(yr) * y_pad_frac
    ylim_all <- yr + c(-pad, pad)
  } else {
    ylim_all <- NULL
  }
  
  op <- par(no.readonly = TRUE); on.exit(par(op), add = TRUE)
  par(mfrow = layout, mar = mar, oma = oma, cex.axis = cex_axis, cex.lab = cex_lab)
  
  for (j in seq_len(J)) {
    idx <- which(gid == j)
    xj  <- xv[idx]; yj <- yv[idx]
    
    xr <- range(xj, na.rm = TRUE); if (!is.finite(diff(xr)) || diff(xr) == 0) xr <- xr + c(-1, 1)
    yr <- if (is.null(ylim_all)) {
      r <- range(yj, na.rm = TRUE); if (!is.finite(diff(r)) || diff(r) == 0) r <- r + c(-1, 1)
      r + c(-1, 1) * diff(r) * y_pad_frac
    } else ylim_all
    
    dx <- diff(xr); xr <- xr + c(-0.04, 0.04) * dx
    
    plot.new(); plot.window(xlim = xr, ylim = yr)
    axis(1); axis(2); box()
    if (!is.null(xlab)) title(xlab = xlab, line = 2)
    if (!is.null(ylab)) title(ylab = ylab, line = 2)
    mtext(bquote(.(group_label)~.(glab[j])), side = 3, line = 0.2, cex = 0.9)
    
    points(xj, yj, pch = point_pch[j], col = point_cols[j], cex = point_cex)
    
    abline(a = B_mcmc[j, 1], b = B_mcmc[j, 2], col = col_mcmc, lty = .fix_lty(lty_mcmc), lwd = lwd_mcmc)
    if (have_vi)  abline(a = B_vi[j, 1],  b = B_vi[j, 2],  col = col_vi,   lty = .fix_lty(lty_vi),   lwd = lwd_vi)
    if (have_svi) abline(a = B_svi[j, 1], b = B_svi[j, 2], col = col_svi,  lty = .fix_lty(lty_svi),  lwd = lwd_svi)
  }
  
  if (!is.null(main) && nzchar(main)) mtext(main, outer = TRUE, cex = 1.05, line = 0.6)
  
  # leyenda global
  labs <- c("MCMC", if (have_vi) "VI", if (have_svi) "SVI")
  cols <- c(col_mcmc, if (have_vi) col_vi, if (have_svi) col_svi)
  ltys <- .fix_lty(c(lty_mcmc, if (have_vi) lty_vi, if (have_svi) lty_svi))
  lwds <- c(lwd_mcmc, if (have_vi) lwd_vi, if (have_svi) lwd_svi)
  
  op2 <- par(no.readonly = TRUE)
  par(fig = c(0, 1, 0, 1), mar = c(0, 0, 0, 0), new = TRUE, xpd = NA)
  plot.new()
  ins <- legend_inset
  if (length(ins) == 1L) ins <- if (legend_side %in% c("top","bottom")) c(0, ins) else c(ins, 0)
  legend(x = legend_side,
         inset = if (legend_outside) ins else c(0, 0),
         horiz = legend_side %in% c("top","bottom"),
         legend = labs, col = cols, lty = ltys, lwd = lwds,
         bty = "n", cex = legend_cex, seg.len = 2.6)
  par(op2)
}


# ===========================================================
# 1) Pares (beta_{j1} vs beta_{j2}) por grupo: MCMC/VI/SVI
#     - nube MCMC (+ contornos opcionales)
#     - medias como puntos
#     - IC marginales como líneas
#     - LEYENDA DINÁMICA: top/bottom/left/right + inset (overlay)
# ===========================================================
plot_hlrm_betas_by_group <- function(
    results_mcmc,              # $chains$Betaj: S x J x P
    hat_values_mcmc,           # $Betaj$mean, $CI_LI, $CI_LS : J x P
    hat_values_vi  = NULL,     # idem (opcional)
    hat_values_svi = NULL,     # idem (opcional)
    
    groups = NULL,             # qué grupos graficar (1..J)
    layout = c(4, 6),
    S = 4000,                  # submuestreo de iteraciones MCMC
    
    # nube y contornos
    point_col = adjustcolor("grey20", 0.12), point_cex = 0.35,
    add_kde = TRUE, kde_n = 80, kde_levels = 6,
    kde_lwd = 0.8, kde_col = adjustcolor("black", 0.35),
    
    # medias (puntos)
    mcmc_col = "red", vi_col = "black", svi_col = "darkgreen",
    pch_mcmc = 16, pch_vi = 15, pch_svi = 17, cex_mean = 1.8,
    
    # IC marginales (líneas)
    lwd_ci = 0.9, lty_mcmc_ci = "13", lty_vi_ci = "42", lty_svi_ci = "44",
    ci_nudge_frac = 0.00,   # para despegar SVI/VI si se tapan
    
    # cruces en la media (opcional)
    draw_mean_cross = FALSE, mean_cross_lwd = 1.2, mean_cross_lty = 1,
    
    # títulos, ejes, márgenes y leyenda
    main = expression(paste("Dispersión ", beta[j,1], " vs ", beta[j,2],
                            " por grupo (MCMC + VI + SVI)")),
    mar = c(3, 3, 1.1, 0.6), oma = c(1.5, 1.5, 2, 12),
    cex_axes = 0.9, cex_lab = 0.95,
    legend_outside = TRUE, legend_side = "right",
    legend_inset = -0.10, legend_cex = 0.95,
    group_label = "Group"
){
  stopifnot(is.list(results_mcmc), is.list(hat_values_mcmc))
  Bm <- results_mcmc$chains$Betaj
  stopifnot(length(dim(Bm)) == 3L)
  S_all <- dim(Bm)[1]; J <- dim(Bm)[2]; P <- dim(Bm)[3]
  if (P < 2L) stop("Se requieren al menos 2 parámetros por grupo (intercepto y pendiente).")
  
  # grupos a mostrar
  if (is.null(groups)) groups <- seq_len(min(J, prod(layout)))
  groups <- intersect(groups, seq_len(J))
  
  # submuestreo de iteraciones para la nube
  idx <- if (S >= S_all) seq_len(S_all) else sample.int(S_all, S)
  
  # resúmenes
  m_mean <- hat_values_mcmc$Betaj$mean
  m_LI   <- hat_values_mcmc$Betaj$CI_LI
  m_LS   <- hat_values_mcmc$Betaj$CI_LS
  
  have_vi  <- !is.null(hat_values_vi)
  if (have_vi) {
    v_mean <- hat_values_vi$Betaj$mean
    v_LI   <- hat_values_vi$Betaj$CI_LI
    v_LS   <- hat_values_vi$Betaj$CI_LS
  }
  have_svi <- !is.null(hat_values_svi)
  if (have_svi) {
    s_mean <- hat_values_svi$Betaj$mean
    s_LI   <- hat_values_svi$Betaj$CI_LI
    s_LS   <- hat_values_svi$Betaj$CI_LS
  }
  
  # helpers etiquetas
  xlabs <- lapply(groups, function(j) bquote(beta[.(j)*","*1]))
  ylabs <- lapply(groups, function(j) bquote(beta[.(j)*","*2]))
  
  # abrir paneles
  old <- par(no.readonly = TRUE); on.exit(par(old))
  par(mfrow = layout, mar = mar, oma = oma, xpd = NA, cex.axis = cex_axes, cex.lab = cex_lab)
  
  for (k in seq_along(groups)) {
    j <- groups[k]
    x <- Bm[idx, j, 1]   # beta_{j,1}
    y <- Bm[idx, j, 2]   # beta_{j,2}
    
    # límites que incluyan IC de todos los métodos
    rx <- range(x, m_LI[j,1], m_LS[j,1],
                if (have_vi)  c(v_LI[j,1], v_LS[j,1]) else NULL,
                if (have_svi) c(s_LI[j,1], s_LS[j,1]) else NULL, na.rm = TRUE)
    ry <- range(y, m_LI[j,2], m_LS[j,2],
                if (have_vi)  c(v_LI[j,2], v_LS[j,2]) else NULL,
                if (have_svi) c(s_LI[j,2], s_LS[j,2]) else NULL, na.rm = TRUE)
    xlim <- rx + c(-0.06, 0.06) * diff(rx)
    ylim <- ry + c(-0.06, 0.06) * diff(ry)
    
    plot.new(); plot.window(xlim = xlim, ylim = ylim)
    axis(1); axis(2); box()
    title(xlab = as.expression(xlabs[[k]]), ylab = as.expression(ylabs[[k]]), line = 1.8)
    mtext(sprintf("%s %d", group_label, j), side = 3, line = 0.2, cex = 0.9)
    
    # nube + contorno
    points(x, y, pch = 16, cex = point_cex, col = point_col)
    if (add_kde && requireNamespace("MASS", quietly = TRUE)) {
      kd <- MASS::kde2d(x, y, n = kde_n, lims = c(xlim, ylim))
      contour(kd, add = TRUE, drawlabels = FALSE, nlevels = kde_levels,
              lwd = kde_lwd, col = kde_col)
    }
    
    # medias (puntos)
    points(m_mean[j,1], m_mean[j,2], pch = pch_mcmc, cex = cex_mean, col = mcmc_col)
    if (have_vi)  points(v_mean[j,1], v_mean[j,2], pch = pch_vi,  cex = cex_mean, col = vi_col)
    if (have_svi) points(s_mean[j,1], s_mean[j,2], pch = pch_svi, cex = cex_mean, col = svi_col)
    
    # cruces opcionales
    if (draw_mean_cross) {
      abline(v = m_mean[j,1], h = m_mean[j,2], col = mcmc_col, lwd = mean_cross_lwd, lty = mean_cross_lty)
      if (have_vi)  abline(v = v_mean[j,1], h = v_mean[j,2], col = vi_col,  lwd = mean_cross_lwd, lty = mean_cross_lty)
      if (have_svi) abline(v = s_mean[j,1], h = s_mean[j,2], col = svi_col, lwd = mean_cross_lwd, lty = mean_cross_lty)
    }
    
    # IC marginales
    segments(m_LI[j,1], par("usr")[3], m_LI[j,1], par("usr")[4], col = mcmc_col, lwd = lwd_ci, lty = lty_mcmc_ci)
    segments(m_LS[j,1], par("usr")[3], m_LS[j,1], par("usr")[4], col = mcmc_col, lwd = lwd_ci, lty = lty_mcmc_ci)
    segments(par("usr")[1], m_LI[j,2], par("usr")[2], m_LI[j,2], col = mcmc_col, lwd = lwd_ci, lty = lty_mcmc_ci)
    segments(par("usr")[1], m_LS[j,2], par("usr")[2], m_LS[j,2], col = mcmc_col, lwd = lwd_ci, lty = lty_mcmc_ci)
    
    if (have_vi) {
      segments(v_LI[j,1], par("usr")[3], v_LI[j,1], par("usr")[4], col = vi_col, lwd = lwd_ci, lty = lty_vi_ci)
      segments(v_LS[j,1], par("usr")[3], v_LS[j,1], par("usr")[4], col = vi_col, lwd = lwd_ci, lty = lty_vi_ci)
      segments(par("usr")[1], v_LI[j,2], par("usr")[2], v_LI[j,2], col = vi_col, lwd = lwd_ci, lty = lty_vi_ci)
      segments(par("usr")[1], v_LS[j,2], par("usr")[2], v_LS[j,2], col = vi_col, lwd = lwd_ci, lty = lty_vi_ci)
    }
    
    if (have_svi) {
      nx <- ci_nudge_frac * diff(xlim); ny <- ci_nudge_frac * diff(ylim)
      segments(s_LI[j,1] - nx, par("usr")[3], s_LI[j,1] - nx, par("usr")[4], col = svi_col, lwd = lwd_ci, lty = lty_svi_ci)
      segments(s_LS[j,1] + nx, par("usr")[3], s_LS[j,1] + nx, par("usr")[4], col = svi_col, lwd = lwd_ci, lty = lty_svi_ci)
      segments(par("usr")[1], s_LI[j,2] - ny, par("usr")[2], s_LI[j,2] - ny, col = svi_col, lwd = lwd_ci, lty = lty_svi_ci)
      segments(par("usr")[1], s_LS[j,2] + ny, par("usr")[2], s_LS[j,2] + ny, col = svi_col, lwd = lwd_ci, lty = lty_svi_ci)
    }
  }
  
  if (!is.null(main)) mtext(main, outer = TRUE, line = 0.2, cex = 1)
  
  # ------ LEYENDA DINÁMICA (overlay como en lines) ------
  leg_lab <- c("MCMC posterior mean",
               if (have_vi)  "VI posterior mean"  else NULL,
               if (have_svi) "SVI posterior mean" else NULL,
               "MCMC 95% CI",
               if (have_vi)  "VI 95% CI"  else NULL,
               if (have_svi) "SVI 95% CI" else NULL)
  
  leg_col <- c(mcmc_col,
               if (have_vi)  vi_col  else NULL,
               if (have_svi) svi_col else NULL,
               mcmc_col,
               if (have_vi)  vi_col  else NULL,
               if (have_svi) svi_col else NULL)
  
  leg_pch <- c(pch_mcmc,
               if (have_vi)  pch_vi  else NULL,
               if (have_svi) pch_svi else NULL,
               rep(NA, 1 + as.integer(have_vi) + as.integer(have_svi)))
  
  leg_lty <- c(rep(NA, 1 + as.integer(have_vi) + as.integer(have_svi)),
               lty_mcmc_ci,
               if (have_vi)  lty_vi_ci  else NULL,
               if (have_svi) lty_svi_ci else NULL)
  
  leg_lwd <- c(rep(NA, 1 + as.integer(have_vi) + as.integer(have_svi)),
               lwd_ci,
               if (have_vi)  lwd_ci else NULL,
               if (have_svi) lwd_ci else NULL)
  
  op2 <- par(no.readonly = TRUE)
  par(fig = c(0, 1, 0, 1), mar = c(0, 0, 0, 0), new = TRUE, xpd = NA)
  plot.new()
  ins <- legend_inset
  if (length(ins) == 1L) ins <- if (legend_side %in% c("top","bottom")) c(0, ins) else c(ins, 0)
  legend(x = legend_side,
         inset = if (legend_outside) ins else c(0, 0),
         horiz = legend_side %in% c("top","bottom"),
         legend = leg_lab, col = leg_col, pch = leg_pch, pt.cex = cex_mean,
         lty = leg_lty, lwd = leg_lwd, bty = "n", cex = legend_cex, seg.len = 2.6)
  par(op2)
}


# ===========================================================
# 2) sigma^2_j por grupo (puntos + IC) para MCMC/VI/SVI
#     - LEYENDA DINÁMICA (overlay) igual que arriba
# ===========================================================
plot_hlrm_sigma2j_by_group <- function(
    hat_values_mcmc,
    hat_values_vi   = NULL,
    hat_values_svi  = NULL,
    
    groups  = NULL,
    layout  = c(4, 6),
    sync_ylim = TRUE,
    ylim      = NULL,
    y_pad_frac = 0.06,
    
    # estética medias
    col_mcmc = "red", col_vi = "black", col_svi = "darkgreen",
    pch_mcmc = 16,     pch_vi = 15,       pch_svi = 17,
    cex_mean = 1.5,
    
    # estética IC
    lwd_ci = 0.9, cap_size = 0.08,
    lty_ci_mcmc = "13", lty_ci_vi = "42", lty_ci_svi = "44",
    
    # baseline (opcional)
    ols_sigma2j = NULL, col_ols = "orangered", lty_ols = 1, lwd_ols = 1,
    
    # ejes / márgenes / títulos
    main = expression(paste(sigma[j]^2, " por grupo: MCMC vs VI vs SVI")),
    xlab = NULL,                 # etiqueta común en X
    ylab_global = NULL,          # si NULL, usa sigma[j]^2 por panel
    mar = c(3, 3, 1.1, 0.6), oma = c(1.5, 1.5, 2, 10),
    cex_axes = 0.9, cex_lab = 0.95,
    
    # leyenda dinámica
    legend_outside = TRUE, legend_side = "right", legend_inset = -0.10,
    legend_cex = 0.95,
    
    group_label = "Grupo"
){
  stopifnot(is.list(hat_values_mcmc), !is.null(hat_values_mcmc$sigma2j))
  m <- hat_values_mcmc$sigma2j
  J <- length(m$mean)
  
  have_vi  <- !is.null(hat_values_vi)
  have_svi <- !is.null(hat_values_svi)
  
  if (is.null(groups)) groups <- seq_len(J)
  groups <- intersect(groups, seq_len(J))
  if (!length(groups)) stop("‘groups’ no tiene intersección con 1:J.")
  
  # límites comunes del eje Y
  get_range <- function(h) range(c(h$CI_LI, h$CI_LS), na.rm = TRUE)
  if (is.null(ylim)) {
    rngs <- list(get_range(m))
    if (have_vi)  rngs <- c(rngs, list(get_range(hat_values_vi$sigma2j)))
    if (have_svi) rngs <- c(rngs, list(get_range(hat_values_svi$sigma2j)))
    yr  <- range(do.call(c, rngs), na.rm = TRUE)
    pad <- diff(yr) * y_pad_frac
    ylim_all <- c(yr[1] - pad, yr[2] + pad)
  } else ylim_all <- ylim
  
  # preparar paneles
  op <- par(no.readonly = TRUE); on.exit(par(op))
  par(mfrow = layout, mar = mar, oma = oma, xpd = NA,
      cex.axis = cex_axes, cex.lab = cex_lab)
  
  xpos <- c(MCMC = 1, VI = 2, SVI = 3)
  
  draw_ci <- function(x, lo, hi, col, lty, lwd, cap){
    if (is.finite(lo) && is.finite(hi)) {
      segments(x, lo, x, hi, col = col, lty = lty, lwd = lwd)
      if (cap > 0) {
        segments(x - cap, lo, x + cap, lo, col = col, lty = lty, lwd = lwd)
        segments(x - cap, hi, x + cap, hi, col = col, lty = lty, lwd = lwd)
      }
    }
  }
  
  for (j in groups) {
    yl <- if (sync_ylim) ylim_all else {
      rr <- range(c(m$CI_LI[j], m$CI_LS[j],
                    if (have_vi)  c(hat_values_vi$sigma2j$CI_LI[j],
                                    hat_values_vi$sigma2j$CI_LS[j])  else NULL,
                    if (have_svi) c(hat_values_svi$sigma2j$CI_LI[j],
                                    hat_values_svi$sigma2j$CI_LS[j]) else NULL),
                  na.rm = TRUE)
      rr + c(-1, 1) * diff(rr) * y_pad_frac
    }
    
    plot.new(); plot.window(xlim = c(0.5, 3.5), ylim = yl)
    axis(2); box()
    axis(1, at = xpos,
         labels = c("MCMC", if (have_vi) "VI" else "", if (have_svi) "SVI" else ""))
    if (!is.null(xlab)) title(xlab = xlab, line = 2)
    ylab_panel <- if (is.null(ylab_global)) bquote(sigma[.(j)]^2) else ylab_global
    title(ylab = ylab_panel, line = 2)
    mtext(sprintf("%s %d", group_label, j), side = 3, line = 0.2, cex = 0.9)
    
    if (!is.null(ols_sigma2j)) {
      val <- if (length(ols_sigma2j) == 1L) ols_sigma2j else ols_sigma2j[j]
      abline(h = val, col = col_ols, lty = lty_ols, lwd = lwd_ols)
    }
    
    # MCMC / VI / SVI
    draw_ci(xpos["MCMC"], m$CI_LI[j], m$CI_LS[j], col_mcmc, lty_ci_mcmc, lwd_ci, cap_size)
    points(xpos["MCMC"], m$mean[j], pch = pch_mcmc, cex = cex_mean, col = col_mcmc)
    
    if (have_vi) {
      v <- hat_values_vi$sigma2j
      draw_ci(xpos["VI"], v$CI_LI[j], v$CI_LS[j], col_vi, lty_ci_vi, lwd_ci, cap_size)
      points(xpos["VI"], v$mean[j], pch = pch_vi, cex = cex_mean, col = col_vi)
    }
    
    if (have_svi) {
      s <- hat_values_svi$sigma2j
      draw_ci(xpos["SVI"], s$CI_LI[j], s$CI_LS[j], col_svi, lty_ci_svi, lwd_ci, cap_size)
      points(xpos["SVI"], s$mean[j], pch = pch_svi, cex = cex_mean, col = col_svi)
    }
  }
  
  if (!is.null(main)) mtext(main, outer = TRUE, cex = 1.05, line = 0.4)
  
  # ------ LEYENDA DINÁMICA (overlay) -------
  leg_lab <- c("MCMC posterior mean",
               if (have_vi)  "VI posterior mean"  else NULL,
               if (have_svi) "SVI posterior mean" else NULL,
               "MCMC 95% CI",
               if (have_vi)  "VI 95% CI"  else NULL,
               if (have_svi) "SVI 95% CI" else NULL)
  
  leg_col <- c(col_mcmc,
               if (have_vi)  col_vi  else NULL,
               if (have_svi) col_svi else NULL,
               col_mcmc,
               if (have_vi)  col_vi  else NULL,
               if (have_svi) col_svi else NULL)
  
  leg_pch <- c(pch_mcmc,
               if (have_vi)  pch_vi  else NULL,
               if (have_svi) pch_svi else NULL,
               rep(NA, 1 + as.integer(have_vi) + as.integer(have_svi)))
  
  leg_lty <- c(rep(NA, 1 + as.integer(have_vi) + as.integer(have_svi)),
               lty_ci_mcmc,
               if (have_vi)  lty_ci_vi  else NULL,
               if (have_svi) lty_ci_svi else NULL)
  
  leg_lwd <- c(rep(NA, 1 + as.integer(have_vi) + as.integer(have_svi)),
               lwd_ci,
               if (have_vi)  lwd_ci else NULL,
               if (have_svi) lwd_ci else NULL)
  
  op2 <- par(no.readonly = TRUE)
  par(fig = c(0, 1, 0, 1), mar = c(0, 0, 0, 0), new = TRUE, xpd = NA)
  plot.new()
  ins <- legend_inset
  if (length(ins) == 1L) ins <- if (legend_side %in% c("top","bottom")) c(0, ins) else c(ins, 0)
  legend(x = legend_side,
         inset = if (legend_outside) ins else c(0, 0),
         horiz = legend_side %in% c("top","bottom"),
         legend = leg_lab, col = leg_col, pch = leg_pch, pt.cex = cex_mean,
         lty = leg_lty, lwd = leg_lwd, bty = "n", cex = legend_cex, seg.len = 2.6)
  par(op2)
}


add_group_labels <- function(data, x, y, group,
                             fun = median,                 # centro por grupo: median/mean
                             cex = 0.85, offset = 0.5,     # tamaño y distancia de la etiqueta
                             col_map = NULL,               # paleta de colores (vector)
                             use_plotrix = TRUE,           # intenta colocar inteligentemente
                             draw_boxes = FALSE,           # dibuja cajita blanca detrás
                             box_bg = "white", box_border = NA) {
  # --- extrae vectores
  xv <- data[[deparse(substitute(x))]]
  yv <- data[[deparse(substitute(y))]]
  gv <- data[[deparse(substitute(group))]]
  if (is.null(xv) || is.null(yv) || is.null(gv)) stop("Revisa x, y y group.")
  
  # --- centros por grupo (robusto: mediana por defecto)
  labs <- aggregate(list(x = xv, y = yv), by = list(g = gv), FUN = fun, na.rm = TRUE)
  
  # --- paleta: si no te pasan, usa la misma de tus puntos si existe; si no, HCL
  if (is.null(col_map)) {
    if ("paletteer_d" %in% getNamespaceExports("paletteer")) {
      col_map <- paletteer::paletteer_d("ggsci::default_aaas")
    } else {
      col_map <- grDevices::hcl.colors(max(length(unique(gv)), 3), "Dark 3")
    }
  }
  # mapea color por grupo (sea numérico o factor)
  g_index <- as.integer(factor(labs$g))              # 1..J en el orden de 'labs'
  cols_lab <- rep(col_map, length.out = max(g_index))[g_index]
  
  # --- calcula posición "inteligente" si hay plotrix, si no, usa pos=3 (arriba)
  pos <- rep(3, nrow(labs)) # por defecto
  if (use_plotrix && requireNamespace("plotrix", quietly = TRUE)) {
    pos <- plotrix::thigmophobe(labs$x, labs$y) # 1:below, 2:left, 3:above, 4:right
  }
  
  # --- dibuja etiquetas
  if (draw_boxes && requireNamespace("plotrix", quietly = TRUE)) {
    # cajita con fondo (no evita solapes, pero mejora legibilidad)
    plotrix::boxed.labels(labs$x, labs$y, labels = as.character(labs$g),
                          bg = box_bg, border = box_border, cex = cex, col = cols_lab)
  } else {
    graphics::text(labs$x, labs$y, labels = as.character(labs$g),
                   pos = pos, offset = offset, cex = cex, col = cols_lab, xpd = NA)
  }
  
  invisible(labs)
}


## =========================
## 1) CORRER TODA LA GRILLA
## =========================
run_vi_hlrm_grid <- function(
    y, X, hypers,
    p1_vals, kappa_vals, batch_vals,        # vectores; usa NA en batch para corrida "batch"
    iterations   = 20000,
    samples_gof  = 1000,
    standardize  = FALSE,
    label_fun = function(p1, kappa, batch, stochastic)
      if (!stochastic) "batch" else sprintf("minibatch: %s", batch),
    seed = NULL,
    verbose = TRUE
){
  # NO asumimos clase de X; solo necesitamos N para validar minibatch
  N <- length(y)
  
  grid <- expand.grid(
    p1    = p1_vals,
    kappa = kappa_vals,
    batch = batch_vals,
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
  )
  
  runs <- vector("list", nrow(grid))
  
  for (i in seq_len(nrow(grid))) {
    p1    <- grid$p1[i]
    kappa <- grid$kappa[i]
    batch <- grid$batch[i]
    
    stochastic <- !is.na(batch)
    if (stochastic) {
      if (batch > N) stop("minibatch_m > N; ajusta 'batch'.")
      tau <- (p1)^(-1/kappa) - 1
      run_settings <- list(
        iterations   = iterations,
        samples_gof  = samples_gof,
        stochastic   = TRUE,
        minibatch_m  = as.integer(batch),
        tau          = tau,
        kappa        = kappa
      )
    } else {
      run_settings <- list(
        iterations   = iterations,
        samples_gof  = samples_gof,
        stochastic   = FALSE
      )
    }
    
    lab <- label_fun(p1, kappa, batch, stochastic)
    if (!is.null(seed)) set.seed(seed + i)
    
    vi <- vi_hlrm$new(y, X, standarize = standardize, hypers, run_settings)
    res <- vi$fit(verbose = verbose)
    
    runs[[i]] <- list(results = res, label = lab, settings = run_settings)
  }
  
  names(runs) <- vapply(runs, `[[`, character(1), "label")
  runs
}


## =============================================
## 2) PROMEDIO MÓVIL + SD (sin repetir corridas)
## =============================================
summarize_elbo_run <- function(results, k = 200, step = k, start_at = k) {
  stopifnot(is.list(results),
            !is.null(results$timing$up_to_iter),
            !is.null(results$params$elbo))
  t <- as.numeric(results$timing$up_to_iter)
  y <- as.numeric(results$params$elbo)
  n <- min(length(t), length(y))
  t <- t[seq_len(n)]; y <- y[seq_len(n)]
  
  cs  <- c(0, cumsum(y))
  cs2 <- c(0, cumsum(y^2))
  i <- seq_len(n)
  left <- pmax(0, i - k)
  w <- i - left
  mean_i <- (cs[i + 1]  - cs[left + 1]) / w
  var_i  <- pmax((cs2[i + 1] - cs2[left + 1]) / w - mean_i^2, 0)
  sd_i   <- sqrt(var_i)
  
  keep <- seq.int(from = max(1L, start_at), to = n, by = max(1L, step))
  data.frame(time = t[keep], mean = mean_i[keep], sd = sd_i[keep])
}

# =================== HELPERS ===================
tau_from_p1 <- function(p1, kappa) (p1)^(-1/kappa) - 1
.near <- function(a,b,tol=1e-8) abs(a-b) <= tol

nearest_tau_in_runs <- function(runs, kappa, tau_target, tol = 1e-8) {
  taus_k <- vapply(runs, function(r) {
    s <- r$settings; if (isTRUE(s$kappa == kappa)) as.numeric(s$tau) else NA_real_
  }, numeric(1))
  taus_k <- taus_k[!is.na(taus_k)]
  if (!length(taus_k)) stop("No hay corridas con ese kappa en 'runs'.")
  tau_closest <- taus_k[which.min(abs(taus_k - tau_target))]
  list(tau = tau_closest, exact = .near(tau_closest, tau_target, tol), dist = abs(tau_closest - tau_target))
}

.get_run <- function(runs, m, kappa, tau, tol = 1e-8) {
  idx <- which(vapply(runs, function(r) {
    s <- r$settings
    isTRUE(s$minibatch_m == m) && isTRUE(s$kappa == kappa) && .near(s$tau, tau, tol)
  }, logical(1)))
  if (!length(idx)) return(NULL)
  runs[[idx[1]]]
}

# media y sd en ventanas deslizantes (salida en el extremo derecho de cada ventana)
roll_mean_sd <- function(x, t = NULL, k = 300, step = 300, start_at = 1L) {
  n <- length(x); if (is.null(t)) t <- seq_len(n)
  if (start_at + k - 1 > n) stop("start_at + k - 1 excede el largo de la serie.")
  starts <- seq.int(from = start_at, to = n - k + 1L, by = step)
  ends   <- starts + k - 1L
  m  <- vapply(seq_along(starts), function(i) mean(x[starts[i]:ends[i]]), numeric(1))
  sd <- vapply(seq_along(starts), function(i)  sd(x[starts[i]:ends[i]]),  numeric(1))
  data.frame(x = t[ends], mean = m, sd = sd, stringsAsFactors = FALSE)
}

.cols_n <- function(n) grDevices::hcl.colors(n, "Dark 3")

# Leyenda flexible: fuera o dentro, con ncol/nrow, labels custom
.draw_legend <- function(labs, cols, ltys, lwds,
                         legend_side = "top", legend_outside = TRUE, legend_inset = -0.10,
                         legend_cex = 0.9, legend_ncol = NULL, legend_nrow = NULL,
                         horiz = NULL, bg = NULL) {
  # calcular ncol si usuario dio nrow
  if (!is.null(legend_nrow) && is.null(legend_ncol)) {
    legend_ncol <- ceiling(length(labs) / legend_nrow)
  }
  if (is.null(legend_ncol)) legend_ncol <- 1L
  
  if (legend_outside) {
    op2 <- par(no.readonly = TRUE); on.exit(par(op2), add = TRUE)
    par(fig = c(0,1,0,1), mar = c(0,0,0,0), new = TRUE, xpd = NA)
    plot.new()
    ins <- legend_inset
    if (length(ins) == 1L) ins <- if (legend_side %in% c("top","bottom")) c(0, ins) else c(ins, 0)
    legend(x = legend_side,
           inset = ins,
           horiz = if (is.null(horiz)) legend_side %in% c("top","bottom") else isTRUE(horiz),
           legend = labs, col = cols, lty = ltys, lwd = lwds, ncol = legend_ncol,
           bty = "n", cex = legend_cex, seg.len = 2.6, bg = bg)
  } else {
    # dentro del panel actual
    legend(x = legend_side,
           legend = labs, col = cols, lty = ltys, lwd = lwds, ncol = legend_ncol,
           bty = "n", cex = legend_cex, seg.len = 2.6, bg = bg)
  }
}

# Dibuja "error bars" verticales: media ± sd con pequeños "caps"
.draw_errbars <- function(x, ymean, ysd, col, lwd = 1.4, lty = 1, cap_frac = 0.015, xlim = NULL) {
  y0 <- ymean - ysd; y1 <- ymean + ysd
  segments(x0 = x, y0 = y0, x1 = x, y1 = y1, col = col, lwd = lwd, lty = lty)
  if (!is.null(xlim) && length(x) > 1L) {
    cap <- diff(range(xlim)) * cap_frac
  } else {
    cap <- if (length(x) > 1L) diff(range(x)) * cap_frac else 1 * cap_frac
  }
  segments(x - cap, y0, x + cap, y0, col = col, lwd = lwd, lty = lty)
  segments(x - cap, y1, x + cap, y1, col = col, lwd = lwd, lty = lty)
}

# =================== 1) ELBO ~ tiempo, comparando m (minibatch) ===================
plot_elbo_by_minibatch <- function(
    runs,
    kappa, p1 = NULL, tau = NULL,                  # define p1 o tau
    minibatch_vals = NULL,                          # por defecto: todos los m disponibles
    k = 300, step = 300, start_at = 1,
    main = "", xlab = "Tiempo acumulado (s)", ylab = "ELBO (media móvil)",
    ylim = NULL, y_pad_frac = 0.06,
    # estética
    cols = NULL, ltys = NULL, lwds = NULL, draw_points = FALSE, pch = 16, cex = 0.7,
    err_lwd = 1.2, err_lty = 1, err_cap_frac = 0.015,   # barras de error con segments
    # márgenes y leyenda
    mar = c(3.2, 3.2, 1.2, 0.8), oma = c(1.5, 1.5, 2, 1),
    legend_labels = NULL,
    legend_side = "top", legend_outside = TRUE, legend_inset = -0.10, legend_cex = 0.9,
    legend_ncol = NULL, legend_nrow = NULL, legend_bg = NULL
){
  # Resolver tau si viene p1
  if (is.null(tau) && !is.null(p1)) {
    tau_target <- tau_from_p1(p1, kappa)
    pick <- nearest_tau_in_runs(runs, kappa, tau_target)
    tau <- pick$tau
    if (!pick$exact) message(sprintf("Aviso: tau objetivo=%.6g; usando tau más cercano en runs=%.6g (Δ=%.3g).",
                                     tau_target, tau, pick$dist))
  } else if (is.null(tau) && is.null(p1)) stop("Pasa 'p1' o 'tau'.")
  
  # m disponibles
  if (is.null(minibatch_vals)) {
    minibatch_vals <- sort(unique(vapply(runs, function(r) r$settings$minibatch_m, numeric(1))))
  }
  
  # Construir series (media móvil y sd)
  series <- list()
  for (m in minibatch_vals) {
    rr <- .get_run(runs, m = m, kappa = kappa, tau = tau)
    if (is.null(rr)) next
    elbo <- rr$results$params$elbo
    tt   <- rr$timing$up_to_iter
    df   <- roll_mean_sd(elbo, t = tt, k = k, step = step, start_at = start_at)
    series[[as.character(m)]] <- df
  }
  if (!length(series)) stop("No hay corridas con (kappa, tau/p1) y minibatches solicitados.")
  
  # Estética
  nm <- length(series)
  if (is.null(cols)) cols <- .cols_n(nm)
  if (length(cols) < nm) cols <- rep(cols, length.out = nm)
  if (is.null(ltys)) ltys <- rep("solid", nm)
  if (is.null(lwds)) lwds <- rep(2, nm)
  if (is.null(legend_labels)) legend_labels <- paste0("m = ", names(series))
  
  # Rangos
  xr <- range(do.call(c, lapply(series, `[[`, "x")), na.rm = TRUE)
  if (is.null(ylim)) {
    yr <- range(do.call(c, lapply(series, function(d) c(d$mean - d$sd, d$mean + d$sd))), na.rm = TRUE)
    pad <- diff(yr) * y_pad_frac
    ylim <- c(yr[1] - pad, yr[2] + pad)
  }
  
  # Plot base
  op <- par(no.readonly = TRUE); on.exit(par(op), add = TRUE)
  par(mar = mar, oma = oma, xpd = NA)
  plot.new(); plot.window(xlim = xr, ylim = ylim)
  axis(1); axis(2); box()
  title(main = main, xlab = xlab, ylab = ylab)
  
  # Dibujar por serie: línea + errorbars (segments) y opcionalmente puntos
  i <- 0
  for (nm_i in names(series)) {
    i <- i + 1
    d <- series[[nm_i]]
    lines(d$x, d$mean, col = cols[i], lty = ltys[i], lwd = lwds[i])
    .draw_errbars(d$x, d$mean, d$sd, col = cols[i], lwd = err_lwd, lty = err_lty,
                  cap_frac = err_cap_frac, xlim = xr)
    if (draw_points) points(d$x, d$mean, pch = pch, cex = cex, col = cols[i])
  }
  
  # Leyenda
  .draw_legend(
    labs = legend_labels, cols = cols, ltys = ltys, lwds = lwds,
    legend_side = legend_side, legend_outside = legend_outside,
    legend_inset = legend_inset, legend_cex = legend_cex,
    legend_ncol = legend_ncol, legend_nrow = legend_nrow, bg = legend_bg
  )
  
  invisible(list(series = series, kappa = kappa, tau = tau, m = minibatch_vals))
}

# =================== 2) ELBO ~ tiempo, comparando (kappa, tau) con m fijo ===================
plot_elbo_by_kappa_tau <- function(
    runs,
    minibatch_m,
    kappas, p1_vals = NULL, tau_vals = NULL,   # define p1_vals (recomendado) o tau_vals
    k = 300, step = 300, start_at = 1,
    main = "", xlab = "Tiempo acumulado (s)", ylab = "ELBO (media móvil)",
    ylim = NULL, y_pad_frac = 0.06,
    # estética
    cols = NULL, ltys = NULL, lwds = NULL, draw_points = FALSE, pch = 16, cex = 0.7,
    err_lwd = 1.2, err_lty = 1, err_cap_frac = 0.015,
    # márgenes y leyenda
    mar = c(3.2, 3.2, 1.2, 0.8), oma = c(1.5, 1.5, 2, 1),
    legend_labels = NULL,  # si NULL se construyen automáticamente
    legend_side = "top", legend_outside = TRUE, legend_inset = -0.10, legend_cex = 0.9,
    legend_ncol = NULL, legend_nrow = NULL, legend_bg = NULL,
    label_fun = function(kappa, p1, tau) {
      if (!is.null(p1)) sprintf("κ=%.2g, p1=%.3g (τ=%.3g)", kappa, p1, tau) else sprintf("κ=%.2g, τ=%.3g", kappa, tau)
    }
){
  K <- length(kappas)
  if (!is.null(p1_vals) && length(p1_vals) != K) stop("'p1_vals' debe tener misma longitud que 'kappas'.")
  if (!is.null(tau_vals) && length(tau_vals) != K) stop("'tau_vals' debe tener misma longitud que 'kappas'.")
  if (is.null(p1_vals) && is.null(tau_vals)) stop("Pasa 'p1_vals' o 'tau_vals'.")
  
  series <- list(); labs <- character(0)
  for (i in seq_len(K)) {
    kap <- kappas[i]
    if (!is.null(p1_vals)) {
      tau_target <- tau_from_p1(p1_vals[i], kap)
      pick <- nearest_tau_in_runs(runs, kap, tau_target)
      tau_use <- pick$tau
      if (!pick$exact) message(sprintf("[κ=%.2g] tau objetivo=%.6g; usando tau en runs=%.6g (Δ=%.3g).",
                                       kap, tau_target, tau_use, pick$dist))
      lab <- label_fun(kap, p1_vals[i], tau_use)
    } else {
      tau_use <- tau_vals[i]
      pick <- nearest_tau_in_runs(runs, kap, tau_use)  # por si no existe el exacto
      tau_use <- pick$tau
      lab <- label_fun(kap, NULL, tau_use)
    }
    
    rr <- .get_run(runs, m = minibatch_m, kappa = kap, tau = tau_use)
    if (is.null(rr)) { message(sprintf("No hay corrida para m=%s, κ=%.2g, τ≈%.3g — se omite.", minibatch_m, kap, tau_use)); next }
    d <- roll_mean_sd(rr$results$params$elbo, t = rr$timing$up_to_iter,
                      k = k, step = step, start_at = start_at)
    series[[paste0("k", kap)]] <- d
    labs <- c(labs, lab)
  }
  if (!length(series)) stop("No hay series para graficar con los parámetros dados.")
  
  ns <- length(series)
  if (is.null(cols)) cols <- .cols_n(ns)
  if (length(cols) < ns) cols <- rep(cols, length.out = ns)
  if (is.null(ltys)) ltys <- rep("solid", ns)
  if (is.null(lwds)) lwds <- rep(2, ns)
  if (is.null(legend_labels)) legend_labels <- labs
  
  xr <- range(do.call(c, lapply(series, `[[`, "x")), na.rm = TRUE)
  if (is.null(ylim)) {
    yr <- range(do.call(c, lapply(series, function(d) c(d$mean - d$sd, d$mean + d$sd))), na.rm = TRUE)
    pad <- diff(yr) * y_pad_frac
    ylim <- c(yr[1] - pad, yr[2] + pad)
  }
  
  op <- par(no.readonly = TRUE); on.exit(par(op), add = TRUE)
  par(mar = mar, oma = oma, xpd = NA)
  plot.new(); plot.window(xlim = xr, ylim = ylim)
  axis(1); axis(2); box()
  title(main = main, xlab = xlab, ylab = ylab)
  
  i <- 0
  for (nm_i in names(series)) {
    i <- i + 1
    d <- series[[nm_i]]
    lines(d$x, d$mean, col = cols[i], lty = ltys[i], lwd = lwds[i])
    .draw_errbars(d$x, d$mean, d$sd, col = cols[i], lwd = err_lwd, lty = err_lty,
                  cap_frac = err_cap_frac, xlim = xr)
    if (draw_points) points(d$x, d$mean, pch = pch, cex = cex, col = cols[i])
  }
  
  .draw_legend(
    labs = legend_labels, cols = cols, ltys = ltys, lwds = lwds,
    legend_side = legend_side, legend_outside = legend_outside,
    legend_inset = legend_inset, legend_cex = legend_cex,
    legend_ncol = legend_ncol, legend_nrow = legend_nrow, bg = legend_bg
  )
  
  invisible(list(series = series, m = minibatch_m))
}


make_latex_goodness_table <- function(
    goodness_mcmc, goodness_vi, goodness_svi,
    methods = c("MCMC", "VI", "SVI"),
    # columnas y redondeo
    cols_order   = c("waic", "dic", "mse", "r2"),
    col_headers  = c("WAIC", "DIC", "MSE", "$R^2$"),
    digits       = c(waic = 0, dic = 0, mse = 3, r2 = 3),
    # estilo LaTeX
    booktabs     = TRUE,        # usa \toprule \midrule \bottomrule
    siunitx      = FALSE,       # usa columnas S de siunitx
    align        = NULL,        # por defecto: l r r r r (o l S S S S si siunitx=TRUE)
    caption      = "Métricas de bondad de ajuste",
    label        = "tab:goodness",
    note         = NULL,        # texto opcional bajo la tabla
    na_print     = "--",
    float_env    = "table",     # "table" o "table*"
    file         = NULL         # si no es NULL, guarda el .tex
){
  `%||%` <- function(a,b) if (!is.null(a)) a else b
  
  # --- helper para extraer seguro ---
  extract_metrics <- function(g) {
    list(
      waic = tryCatch(as.numeric(g$waic$WAIC),  error = function(e) NA_real_),
      dic  = tryCatch(as.numeric(g$dic$DIC),    error = function(e) NA_real_),
      mse  = tryCatch(as.numeric(g$mse),        error = function(e) NA_real_),
      r2   = tryCatch(as.numeric(g$r2),         error = function(e) NA_real_)
    )
  }
  
  fmt_num <- function(x, d) {
    ifelse(is.na(x), na_print, formatC(x, format = "f", digits = d))
  }
  
  # normaliza digits
  d_waic <- digits[["waic"]] %||% 0
  d_dic  <- digits[["dic"]]  %||% 0
  d_mse  <- digits[["mse"]]  %||% 3
  d_r2   <- digits[["r2"]]   %||% 3
  d_vec  <- c(waic = d_waic, dic = d_dic, mse = d_mse, r2 = d_r2)
  
  # arma data.frame método x métricas
  goods <- list(goodness_mcmc, goodness_vi, goodness_svi)
  vals  <- lapply(goods, extract_metrics)
  df    <- do.call(rbind, lapply(vals, function(v) unlist(v)[cols_order]))
  rownames(df) <- methods
  
  # formateo (siunitx no necesita formateo como texto, pero igual dejamos fijos)
  df_fmt <- as.data.frame(df, stringsAsFactors = FALSE)
  for (nm in cols_order) df_fmt[[nm]] <- fmt_num(df_fmt[[nm]], d_vec[[nm]])
  
  # columnas & alineación
  if (isTRUE(siunitx)) {
    # requiere \usepackage{siunitx} en el preámbulo
    colspec <- paste(c("l", rep("S", length(cols_order))), collapse = " ")
  } else {
    colspec <- paste(c("l", rep("r", length(cols_order))), collapse = " ")
  }
  if (!is.null(align)) colspec <- align
  
  # líneas LaTeX
  ln <- c()
  ln <- c(ln, sprintf("\\begin{%s}[!ht]", float_env))
  ln <- c(ln, "  \\centering")
  if (!is.null(caption)) ln <- c(ln, sprintf("  \\caption{%s}", caption))
  if (!is.null(label))   ln <- c(ln, sprintf("  \\label{%s}",   label))
  ln <- c(ln, sprintf("  \\begin{tabular}{%s}", colspec))
  if (booktabs) {
    ln <- c(ln, "    \\toprule")
  } else {
    ln <- c(ln, "    \\hline")
  }
  
  # encabezados
  header <- paste(c("", col_headers), collapse = " & ")
  ln <- c(ln, paste0("    ", header, " \\\\"))
  if (booktabs) {
    ln <- c(ln, "    \\midrule")
  } else {
    ln <- c(ln, "    \\hline")
  }
  
  # filas
  for (i in seq_len(nrow(df_fmt))) {
    row <- paste(c(rownames(df_fmt)[i], unname(unlist(df_fmt[i, ]))), collapse = " & ")
    ln  <- c(ln, paste0("    ", row, " \\\\"))
  }
  
  if (booktabs) {
    ln <- c(ln, "    \\bottomrule")
  } else {
    ln <- c(ln, "    \\hline")
  }
  ln <- c(ln, "  \\end{tabular}")
  
  if (!is.null(note) && nzchar(note)) {
    ln <- c(ln, sprintf("  \\\\[-0.6em]\\footnotesize\\emph{%s}", note))
  }
  
  ln <- c(ln, sprintf("\\end{%s}", float_env))
  
  out <- paste(ln, collapse = "\n")
  if (!is.null(file)) {
    writeLines(out, con = file)
    message(sprintf("Tabla LaTeX guardada en: %s", normalizePath(file)))
  }
  return(out)
}



## =========================================================
## Helper: extraer coeficiente j respetando orientación KxP o PxK
## =========================================================
.pick_coef_vec <- function(A, K, j) {
  stopifnot(is.matrix(A))
  if (nrow(A) == K) A[, j] else if (ncol(A) == K) A[j, ] else
    stop("Dimensión de Betak no concuerda con K.")
}

## =========================================================
## 1) β_k: medias + IC, usando maps (sin tocar los objetos)
## =========================================================
plot_betak_ci_maps <- function(
    hat_mcmc, map_mcmc,
    hat_vi   = NULL, map_vi  = NULL,
    hat_svi  = NULL, map_svi = NULL,
    beta_index = 1,                  # cuál coeficiente j
    groups     = NULL,               # subconjunto de grupos canónicos (1..K)
    layout     = c(1, 3),
    
    # NUEVO:
    group_labels = NULL,             # etiquetas por grupo canónico (vector de largo K)
    xlab_common  = NULL,             # etiqueta común en el borde inferior
    
    # ejes
    sync_ylim  = FALSE,              # por defecto NO alinea eje vertical
    y_pad_frac = 0.06,
    
    # estética
    cols  = c(MCMC="red", VI="black", SVI="darkgreen"),
    pch   = c(MCMC=16, VI=15, SVI=17),
    lty   = c(MCMC="13", VI="42", SVI="44"),
    lwd_ci = 1.0, cap_size = 0.08, cex_mean = 1.5,
    
    # títulos / leyenda
    main = expression(paste(beta[j], " por grupo: MCMC vs VI vs SVI")),
    zero_line = TRUE,
    legend_side  = "right", legend_inset = -0.12, legend_cex = 0.95,
    
    # márgenes / cex
    mar = c(3,3,1.1,0.6), oma = c(1.8,1.5,2,8), cex_axes = 0.9, cex_lab = 0.95
){
  stopifnot(is.list(hat_mcmc), is.numeric(map_mcmc))
  K <- length(map_mcmc)
  
  # Detecta orientación Betak en MCMC
  Bm <- hat_mcmc$Betak$mean
  if (!(is.matrix(Bm) && (nrow(Bm) == K || ncol(Bm) == K)))
    stop("Betak (MCMC) no concuerda con length(map_mcmc).")
  
  if (is.null(groups)) groups <- seq_len(K)
  
  # validar/usar group_labels
  if (!is.null(group_labels)) {
    if (length(group_labels) < max(groups))
      stop("length(group_labels) debe cubrir el índice máximo en 'groups'.")
    group_labels <- as.character(group_labels)
  }
  
  # Construye lista de métodos presentes
  methods <- list(MCMC = list(h=hat_mcmc, map=map_mcmc))
  if (!is.null(hat_vi)  && !is.null(map_vi))  methods$VI  <- list(h=hat_vi,  map=map_vi)
  if (!is.null(hat_svi) && !is.null(map_svi)) methods$SVI <- list(h=hat_svi, map=map_svi)
  
  # Tripletas (li, mean, ls) ya reordenadas por map, por método
  get_triplet_mapped <- function(h, map, j){
    list(
      li   = .pick_coef_vec(if (is.null(h$Betak$CI_LI)) stop("Falta Betak$CI_LI") else h$Betak$CI_LI, K, j)[map],
      mean = .pick_coef_vec(if (is.null(h$Betak$mean))  stop("Falta Betak$mean")  else h$Betak$mean,  K, j)[map],
      ls   = .pick_coef_vec(if (is.null(h$Betak$CI_LS)) stop("Falta Betak$CI_LS") else h$Betak$CI_LS, K, j)[map]
    )
  }
  trips <- lapply(methods, function(m) get_triplet_mapped(m$h, m$map, beta_index))
  
  # Rango Y global (si sync_ylim=TRUE)
  if (isTRUE(sync_ylim)) {
    yr  <- range(unlist(lapply(trips, function(t) c(t$li[groups], t$ls[groups]))), na.rm = TRUE)
    pad <- diff(yr) * y_pad_frac
    ylim_all <- c(yr[1]-pad, yr[2]+pad)
  }
  
  # Plot
  op <- par(no.readonly = TRUE); on.exit(par(op))
  par(mfrow = layout, mar = mar, oma = oma, xpd = NA, cex.axis = cex_axes, cex.lab = cex_lab)
  
  mnames <- names(methods)
  xpos   <- setNames(seq_along(mnames), mnames)
  
  draw_ci <- function(x, lo, hi, col, lty, lwd, cap) {
    segments(x, lo, x, hi, col = col, lty = lty, lwd = lwd)
    if (cap > 0) {
      segments(x-cap, lo, x+cap, lo, col = col, lty = lty, lwd = lwd)
      segments(x-cap, hi, x+cap, hi, col = col, lty = lty, lwd = lwd)
    }
  }
  
  for (k in groups) {
    if (isTRUE(sync_ylim)) {
      yl <- ylim_all
    } else {
      vals_k <- unlist(lapply(trips, function(t) c(t$li[k], t$ls[k])))
      rk  <- range(vals_k, na.rm = TRUE)
      pad <- diff(rk) * y_pad_frac
      yl  <- c(rk[1]-pad, rk[2]+pad)
    }
    
    plot.new(); plot.window(xlim = c(0.5, length(mnames)+0.5), ylim = yl)
    axis(2); box(); axis(1, at = xpos, labels = mnames)
    if (zero_line) abline(h = 0, lty = 2, col = "gray80")
    
    panel_lab <- if (is.null(group_labels)) sprintf("Grupo %d", k) else group_labels[k]
    mtext(panel_lab, side = 3, line = 0.2, cex = 0.9)
    title(ylab = bquote(beta[.(beta_index)]), line = 2)
    
    for (nm in mnames) {
      t  <- trips[[nm]]
      xx <- xpos[[nm]]
      draw_ci(xx, t$li[k], t$ls[k], col = cols[[nm]], lty = lty[[nm]], lwd = lwd_ci, cap = cap_size)
      points(xx, t$mean[k], pch = pch[[nm]], col = cols[[nm]], cex = cex_mean)
    }
  }
  
  if (!is.null(main)) mtext(main, outer = TRUE, cex = 1.05, line = 0.4)
  if (!is.null(xlab_common)) mtext(xlab_common, side = 1, outer = TRUE, line = 0.7)
  
  # Leyenda fuera
  par(fig = c(0,1,0,1), mar = c(0,0,0,0), new = TRUE, xpd = NA); plot.new()
  legend(legend_side, inset = c(legend_inset, 0), bty = "n", cex = legend_cex,
         legend = c(paste(mnames, "posterior mean"), paste(mnames, "95% CI")),
         col = c(unname(cols[mnames]), unname(cols[mnames])),
         pch = c(unname(pch[mnames]), rep(NA, length(mnames))),
         lty = c(rep(NA, length(mnames)), unname(lty[mnames])),
         lwd = c(rep(NA, length(mnames)), rep(lwd_ci, length(mnames))))
}

## =========================================================
## 2) σ_k^2: medias + IC, usando maps (sin tocar objetos)
## =========================================================
plot_sigma2k_ci_maps <- function(
    hat_mcmc, map_mcmc,
    hat_vi   = NULL, map_vi  = NULL,
    hat_svi  = NULL, map_svi = NULL,
    groups     = NULL,
    layout     = c(1, 3),
    
    # NUEVO:
    group_labels = NULL,             # etiquetas por grupo canónico (vector de largo K)
    xlab_common  = NULL,             # etiqueta común en el borde inferior
    
    # ejes
    sync_ylim  = FALSE,              # por defecto NO alinea
    y_pad_frac = 0.06,
    
    # estética
    cols  = c(MCMC="red", VI="black", SVI="darkgreen"),
    pch   = c(MCMC=16, VI=15, SVI=17),
    lty   = c(MCMC="13", VI="42", SVI="44"),
    lwd_ci = 1.0, cap_size = 0.08, cex_mean = 1.5,
    
    # títulos / leyenda
    main = expression(paste(sigma[k]^2, " por grupo: MCMC vs VI vs SVI")),
    legend_side  = "right", legend_inset = -0.12, legend_cex = 0.95,
    
    # márgenes / cex
    mar = c(3,3,1.1,0.6), oma = c(1.8,1.5,2,8), cex_axes = 0.9, cex_lab = 0.95
){
  # detectar sigma2k/sigma2j en MCMC y K
  s2name <- if (!is.null(hat_mcmc$sigma2k)) "sigma2k" else if (!is.null(hat_mcmc$sigma2j)) "sigma2j" else
    stop("El objeto MCMC no tiene sigma2k/sigma2j.")
  K <- length(hat_mcmc[[s2name]]$mean)
  if (is.null(groups)) groups <- seq_len(K)
  
  # validar/usar group_labels
  if (!is.null(group_labels)) {
    if (length(group_labels) < max(groups))
      stop("length(group_labels) debe cubrir el índice máximo en 'groups'.")
    group_labels <- as.character(group_labels)
  }
  
  methods <- list(MCMC = list(h=hat_mcmc, map=map_mcmc))
  if (!is.null(hat_vi)  && !is.null(map_vi))  methods$VI  <- list(h=hat_vi,  map=map_vi)
  if (!is.null(hat_svi) && !is.null(map_svi)) methods$SVI <- list(h=hat_svi, map=map_svi)
  
  # Vector (li, mean, ls) ya reordenado por map, por método
  get_triplet_mapped <- function(h, map, s2name){
    s2 <- h[[s2name]]
    list(
      li   = as.numeric(s2$CI_LI)[map],
      mean = as.numeric(s2$mean)[map],
      ls   = as.numeric(s2$CI_LS)[map]
    )
  }
  trips <- lapply(methods, function(m) get_triplet_mapped(m$h, m$map, s2name))
  
  # Rango Y global si sync=TRUE
  if (isTRUE(sync_ylim)) {
    yr  <- range(unlist(lapply(trips, function(t) c(t$li[groups], t$ls[groups]))), na.rm = TRUE)
    pad <- diff(yr) * y_pad_frac
    ylim_all <- c(yr[1]-pad, yr[2]+pad)
  }
  
  op <- par(no.readonly = TRUE); on.exit(par(op))
  par(mfrow = layout, mar = mar, oma = oma, xpd = NA, cex.axis = cex_axes, cex.lab = cex_lab)
  
  mnames <- names(methods)
  xpos   <- setNames(seq_along(mnames), mnames)
  
  draw_ci <- function(x, lo, hi, col, lty, lwd, cap) {
    segments(x, lo, x, hi, col = col, lty = lty, lwd = lwd)
    if (cap > 0) {
      segments(x-cap, lo, x+cap, lo, col = col, lty = lty, lwd = lwd)
      segments(x-cap, hi, x+cap, hi, col = col, lty = lty, lwd = lwd)
    }
  }
  
  for (k in groups) {
    if (isTRUE(sync_ylim)) {
      yl <- ylim_all
    } else {
      vals_k <- unlist(lapply(trips, function(t) c(t$li[k], t$ls[k])))
      rk  <- range(vals_k, na.rm = TRUE)
      pad <- diff(rk) * y_pad_frac
      yl  <- c(rk[1]-pad, rk[2]+pad)
    }
    
    plot.new(); plot.window(xlim = c(0.5, length(mnames)+0.5), ylim = yl)
    axis(2); box(); axis(1, at = xpos, labels = mnames)
    
    panel_lab <- if (is.null(group_labels)) sprintf("Grupo %d", k) else group_labels[k]
    mtext(panel_lab, side = 3, line = 0.2, cex = 0.9)
    title(ylab = bquote(sigma[.(k)]^2), line = 2)
    
    for (nm in mnames) {
      t  <- trips[[nm]]
      xx <- xpos[[nm]]
      draw_ci(xx, t$li[k], t$ls[k], col = cols[[nm]], lty = lty[[nm]], lwd = lwd_ci, cap = cap_size)
      points(xx, t$mean[k], pch = pch[[nm]], col = cols[[nm]], cex = cex_mean)
    }
  }
  
  if (!is.null(main)) mtext(main, outer = TRUE, cex = 1.05, line = 0.4)
  if (!is.null(xlab_common)) mtext(xlab_common, side = 1, outer = TRUE, line = 0.7)
  
  par(fig = c(0,1,0,1), mar = c(0,0,0,0), new = TRUE, xpd = NA); plot.new()
  legend(legend_side, inset = c(legend_inset, 0), bty = "n", cex = legend_cex,
         legend = c(paste(mnames, "posterior mean"), paste(mnames, "95% CI")),
         col = c(unname(cols[mnames]), unname(cols[mnames])),
         pch = c(unname(pch[mnames]), rep(NA, length(mnames))),
         lty = c(rep(NA, length(mnames)), unname(lty[mnames])),
         lwd = c(rep(NA, length(mnames)), rep(lwd_ci, length(mnames))))
}

plot_3d_by_group <- function(
    X, y,
    # --- leyendas y grupos ---
    group_ids = NULL,              # opcional: vector con id de grupo por fila de Xall; si NULL, se infiere por lista X
    group_levels = NULL,           # opcional: nombres de grupos; por defecto 1..G
    col_groups = NULL,             # vector de colores por grupo (longitud G)
    pch_groups = NULL,             # vector de pch por grupo (longitud G)
    # --- layout de leyenda ---
    legend_title = "Grupo",
    legend_pos   = "topright",     # "topright","bottomleft", coordenadas c(x,y), o "right" con inset negativo
    legend_inset = 0.02,           # puede ser negativo para poner fuera del plot si xpd=NA
    legend_bty   = "n",
    legend_cex   = 0.95,
    legend_ncol  = NULL,           # número de columnas en la leyenda
    legend_nrow  = NULL,           # número de filas deseadas (se ajusta ncol en consecuencia)
    legend_byrow = FALSE,          # TRUE: llenar por filas (reordenamos entradas)
    # --- puntos/graf ---
    pt_cex = 0.9,                  # tamaño de puntos
    main   = "Scatter 3D por grupo",
    xlab   = "x1", ylab = "x2", zlab = "y",
    grid   = TRUE, box = TRUE,
    xpd_outside_legend = TRUE      # permite dibujar leyenda fuera del gráfico
){
  stopifnot(is.list(X))
  
  # --- Unir datos ---
  Xall <- do.call(rbind, X)
  yall <- if (is.list(y)) unlist(y) else as.numeric(y)
  stopifnot(nrow(Xall) == length(yall))
  
  # --- IDs de grupo por punto ---
  if (is.null(group_ids)) {
    gid <- rep(seq_along(X), vapply(X, nrow, integer(1)))
  } else {
    gid <- as.integer(group_ids)
    stopifnot(length(gid) == nrow(Xall))
  }
  gid_fac <- factor(gid)
  G <- nlevels(gid_fac)
  
  # --- Nombres de grupos ---
  if (is.null(group_levels)) {
    group_levels <- levels(gid_fac)
  } else {
    stopifnot(length(group_levels) == G)
  }
  
  # --- Colores/pch por grupo ---
  if (is.null(col_groups)) col_groups <- grDevices::rainbow(G)
  if (length(col_groups) < G) col_groups <- rep(col_groups, length.out = G)
  
  if (is.null(pch_groups)) pch_groups <- rep(16, G) # círculos sólidos por defecto
  if (length(pch_groups) < G) pch_groups <- rep(pch_groups, length.out = G)
  
  # mapear por punto
  cols_pts <- col_groups[gid_fac]
  pch_pts  <- pch_groups[gid_fac]
  
  # --- Dibujar scatterplot3d ---
  s3d <- scatterplot3d::scatterplot3d(
    x = Xall[, 1], y = Xall[, 2], z = yall,
    color = cols_pts, pch = pch_pts, cex.symbols = pt_cex,
    main = main, xlab = xlab, ylab = ylab, zlab = zlab,
    grid = grid, box = box
  )
  
  # --- Preparar entradas de leyenda ---
  leg_cols <- col_groups
  leg_pchs <- pch_groups
  leg_labs <- group_levels
  
  # Resolver ncol/nrow
  n <- length(leg_labs)
  if (!is.null(legend_nrow) && is.null(legend_ncol)) {
    legend_ncol <- ceiling(n / legend_nrow)
  }
  if (is.null(legend_ncol)) legend_ncol <- 1L
  
  # Reordenar si se quiere llenar por filas
  if (isTRUE(legend_byrow)) {
    nrow_eff <- ceiling(n / legend_ncol)
    M <- matrix(seq_len(n), nrow = nrow_eff, ncol = legend_ncol, byrow = TRUE)
    ord <- as.vector(M)
    ord <- ord[ord <= n]
    leg_cols <- leg_cols[ord]
    leg_pchs <- leg_pchs[ord]
    leg_labs <- leg_labs[ord]
  }
  
  # --- Dibujar leyenda (posible fuera del área) ---
  op <- par(no.readonly = TRUE); on.exit(par(op))
  if (isTRUE(xpd_outside_legend)) par(xpd = NA)
  
  legend(x = legend_pos,
         legend = leg_labs,
         title  = legend_title,
         col    = leg_cols,
         pch    = leg_pchs,
         ncol   = legend_ncol,
         inset  = legend_inset,
         bty    = legend_bty,
         cex    = legend_cex,
         pt.cex = pt_cex)
  
  invisible(s3d)
}

# Graficar el último ELBO vs K a partir de una lista de resultados tipo run_vi_over_K(...)
plot_elbo_last_by_K <- function(
    results_list,
    Ks = NULL,                    # subset opcional de K a graficar; por defecto, todos los presentes
    na_rm = TRUE,                 # remover K sin ELBO válido
    # Estética del gráfico
    type = "b", pch = 16, col = "#2c7bb6", lwd = 2, lty = 1, cex = 1.1,
    xlab = "K (número de clusters)",
    ylab = "ELBO (valor final)",
    main = "ELBO final por K",
    xlim = NULL, ylim = NULL,
    grid = TRUE,
    # Marcar el mejor K (máximo ELBO)
    mark_best = TRUE,
    best_point_col = "firebrick",
    best_point_pch = 19,
    best_point_cex = 1.3,
    best_label = TRUE,
    best_label_fmt = "Mejor K = %d",
    best_label_pos = NULL        # NULL => al lado del punto; o usa c(x,y) en coords del plot
) {
  stopifnot(is.list(results_list))
  
  # --- helper: extraer K de cada entrada ---
  extract_K <- function(name_i, obj_i) {
    if (!is.null(obj_i$K)) return(as.integer(obj_i$K))
    # intentar parsear del nombre "K=14"
    if (!is.null(name_i) && grepl("^K=", name_i)) {
      val <- suppressWarnings(as.integer(sub("^K=", "", name_i)))
      if (!is.na(val)) return(val)
    }
    NA_integer_
  }
  
  # --- helper: último ELBO ---
  extract_last_elbo <- function(obj_i) {
    v <- NULL
    if (!is.null(obj_i$params) && !is.null(obj_i$params$elbo)) v <- obj_i$params$elbo
    else if (!is.null(obj_i$elbo)) v <- obj_i$elbo
    else return(NA_real_)
    v <- as.numeric(v)
    v <- v[is.finite(v)]
    if (length(v) == 0) return(NA_real_)
    tail(v, 1)
  }
  
  # Recolectar K y ELBO
  n <- length(results_list)
  Ks_vec   <- integer(n)
  elbo_vec <- numeric(n)
  
  for (i in seq_len(n)) {
    nm <- names(results_list)[i]
    obj <- results_list[[i]]
    Ks_vec[i]   <- extract_K(nm, obj)
    elbo_vec[i] <- extract_last_elbo(obj)
  }
  
  # Armar data.frame base
  df <- data.frame(K = Ks_vec, ELBO_last = elbo_vec, stringsAsFactors = FALSE)
  
  # Filtrar por Ks si se pide
  if (!is.null(Ks)) df <- df[df$K %in% Ks, , drop = FALSE]
  
  # Ordenar por K creciente
  df <- df[order(df$K), , drop = FALSE]
  
  # Opcional: remover NA
  if (na_rm) df <- df[is.finite(df$ELBO_last) & !is.na(df$K), , drop = FALSE]
  
  if (nrow(df) == 0L) {
    stop("No hay valores de ELBO válidos para graficar (revisa 'results_list' o desactiva na_rm).")
  }
  
  # Rango de ejes
  if (is.null(xlim)) xlim <- range(df$K, na.rm = TRUE)
  if (is.null(ylim)) {
    yr <- range(df$ELBO_last, na.rm = TRUE)
    pad <- diff(yr) * 0.05
    ylim <- c(yr[1] - pad, yr[2] + pad)
  }
  
  # Plot base
  plot(df$K, df$ELBO_last, type = type, pch = pch, col = col, lwd = lwd, lty = lty, cex = cex,
       xlab = xlab, ylab = ylab, main = main, xlim = xlim, ylim = ylim)
  
  if (isTRUE(grid)) {
    abline(h = pretty(ylim), v = pretty(xlim), col = "gray90", lwd = 1, lty = "dotted")
  }
  
  # Marcar mejor K (máximo ELBO)
  if (isTRUE(mark_best)) {
    idx_best <- which.max(df$ELBO_last)
    points(df$K[idx_best], df$ELBO_last[idx_best],
           pch = best_point_pch, col = best_point_col, cex = best_point_cex)
    if (isTRUE(best_label)) {
      if (is.null(best_label_pos)) {
        # etiqueta pegada al punto
        text(df$K[idx_best], df$ELBO_last[idx_best],
             labels = sprintf(best_label_fmt, df$K[idx_best]),
             pos = 4, cex = 0.9, col = best_point_col, xpd = NA)
      } else {
        # coordenadas absolutas (c(x,y))
        text(best_label_pos[1], best_label_pos[2],
             labels = sprintf(best_label_fmt, df$K[idx_best]),
             cex = 0.9, col = best_point_col, xpd = NA)
      }
    }
  }
  
  invisible(df)
}

# install.packages(c("knitr", "kableExtra"))  # si no los tienes
library(knitr)
library(kableExtra)

# ------------------------------------------------------------
# Construye dos tablas LaTeX (ajustables) desde 3 objetos goodness_* y tiempos
# ------------------------------------------------------------
# install.packages(c("knitr", "kableExtra"))  # si hace falta
library(knitr)
library(kableExtra)

make_latex_goodness_tables <- function(
    goodness_mcmc, goodness_vi, goodness_svi,
    times = c(MCMC = NA, VI = NA, SVI = NA),
    
    method_names = c("MCMC", "VI", "SVI"),
    col_method   = "Método",
    
    caption_fit = "Criterios de ajuste por método",
    label_fit   = "tab:fit",
    digits_fit  = c(WAIC = 0, DIC = 0, MSE = 2, R2 = 3, Tiempo = 1),
    align_fit   = c("l", "r", "r", "r", "r", "r"),
    
    caption_ppp = "Posterior predictive p-values por método",
    label_ppp   = "tab:ppp",
    ppp_order   = c("min", "max", "IQR", "mean", "med", "sd"),
    ppp_digits  = 2,
    align_ppp   = NULL,
    
    latex_hold_position = TRUE,
    booktabs = TRUE,
    save_tex_prefix = NULL
) {
  gl <- setNames(list(goodness_mcmc, goodness_vi, goodness_svi), method_names)
  
  .ext_fit <- function(g) {
    c(
      WAIC = suppressWarnings(as.numeric(g$waic$WAIC)),
      DIC  = suppressWarnings(as.numeric(g$dic$DIC)),
      MSE  = suppressWarnings(as.numeric(g$mse)),
      R2   = suppressWarnings(as.numeric(g$r2))
    )
  }
  .ext_ppp <- function(g, order) {
    vals <- setNames(rep(NA_real_, length(order)), order)
    if (!is.null(g$ppval) && !is.null(g$ppval$ppp)) {
      pp <- g$ppval$ppp
      inter <- intersect(names(pp), order)
      vals[inter] <- as.numeric(pp[inter])
    }
    vals
  }
  
  ## -------- TABLA 1 --------
  fit_mat <- t(vapply(gl, .ext_fit, numeric(4)))
  # tiempos alineados por nombre; si no tiene nombre, asume orden de method_names
  times_vec <- if (length(names(times))) as.numeric(times[rownames(fit_mat)]) else as.numeric(times)[seq_len(nrow(fit_mat))]
  
  df_fit <- data.frame(
    rownames(fit_mat),            # nombre del método (temporal)
    as.data.frame(fit_mat),
    Tiempo = times_vec,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  names(df_fit)[1] <- col_method    # renombra 1ra columna
  
  digits_vec_fit <- c(NA, digits_fit[c("WAIC","DIC","MSE","R2","Tiempo")])
  if (any(is.na(digits_vec_fit[-1]))) stop("Revisa 'digits_fit' (WAIC, DIC, MSE, R2, Tiempo).")
  
  k_fit <- kable(
    df_fit, format = "latex",
    caption = caption_fit, label = label_fit,
    booktabs = booktabs,
    digits = digits_vec_fit,
    align = align_fit,
    escape = TRUE
  )
  if (isTRUE(latex_hold_position)) k_fit <- k_fit %>% kable_styling(latex_options = "hold_position")
  
  ## -------- TABLA 2 --------
  ppp_mat <- t(vapply(gl, .ext_ppp, numeric(length(ppp_order)), order = ppp_order))
  df_ppp <- data.frame(
    rownames(ppp_mat),
    as.data.frame(ppp_mat, stringsAsFactors = FALSE),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  names(df_ppp)[1] <- col_method
  
  digits_vec_ppp <- c(NA, rep(ppp_digits, length(ppp_order)))
  if (is.null(align_ppp)) align_ppp <- c("l", rep("r", length(ppp_order)))
  
  k_ppp <- kable(
    df_ppp, format = "latex",
    caption = caption_ppp, label = label_ppp,
    booktabs = booktabs,
    digits = digits_vec_ppp,
    align = align_ppp,
    escape = TRUE
  )
  if (isTRUE(latex_hold_position)) k_ppp <- k_ppp %>% kable_styling(latex_options = "hold_position")
  
  if (!is.null(save_tex_prefix)) {
    save_kable(k_fit, paste0(save_tex_prefix, "_fit.tex"))
    save_kable(k_ppp, paste0(save_tex_prefix, "_ppp.tex"))
  }
  
  list(
    df_fit     = df_fit,
    df_ppp     = df_ppp,
    latex_fit  = k_fit,
    latex_ppp  = k_ppp
  )
}

# ------------------------------------------------------------
# Ejecuta VI-CHLRM variando K (por defecto 1:15) permitiendo
# que 'hypers' dependa de K (vía lista o función).
# ------------------------------------------------------------
run_vi_over_K <- function(
    y, X,
    Ks = 1:15,
    
    # --- HYPERS dependientes de K (elige UNA de estas o combina) ---
    hypers_byK = NULL,      # lista de hypers por K (ver notas; prioridad más alta)
    hypers_fn  = NULL,      # función K -> hypers  (o function(K,y,X) -> hypers)
    hypers     = NULL,      # hypers fijo (si no das ni hypers_byK ni hypers_fn)
    
    # --- settings de VI ---
    iterations = 20000,
    samples_gof = 1000,
    epsilon = NULL,
    stochastic = NULL,
    minibatch_m = NULL,
    tau = NULL,
    kappa = NULL,
    
    # --- opciones varias ---
    standarize = FALSE,   # (sic) se respeta el nombre de tu constructor
    verbose = TRUE,
    seeds = NULL          # opcional: vector de seeds por K
) {
  stopifnot(is.list(X) || is.matrix(X), length(Ks) >= 1)
  
  # Construir run_settings_vi con tus nombres de campos
  run_settings_vi <- list(
    iterations  = iterations,
    samples_gof = samples_gof
  )
  if (!is.null(epsilon))     run_settings_vi$epsilon     <- epsilon
  if (!is.null(stochastic))  run_settings_vi$stochastic  <- stochastic
  if (!is.null(minibatch_m)) run_settings_vi$minibatch_m <- minibatch_m
  if (!is.null(tau))         run_settings_vi$tau         <- tau
  if (!is.null(kappa))       run_settings_vi$kappa       <- kappa
  
  # Helper: extraer hypers para un K dado respetando prioridad
  fetch_hypers_for_K <- function(K) {
    # 1) hypers_byK (permite nombres "K=3" o índices 3)
    if (!is.null(hypers_byK)) {
      # intenta por nombre "K=3"
      nm <- paste0("K=", K)
      if (is.list(hypers_byK) && !is.null(hypers_byK[[nm]])) return(hypers_byK[[nm]])
      # intenta por índice K (si la lista está en el mismo orden que Ks)
      if (length(hypers_byK) >= K && !is.null(hypers_byK[[K]])) return(hypers_byK[[K]])
      stop(sprintf("No encontré hypers_byK para K=%d (ni nombre '%s' ni índice %d).", K, nm, K))
    }
    # 2) hypers_fn
    if (is.function(hypers_fn)) {
      # permite firmas hypers_fn(K), hypers_fn(K,y,X)
      out <- tryCatch({
        if (length(formals(hypers_fn)) == 1L) hypers_fn(K) else hypers_fn(K, y, X)
      }, error = function(e) e)
      if (inherits(out, "error")) stop(sprintf("Error al construir hypers con hypers_fn(K=%d): %s", K, out$message))
      return(out)
    }
    # 3) hypers fijo
    if (!is.null(hypers)) return(hypers)
    stop("Debes proveer 'hypers_byK' o 'hypers_fn' o un 'hypers' fijo.")
  }
  
  out_list <- vector("list", length(Ks))
  names(out_list) <- paste0("K=", Ks)
  
  # barra de progreso
  if (interactive()) {
    pb <- txtProgressBar(min = 0, max = length(Ks), initial = 0, style = 3)
    on.exit(try(close(pb), silent = TRUE), add = TRUE)
  }
  
  for (i in seq_along(Ks)) {
    K <- Ks[i]
    if (!is.null(seeds)) {
      s <- seeds[min(i, length(seeds))]
      if (!is.na(s)) set.seed(s)
    }
    
    # construir hypers específico para este K
    hypers_K <- fetch_hypers_for_K(K)
    
    warn_msgs <- character()
    res_i <- withCallingHandlers(
      expr = {
        vi <- vi_chlrm$new(y, X, K, standarize = standarize, hypers_K, run_settings_vi)
        results_vi <- vi$fit(verbose = verbose)
        
        list(
          K            = K,
          run_settings = run_settings_vi,
          hypers_used  = hypers_K,                               # << guardamos los hypers usados
          params       = if (!is.null(results_vi$params)) results_vi$params else NULL,
          timing_total = tryCatch(results_vi$timing$total, error = function(e) NA_real_),
          inferences   = tryCatch(vi$inferences(), error = function(e) NULL),
          error        = NULL,
          warnings     = NULL
        )
      },
      warning = function(w) { warn_msgs <<- c(warn_msgs, conditionMessage(w)); invokeRestart("muffleWarning") }
    )
    
    # si algo explotó dentro:
    if (inherits(res_i, "try-error")) {
      out_list[[i]] <- list(
        K            = K,
        run_settings = run_settings_vi,
        hypers_used  = hypers_K,
        params       = NULL,
        timing_total = NA_real_,
        inferences   = NULL,
        error        = as.character(res_i),
        warnings     = warn_msgs
      )
    } else {
      res_i$warnings <- warn_msgs
      out_list[[i]] <- res_i
    }
    
    if (interactive()) setTxtProgressBar(pb, i)
    if (isTRUE(verbose)) message(sprintf("[VI] K=%d terminado.", K))
  }
  
  class(out_list) <- c("vi_grid_results", class(out_list))
  out_list
}
