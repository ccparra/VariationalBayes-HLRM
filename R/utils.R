n_samples_mcmc <- function(run_settings){
  if (!is.null(run_settings$thin)) {
    B <- run_settings$samples*run_settings$thin + run_settings$burn_in
  } else {
    B <- run_settings$samples + run_settings$burn_in
  }
  return(B)
}

thin_chain_mcmc <- function(b,run_settings){
  burn <- if (is.null(run_settings$burn_in)) 0L else as.integer(run_settings$burn_in)
  thin <- if (is.null(run_settings$thin)) 1L else as.integer(run_settings$thin)
  if (thin < 1L) stop("run_settings$thin must be a positive integer")

  post_burn <- as.integer(b) - burn
  if (post_burn < 1L || post_burn %% thin != 0L) return(NULL)
  post_burn %/% thin
}

logmeanexp <- function(x) {
  if (!length(x)) return(-Inf)
  m <- max(x)
  if (!is.finite(m)) return(m)
  m + log(mean(exp(x - m)))
}

logsumexp <- function(x) {
  if (!length(x)) return(-Inf)
  x_max <- max(x)
  if (!is.finite(x_max)) return(x_max)
  x_max + log(sum(exp(x - x_max)))
}

mse <- function(y,yhat){
  return(mean((y-yhat)^2))
}

step_size <- function(t, tau = 10, kappa = 0.7) {
  if (any(!is.finite(c(t, tau, kappa))) || any(t < 1) || tau < 0 ||
      kappa <= 0.5 || kappa > 1) {
    stop("Require t >= 1, tau >= 0 and 0.5 < kappa <= 1")
  }
  (t + tau)^(-kappa)
}

log_multigamma <- function(a, p) (p * (p - 1) / 4) * log(pi) + sum(lgamma(a + (1 - seq_len(p)) / 2))

expected_logdet_inverse_wishart <- function(nu, S) {
  p <- nrow(S)
  return(-sum(digamma((nu + 1 - seq_len(p)) / 2)) - p * log(2) + logdet_spd(S))
}

.log_2pi <- log(2 * pi)

expected_log_inverse_gamma <- function(shape, scale) log(scale) - digamma(shape)
expected_inverse_inverse_gamma <- function(shape, scale) shape / scale

expected_log_dirichlet <- function(alpha) digamma(alpha) - digamma(sum(alpha))


r2 <- function(y, yhat) {
  y    <- as.numeric(y)
  yhat <- as.numeric(yhat)
  sse <- sum((y - yhat)^2)
  sst <- sum((y - mean(y))^2)
  if (sst == 0) return(NA_real_)
 return(1 - sse/sst)
}

project_root <- function(start = getwd()) {
  configured <- Sys.getenv("VB_HLRM_ROOT", unset = "")
  if (nzchar(configured)) return(normalizePath(configured, mustWork = TRUE))

  current <- normalizePath(start, mustWork = TRUE)
  repeat {
    if (dir.exists(file.path(current, ".git")) &&
        dir.exists(file.path(current, "R"))) return(current)
    parent <- dirname(current)
    if (identical(parent, current)) {
      stop("Cannot locate the project root; set VB_HLRM_ROOT")
    }
    current <- parent
  }
}

project_path <- function(...) file.path(project_root(), ...)

project_dir <- function(...) {
  path <- project_path(...)
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  normalizePath(path, mustWork = TRUE)
}

validate_grouped_data <- function(X, y, K = NULL) {
  if (!is.list(X) || !is.list(y) || !length(X) || length(X) != length(y)) {
    stop("X and y must be non-empty lists with the same number of groups")
  }
  p <- ncol(as.matrix(X[[1L]]))
  if (is.null(p) || p < 1L) stop("Each X[[j]] must be a matrix with at least one column")
  for (j in seq_along(X)) {
    Xj <- as.matrix(X[[j]])
    yj <- as.numeric(y[[j]])
    if (ncol(Xj) != p || nrow(Xj) != length(yj) || !length(yj)) {
      stop(sprintf("Incompatible dimensions in group %d", j))
    }
    if (any(!is.finite(Xj)) || any(!is.finite(yj))) {
      stop(sprintf("Non-finite data in group %d", j))
    }
  }
  if (!is.null(K) && (length(K) != 1L || !is.finite(K) || K < 1L || K != as.integer(K))) {
    stop("K must be a positive integer")
  }
  invisible(list(m = length(X), p = p, nj = lengths(y)))
}

validate_mcmc_settings <- function(run_settings) {
  required <- c("samples", "burn_in")
  missing <- setdiff(required, names(run_settings))
  if (length(missing)) stop("Missing MCMC settings: ", paste(missing, collapse = ", "))
  vals <- unlist(run_settings[c("samples", "burn_in", "thin")], use.names = FALSE)
  vals <- vals[!is.na(vals)]
  if (any(!is.finite(vals)) || run_settings$samples < 1 || run_settings$burn_in < 0 ||
      run_settings$samples != as.integer(run_settings$samples) ||
      run_settings$burn_in != as.integer(run_settings$burn_in) ||
      (!is.null(run_settings$thin) &&
       (run_settings$thin < 1 || run_settings$thin != as.integer(run_settings$thin)))) {
    stop("samples/thin must be positive integers and burn_in a non-negative integer")
  }
  invisible(TRUE)
}

# Vehtari et al.'s variance form of WAIC. Rows are predictive units and
# columns posterior draws; grouping observations changes the estimand.
waic_from_loglik <- function(log_lik) {
  contract <- if (inherits(log_lik, "vb_log_likelihood")) log_lik else NULL
  values <- as.matrix(log_lik)
  if (!is.numeric(values) || !nrow(values) || any(!is.finite(values)))
    stop("WAIC requires a nonempty finite numeric log-likelihood matrix")
  if (ncol(values) < 2L) stop("WAIC requires at least two posterior draws")
  lppd_i <- apply(values, 1L, logmeanexp)
  p_i <- apply(values, 1L, stats::var)
  out <- summarize_waic_pointwise(data.frame(lppd = lppd_i, p_waic = p_i))
  if (!is.null(contract)) {
    out$estimand <- log_likelihood_estimand(contract)
    out$predictive_units <- contract$predictive_units
  }
  out
}
