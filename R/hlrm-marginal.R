# Integrate beta_j analytically, then sigma2_j by deterministic quadrature.
# beta_j | Beta,Sigma ~ N(Beta,Sigma), sigma2_j | xi2 ~ IG(nu0/2,nu0*xi2/2),
# independently. Their independence means the marginal is NOT generally t.

gamma_quadrature_rule <- function(n, shape) {
  # Golub-Welsch for Gamma(shape, rate=1), with normalized probability weights.
  i <- seq_len(n)
  J <- diag(2 * i + shape - 2, n)
  if (n > 1) {
    off <- sqrt(seq_len(n - 1) * (seq_len(n - 1) + shape - 1))
    J[cbind(seq_len(n - 1), 2:n)] <- off
    J[cbind(2:n, seq_len(n - 1))] <- off
  }
  e <- eigen(J, symmetric = TRUE)
  keep <- e$values > 0 & e$vectors[1, ] != 0
  list(nodes = e$values[keep], log_weights = 2 * log(abs(e$vectors[1, keep])))
}

loglik_group_marginal_hlrm <- function(y_list, X_list, draws, nu0,
                                        tolerance = 1e-6, max_nodes = 256L,
                                        max_evaluations = 1e6, max_output_mb = 256) {
  info <- validate_grouped_data(X_list, y_list)
  S <- nrow(draws$Beta); p <- info$p
  if (!is.matrix(draws$Beta) || !is.numeric(draws$Beta) || is.null(S) || S < 1 ||
      ncol(draws$Beta) != p || any(!is.finite(draws$Beta))) stop("Beta must be finite [S,p]")
  if (!is.numeric(draws$Sigma) || !identical(dim(draws$Sigma), as.integer(c(S, p, p))) ||
      any(!is.finite(draws$Sigma))) stop("Sigma must be finite [S,p,p]")
  if (!is.numeric(draws$xi2) || length(draws$xi2) != S || any(!is.finite(draws$xi2)) ||
      any(draws$xi2 <= 0)) stop("xi2 must be a positive finite S-vector")
  for (item in list(nu0 = nu0, tolerance = tolerance, max_evaluations = max_evaluations))
    if (!is.numeric(item) || length(item) != 1L || !is.finite(item) || item <= 0)
      stop("nu0, tolerance and max_evaluations must be positive finite scalars")
  if (!is.numeric(max_nodes) || length(max_nodes) != 1L || !is.finite(max_nodes) ||
      max_nodes < 32 || max_nodes != floor(max_nodes) || max_nodes > 1024)
    stop("max_nodes must be an integer between 32 and 1024")
  # Includes eigensystem work arrays, output and rank-p local decompositions.
  validate_evaluation_size(max(2L, S), max_output_mb,
    as.double(info$m) * S + 4 * max_nodes^2 + 3 * max(info$nj) * p)
  levels <- unique(c(2^(4:floor(log2(max_nodes))), max_nodes))
  if (as.double(info$m) * S * sum(levels) > max_evaluations)
    stop("Marginal quadrature exceeds max_evaluations; reduce draws/groups or raise the explicit work budget")
  cache <- new.env(parent = emptyenv())
  rule <- function(nodes, shape) {
    key <- paste(nodes, format(shape, digits = 17), sep = ":")
    if (!exists(key, cache, inherits = FALSE)) assign(key, gamma_quadrature_rule(nodes, shape), cache)
    get(key, cache, inherits = FALSE)
  }
  # Compute one posterior draw at a time. No n_j x n_j covariance is formed.
  out <- matrix(NA_real_, info$m, S)
  used <- matrix(0L, info$m, S); change <- matrix(NA_real_, info$m, S)
  a <- nu0 / 2
  for (s in seq_len(S)) {
    Sigma <- matrix(draws$Sigma[s, , ], p, p)
    if (!isSymmetric(Sigma)) stop("Sigma must be symmetric")
    # PSD includes Sigma=0 solely as a useful limiting-case validation.
    eig <- eigen(Sigma, symmetric = TRUE)
    if (min(eig$values) < -1e-10 * max(1, max(abs(eig$values)))) stop("Sigma must be positive semidefinite")
    root <- sweep(eig$vectors, 2, sqrt(pmax(eig$values, 0)), `*`)
    b <- nu0 * draws$xi2[s] / 2
    for (j in seq_len(info$m)) {
      r <- y_list[[j]] - drop(X_list[[j]] %*% draws$Beta[s, ])
      B <- X_list[[j]] %*% root
      dec <- svd(B, nu = min(dim(B)), nv = 0)
      active <- which(dec$d > 0)
      U <- dec$u[, active, drop = FALSE]
      projection <- drop(crossprod(U, r))
      residual <- r - drop(U %*% projection)
      lambda <- dec$d[active]^2
      A <- a + info$nj[j] / 2
      rate <- b + sum(residual^2) / 2
      constant <- -info$nj[j] / 2 * log(2 * pi) + a * log(b) - lgamma(a) + lgamma(A) - A * log(rate)
      previous <- NULL; converged <- FALSE
      for (nodes in levels) {
        q <- rule(nodes, A); u <- q$nodes / rate
        terms <- vapply(u, function(v) -.5 * sum(log1p(lambda * v)) -
          .5 * sum(projection^2 * v / (1 + lambda * v)), numeric(1))
        value <- constant + logsumexp(q$log_weights + terms)
        if (!is.finite(value)) stop("Non-finite marginal likelihood at group ", j, ", draw ", s)
        if (!is.null(previous)) {
          change[j, s] <- abs(value - previous)
          if (change[j, s] <= tolerance) { converged <- TRUE; break }
        }
        previous <- value
      }
      if (!converged) stop("Marginal quadrature did not converge at group ", j, ", draw ", s,
                           "; increase max_nodes/work budget or inspect this case")
      out[j, s] <- value; used[j, s] <- nodes
    }
  }
  attr(out, "quadrature") <- list(nodes_used = used, log_change = change,
    tolerance = tolerance, max_nodes = max_nodes,
    note = "Successive-rule agreement is a numerical diagnostic, not a rigorous error bound")
  out
}

vi_global_draws_hlrm <- function(results, S) {
  list(Beta = mvtnorm::rmvnorm(S, drop(results$MuBeta), results$SigmaBeta),
       Sigma = array(unlist(lapply(seq_len(S), function(s)
         as.matrix(LaplacesDemon::rinvwishart(results$nuS, results$SS)))),
         dim = c(nrow(results$SS), nrow(results$SS), S)) |> aperm(c(3, 1, 2)),
       xi2 = rgamma(S, results$axi2, rate = results$bxi2))
}

predictive_criteria_hlrm_marginal <- function(y_list, X_list, draws, nu0,
                                               method,
                                               control = list(),
                                               return_log_likelihood = FALSE,
                                               max_output_mb = 256) {
  method <- match.arg(method, c("MCMC", "CAVI", "SVI"))
  if (!is.logical(return_log_likelihood) ||
      length(return_log_likelihood) != 1L || is.na(return_log_likelihood))
    stop("return_log_likelihood must be TRUE or FALSE")
  if (!is.list(control) || (length(control) &&
      (is.null(names(control)) || anyDuplicated(names(control)) ||
       any(!names(control) %in% c("tolerance", "max_nodes", "max_evaluations")))))
    stop("Unknown or duplicated marginal_control fields")
  ll <- log_likelihood(
    "HLRM", y_list, X_list, draws, unit = "group",
    conditioning = "group_marginal", nu0 = nu0,
    marginal_control = control, max_output_mb = max_output_mb)
  new_predictive_criteria_result(
    "HLRM", method, "new_groups",
    criteria = list(waic_group_marginal = waic_from_loglik(ll)),
    metadata = list(
      draw_source = if (method == "MCMC") "posterior_mcmc" else "variational_q",
      likelihood = "group_marginal_integrating_beta_j_and_sigma2_j",
      log_likelihood_estimand = log_likelihood_estimand(ll),
      marginal_quadrature = ll$metadata$quadrature,
      interpretation = "new-group predictive criterion; LOGO refitting is still required for out-of-sample assessment"),
    log_likelihood = if (return_log_likelihood) list(group = ll) else NULL)
}
