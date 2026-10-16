# Shared component-draw contract; dimensions are never inferred from group arrays.
validate_component_draws_chlrm <- function(draws, p = NULL) {
  beta <- draws$Betak
  dims <- dim(beta)
  if (!is.numeric(beta) || length(dims) != 3L || any(dims < 1L) ||
      (!is.null(p) && dims[3] != p) || any(!is.finite(beta)))
    stop("draws$Betak must be a finite [S,K,p] array")
  S <- dims[1]; K <- dims[2]
  for (nm in c("sigma2k", "omega")) {
    a <- draws[[nm]]
    if (!is.matrix(a) || !is.numeric(a) || !identical(dim(a), c(S, K)) ||
        any(!is.finite(a))) stop("draws$", nm, " must be a finite [S,K] matrix")
  }
  if (any(draws$sigma2k <= 0)) stop("sigma2k must be strictly positive")
  if (any(draws$omega < 0) || any(abs(rowSums(draws$omega) - 1) > 1e-10))
    stop("omega rows must be nonnegative and sum to one")
  invisible(list(S = S, K = K, p = dims[3]))
}

#' Observed-data log likelihood of complete groups in CHLRM
#'
#' Integrates the one allocation shared by all observations in each group:
#' log sum_k omega_k prod_i Normal(y_ij | x_ij beta_k, sigma2_k).
#' Does not use posterior responsibilities or sampled allocations.
#' @param y_list Response vectors, one per group.
#' @param X_list Design matrices, including an intercept if needed.
#' @param draws List containing Betak [S,K,p], sigma2k [S,K], omega [S,K].
#' @param max_output_mb Maximum size in MiB of the returned matrix (double data).
#'   This allocation guard is not a bound on total process memory.
#' @return Numeric matrix [groups, draws], suitable for waic_from_loglik().
loglik_group_marginal_chlrm <- function(y_list, X_list, draws,
                                       max_output_mb = 256) {
  info <- validate_grouped_data(X_list, y_list)
  dims <- validate_component_draws_chlrm(draws, info$p)
  S <- dims$S; K <- dims$K; beta <- draws$Betak
  if (!is.numeric(max_output_mb) || length(max_output_mb) != 1L ||
      !is.finite(max_output_mb) || max_output_mb <= 0)
    stop("max_output_mb must be a positive finite scalar")
  required_mb <- 8 * as.double(info$m) * S / 1024^2
  if (required_mb > max_output_mb)
    stop("Output exceeds max_output_mb; evaluate subsets of groups or increase the explicit budget")
  out <- matrix(NA_real_, info$m, S)
  # Only one group's residual vector is allocated per component/draw. This
  # avoids n x S and n x K x S arrays, and avoids SSE cancellation from y'y.
  for (s in seq_len(S)) {
    for (j in seq_len(info$m)) {
      by_component <- vapply(seq_len(K), function(k) {
        if (draws$omega[s, k] == 0) return(-Inf)
        mu <- drop(X_list[[j]] %*% beta[s, k, ])
        log(draws$omega[s, k]) + sum(stats::dnorm(
          y_list[[j]], mu, sqrt(draws$sigma2k[s, k]), log = TRUE))
      }, numeric(1L))
      out[j, s] <- logsumexp(by_component)
    }
  }
  out
}
