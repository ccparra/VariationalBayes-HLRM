vi_sample_draws_chlrm <- function(chlrm_vi_object, results, S, seed = NULL,
                                  include_allocations = TRUE) {
  if (!is.numeric(S) || length(S) != 1L || !is.finite(S) ||
      S < 1 || S != floor(S) || S > .Machine$integer.max)
    stop("S must be a positive integer")
  if (!is.logical(include_allocations) || length(include_allocations) != 1L ||
      is.na(include_allocations)) stop("include_allocations must be TRUE or FALSE")
  if (!is.null(seed)) set.seed(seed)

  # Dimensions.
  m <- chlrm_vi_object$stats$m
  p <- chlrm_vi_object$stats$p
  K <- chlrm_vi_object$stats$K

  # Accept either list(params = ...) or a flat parameter list.
  params <- if (!is.null(results$params)) results$params else results
  if (!is.numeric(params$alpha) || length(params$alpha) != K ||
      any(!is.finite(params$alpha)) || any(params$alpha <= 0))
    stop("params$alpha must be a positive finite K-vector for q(omega)")

  # Required variational parameters.
  MuBetak     <- params$MuBetak                 # K x p
  SigmaBetak  <- lapply(params$SigmaBetak, as.matrix)  # K matrices of p x p
  ak          <- as.numeric(params$ak)          # K
  bk          <- as.numeric(params$bk)          # K
  rho         <- params$rho                     # m x K group responsibilities

  # Local sampling helpers.
  .rcat <- function(prob) sample.int(length(prob), size = 1L, prob = prob)

  # Preallocation.
  Betak_arr <- array(NA_real_, dim = c(S, K, p))  # S×K×p
  s2k_mat   <- matrix(NA_real_, nrow = S, ncol = K)
  if (include_allocations) {
    z_mat     <- matrix(NA_integer_, nrow = m, ncol = S)  # m×S
    Betaj     <- array(NA_real_, dim = c(S, m, p))        # S×m×p
    s2j_mat   <- matrix(NA_real_, nrow = S, ncol = m)     # S×m
  }

  # Component sampling, independent across draws.
  for (k in seq_len(K)) {
    Betak_arr[, k, ] <- mvtnorm::rmvnorm(
      n      = S,
      mean   = MuBetak[k, ],
      sigma  = SigmaBetak[[k]],
      method = "chol"
    )
    # σ_k^2 ~ IG(a_k, b_k)  (shape–rate) ⇒ σ^2 = 1 / Gamma(a_k, rate=b_k)
    s2k_mat[, k] <- 1 / stats::rgamma(n = S, shape = ak[k], rate = bk[k])
  }

  # Group allocations from rho_j and their induced local parameters.
  if (include_allocations) for (s in seq_len(S)) {
    for (j in seq_len(m)) {
      z_js <- .rcat(rho[j, ])       # z_j | y ~ Cat(rho_j)
      z_mat[j, s] <- z_js
      Betaj[s, j, ] <- Betak_arr[s, z_js, ]
      s2j_mat[s, j] <- s2k_mat[s, z_js]
    }
  }

  # Sample q(omega), not rho_j: these are population mixture weights.
  # Draw mixture weights after component/allocation draws to preserve RNG order.
  out <- list(Betak = Betak_arr, sigma2k = s2k_mat,
              omega = gtools::rdirichlet(S, params$alpha))
  if (include_allocations) {
    out <- c(list(Betaj = Betaj, sigma2j = s2j_mat, z = z_mat), out)
  }
  out
}
