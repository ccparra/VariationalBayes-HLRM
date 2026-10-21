# Label alignment for finite-mixture CHLRM output.  The reference is learned
# from component locations eta_k = (beta_k, log(sigma2_k)) and can then be
# reused unchanged for MCMC, CAVI and SVI fits.

.chlrm_component_locations <- function(draws) {
  dims <- validate_component_draws_chlrm(draws)
  eta <- array(NA_real_, c(dims$S, dims$K, dims$p + 1L))
  eta[, , seq_len(dims$p)] <- draws$Betak
  eta[, , dims$p + 1L] <- log(draws$sigma2k)
  eta
}

.chlrm_assignment <- function(current, target, scale) {
  K <- nrow(target)
  cost <- outer(seq_len(K), seq_len(K), Vectorize(function(k, r)
    sum(((current[k, ] - target[r, ]) / scale)^2)))
  as.integer(clue::solve_LSAP(cost)) # current component -> reference component
}

.validate_chlrm_label_reference <- function(reference, K = NULL, p = NULL) {
  if (!inherits(reference, "chlrm_label_reference") ||
      !is.matrix(reference$eta) || !is.numeric(reference$eta) ||
      any(!is.finite(reference$eta)) || !is.numeric(reference$scale) ||
      length(reference$scale) != ncol(reference$eta) ||
      any(!is.finite(reference$scale)) || any(reference$scale <= 0))
    stop("reference must be a valid chlrm_label_reference object")
  if (!is.null(K) && nrow(reference$eta) != K)
    stop("reference and object have different numbers of components")
  if (!is.null(p) && ncol(reference$eta) != p + 1L)
    stop("reference and object have different coefficient dimensions")
  invisible(reference)
}

# Fit the iterative reference specified by the computational protocol:
# initial draw (normally maximum log joint), pooled coordinate scales, Hungarian
# assignments, then reference updates until all assignments stop changing.
fit_label_reference_chlrm <- function(draws, log_joint = NULL,
                                      reference_draw = NULL, max_iter = 100L) {
  if (!requireNamespace("clue", quietly = TRUE))
    stop("Package 'clue' is required for Hungarian relabeling")
  dims <- validate_component_draws_chlrm(draws)
  if (!is.numeric(max_iter) || length(max_iter) != 1L || !is.finite(max_iter) ||
      max_iter < 2L || max_iter != floor(max_iter))
    stop("max_iter must be an integer >= 2")
  if (is.null(log_joint)) log_joint <- draws$logjoint
  if (is.null(reference_draw)) {
    if (!is.numeric(log_joint) || length(log_joint) != dims$S ||
        any(!is.finite(log_joint)))
      stop("Provide finite log_joint or an explicit reference_draw")
    reference_draw <- which.max(log_joint)
  }
  if (!is.numeric(reference_draw) || length(reference_draw) != 1L ||
      !is.finite(reference_draw) || reference_draw != floor(reference_draw) ||
      reference_draw < 1L || reference_draw > dims$S)
    stop("reference_draw must identify one stored draw")

  eta <- .chlrm_component_locations(draws)
  pooled <- matrix(eta, nrow = dims$S * dims$K, ncol = dims$p + 1L)
  scale <- apply(pooled, 2L, stats::sd)
  scale[!is.finite(scale) | scale < sqrt(.Machine$double.eps)] <- 1
  target <- matrix(eta[reference_draw, , ], dims$K, dims$p + 1L)
  previous <- NULL
  converged <- FALSE

  for (iteration in seq_len(as.integer(max_iter))) {
    permutations <- matrix(NA_integer_, dims$S, dims$K)
    aligned_eta <- array(NA_real_, dim(eta))
    for (s in seq_len(dims$S)) {
      current <- matrix(eta[s, , ], dims$K, dims$p + 1L)
      current_to_reference <- .chlrm_assignment(current, target, scale)
      permutations[s, ] <- current_to_reference
      reference_to_current <- match(seq_len(dims$K), current_to_reference)
      aligned_eta[s, , ] <- current[reference_to_current, , drop = FALSE]
    }
    if (!is.null(previous) && identical(permutations, previous)) {
      converged <- TRUE
      break
    }
    previous <- permutations
    target <- matrix(apply(aligned_eta, c(2L, 3L), mean),
                     dims$K, dims$p + 1L)
  }
  if (!converged)
    stop("CHLRM label alignment did not stabilize within max_iter")

  structure(list(
    eta = target,
    scale = scale,
    K = dims$K,
    p = dims$p,
    initial_reference_draw = as.integer(reference_draw),
    iterations = iteration,
    converged = TRUE,
    permutations = permutations,
    algorithm = "iterative Hungarian; pooled coordinate scales"
  ), class = "chlrm_label_reference")
}

# Bound temporary logical/double vectors while checking every retained label.
.valid_chlrm_allocation_labels <- function(z, S, K) {
  if (!is.matrix(z) || ncol(z) != S) return(FALSE)
  if (length(z)) for (first in seq.int(1, length(z), by = 65536L)) {
    values <- z[seq.int(first, min(length(z), first + 65535L))]
    if (any(!is.finite(values)) || any(values != floor(values)) ||
        any(values < 1L | values > K)) return(FALSE)
  }
  TRUE
}

.align_component_draws_chlrm <- function(draws, reference) {
  dims <- validate_component_draws_chlrm(draws)
  .validate_chlrm_label_reference(reference, dims$K, dims$p)
  if (!is.null(draws$z) && !.valid_chlrm_allocation_labels(draws$z, dims$S, dims$K))
    stop("draws$z must be an integer [groups,S] matrix with labels in 1:K")
  eta <- .chlrm_component_locations(draws)
  out <- draws
  permutations <- matrix(NA_integer_, dims$S, dims$K)
  for (s in seq_len(dims$S)) {
    current <- matrix(eta[s, , ], dims$K, dims$p + 1L)
    current_to_reference <- .chlrm_assignment(current, reference$eta, reference$scale)
    permutations[s, ] <- current_to_reference
    reference_to_current <- match(seq_len(dims$K), current_to_reference)
    out$Betak[s, , ] <- draws$Betak[s, reference_to_current, ]
    out$sigma2k[s, ] <- draws$sigma2k[s, reference_to_current]
    out$omega[s, ] <- draws$omega[s, reference_to_current]
    if (!is.null(draws$z)) out$z[, s] <- current_to_reference[draws$z[, s]]
  }
  attr(out, "label_reference") <- reference
  attr(out, "permutations") <- permutations
  out
}

.align_variational_result_chlrm <- function(result, reference) {
  required <- c("MuBetak", "SigmaBetak", "ak", "bk", "alpha", "rho")
  if (!is.list(result) || any(!required %in% names(result)))
    stop("Variational CHLRM result lacks component parameters")
  if (!is.matrix(result$MuBetak) || !is.numeric(result$MuBetak) ||
      any(dim(result$MuBetak) < 1L) || any(!is.finite(result$MuBetak)))
    stop("MuBetak must be a finite [K,p] matrix")
  K <- nrow(result$MuBetak); p <- ncol(result$MuBetak)
  valid_covariance <- function(x) is.matrix(x) && is.numeric(x) &&
    identical(dim(x), as.integer(c(p, p))) && all(is.finite(x)) && isSymmetric(x)
  if (length(result$SigmaBetak) != K ||
      !all(vapply(result$SigmaBetak, valid_covariance, logical(1))) ||
      length(result$ak) != K || length(result$bk) != K ||
      length(result$alpha) != K || !is.matrix(result$rho) || ncol(result$rho) != K ||
      any(!is.finite(result$ak)) || any(!is.finite(result$bk)) ||
      any(result$ak <= 0) || any(result$bk <= 0) ||
      any(!is.finite(result$alpha)) || any(result$alpha <= 0) ||
      any(!is.finite(result$rho)) || any(result$rho < 0) ||
      any(abs(rowSums(result$rho) - 1) > 1e-8))
    stop("Invalid variational CHLRM component parameters")
  .validate_chlrm_label_reference(reference, K, p)
  # E_q[log sigma_k^2] exists for every positive inverse-gamma shape/rate,
  # including empty components whose arithmetic variance mean may not exist.
  eta <- cbind(result$MuBetak, log(result$bk) - digamma(result$ak))
  current_to_reference <- .chlrm_assignment(eta, reference$eta, reference$scale)
  reference_to_current <- match(seq_len(K), current_to_reference)
  out <- result
  out$MuBetak <- result$MuBetak[reference_to_current, , drop = FALSE]
  out$SigmaBetak <- result$SigmaBetak[reference_to_current]
  out$ak <- result$ak[reference_to_current]
  out$bk <- result$bk[reference_to_current]
  out$alpha <- result$alpha[reference_to_current]
  out$rho <- result$rho[, reference_to_current, drop = FALSE]
  attr(out, "label_reference") <- reference
  attr(out, "permutation") <- current_to_reference
  out
}

align_to_reference_chlrm <- function(object, reference) {
  if (!is.list(object)) stop("object must be a CHLRM draw or variational-result list")
  if (!is.null(object$Betak)) return(.align_component_draws_chlrm(object, reference))
  if (!is.null(object$MuBetak)) return(.align_variational_result_chlrm(object, reference))
  stop("Cannot identify CHLRM component fields in object")
}
