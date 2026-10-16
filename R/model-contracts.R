# Model-level mathematical contracts and SPD linear algebra.  These helpers do
# not repair invalid matrices silently: a covariance/scale/precision that is not
# symmetric positive definite is a model error, not a numerical warning.

chol_spd <- function(matrix, name = deparse(substitute(matrix)),
                     symmetry_tolerance = sqrt(.Machine$double.eps)) {
  matrix <- as.matrix(matrix)
  if (!is.numeric(matrix) || !nrow(matrix) || nrow(matrix) != ncol(matrix) ||
      any(!is.finite(matrix)))
    stop(name, " must be a finite numeric square matrix")
  if (!isSymmetric(matrix, tol = symmetry_tolerance))
    stop(name, " must be symmetric")
  tryCatch(chol(matrix), error = function(e)
    stop(name, " must be positive definite", call. = FALSE))
}

symmetrize_numeric <- function(matrix, name = deparse(substitute(matrix)),
                               tolerance = 1e-8) {
  matrix <- as.matrix(matrix)
  if (!is.numeric(matrix) || !nrow(matrix) || nrow(matrix) != ncol(matrix) ||
      any(!is.finite(matrix)))
    stop(name, " must be a finite numeric square matrix")
  discrepancy <- max(abs(matrix - t(matrix)))
  scale <- 1 + max(abs(matrix))
  if (discrepancy > tolerance * scale)
    stop(name, " is materially asymmetric and cannot be repaired by symmetrization")
  0.5 * (matrix + t(matrix))
}

solve_spd <- function(matrix, rhs = NULL,
                      name = deparse(substitute(matrix))) {
  factor <- chol_spd(matrix, name)
  if (is.null(rhs)) {
    inverse <- chol2inv(factor)
    return(symmetrize_numeric(inverse, paste0("inverse of ", name)))
  }
  rhs <- as.matrix(rhs)
  if (!is.numeric(rhs) || nrow(rhs) != nrow(factor) || any(!is.finite(rhs)))
    stop("rhs must be finite with one row per matrix dimension")
  answer <- backsolve(factor, forwardsolve(t(factor), rhs))
  if (ncol(answer) == 1L) drop(answer) else answer
}

logdet_spd <- function(matrix, name = deparse(substitute(matrix))) {
  factor <- chol_spd(matrix, name)
  2 * sum(log(diag(factor)))
}

.positive_scalar <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x <= 0)
    stop(name, " must be a positive finite scalar")
  as.numeric(x)
}

.finite_vector <- function(x, length_required, name) {
  if (!is.numeric(x) || length(x) != length_required || any(!is.finite(x)))
    stop(name, " must be a finite numeric vector of length ", length_required)
  as.numeric(x)
}

.spd_matrix <- function(x, dimension, name) {
  if (!is.matrix(x) || !is.numeric(x) ||
      !identical(dim(x), as.integer(c(dimension, dimension))))
    stop(name, " must be a numeric ", dimension, " by ", dimension, " matrix")
  chol_spd(x, name)
  unname(x)
}

validate_lrm_hyperparams <- function(hyperparams, p) {
  if (!is.list(hyperparams)) stop("hyperparams must be a list")
  required <- c("Beta0", "Sigma0", "nu0", "sigma20")
  missing <- setdiff(required, names(hyperparams))
  if (length(missing)) stop("Missing LRM hyperparameters: ", paste(missing, collapse = ", "))
  out <- hyperparams
  out$Beta0 <- .finite_vector(out$Beta0, p, "Beta0")
  out$Sigma0 <- .spd_matrix(out$Sigma0, p, "Sigma0")
  out$nu0 <- .positive_scalar(out$nu0, "nu0")
  out$sigma20 <- .positive_scalar(out$sigma20, "sigma20")
  inverse <- solve_spd(out$Sigma0, name = "Sigma0")
  if (!is.null(out$Sigma0inv)) {
    supplied <- as.matrix(out$Sigma0inv)
    if (!identical(dim(supplied), as.integer(c(p, p))) ||
        any(!is.finite(supplied)) ||
        max(abs(supplied - inverse)) > 1e-7 * (1 + max(abs(inverse))))
      stop("Sigma0inv is inconsistent with Sigma0")
  }
  out$Sigma0inv <- inverse
  out
}

validate_hlrm_hyperparams <- function(hyperparams, p, m,
                                      model = c("HLRM", "CHLRM"), K = NULL) {
  model <- match.arg(model)
  if (!is.list(hyperparams)) stop("hyperparams must be a list")
  required <- c("Mu0", "Lambda0", "n0", "S0", "nu0", "a0", "b0")
  if (model == "CHLRM") required <- c(required, "alpha0")
  missing <- setdiff(required, names(hyperparams))
  if (length(missing)) stop("Missing ", model, " hyperparameters: ",
                            paste(missing, collapse = ", "))
  out <- hyperparams
  out$Mu0 <- .finite_vector(out$Mu0, p, "Mu0")
  out$Lambda0 <- .spd_matrix(out$Lambda0, p, "Lambda0")
  out$S0 <- .spd_matrix(out$S0, p, "S0")
  out$n0 <- .positive_scalar(out$n0, "n0")
  if (out$n0 <= p - 1) stop("n0 must be greater than p - 1 for a proper inverse-Wishart prior")
  out$nu0 <- .positive_scalar(out$nu0, "nu0")
  out$a0 <- .positive_scalar(out$a0, "a0")
  out$b0 <- .positive_scalar(out$b0, "b0")
  out$Lambda0inv <- solve_spd(out$Lambda0, name = "Lambda0")
  if (!is.null(out$Betaj_ols)) {
    B <- as.matrix(out$Betaj_ols)
    if (!identical(dim(B), as.integer(c(m, p))) || any(!is.finite(B)))
      stop("Betaj_ols must be a finite m by p matrix when supplied")
    out$Betaj_ols <- B
  }
  if (model == "CHLRM") {
    if (is.null(K)) stop("K is required to validate CHLRM hyperparameters")
    out$alpha0 <- .finite_vector(out$alpha0, K, "alpha0")
    if (any(out$alpha0 <= 0)) stop("alpha0 entries must be positive")
  }
  out
}
