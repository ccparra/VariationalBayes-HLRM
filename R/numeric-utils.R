logmeanexp <- function(x) {
  if (!length(x)) return(-Inf)
  maximum <- max(x)
  if (!is.finite(maximum)) return(maximum)
  maximum + log(mean(exp(x - maximum)))
}

logsumexp <- function(x) {
  if (!length(x)) return(-Inf)
  maximum <- max(x)
  if (!is.finite(maximum)) return(maximum)
  maximum + log(sum(exp(x - maximum)))
}

log_multigamma <- function(a, p) {
  (p * (p - 1) / 4) * log(pi) +
    sum(lgamma(a + (1 - seq_len(p)) / 2))
}

expected_logdet_inverse_wishart <- function(nu, scale) {
  dimension <- nrow(scale)
  -sum(digamma((nu + 1 - seq_len(dimension)) / 2)) -
    dimension * log(2) + logdet_spd(scale)
}

.log_2pi <- log(2 * pi)

expected_log_inverse_gamma <- function(shape, scale) {
  log(scale) - digamma(shape)
}

expected_inverse_inverse_gamma <- function(shape, scale) shape / scale

expected_log_dirichlet <- function(alpha) {
  digamma(alpha) - digamma(sum(alpha))
}
