test_that("Gamma quadrature probability weights reproduce moments", {
  for (shape in c(.5, 2, 10)) {
    q <- gamma_quadrature_rule(32, shape)
    expect_equal(sum(exp(q$log_weights)), 1, tolerance = 1e-12)
    expect_equal(sum(exp(q$log_weights) * q$nodes), shape, tolerance = 1e-12)
    expect_equal(sum(exp(q$log_weights) * q$nodes^2), shape * (shape + 1), tolerance = 1e-12)
  }
})

test_that("HLRM marginal matches its analytic zero-covariance limit", {
  X <- list(matrix(1, 3, 1)); y <- list(c(1, 2, 3))
  draws <- list(Beta = matrix(c(1, 2), 2, 1), Sigma = array(0, c(2, 1, 1)), xi2 = c(.7, 1.2))
  ll <- loglik_group_marginal_hlrm(y, X, draws, nu0 = 4)
  expected <- vapply(1:2, function(s) {
    a <- 2; b <- 2 * draws$xi2[s]; rss <- sum((y[[1]] - draws$Beta[s, 1])^2)
    lgamma(a + 1.5) - lgamma(a) + a * log(b) - (a + 1.5) * log(b + rss / 2) - 1.5 * log(2 * pi)
  }, numeric(1))
  expect_equal(as.numeric(ll), expected, tolerance = 1e-10)
  expect_equal(dim(ll), c(1, 2))
  expect_error(loglik_group_marginal_hlrm(y, X, draws, 4, max_evaluations = 1), "max_evaluations")
  expect_error(loglik_group_marginal_hlrm(y, X, draws, 4, max_output_mb = 1e-8), "exceed")
})

test_that("HLRM marginal agrees with independent adaptive integration over precision", {
  X <- list(matrix(c(1, 1), 2, 1)); y <- list(c(.2, 1))
  draws <- list(Beta = matrix(.4), Sigma = array(.6, c(1, 1, 1)), xi2 = 1.2)
  ll <- loglik_group_marginal_hlrm(y, X, draws, nu0 = 4, tolerance = 1e-8)
  independent <- integrate(function(u) vapply(u, function(v) {
    if (v == 0) return(0)
    covariance <- .6 * tcrossprod(X[[1]]) + diag(1 / v, 2)
    exp(mvtnorm::dmvnorm(y[[1]], rep(.4, 2), covariance, log = TRUE) +
          dgamma(v, shape = 2, rate = 2.4, log = TRUE))
  }, numeric(1)), 0, Inf, rel.tol = 1e-10)
  expect_equal(as.numeric(ll), log(independent$value), tolerance = 1e-7)
  expect_lte(max(attr(ll, "quadrature")$log_change), 1e-8)
  expect_error(loglik_group_marginal_hlrm(y, X, draws, 4, tolerance = 1e-15, max_nodes = 32),
               "did not converge")
})

test_that("global q draws have the required dimensions and moments", {
  set.seed(44)
  params <- list(MuBeta = .3, SigmaBeta = matrix(.2), nuS = 10, SS = matrix(2), axi2 = 3, bxi2 = 2)
  draws <- vi_global_draws_hlrm(params, 2000)
  expect_equal(dim(draws$Beta), c(2000, 1))
  expect_equal(dim(draws$Sigma), c(2000, 1, 1))
  expect_equal(mean(draws$Beta), .3, tolerance = .04)
  expect_equal(mean(draws$Sigma), 2 / 8, tolerance = .03)
  expect_equal(mean(draws$xi2), 1.5, tolerance = .06)
})

test_that("multivariate random coefficients match independent integration and row permutation", {
  X <- list(cbind(1, c(-1, .2, 1.1))); y <- list(c(-.2, .5, 1.4))
  Sigma <- matrix(c(.4, .1, .1, .3), 2)
  draws <- list(Beta = matrix(c(.3, .5), 1), Sigma = array(Sigma, c(1, 2, 2)), xi2 = .8)
  ll <- loglik_group_marginal_hlrm(y, X, draws, 5, tolerance = 1e-8)
  truth <- integrate(function(u) vapply(u, function(v) {
    if (v == 0) return(0)
    covariance <- X[[1]] %*% Sigma %*% t(X[[1]]) + diag(1 / v, 3)
    exp(mvtnorm::dmvnorm(y[[1]], drop(X[[1]] %*% t(draws$Beta)), covariance, log = TRUE) +
          dgamma(v, 2.5, rate = 2, log = TRUE))
  }, numeric(1)), 0, Inf, rel.tol = 1e-10)
  expect_equal(as.numeric(ll), log(truth$value), tolerance = 1e-7)
  perm <- c(3, 1, 2)
  permuted <- loglik_group_marginal_hlrm(list(y[[1]][perm]), list(X[[1]][perm, ]), draws, 5)
  expect_equal(as.numeric(permuted), as.numeric(ll), tolerance = 1e-6)
})
