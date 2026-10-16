test_that("SPD algebra agrees with base reference calculations", {
  A <- matrix(c(2, .3, .3, 1.5), 2)
  b <- c(1, -2)

  expect_equal(solve_spd(A), solve(A), tolerance = 1e-12)
  expect_equal(solve_spd(A, b), drop(solve(A, b)), tolerance = 1e-12)
  expect_equal(logdet_spd(A), as.numeric(determinant(A, logarithm = TRUE)$modulus),
               tolerance = 1e-12)
})

test_that("numerical symmetrization repairs roundoff but rejects model errors", {
  almost <- matrix(c(2, .4 + 1e-12, .4, 1), 2)
  repaired <- symmetrize_numeric(almost, "almost")
  expect_true(isSymmetric(repaired))
  expect_equal(repaired[1, 2], mean(c(almost[1, 2], almost[2, 1])))

  material <- matrix(c(2, .7, .2, 1), 2)
  expect_error(symmetrize_numeric(material, "material"), "materially asymmetric")
  expect_error(chol_spd(material, "material"), "symmetric")
  expect_error(chol_spd(matrix(c(1, 2, 2, 1), 2), "indefinite"),
               "positive definite")
})

test_that("LRM hyperparameter validation canonicalizes and checks the inverse", {
  prior <- list(Beta0 = c(0, 0), Sigma0 = diag(2), nu0 = 3, sigma20 = 2)
  checked <- validate_lrm_hyperparams(prior, 2)
  expect_equal(checked$Sigma0inv, diag(2))

  prior$Sigma0inv <- 2 * diag(2)
  expect_error(validate_lrm_hyperparams(prior, 2), "inconsistent")
  prior$Sigma0inv <- NULL
  prior$sigma20 <- 0
  expect_error(validate_lrm_hyperparams(prior, 2), "positive finite scalar")
})

test_that("hierarchical priors enforce propriety without inventing a mean", {
  p <- 3L
  prior <- list(Mu0 = numeric(p), Lambda0 = diag(p), n0 = p,
                S0 = diag(p), nu0 = 3, a0 = 2, b0 = 1)
  checked <- validate_hlrm_hyperparams(prior, p, m = 4, model = "HLRM")
  expect_equal(checked$n0, p)
  expect_equal(checked$Lambda0inv, diag(p))

  improper <- prior
  improper$n0 <- p - 1
  expect_error(validate_hlrm_hyperparams(improper, p, 4, "HLRM"),
               "greater than p - 1")

  mixture <- c(prior, list(alpha0 = c(1, 0)))
  expect_error(validate_hlrm_hyperparams(mixture, p, 4, "CHLRM", K = 2),
               "must be positive")
})

test_that("OLS empirical-Bayes prior stays proper under a perfect fit", {
  X <- cbind(1, -2:2)
  y <- 1 + 3 * X[, 2]
  prior <- hyperparams_lrm(list(X = X, y = y), "OLS")

  expect_gt(prior$sigma20, 0)
  expect_silent(chol_spd(prior$Sigma0, "Sigma0"))
})
