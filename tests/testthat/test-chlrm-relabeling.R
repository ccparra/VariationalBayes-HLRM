relabel_fixture <- function() {
  S <- 4L; K <- 2L
  draws <- list(
    Betak = array(NA_real_, c(S, K, 1)),
    sigma2k = matrix(NA_real_, S, K),
    omega = matrix(NA_real_, S, K),
    z = matrix(NA_integer_, 3, S),
    logjoint = c(10, 4, 3, 2)
  )
  for (s in seq_len(S)) {
    if (s %% 2L) {
      draws$Betak[s, , 1] <- c(-2, 2)
      draws$sigma2k[s, ] <- c(1, 4)
      draws$omega[s, ] <- c(.3, .7)
      draws$z[, s] <- c(1, 2, 1)
    } else {
      draws$Betak[s, , 1] <- c(2, -2)
      draws$sigma2k[s, ] <- c(4, 1)
      draws$omega[s, ] <- c(.7, .3)
      draws$z[, s] <- c(2, 1, 2)
    }
  }
  draws
}

test_that("iterative CHLRM reference uses max log joint and pooled scales", {
  draws <- relabel_fixture()
  reference <- fit_label_reference_chlrm(draws)
  expect_s3_class(reference, "chlrm_label_reference")
  expect_equal(reference$initial_reference_draw, 1L)
  expect_true(reference$converged)
  eta <- cbind(as.vector(draws$Betak[, , 1]), as.vector(log(draws$sigma2k)))
  expect_equal(reference$scale, apply(eta, 2, sd))

  aligned <- align_to_reference_chlrm(draws, reference)
  expect_equal(drop(aligned$Betak[, 1, ]), rep(-2, 4))
  expect_equal(drop(aligned$Betak[, 2, ]), rep(2, 4))
  expect_equal(aligned$sigma2k[, 1], rep(1, 4))
  expect_equal(aligned$z[, 1], aligned$z[, 2])
  expect_equal(incidence_matrix_mcmc_chlrm(t(draws$z)),
               incidence_matrix_mcmc_chlrm(t(aligned$z)))
  expect_equal(align_to_reference_chlrm(aligned, reference)$Betak, aligned$Betak)
})

test_that("CAVI and SVI result objects align to an MCMC reference", {
  reference <- fit_label_reference_chlrm(relabel_fixture())
  shape <- 3
  result <- list(
    MuBetak = matrix(c(2, -2), 2, 1),
    SigmaBetak = list(matrix(.2), matrix(.1)),
    ak = c(shape, shape),
    bk = exp(c(log(4), log(1)) + digamma(shape)),
    alpha = c(7, 3),
    rho = rbind(c(.9, .1), c(.2, .8)),
    MuBeta = 0, SigmaBeta = matrix(1), nuS = 4, SS = matrix(1),
    axi2 = 2, bxi2 = 1
  )
  aligned <- align_to_reference_chlrm(result, reference)
  expect_equal(drop(aligned$MuBetak), c(-2, 2))
  expect_equal(aligned$alpha, c(3, 7))
  expect_equal(aligned$rho, result$rho[, 2:1])
  expect_equal(aligned$SigmaBetak, result$SigmaBetak[2:1])
  expect_equal(attr(aligned, "permutation"), c(2L, 1L))
})

test_that("label-reference contracts reject ambiguous or incompatible input", {
  draws <- relabel_fixture()
  compact <- draws; compact$logjoint <- NULL
  expect_error(fit_label_reference_chlrm(compact), "log_joint")
  expect_error(fit_label_reference_chlrm(draws, max_iter = 1), "integer >= 2")
  reference <- fit_label_reference_chlrm(draws)
  bad <- draws; bad$Betak <- array(1, c(4, 2, 2))
  expect_error(align_to_reference_chlrm(bad, reference), "coefficient dimensions")
})

test_that("MCMC and variational classes share one fitted label reference", {
  d <- .test_chlrm_simulation(4, 6, 2, 5, seed = 510)
  prior <- .test_chlrm_native_prior(5, 2)
  mcmc <- mcmc_chlrm$new(d$train$y, d$train$X, 2, prior,
    list(samples = 5, burn_in = 2, thin = 1, seed = 511))
  mcmc$fit(FALSE)
  mcmc_hat <- mcmc$summarize_posterior()
  expect_s3_class(mcmc$label_reference, "chlrm_label_reference")
  expect_true(mcmc_hat$relabeling$iterations >= 2)

  vi <- vi_chlrm$new(d$train$y, d$train$X, 2, prior,
    list(iterations = 6, epsilon = NULL, stochastic = FALSE, seed = 512))
  vi$fit(FALSE)
  vi_hat <- vi$summarize_posterior(include_incidence = FALSE,
                          label_reference = mcmc$label_reference)
  expect_true(vi_hat$relabeling$applied)
  expect_identical(vi$label_reference, mcmc$label_reference)
  expect_equal(dim(vi$aligned_results$MuBetak), c(2, 5))
})
