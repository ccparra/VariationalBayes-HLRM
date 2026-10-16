test_that("LRM separates latent-mean and response prediction and reuses y_rep", {
  X <- cbind(1, c(-1, 1))
  draws <- list(Beta = rbind(c(1, 2), c(3, 4), c(5, 6)), sigma2 = c(1, 2, 3))
  mean_pred <- posterior_prediction_lrm(X, draws, "MCMC", return_draws = TRUE)
  expect_equal(mean_pred$draws[[1]], X %*% t(draws$Beta))
  expect_equal(mean_pred$target, "mean")
  expect_equal(unique(mean_pred$summary$group_status), "not_applicable")

  set.seed(20); before <- rng_snapshot()
  response_pred <- local({
    scope <- inference_rng_scope(21); on.exit(scope$restore())
    posterior_prediction_lrm(X, draws, "MCMC", include_noise = TRUE,
                             return_draws = TRUE)
  })
  expect_identical(before, rng_snapshot())
  expect_equal(response_pred$target, "response")
  expect_false(isTRUE(all.equal(response_pred$draws, mean_pred$draws)))
  ppc <- pp_check_from_prediction(c(0, 1), response_pred, return_yrep = TRUE)
  expect_equal(ppc$repdata$mean, colMeans(ppc$y_rep))
  expect_equal(ppc$prediction_metadata$model, "LRM")
  expect_error(pp_check_from_prediction(c(0, 1), mean_pred), "response")
})

test_that("HLRM prediction distinguishes observed and genuinely new groups", {
  S <- 4L
  local <- list(Betaj = array(c(1, 2, 3, 4), c(S, 1, 1)),
                sigma2j = matrix(1, S, 1))
  global <- list(Beta = matrix(0, S, 1), Sigma = array(1, c(S, 1, 1)),
                 xi2 = rep(1, S))
  X <- list(matrix(1, 2, 1), matrix(1, 2, 1))
  set.seed(30)
  pred <- posterior_prediction_hlrm(X, local, global, c(1, NA), 4, "MCMC",
                                    return_draws = TRUE)
  expect_equal(pred$draws[[1]][1, ], 1:4)
  # A new group's beta is drawn once per posterior draw and shared by its rows.
  expect_equal(pred$draws[[2]][1, ], pred$draws[[2]][2, ])
  expect_equal(pred$summary$group_status, rep(c("observed", "new"), each = 2))
  only_new <- posterior_prediction_hlrm(X[1], local, global, NA, 4, "MCMC")
  expect_equal(unique(only_new$summary$group_status), "new")
  expect_error(posterior_prediction_hlrm(X, local, global, NULL, 4, "MCMC"),
               "group_ids is required")
})

test_that("CHLRM prediction integrates allocations and is label invariant draw by draw", {
  S <- 4L; K <- 2L
  draws <- list(
    Betak = array(c(-2, -1, 2, 1, -3, -2, 3, 2), c(S, K, 1)),
    sigma2k = matrix(1, S, K),
    omega = matrix(rep(c(.25, .75), S), S, K, byrow = TRUE))
  z <- matrix(c(1, 2, 1, 2), 1, S)
  X <- list(matrix(1, 2, 1), matrix(1, 2, 1))
  set.seed(40)
  original <- posterior_prediction_chlrm(X, draws, c(1, NA), "MCMC", z = z,
                                         return_draws = TRUE)
  expect_equal(original$draws[[1]][1, ], original$draws[[1]][2, ])
  expect_equal(rowSums(matrix(original$allocation$probability, 2, K, byrow = TRUE)), c(1, 1))

  switched <- draws
  switched$Betak[, 1, ] <- draws$Betak[, 2, ]
  switched$Betak[, 2, ] <- draws$Betak[, 1, ]
  switched$sigma2k <- draws$sigma2k[, 2:1, drop = FALSE]
  switched$omega <- draws$omega[, 2:1, drop = FALSE]
  set.seed(40)
  relabeled <- posterior_prediction_chlrm(X, switched, c(1, NA), "MCMC",
                                          z = 3L - z, return_draws = TRUE)
  expect_equal(relabeled$draws, original$draws)
})

test_that("all fitted class families expose the same prediction surface", {
  generators <- list(mcmc_lrm, vi_lrm, mcmc_hlrm, vi_hlrm, mcmc_chlrm, vi_chlrm)
  for (generator in generators) {
    expect_true(all(c("predict", "posterior_predict") %in% names(generator$public_methods)))
  }
})
