test_that("LRM log likelihood preserves observation and supplied-group units", {
  X <- cbind(1, c(-1, 0, 1))
  y <- c(-.5, 1, 2.5)
  draws <- list(Beta = rbind(c(1, 1), c(.5, 1.5)), sigma2 = c(1, 2))
  observation <- log_likelihood(
    "LRM", y, X, draws, unit = "observation",
    conditioning = "parameter_conditional")
  manual <- vapply(1:2, function(s)
    dnorm(y, drop(X %*% draws$Beta[s, ]), sqrt(draws$sigma2[s]), log = TRUE),
    numeric(3))
  expect_s3_class(observation, "vb_log_likelihood")
  expect_identical(observation$schema_version, "log-likelihood-1")
  expect_equal(as.matrix(observation), manual)
  expect_identical(observation$unit, "observation")
  expect_identical(observation$predictive_units$observation, 1:3)

  grouped <- log_likelihood(
    "LRM", y, X, draws, unit = "group",
    conditioning = "parameter_conditional", groups = c("a", "b", "b"))
  expect_equal(unname(as.matrix(grouped)),
               unname(rowsum(manual, c(1, 2, 2), reorder = FALSE)))
  expect_identical(grouped$predictive_units$group, c("a", "b"))
  waic <- waic_from_loglik(grouped)
  expect_identical(waic$estimand$unit, "group")
  expect_identical(waic$predictive_units, grouped$predictive_units)
})

test_that("HLRM conditional likelihood uses the same matrix for both units", {
  X <- list(matrix(1, 2, 1), matrix(1, 1, 1))
  y <- list(c(-1, 1), 2)
  draws <- list(
    Betaj = array(c(0, 2, 1, 3), dim = c(2, 2, 1)),
    sigma2j = matrix(c(1, 2, .5, 1.5), 2, 2))
  observation <- log_likelihood(
    "HLRM", y, X, draws, "observation", "parameter_conditional")
  grouped <- log_likelihood(
    "HLRM", y, X, draws, "group", "parameter_conditional")
  expect_equal(unname(as.matrix(grouped)),
    unname(rowsum(as.matrix(observation), c(1, 1, 2), reorder = FALSE)))
  expect_identical(grouped$metadata$conditioned_on,
                   c("group_coefficients", "group_variances"))
  expect_length(grouped$integrated_over, 0)
})

test_that("HLRM and CHLRM marginal routes retain integrated variables", {
  X <- list(matrix(1, 2, 1))
  y <- list(c(-.5, .5))
  hlrm_draws <- list(
    Beta = matrix(c(0, .2), 2, 1),
    Sigma = array(0, c(2, 1, 1)), xi2 = c(1, 1.2))
  hlrm <- log_likelihood(
    "HLRM", y, X, hlrm_draws, "group", "group_marginal", nu0 = 4)
  expect_equal(as.matrix(hlrm),
    loglik_group_marginal_hlrm(y, X, hlrm_draws, 4))
  expect_identical(hlrm$integrated_over,
                   c("group_coefficients", "group_variances"))
  expect_true(is.list(hlrm$metadata$quadrature))

  chlrm_draws <- list(
    Betak = array(c(-1, 1, -2, 2), c(2, 2, 1)),
    sigma2k = matrix(c(1, 2, 1.5, .8), 2, 2),
    omega = matrix(c(.4, .6, .7, .3), 2, 2, byrow = TRUE))
  chlrm <- log_likelihood(
    "CHLRM", y, X, chlrm_draws, "group", "group_marginal")
  expect_equal(as.matrix(chlrm),
               loglik_group_marginal_chlrm(y, X, chlrm_draws))
  expect_identical(chlrm$integrated_over, "component_allocation")
})

test_that("invalid estimands fail rather than silently changing likelihood", {
  X <- matrix(1, 2, 1); y <- c(0, 1)
  lrm_draws <- list(Beta = matrix(c(0, 1), 2, 1), sigma2 = c(1, 1))
  expect_error(log_likelihood("LRM", y, X, lrm_draws, "group",
                              "parameter_conditional"), "groups")
  expect_error(log_likelihood("LRM", y, X, lrm_draws, "group",
                              "group_marginal", groups = c(1, 1)),
               "no group latent")
  expect_error(log_likelihood("CHLRM", list(y), list(X), list(),
                              "observation", "group_marginal"),
               "requires unit")
  expect_error(log_likelihood("HLRM", list(y), list(X), list(),
                              "group", "group_marginal"), "nu0")
})
