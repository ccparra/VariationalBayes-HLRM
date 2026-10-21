.report_view_fixture <- function(chains = 2L, scaling = NULL) {
  X <- cbind(1, seq(-1, 1, length.out = 6))
  data <- list(X = rep(list(X), 4), y = lapply(1:4, function(j)
    drop(X %*% c(j / 3, (-1)^j)) + sin(1:6) / 10))
  fit <- fit_model(chlrm_model(2, 2), data,
    chlrm_prior(c(0, 0), diag(2), 4, diag(2), 4, 2, 2, c(1, 1)),
    "mcmc", initialization_kmeans(),
    mcmc_control(warmup = 2L, post_warmup_draws = 8L, chains = chains,
      stopping = if (chains > 1L) mcmc_stopping(min_post_warmup_draws = 8L,
        check_every = 8L, rhat_max = 100, ess_bulk_min = .1, ess_tail_min = .1) else NULL),
    scaling = scaling, seed = 38401, verbose = FALSE)
  list(fit = fit, data = data)
}

test_that("compact CHLRM views preserve public reports and caller/source state", {
  for (chains in c(1L, 2L)) for (scaled in c(FALSE, TRUE)) {
    x <- .report_view_fixture(chains, if (scaled) data_scaling() else NULL)
    fit <- x$fit; before <- .typed_mcmc_checkpoint(fit)
    set.seed(38402); rng <- .Random.seed
    view <- mcmc_reporting_view(fit)
    expect_s3_class(view, "vb_mcmc_reporting_view")
    expect_identical(mcmc_reporting_view(view), view)
    expect_identical(view$parameters, fit$parameters)
    expect_identical(view$latent_variables, fit$latent_variables)
    expect_identical(view$objective, fit$objective)
    expect_identical(.Random.seed, rng)
    operations <- list(
      summary = function(z) summarize_posterior(z),
      components = components, population = population_parameters,
      probabilities = responsibilities, classification = classify_groups,
      uncertainty = classification_uncertainty, clustering = co_clustering,
      draws = function(z) posterior_draws(z, n_draws = 4, seed = 38403, include_latent = TRUE),
      prediction = function(z) predict(z, x$data$X, n_draws = 4, seed = 38403),
      posterior_prediction = function(z) posterior_predict(z, x$data$X, n_draws = 4, seed = 38403),
      score = function(z) predictive_performance(z, x$data,
        predictive_target = "observed_groups", predictive_unit = "observation",
        group_ids = 1:4, n_draws = 4, seed = 38403),
      likelihood = function(z) predictive_log_likelihood(z, x$data,
        predictive_target = "observed_groups", predictive_unit = "observation",
        group_ids = 1:4, n_draws = 4, seed = 38403))
    if (chains > 1L) operations$orchestration <- fit_orchestration
    for (name in names(operations)) expect_identical(operations[[name]](view),
      operations[[name]](fit), info = name)
    expect_identical(.typed_mcmc_checkpoint(fit), before)
    expect_identical(.Random.seed, rng)
    backends <- c(list(.fit_backend(view)), attr(view, "backends", exact = TRUE))
    for (backend in backends) {
      expect_null(backend$chains$Betaj); expect_null(backend$chains$sigma2j)
      expect_null(backend$relabeled_chains$Betaj); expect_null(backend$relabeled_chains$sigma2j)
    }
    orchestration <- attr(view, "orchestration", exact = TRUE)
    if (!is.null(orchestration)) for (draws in c(orchestration$chains, orchestration$aligned_chains)) {
      expect_null(draws$Betaj); expect_null(draws$sigma2j)
    }
  }
})

test_that("view backends are independent and cannot publish incomplete checkpoints", {
  fit <- .report_view_fixture()$fit
  view <- mcmc_reporting_view(fit)
  original <- .typed_mcmc_checkpoint(fit)
  backend <- .fit_backend(view)
  value <- .fit_backend(fit)$chains$Betak[1L]
  backend$chains$Betak[1L] <- value + 100
  expect_identical(.fit_backend(fit)$chains$Betak[1L], value)
  expect_identical(attr(view, "backends")[[1L]]$chains$Betak,
    attr(fit, "backends")[[1L]]$chains$Betak)
  expect_identical(.typed_mcmc_checkpoint(fit), original)
  for (candidate in list(view, unserialize(serialize(view, NULL)))) {
    expect_error(resume_mcmc(candidate, 2L), "reporting view")
    for (store in list(NULL, tempfile("view-store-"))) {
      path <- tempfile("view-checkpoint-")
      expect_error(save_mcmc_checkpoint(candidate, path, store = store), "reporting view")
      expect_false(file.exists(path))
      if (!is.null(store)) expect_false(dir.exists(store))
    }
    expect_error(.typed_mcmc_checkpoint(candidate), "reporting view")
  }
  view <- mcmc_reporting_view(fit) # Discard the deliberately mutated test view.
  reference <- fit_orchestration(fit)$label_reference
  reference$eta <- reference$eta[2:1, , drop = FALSE]
  aligned <- relabel_fit(view, reference)
  expect_s3_class(aligned, "vb_mcmc_reporting_view")
  expect_identical(aligned$parameters, relabel_fit(fit, reference)$parameters)
  expect_error(resume_mcmc(aligned, 2L), "reporting view")
  path <- tempfile(); save_mcmc_checkpoint(fit, path)
  restored <- read_mcmc_checkpoint(path)
  expect_identical(.typed_mcmc_checkpoint(restored), original)
  resume <- function(z) withCallingHandlers(resume_mcmc(z, 2L), warning = function(w) {
    if (identical(trimws(conditionMessage(w)), "The ESS has been capped to avoid unstable estimates."))
      invokeRestart("muffleWarning")
  })
  a <- resume(fit); b <- resume(restored)
  expect_identical(a$parameters, b$parameters)
  expect_identical(a$objective, b$objective)
  expect_identical(.fit_backend(a)$rng_state, .fit_backend(b)$rng_state)
  expect_identical(.fit_backend(a)$run_settings, .fit_backend(b)$run_settings)
})

test_that("view construction rejects unsupported fits and mutable metadata", {
  expect_error(mcmc_reporting_view(list()), "CHLRM MCMC")
  fit <- .report_view_fixture(1L)$fit
  other <- fit; other$model <- "HLRM"
  expect_error(mcmc_reporting_view(other), "CHLRM MCMC")
  other <- fit; other$method <- "CAVI"
  expect_error(mcmc_reporting_view(other), "CHLRM MCMC")
  backend <- .fit_backend(fit)
  backend$fit_metadata$mutable <- new.env(parent = emptyenv())
  expect_error(mcmc_reporting_view(fit), "plain atomic/list")
})
