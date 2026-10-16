test_that("LRM multistart is reproducible and selects only eligible ELBOs", {
  d <- .test_lrm_simulation(20, c(1, -0.5), 1, matrix(c(-1, 1), 1, 2), seed = 801)
  h <- hyperparams_lrm(d$train, "OLS")
  args <- list(model = "LRM", method = "CAVI", y = d$train$y,
               X = d$train$X,
               hyperparams = h, run_settings = list(iterations = 10, epsilon = 1e10),
               n_starts = 3, seed = 802)
  before <- rng_snapshot()
  a <- do.call(run_vi_multistart, args)
  expect_identical(rng_snapshot(), before)
  b <- do.call(run_vi_multistart, args)

  expect_s3_class(a, "vi_multistart_fit")
  expect_identical(lapply(a$attempts, `[[`, "results"),
                   lapply(b$attempts, `[[`, "results"))
  expect_identical(a$summary[, setdiff(names(a$summary), "elapsed_seconds")],
                   b$summary[, setdiff(names(b$summary), "elapsed_seconds")])
  eligible <- a$summary$status == "complete" & a$summary$converged
  expect_equal(a$selected_start,
               which(eligible)[which.max(a$summary$final_elbo[eligible])])
  expect_identical(selected_vi_model(a)$fit_status, "fitted")
  states <- make_vi_initial_states(
    "LRM", d$train$y, d$train$X, h, args$run_settings,
    n_starts = 2, seed = 803)
  svi <- run_vi_multistart(
    "LRM", "SVI", d$train$y, d$train$X, h,
    list(iterations = 4, stochastic = TRUE, minibatch_m = 5,
         epsilon = NULL, elbo_every = 2, min_iter = 1),
    n_starts = 2, seed = 804, initial_states = states,
    require_converged = FALSE)
  expect_identical(svi$status, "selected")
  expect_true(all(svi$summary$status == "complete"))
  expect_identical(selected_vi_model(svi)$run_settings$stochastic, TRUE)
})

test_that("no unconverged start is selected unless explicitly allowed", {
  d <- .test_lrm_simulation(15, c(0, 1), 1, matrix(c(-1, 1), 1, 2), seed = 810)
  h <- hyperparams_lrm(d$train, "OLS")
  settings <- list(iterations = 2, epsilon = NULL)
  strict <- run_vi_multistart("LRM", "CAVI", d$train$y, d$train$X, h, settings,
                              n_starts = 2, seed = 811)
  expect_identical(strict$status, "no_valid_start")
  expect_true(is.na(strict$selected_start))
  expect_error(selected_vi_model(strict), "No start")

  exploratory <- run_vi_multistart("LRM", "CAVI", d$train$y, d$train$X, h,
                                   settings,
                                   n_starts = 2, seed = 811,
                                   require_converged = FALSE)
  expect_identical(exploratory$status, "selected")
  expect_false(exploratory$summary$converged[exploratory$selected_start])
})

test_that("HLRM CAVI and SVI reuse exactly the same serializable initial states", {
  d <- .test_hlrm_hierarchy(rep(5, 4), c(0, 0.5), diag(2), 1,
                               seed = 820)
  h <- hyperparams_hlrm(d$train, "OLS")
  base <- list(iterations = 10, epsilon = 1e10)
  states <- make_vi_initial_states("HLRM", d$train$y, d$train$X, h, base,
                                   n_starts = 2, seed = 821)
  restored <- unserialize(serialize(states, NULL))
  expect_identical(states, restored)

  vi <- run_vi_multistart("HLRM", "CAVI", d$train$y, d$train$X, h, base,
    n_starts = 2, seed = 822, initial_states = states)
  svi <- run_vi_multistart("HLRM", "SVI", d$train$y, d$train$X, h,
    list(iterations = 12, epsilon = 1e10, minibatch_m = 2,
         elbo_every = 2, elbo_window = 2, min_iter = 6, patience = 2),
    n_starts = 2, seed = 823, initial_states = states)
  expect_identical(vi$initial_states$states, svi$initial_states$states)
  expect_true(all(vi$summary$status == "complete"))
  expect_true(all(svi$summary$status == "complete"))
  expect_true(all(c("initialization_seed", "fit_seed", "iterations",
                    "final_elbo", "elapsed_seconds", "selected") %in%
                    names(vi$summary)))
})

test_that("CHLRM aligns every completed start to one reusable reference", {
  d <- .test_chlrm_simulation(4, 5, 2, 5, seed = 830)
  h <- .test_chlrm_native_prior(5, 2)
  settings <- list(iterations = 10, epsilon = 1e10)
  states <- make_vi_initial_states("CHLRM", d$train$y, d$train$X, h,
    settings, K = 2, n_starts = 2, seed = 831)
  vi <- run_vi_multistart("CHLRM", "CAVI", d$train$y, d$train$X, h,
    settings, K = 2, n_starts = 2, seed = 832, initial_states = states)

  expect_s3_class(vi$label_reference, "chlrm_label_reference")
  expect_true(all(vapply(vi$attempts, function(x)
    !is.null(x$aligned_results), logical(1))))
  expect_true(all(vapply(vi$attempts, function(x)
    identical(attr(x$aligned_results, "label_reference"), vi$label_reference),
    logical(1))))
  selected <- selected_vi_model(vi)
  expect_identical(selected$label_reference, vi$label_reference)
  selected$summarize_posterior(include_incidence = FALSE)
  expect_identical(selected$label_reference, vi$label_reference)

  svi <- run_vi_multistart("CHLRM", "SVI", d$train$y, d$train$X, h,
    list(iterations = 12, epsilon = 1e10, minibatch_m = 2,
         elbo_every = 2, elbo_window = 2, min_iter = 6, patience = 2),
    K = 2, n_starts = 2, seed = 833, initial_states = states,
    label_reference = vi$label_reference)
  expect_identical(svi$label_reference, vi$label_reference)
  expect_identical(svi$initial_states$states, vi$initial_states$states)
})

test_that("multistart validates state dimensions and memory before fitting", {
  d <- .test_lrm_simulation(12, c(0, 1), 1, matrix(c(-1, 1), 1, 2), seed = 840)
  h <- hyperparams_lrm(d$train, "OLS")
  settings <- list(iterations = 4, epsilon = NULL)
  states <- make_vi_initial_states("LRM", d$train$y, d$train$X, h, settings,
                                   n_starts = 2, seed = 841)
  expect_error(run_vi_multistart("LRM", "CAVI", d$train$y, d$train$X, h, settings,
    n_starts = 1, seed = 842, initial_states = states), "n_starts")
  expect_error(run_vi_multistart("LRM", "CAVI", d$train$y, d$train$X, h, settings,
    n_starts = 2, seed = 842, initial_states = states,
    max_memory_mb = 1e-12), "exceeds")
})

test_that("LRM and HLRM CAVI enforce fitted-state contracts", {
  dl <- .test_lrm_simulation(12, c(0, 1), 1, matrix(c(-1, 1), 1, 2), seed = 850)
  lrm <- vi_lrm$new(dl$train$y, dl$train$X,
                    hyperparams_lrm(dl$train, "OLS"),
                    list(iterations = 3, epsilon = NULL, seed = 851))
  expect_error(lrm$summarize_posterior(), "Call fit")
  lrm$fit(FALSE)
  expect_identical(lrm$fit_status, "fitted")
  expect_error(lrm$fit(FALSE), "already been fitted")

  dh <- .test_hlrm_hierarchy(rep(4, 3), c(0, 1), diag(2), 1,
                                seed = 852)
  hlrm <- vi_hlrm$new(dh$train$y, dh$train$X,
                      hyperparams_hlrm(dh$train, "OLS"),
                      list(iterations = 3, epsilon = NULL, seed = 853))
  expect_error(hlrm$predictive_criteria(), "Call fit")
  hlrm$fit(FALSE)
  expect_identical(hlrm$fit_status, "fitted")
  expect_error(hlrm$fit(FALSE), "already been fitted")
})
