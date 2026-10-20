.ordered_draw_fits <- local({
  cache <- new.env(parent = emptyenv())
  function(model = "CHLRM") {
    if (exists(model, cache, inherits = FALSE)) return(cache[[model]])
    X <- cbind(1, seq(-1, 1, length.out = 8))
    data <- list(X = rep(list(X), 6), y = lapply(seq_len(6), function(j)
      drop(X %*% c(j / 3, (-1)^j)) + sin(seq_len(8)) / 10))
    spec <- switch(model, LRM = lrm_model(2), HLRM = hlrm_model(2),
                   CHLRM = chlrm_model(2, 2))
    prior <- switch(model,
      LRM = lrm_prior(c(0, 0), diag(2), 4, 1),
      HLRM = hlrm_prior(c(0, 0), diag(2), 4, diag(2), 4, 2, 2),
      CHLRM = chlrm_prior(c(0, 0), diag(2), 4, diag(2), 4, 2, 2, c(1, 1)))
    if (model == "LRM") data <- list(X = X, y = data$y[[1]])
    cache[[model]] <- fit_model(
      spec, data, prior, "mcmc",
      if (model == "CHLRM") initialization_kmeans() else initialization_ols(),
      mcmc_control(warmup = 1, post_warmup_draws = 8, chains = 2,
        stopping = mcmc_stopping(min_post_warmup_draws = 8, check_every = 8,
          rhat_max = 100, ess_bulk_min = .1, ess_tail_min = .1),
        max_memory_mb = 50), seed = 8021, verbose = FALSE)
    cache[[model]]
  }
})

test_that("ordered chunks recover exact tensors and original chain provenance", {
  for (model in c("LRM", "HLRM", "CHLRM")) {
    fit <- .ordered_draw_fits(model)
    before <- serialize(fit, NULL)
    full <- posterior_draws(fit, include_latent = model == "CHLRM")
    n <- nrow(full$draw_index)
    # The middle block crosses the boundary between stored chains.
    ranges <- list(1:3, 4:(n - 2L), (n - 1L):n)
    chunks <- lapply(ranges, function(indices) posterior_draws(
      fit, draw_indices = indices, include_latent = model == "CHLRM"))
    expect_identical(do.call(rbind, lapply(chunks, function(x)
      x$draw_index[, c("stored_draw", "chain", "iteration")])),
      full$draw_index[, c("stored_draw", "chain", "iteration")])
    for (name in names(full$parameters)) {
      combined <- do.call(rbind, lapply(chunks, function(x)
        matrix(x$parameters[[name]], nrow = x$metadata$draw_count)))
      expect_identical(unname(combined),
        unname(matrix(full$parameters[[name]], nrow = n)))
      expect_identical(full$parameters[[name]], fit$parameters[[name]]$draws)
    }
    if (model == "CHLRM")
      expect_identical(do.call(rbind, lapply(chunks, function(x)
        x$latent_variables$group_allocations)),
        t(fit$latent_variables$group_allocations))
    expect_identical(serialize(fit, NULL), before)
  }
})

test_that("parameter selection preserves axes and budgets only selected output", {
  for (model in c("LRM", "HLRM", "CHLRM")) {
    fit <- .ordered_draw_fits(model)
    names <- parameter_names_for_model(model)
    selected <- posterior_draws(fit, draw_indices = 2,
      parameters = rev(names[c(1, length(names))]))
    expect_named(selected$parameters, names[c(1, length(names))])
    expect_identical(selected$schema_version,
      if (model == "LRM") "posterior-draws-1" else "posterior-draws-2")
    if (model != "LRM")
      expect_identical(selected$metadata$parameter_selection,
                       names[c(1, length(names))])
    for (name in names(selected$parameters))
      expect_identical(selected$parameters[[name]],
        .draw_first_subset(fit$parameters[[name]]$draws, 2L))
    full <- posterior_draws(fit, draw_indices = 1:8)
    subset <- posterior_draws(fit, draw_indices = 1:8, parameters = names[1])
    budget <- mean(c(full$metadata$estimated_storage_mb,
                     subset$metadata$estimated_storage_mb))
    expect_lt(subset$metadata$estimated_storage_mb,
               full$metadata$estimated_storage_mb)
    expect_error(posterior_draws(fit, draw_indices = 1:8,
      max_output_mb = budget), "exceed max_output_mb")
    expect_identical(posterior_draws(fit, draw_indices = 1:8,
      parameters = names[1], max_output_mb = budget), subset)
    expect_error(posterior_draws(fit, draw_indices = 1:8,
      parameters = names[1], max_output_mb =
        subset$metadata$estimated_storage_mb / 2), "exceed max_output_mb")
  }
})

test_that("latent-only extraction and conversion retain exact allocations", {
  fit <- .ordered_draw_fits()
  draws <- posterior_draws(fit, draw_indices = 7:10,
    parameters = character(), include_latent = TRUE)
  expect_length(draws$parameters, 0)
  expect_identical(draws$metadata$parameter_selection, character())
  expect_equal(draws$latent_variables$group_allocations,
    t(fit$latent_variables$group_allocations[, 7:10, drop = FALSE]))
  converted <- posterior::as_draws(draws, include_latent = TRUE)
  expect_equal(converted$`group_allocations[1]`,
    draws$latent_variables$group_allocations[, 1])
  expect_equal(converted$source_stored_draw, 7:10)
  expect_equal(converted$source_iteration, draws$draw_index$iteration)
  expect_equal(converted$.chain, draws$draw_index$chain)
  expect_equal(as.data.frame(posterior::as_draws(fit, draw_indices = 7:10,
    parameters = character(), include_latent = TRUE)), as.data.frame(converted))
  expect_error(posterior::as_draws(draws), "No parameter draws selected")
  expect_error(plot(draws), "parameter must contain")
})

test_that("unselected parameter tensors are never extracted", {
  fit <- .ordered_draw_fits()
  expected <- fit$parameters$component_variances$draws[1:3, , drop = FALSE]
  # A non-subsettable sentinel would fail if all blocks were copied first.
  fit$parameters$population_covariance$draws <- new.env(parent = emptyenv())
  draws <- posterior_draws(fit, draw_indices = 1:3,
    parameters = "component_variances")
  expect_identical(draws$parameters$component_variances, expected)
  expect_named(draws$parameters, "component_variances")
})

test_that("explicit indices neither sample nor change caller RNG state", {
  fit <- .ordered_draw_fits()
  scope <- inference_rng_scope(8022)
  on.exit(scope$restore())
  before <- .Random.seed
  explicit <- posterior_draws(fit, draw_indices = c(1, 8, 9, 16))
  expect_identical(.Random.seed, before)
  expect_identical(explicit$draw_index$stored_draw, c(1L, 8L, 9L, 16L))
  # Legacy n_draws requests still use the same sorted random selection.
  set.seed(8023)
  expected <- sort(sample.int(16, 5))
  before <- .Random.seed
  legacy <- posterior_draws(fit, n_draws = 5, seed = 8023)
  expect_identical(legacy$draw_index$stored_draw, expected)
  expect_identical(.Random.seed, before)
  all <- posterior_draws(fit)
  expect_identical(all$schema_version, "posterior-draws-1")
  expect_false("parameter_selection" %in% names(all$metadata))
})

test_that("invalid selectors fail without silently reordering or thinning", {
  fit <- .ordered_draw_fits()
  for (indices in list(integer(), c(2, 1), c(1, 1), 0, -1, 17, 1.5,
    NA_real_, Inf, "1", TRUE, matrix(1)))
    expect_error(posterior_draws(fit, draw_indices = indices), "draw_indices")
  expect_error(posterior_draws(fit, n_draws = 2, draw_indices = 1:2),
    "mutually exclusive")
  for (parameters in list("unknown", c("mixture_weights", "mixture_weights"),
    NA_character_, 1, TRUE))
    expect_error(posterior_draws(fit, parameters = parameters), "parameters")
  expect_error(posterior_draws(fit, parameters = character()),
    "at least one parameter")
  expect_error(posterior_draws(.ordered_draw_fits("LRM"),
    include_latent = TRUE, parameters = character()), "only to CHLRM")
})

test_that("variational draws keep their default contract and reject selectors", {
  X <- cbind(1, seq(-1, 1, length.out = 8))
  for (method in c("cavi", "svi")) {
    fit <- fit_model(lrm_model(2), list(X = X, y = sin(seq_len(8))),
      lrm_prior(c(0, 0), diag(2), 4, 1), method, initialization_ols(),
      if (method == "cavi") cavi_control(3, evaluation_draws = 4) else
        svi_control(3, 4, evaluation_draws = 4,
          stopping = svi_stopping(min_iterations = 1)),
      seed = 8024, verbose = FALSE)
    draws <- posterior_draws(fit, n_draws = 4, seed = 8025)
    scope <- inference_rng_scope(8025)
    expected <- .sample_variational_parameters(fit, 4, FALSE)
    scope$restore()
    expect_identical(draws$parameters, expected$parameters)
    expect_identical(draws$schema_version, "posterior-draws-1")
    expect_error(posterior_draws(fit, draw_indices = 1:2), "only to stored MCMC")
    expect_error(posterior_draws(fit, parameters = "coefficients"),
      "only to stored MCMC")
  }
})
