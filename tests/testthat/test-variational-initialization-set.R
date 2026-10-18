test_that("one public initialization set is reusable across CAVI and SVI", {
  simulation <- .test_lrm_simulation(
    30, c(1, 2), 0.5, c(-1, 1), seed = 1811)
  model <- lrm_model(2)
  prior <- lrm_prior(c(0, 0), diag(10, 2), 4, 1)
  initial_states <- variational_initialization_set(
    model, simulation$train, prior, initialization_from_prior(),
    starts = 2, seed = 1812)

  expect_s3_class(initial_states, "vb_variational_initialization_set")
  expect_identical(initial_states$schema_version,
                   "variational-initialization-set-1")
  expect_identical(unserialize(serialize(initial_states, NULL)),
                   initial_states)
  expect_equal(length(unique(initial_states$metadata$initialization_seed)), 2)
  expect_false(any(c("MuBeta", "SigmaBeta", "rho") %in%
                   names(initial_states)))

  cavi <- fit_model(
    model, simulation$train, prior, "cavi", initial_states,
    cavi_control(
      4, 3, cavi_stopping(1e20, 1e20, min_iterations = 2, patience = 1),
      starts = 2, require_converged = FALSE),
    seed = 1813, verbose = FALSE)
  svi <- fit_model(
    model, simulation$train, prior, "svi", initial_states,
    svi_control(
      4, 5, evaluation_draws = 3,
      stopping = svi_stopping(
        1e20, 1e20, min_iterations = 2, patience = 1,
        check_every = 1, smoothing_window = 1),
      starts = 2, require_converged = FALSE),
    seed = 1814, verbose = FALSE)

  expect_identical(cavi$specifications$initialization,
                   svi$specifications$initialization)
  cavi_runs <- fit_orchestration(cavi)$summary
  svi_runs <- fit_orchestration(svi)$summary
  expect_identical(cavi_runs$initialization_seed,
                   initial_states$metadata$initialization_seed)
  expect_identical(svi_runs$initialization_seed,
                   initial_states$metadata$initialization_seed)
  expect_false(identical(cavi_runs$fit_seed, svi_runs$fit_seed))
})

test_that("initialization sets validate dimensions, counts and scaling", {
  simulation <- .test_lrm_simulation(
    24, c(1, 2), 0.5, c(-1, 1), seed = 1821)
  model <- lrm_model(2)
  prior <- lrm_prior(c(0, 0), diag(10, 2), 4, 1)
  scaling <- data_scaling()
  initial_states <- variational_initialization_set(
    model, simulation$train, prior, initialization_ols(),
    starts = 1, scaling = scaling, seed = 1822)
  control <- cavi_control(3, 2)

  fit <- fit_model(
    model, simulation$train, prior, "cavi", initial_states, control,
    scaling = scaling, seed = 1823, verbose = FALSE)
  expect_s3_class(fit, "vb_model_fit")

  expect_error(fit_model(
    model, simulation$train, prior, "cavi", initial_states,
    cavi_control(3, 2, starts = 2), scaling = scaling,
    seed = 1824, verbose = FALSE), "starts")
  expect_error(fit_model(
    model, simulation$train, prior, "cavi", initial_states, control,
    seed = 1825, verbose = FALSE), "scaling convention")

  incompatible <- .test_lrm_simulation(
    25, c(1, 2), 0.5, c(-1, 1), seed = 1826)
  expect_error(fit_model(
    model, incompatible$train, prior, "cavi", initial_states, control,
    scaling = scaling, seed = 1827, verbose = FALSE),
    "data dimensions|scaling map")
  expect_error(fit_model(
    model, simulation$train, prior, "mcmc", initial_states,
    mcmc_control(1, 3), seed = 1828, verbose = FALSE),
    "only to CAVI and SVI")
})

test_that("variational initialization sets support every model family", {
  simulations <- list(
    LRM = .test_lrm_simulation(
      24, c(1, 2), 0.5, c(-1, 1), seed = 1831),
    HLRM = .test_hlrm_simulation(
      rep(6, 4), matrix(rep(c(1, 0.5), 4), 4, 2, byrow = TRUE),
      rep(0.5, 4), seed = 1832),
    CHLRM = .test_chlrm_simulation(
      4, rep(6, 4), 2, 2,
      component_coefficients = rbind(c(0, 1), c(2, -1)),
      component_variances = c(0.5, 0.5),
      allocations = c(1, 1, 2, 2), seed = 1833))
  models <- list(LRM = lrm_model(2), HLRM = hlrm_model(2),
                 CHLRM = chlrm_model(2, 2))
  priors <- list(
    LRM = lrm_prior(c(0, 0), diag(10, 2), 4, 1),
    HLRM = hlrm_prior(c(0, 0), diag(10, 2), 4, diag(2), 4, 2, 2),
    CHLRM = chlrm_prior(
      c(0, 0), diag(10, 2), 4, diag(2), 4, 2, 2, c(1, 1)))
  initializations <- list(
    LRM = initialization_ols(), HLRM = initialization_ols(),
    CHLRM = initialization_kmeans())
  controls <- list(
    cavi = cavi_control(3, 2),
    svi = svi_control(
      3, 2, evaluation_draws = 2,
      stopping = svi_stopping(min_iterations = 1, check_every = 1,
                              smoothing_window = 1)))

  index <- 0L
  for (model_name in names(models)) {
    initial_states <- variational_initialization_set(
      models[[model_name]], simulations[[model_name]]$train,
      priors[[model_name]], initializations[[model_name]],
      starts = 1, seed = 1840 + index)
    for (engine in names(controls)) {
      index <- index + 1L
      fit <- fit_model(
        models[[model_name]], simulations[[model_name]]$train,
        priors[[model_name]], engine, initial_states, controls[[engine]],
        seed = 1850 + index, verbose = FALSE)
      expect_identical(fit$model, model_name)
      expect_identical(fit$method, toupper(engine))
      expect_identical("latent_variables" %in% names(fit),
                       identical(model_name, "CHLRM"))
    }
  }
})
