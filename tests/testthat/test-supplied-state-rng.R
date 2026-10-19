test_that("single-start supplied states isolate optimization RNG from constructor strategy", {
  simulation <- .test_chlrm_simulation(8, rep(12, 8), 2, 3, seed = 2901)
  model <- chlrm_model(2, 3)
  prior <- chlrm_prior(rep(0, 3), diag(10, 3), 5, diag(3), 2, 2, 2, c(1, 1))
  initial <- variational_initialization_set(model, simulation$train, prior,
    starts = 1, seed = 2902)
  same_state <- initial
  same_state$initialization <- initialization_kmeans()
  expect_identical(attr(initial, "serialized_states"), attr(same_state, "serialized_states"))
  control <- svi_control(8, 3, evaluation_draws = 2,
    stopping = svi_stopping(0, min_iterations = 1, check_every = 1, smoothing_window = 1))
  set.seed(2903); before <- rng_snapshot(); prior_before <- serialize(prior, NULL)
  run <- function(state, seed = 2904) fit_model(model, simulation$train, prior,
    "svi", state, control, seed = seed, verbose = FALSE)
  a <- run(initial); b <- run(same_state)
  expect_identical(a$parameters, b$parameters)
  expect_identical(a$diagnostics$history, b$diagnostics$history)
  expect_identical(.fit_backend(a)$rng_state, .fit_backend(b)$rng_state)
  expect_identical(rng_snapshot(), before)
  expect_identical(serialize(prior, NULL), prior_before)
  expect_false(identical(a$parameters, run(initial, 2905)$parameters))
  # The shared start must agree with the existing explicitly seeded native
  # multistart worker for the same single optimization seed.
  data <- .validate_fit_data(model, simulation$train)
  native_prior <- .prior_to_native(prior, model, data$dimensions, NULL)
  native_control <- .control_to_native(control, 2904)
  native_control$initialization <- "random_responsibilities"
  native <- .run_one_vi_start("CHLRM", data$y, data$X, 2, native_prior,
    native_control, unserialize(attr(initial, "serialized_states")[[1]]), 2904, FALSE)
  expect_identical(.fit_backend(a)$chlrm_vi_object$params, native$model$chlrm_vi_object$params)
})
