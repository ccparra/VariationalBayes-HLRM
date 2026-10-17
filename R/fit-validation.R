# Public typed construction facade. The existing R6 classes are execution
# backends; this is the only translation layer from semantic public names to
# their current equation-oriented notation.

.validate_fit_data <- function(model, data) {
  if (!inherits(model, "vb_model_spec"))
    stop("model must be created by lrm_model(), hlrm_model() or chlrm_model()")
  if (!is.list(data) || !all(c("X", "y") %in% names(data)))
    stop("data must be a list containing X and y")
  if (model$model == "LRM") {
    X <- as.matrix(data$X)
    y <- as.numeric(data$y)
    if (!nrow(X) || !ncol(X) || nrow(X) != length(y) ||
        any(!is.finite(X)) || any(!is.finite(y)))
      stop("LRM requires finite y and X with compatible nonzero dimensions")
    dimensions <- list(n = length(y), p = ncol(X))
  } else {
    X <- lapply(data$X, as.matrix)
    y <- lapply(data$y, as.numeric)
    dimensions <- validate_grouped_data(
      X, y, if (model$model == "CHLRM") model$component_count else NULL)
  }
  if (!is.null(model$coefficient_count) &&
      model$coefficient_count != dimensions$p)
    stop("model coefficient_count does not match the data")
  list(X = X, y = y, dimensions = dimensions)
}

.validate_prior_for_model <- function(prior, model, dimensions) {
  if (!inherits(prior, "vb_prior_spec"))
    stop("prior must be created by a model-specific prior constructor")
  if (!identical(prior$model, model$model))
    stop("prior and model must describe the same model family")
  if (length(prior$coefficient_mean) != dimensions$p)
    stop("prior coefficient dimension does not match the data")
  if (model$model == "CHLRM" &&
      length(prior$mixture_concentration) != model$component_count)
    stop("mixture_concentration must have one entry per component")
  invisible(prior)
}

.group_ols_initialization <- function(data, jitter) {
  m <- length(data$X)
  p <- ncol(data$X[[1L]])
  coefficients <- matrix(NA_real_, m, p)
  variances <- numeric(m)
  for (j in seq_len(m)) {
    precision <- crossprod(data$X[[j]]) + jitter * diag(p)
    coefficients[j, ] <- solve_spd(
      precision, crossprod(data$X[[j]], data$y[[j]]),
      paste0("regularized group XtX ", j))
    residual <- data$y[[j]] - drop(data$X[[j]] %*% coefficients[j, ])
    variances[j] <- max(mean(residual^2), jitter)
  }
  list(coefficients = coefficients, variances = variances)
}

.lrm_ols_initialization <- function(data, jitter) {
  p <- ncol(data$X)
  coefficient <- solve_spd(
    crossprod(data$X) + jitter * diag(p),
    crossprod(data$X, data$y), "regularized XtX")
  residual <- data$y - drop(data$X %*% coefficient)
  list(coefficients = as.numeric(coefficient),
       residual_variance = max(mean(residual^2), jitter))
}

.prior_to_native <- function(prior, model, dimensions,
                             grouped_ols = NULL) {
  if (model$model == "LRM") {
    return(validate_lrm_hyperparams(list(
      Beta0 = prior$coefficient_mean,
      Sigma0 = prior$coefficient_covariance,
      nu0 = prior$residual_degrees_of_freedom,
      sigma20 = prior$residual_variance_scale
    ), dimensions$p))
  }
  native <- list(
    Mu0 = prior$coefficient_mean,
    Lambda0 = prior$coefficient_covariance,
    n0 = prior$coefficient_dispersion_df,
    S0 = prior$coefficient_dispersion_scale,
    nu0 = prior$residual_degrees_of_freedom,
    a0 = prior$residual_variance_scale_shape,
    b0 = prior$residual_variance_scale_rate
  )
  if (!is.null(grouped_ols)) native$Betaj_ols <- grouped_ols$coefficients
  if (model$model == "CHLRM")
    native$alpha0 <- prior$mixture_concentration
  validate_hlrm_hyperparams(
    native, dimensions$p, dimensions$m, model$model,
    if (model$model == "CHLRM") model$component_count else NULL)
}

.control_to_native <- function(control, seed) {
  validate_rng_seed(seed)
  if (control$method == "MCMC") {
    return(list(samples = control$post_warmup_draws,
                burn_in = control$warmup,
                thin = control$thinning,
                samples_gof = control$evaluation_draws,
                seed = seed))
  }
  stopping <- control$stopping
  out <- list(
    iterations = control$max_iterations,
    epsilon = stopping$objective_tolerance,
    parameter_tolerance = stopping$parameter_tolerance,
    patience = stopping$patience,
    min_iter = stopping$min_iterations,
    elbo_every = stopping$check_every,
    elbo_window = stopping$smoothing_window,
    samples_gof = control$evaluation_draws,
    stochastic = identical(control$method, "SVI"),
    seed = seed)
  if (identical(control$method, "SVI")) {
    out$minibatch_m <- control$minibatch_size
    out$tau <- control$learning_rate$offset
    out$kappa <- control$learning_rate$exponent
  }
  out
}

.new_backend <- function(model, method, data, native_prior, native_control) {
  switch(model$model,
    LRM = if (method == "mcmc")
        mcmc_lrm$new(data$y, data$X, native_prior, native_control) else
        vi_lrm$new(data$y, data$X, native_prior, native_control),
    HLRM = if (method == "mcmc")
      mcmc_hlrm$new(data$y, data$X, native_prior, native_control) else
      vi_hlrm$new(data$y, data$X, native_prior, native_control),
    CHLRM = if (method == "mcmc")
      mcmc_chlrm$new(data$y, data$X, model$component_count,
                     native_prior, native_control) else
      vi_chlrm$new(data$y, data$X, model$component_count,
                   native_prior, native_control))
}

