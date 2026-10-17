# Uniform result schemas shared by LRM, HLRM and CHLRM. Each model exposes only
# mathematically applicable quantities; nonapplicable fields are absent rather
# than represented by NULL or fabricated estimates.

inference_method <- function(run_settings, family = c("MCMC", "CAVI")) {
  family <- match.arg(family)
  if (family == "MCMC") return("MCMC")
  if (isTRUE(run_settings$stochastic)) "SVI" else "CAVI"
}

new_parameter_summary <- function(mean, sd, lower, upper, level = 0.95,
                                  distribution, mean_exists = NULL) {
  if (!is.numeric(level) || length(level) != 1L || !is.finite(level) ||
      level <= 0 || level >= 1) stop("level must be strictly between zero and one")
  if (!is.character(distribution) || length(distribution) != 1L ||
      !nzchar(distribution)) stop("distribution must be a nonempty string")
  template <- mean
  fill_like <- function(value) {
    if (length(value) == length(template)) {
      dim(value) <- dim(template)
      dimnames(value) <- dimnames(template)
      return(value)
    }
    stop("summary fields must have compatible dimensions")
  }
  mean <- as.numeric(mean); dim(mean) <- dim(template); dimnames(mean) <- dimnames(template)
  sd <- fill_like(as.numeric(sd))
  lower <- fill_like(as.numeric(lower))
  upper <- fill_like(as.numeric(upper))
  if (is.null(mean_exists)) mean_exists <- is.finite(mean)
  mean_exists <- fill_like(as.logical(mean_exists))
  list(mean = mean, sd = sd, lower = lower, upper = upper,
       level = level, distribution = distribution,
       mean_exists = mean_exists,
       contains_zero = lower <= 0 & 0 <= upper)
}

summarize_draws_parameter <- function(draws, margin = NULL, level = 0.95,
                                      distribution = "posterior draws") {
  if (!is.numeric(draws) || !length(draws) || any(!is.finite(draws)))
    stop("draws must be a nonempty finite numeric object")
  alpha <- (1 - level) / 2
  summarize <- function(x) c(
    mean = mean(x), sd = stats::sd(x),
    lower = stats::quantile(x, alpha, names = FALSE),
    upper = stats::quantile(x, 1 - alpha, names = FALSE))
  values <- if (is.null(margin)) summarize(draws) else apply(draws, margin, summarize)
  take <- function(name) {
    if (is.null(margin)) return(unname(values[name]))
    subs <- c(list(values, name), rep(list(TRUE), length(dim(values)) - 1L),
              list(drop = TRUE))
    do.call(`[`, subs)
  }
  new_parameter_summary(take("mean"), take("sd"), take("lower"), take("upper"),
                        level, distribution)
}

summarize_normal_parameter <- function(mean, covariance, level = 0.95) {
  mean <- as.numeric(mean)
  covariance <- as.matrix(covariance)
  if (!identical(dim(covariance), as.integer(c(length(mean), length(mean)))))
    stop("covariance dimension must match mean")
  chol_spd(covariance, "normal covariance")
  standard_deviation <- sqrt(diag(covariance))
  critical <- stats::qnorm(1 - (1 - level) / 2)
  new_parameter_summary(mean, standard_deviation,
                        mean - critical * standard_deviation,
                        mean + critical * standard_deviation,
                        level, "normal variational factor")
}

summarize_independent_normals <- function(mean, covariances, level = 0.95) {
  mean <- as.matrix(mean)
  if (!is.list(covariances) || length(covariances) != nrow(mean))
    stop("covariances must contain one matrix per row of mean")
  standard_deviation <- t(vapply(seq_len(nrow(mean)), function(i) {
    covariance <- as.matrix(covariances[[i]])
    chol_spd(covariance, paste0("normal covariance ", i))
    sqrt(diag(covariance))
  }, numeric(ncol(mean))))
  critical <- stats::qnorm(1 - (1 - level) / 2)
  new_parameter_summary(mean, standard_deviation,
                        mean - critical * standard_deviation,
                        mean + critical * standard_deviation,
                        level, "independent normal variational factors")
}

summarize_inverse_gamma_parameter <- function(shape, scale, level = 0.95) {
  shape <- as.numeric(shape); scale <- as.numeric(scale)
  if (length(shape) != length(scale) || any(!is.finite(shape)) ||
      any(!is.finite(scale)) || any(shape <= 0) || any(scale <= 0))
    stop("shape and scale must be positive finite vectors of equal length")
  mean_exists <- shape > 1
  mean <- rep(NA_real_, length(shape)); mean[mean_exists] <- scale[mean_exists] / (shape[mean_exists] - 1)
  sd <- rep(NA_real_, length(shape)); variance_exists <- shape > 2
  sd[variance_exists] <- scale[variance_exists] /
    ((shape[variance_exists] - 1) * sqrt(shape[variance_exists] - 2))
  alpha <- (1 - level) / 2
  lower <- 1 / stats::qgamma(1 - alpha, shape = shape, rate = scale)
  upper <- 1 / stats::qgamma(alpha, shape = shape, rate = scale)
  new_parameter_summary(mean, sd, lower, upper, level,
                        "inverse-gamma variational factor", mean_exists)
}

summarize_gamma_parameter <- function(shape, rate, level = 0.95) {
  shape <- as.numeric(shape); rate <- as.numeric(rate)
  if (length(shape) != length(rate) || any(!is.finite(shape)) ||
      any(!is.finite(rate)) || any(shape <= 0) || any(rate <= 0))
    stop("shape and rate must be positive finite vectors of equal length")
  alpha <- (1 - level) / 2
  new_parameter_summary(shape / rate, sqrt(shape) / rate,
                        stats::qgamma(alpha, shape, rate = rate),
                        stats::qgamma(1 - alpha, shape, rate = rate),
                        level, "gamma variational factor")
}

summarize_dirichlet_parameter <- function(alpha, level = 0.95) {
  alpha <- as.numeric(alpha)
  if (!length(alpha) || any(!is.finite(alpha)) || any(alpha <= 0))
    stop("alpha must contain positive finite values")
  total <- sum(alpha)
  if (length(alpha) == 1L) {
    return(new_parameter_summary(1, 0, 1, 1, level,
                                 "degenerate one-component simplex"))
  }
  mean <- alpha / total
  sd <- sqrt(alpha * (total - alpha) / (total^2 * (total + 1)))
  tail <- (1 - level) / 2
  lower <- vapply(seq_along(alpha), function(k)
    stats::qbeta(tail, alpha[k], total - alpha[k]), numeric(1))
  upper <- vapply(seq_along(alpha), function(k)
    stats::qbeta(1 - tail, alpha[k], total - alpha[k]), numeric(1))
  new_parameter_summary(mean, sd, lower, upper, level,
                        "Dirichlet variational factor")
}

summarize_inverse_wishart_mean <- function(df, scale, level = 0.95) {
  scale <- as.matrix(scale); p <- nrow(scale)
  chol_spd(scale, "inverse-Wishart scale")
  mean_exists <- df > p + 1
  mean <- if (mean_exists) scale / (df - p - 1) else matrix(NA_real_, p, p)
  missing_moment <- matrix(NA_real_, p, p)
  new_parameter_summary(mean, missing_moment, missing_moment, missing_moment,
                        level, "inverse-Wishart variational factor",
                        matrix(mean_exists, p, p))
}

parameter_names_for_model <- function(model) {
  switch(model,
    LRM = c("coefficients", "residual_variance"),
    HLRM = c("population_coefficients", "population_covariance",
             "variance_scale", "group_coefficients", "group_variances"),
    CHLRM = c("population_coefficients", "population_covariance",
              "variance_scale", "component_coefficients",
              "component_variances", "mixture_weights"),
    stop("unknown model"))
}

new_parameter_slots <- function(model) {
  stats::setNames(vector("list", length(parameter_names_for_model(model))),
                  parameter_names_for_model(model))
}

posterior_draw_parameter <- function(draws) {
  list(distribution = "posterior draws", draws = draws)
}

normal_variational_parameter <- function(mean, covariance,
                                         independent = FALSE) {
  list(distribution = if (independent) "independent normal variational factors" else
         "normal variational factor",
       mean = mean, covariance = covariance)
}

inverse_gamma_variational_parameter <- function(shape, scale) {
  list(distribution = "inverse-gamma variational factor",
       shape = shape, scale = scale)
}

canonical_fit_parameters <- function(model, method, native_parameters) {
  parameters <- new_parameter_slots(model)
  if (method == "MCMC") {
    if (model == "LRM") {
      parameters$coefficients <- posterior_draw_parameter(native_parameters$Beta)
      parameters$residual_variance <- posterior_draw_parameter(native_parameters$sigma2)
    } else {
      parameters$population_coefficients <- posterior_draw_parameter(native_parameters$Beta)
      parameters$population_covariance <- posterior_draw_parameter(native_parameters$Sigma)
      parameters$variance_scale <- posterior_draw_parameter(native_parameters$xi2)
      if (model == "HLRM") {
        parameters$group_coefficients <- posterior_draw_parameter(native_parameters$Betaj)
        parameters$group_variances <- posterior_draw_parameter(native_parameters$sigma2j)
      } else {
        parameters$component_coefficients <- posterior_draw_parameter(native_parameters$Betak)
        parameters$component_variances <- posterior_draw_parameter(native_parameters$sigma2k)
        parameters$mixture_weights <- posterior_draw_parameter(native_parameters$omega)
      }
    }
    return(parameters)
  }
  if (model == "LRM") {
    parameters$coefficients <- normal_variational_parameter(
      native_parameters$MuBeta, native_parameters$SigmaBeta)
    parameters$residual_variance <- inverse_gamma_variational_parameter(
      native_parameters$asigma2, native_parameters$bsigma2)
  } else {
    parameters$population_coefficients <- normal_variational_parameter(
      native_parameters$MuBeta, native_parameters$SigmaBeta)
    parameters$population_covariance <- list(
      distribution = "inverse-Wishart variational factor",
      df = native_parameters$nuS, scale = native_parameters$SS)
    parameters$variance_scale <- list(
      distribution = "gamma variational factor",
      shape = native_parameters$axi2, rate = native_parameters$bxi2)
    if (model == "HLRM") {
      parameters$group_coefficients <- normal_variational_parameter(
        native_parameters$MuBetaj, native_parameters$SigmaBetaj,
        independent = TRUE)
      parameters$group_variances <- inverse_gamma_variational_parameter(
        native_parameters$aj, native_parameters$bj)
    } else {
      parameters$component_coefficients <- normal_variational_parameter(
        native_parameters$MuBetak, native_parameters$SigmaBetak,
        independent = TRUE)
      parameters$component_variances <- inverse_gamma_variational_parameter(
        native_parameters$ak, native_parameters$bk)
      parameters$mixture_weights <- list(
        distribution = "Dirichlet variational factor",
        concentration = native_parameters$alpha)
    }
  }
  parameters
}

canonical_fit_latent_variables <- function(model, method, native_parameters) {
  if (model != "CHLRM") return(NULL)
  if (method == "MCMC") list(group_allocations = native_parameters$z) else
    list(allocation_probabilities = native_parameters$rho)
}

canonical_fit_objective <- function(method, native_parameters) {
  if (method == "MCMC") {
    return(list(log_likelihood = native_parameters$loglik,
                log_joint = native_parameters$logjoint))
  }
  list(elbo = native_parameters$elbo,
       evaluation_iterations = native_parameters$elbo_iter)
}

new_posterior_summary_result <- function(model, method, parameters,
                                         allocations = NULL,
                                         relabeling = NULL, level = 0.95) {
  models <- c("LRM", "HLRM", "CHLRM")
  methods <- c("MCMC", "CAVI", "SVI")
  if (!model %in% models) stop("unknown model")
  if (!method %in% methods) stop("unknown inference method")
  required <- parameter_names_for_model(model)
  if (!is.list(parameters) || !identical(names(parameters), required))
    stop("parameters must use the canonical slots for the selected model")
  if (any(vapply(parameters, is.null, logical(1L))))
    stop("applicable parameters cannot be NULL")
  out <- list(
    schema_version = "posterior-summary-1",
    model = model,
    method = method,
    interval = list(level = level, type = "pointwise credible"),
    parameters = parameters)
  if (!is.null(allocations)) out$allocations <- allocations
  if (!is.null(relabeling)) out$relabeling <- relabeling
  structure(out, class = c("vb_posterior_summary", "list"))
}

new_posterior_draws_result <- function(fit, sampled, draw_index, seed,
                                       estimated_storage_mb) {
  if (!inherits(fit, "vb_model_fit"))
    stop("fit must be produced by fit_model()")
  if (!is.list(sampled) || !is.list(sampled$parameters) ||
      !identical(names(sampled$parameters), parameter_names_for_model(fit$model)))
    stop("posterior draws must match the parameters of the selected model")
  if (!is.data.frame(draw_index) || !nrow(draw_index))
    stop("draw_index must be a nonempty data frame")
  out <- list(
    schema_version = "posterior-draws-1",
    model = fit$model,
    method = fit$method,
    distribution = if (fit$method == "MCMC") "posterior" else
      "variational_approximation",
    parameters = sampled$parameters,
    draw_index = draw_index,
    metadata = list(
      draw_count = nrow(draw_index), draw_axis = 1L, seed = seed,
      estimated_storage_mb = estimated_storage_mb))
  if (!is.null(sampled$latent_variables))
    out$latent_variables <- sampled$latent_variables
  structure(out, class = c("vb_posterior_draws", "list"))
}

