# Small execution contracts, independent of model formulas and file I/O.
validate_integer_scalar <- function(x, name, minimum = 1L) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
      x < minimum || x != floor(x) || x > .Machine$integer.max)
    stop(name, " must be an integer >= ", minimum)
  as.integer(x)
}

validate_integer_vector <- function(x, name, length_required = NULL,
                                    minimum = 1L) {
  if (!is.numeric(x) || !length(x) || any(!is.finite(x)) ||
      any(x < minimum) || any(x != floor(x)) ||
      any(x > .Machine$integer.max) ||
      (!is.null(length_required) && length(x) != length_required)) {
    suffix <- if (is.null(length_required)) "" else
      paste0(" with length ", length_required)
    stop(name, " must contain integers >= ", minimum, suffix)
  }
  as.integer(x)
}

validate_rng_seed <- function(seed) {
  if (!is.null(seed) && (!is.numeric(seed) || length(seed) != 1L ||
      !is.finite(seed) || seed < 0 || seed != floor(seed) ||
      seed > .Machine$integer.max)) stop("seed must be a nonnegative integer")
  invisible(seed)
}

derive_independent_seeds <- function(n, seed = NULL, seeds = NULL) {
  n <- validate_integer_scalar(n, "n", 1L)
  if (!is.null(seeds)) {
    if (!is.numeric(seeds) || length(seeds) != n ||
        any(!is.finite(seeds)) || any(seeds < 0) ||
        any(seeds != floor(seeds)) || any(seeds > .Machine$integer.max) ||
        anyDuplicated(seeds))
      stop("seeds must contain one distinct nonnegative integer per run")
    return(as.integer(seeds))
  }
  validate_rng_seed(seed)
  if (is.null(seed)) stop("Provide seed or explicit seeds for reproducible runs")
  scope <- inference_rng_scope(seed)
  on.exit(scope$restore(), add = TRUE)
  as.integer(sample.int(.Machine$integer.max, n, replace = FALSE))
}

rng_snapshot <- function() {
  list(kind = RNGkind(), seed = if (exists(".Random.seed", .GlobalEnv, inherits = FALSE))
    get(".Random.seed", .GlobalEnv) else NULL)
}

restore_rng_snapshot <- function(state) {
  do.call(RNGkind, as.list(state$kind))
  if (is.null(state$seed)) {
    if (exists(".Random.seed", .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv)
  } else assign(".Random.seed", state$seed, envir = .GlobalEnv)
  invisible(NULL)
}

# Caller must register on.exit(scope$restore(), add=TRUE) immediately.
# With neither seed nor state, use R's ambient RNG state.
inference_rng_scope <- function(seed = NULL, state = NULL) {
  validate_rng_seed(seed)
  isolated <- !is.null(seed) || !is.null(state)
  previous <- rng_snapshot()
  if (!is.null(seed)) set.seed(seed)
  else if (!is.null(state)) restore_rng_snapshot(state)
  list(isolated = isolated, restore = function() {
    if (isolated) restore_rng_snapshot(previous)
    invisible(NULL)
  })
}

step_size <- function(iteration, offset = 10, exponent = 0.7) {
  if (any(!is.finite(c(iteration, offset, exponent))) ||
      any(iteration < 1) || offset < 0 || exponent <= 0.5 || exponent > 1) {
    stop("Require iteration >= 1, offset >= 0 and 0.5 < exponent <= 1")
  }
  (iteration + offset)^(-exponent)
}

vi_stopping_control <- function(settings, m) {
  fallback <- function(x, default) if (is.null(x)) default else x
  if (!is.null(settings$stochastic) && (!is.logical(settings$stochastic) ||
      length(settings$stochastic) != 1L || is.na(settings$stochastic)))
    stop("stochastic must be TRUE or FALSE")
  stochastic <- isTRUE(settings$stochastic)
  iterations <- validate_integer_scalar(settings$iterations, "iterations")
  epsilon <- if ("epsilon" %in% names(settings)) settings$epsilon else
    if (stochastic) 1e-4 else 1e-8
  if (!is.null(epsilon) && (!is.numeric(epsilon) || length(epsilon) != 1L ||
      !is.finite(epsilon) || epsilon < 0)) stop("epsilon must be NULL or nonnegative and finite")
  parameter_tolerance <- settings$parameter_tolerance
  if (!is.null(parameter_tolerance) &&
      (!is.numeric(parameter_tolerance) || length(parameter_tolerance) != 1L ||
       !is.finite(parameter_tolerance) || parameter_tolerance < 0))
    stop("parameter_tolerance must be NULL or nonnegative and finite")
  every <- validate_integer_scalar(fallback(settings$elbo_every, if (stochastic) 50L else 1L), "elbo_every")
  window <- validate_integer_scalar(fallback(settings$elbo_window, if (stochastic) 5L else 1L), "elbo_window")
  if (!stochastic && (every != 1L || window != 1L))
    stop("CAVI requires elbo_every = 1 and elbo_window = 1")
  batch <- validate_integer_scalar(fallback(settings$minibatch_m, m), "minibatch_m")
  if (batch > m) stop("minibatch_m cannot exceed the number of groups")
  tau <- fallback(settings$tau, 100)
  kappa <- fallback(settings$kappa, .7)
  if (!is.numeric(tau) || length(tau) != 1L || !is.finite(tau) || tau < 0)
    stop("tau must be a nonnegative finite scalar")
  if (!is.numeric(kappa) || length(kappa) != 1L || !is.finite(kappa) ||
      kappa <= .5 || kappa > 1) stop("kappa must satisfy 0.5 < kappa <= 1")
  list(iterations = iterations, stochastic = stochastic, epsilon = epsilon,
       parameter_tolerance = parameter_tolerance,
       every = every, window = window,
       patience = validate_integer_scalar(fallback(settings$patience, 5L), "patience"),
       min_iter = validate_integer_scalar(fallback(settings$min_iter, if (stochastic) 1000L else 1L), "min_iter"),
       minibatch_m = batch, tau = tau, kappa = kappa)
}

new_vi_monitor <- function() {
  list(history = numeric(), previous = NULL, stable = 0L,
       previous_parameters = NULL, relative_change = NA_real_,
       parameter_change = NA_real_, smoothed_elbo = NA_real_,
       converged = FALSE,
       records = list(
         iteration = integer(), objective = numeric(),
         smoothed_objective = numeric(),
         objective_relative_change = numeric(),
         parameter_relative_change = numeric(), stable_checks = integer()))
}

.flatten_vi_parameters <- function(parameters) {
  values <- unlist(parameters, recursive = TRUE, use.names = FALSE)
  if (!is.numeric(values) || !length(values) || any(!is.finite(values)))
    stop("Variational parameters must be finite numeric values")
  as.numeric(values)
}

.normal_natural_parameters <- function(mean, covariance, name) {
  precision <- solve_spd(covariance, name = paste(name, "covariance"))
  c(drop(precision %*% mean), -0.5 * as.numeric(precision))
}

.independent_normal_natural_parameters <- function(means, covariances, name) {
  unlist(lapply(seq_len(nrow(means)), function(index) {
    .normal_natural_parameters(
      means[index, ], covariances[[index]], paste(name, index))
  }), use.names = FALSE)
}

# Canonical natural coordinates make parameter stopping comparable across
# CAVI and SVI. They are diagnostics only and never replace the stored factors.
vi_convergence_parameters <- function(parameters,
                                      model = c("LRM", "HLRM", "CHLRM")) {
  model <- match.arg(model)
  population <- .normal_natural_parameters(
    parameters$MuBeta, parameters$SigmaBeta, paste(model, "population"))
  if (model == "LRM") {
    return(c(population, -parameters$asigma2 - 1, -parameters$bsigma2))
  }
  local <- if (model == "HLRM")
    .independent_normal_natural_parameters(
      parameters$MuBetaj, parameters$SigmaBetaj, "HLRM group") else
    .independent_normal_natural_parameters(
      parameters$MuBetak, parameters$SigmaBetak, "CHLRM component")
  covariance_dimension <- nrow(parameters$SS)
  covariance <- c(
    -(parameters$nuS + covariance_dimension + 1) / 2,
    -0.5 * as.numeric(parameters$SS))
  residual_variance <- if (model == "HLRM")
    c(rbind(-parameters$aj - 1, -parameters$bj)) else
    c(rbind(-parameters$ak - 1, -parameters$bk))
  variance_scale <- c(parameters$axi2 - 1, -parameters$bxi2)
  if (model == "HLRM")
    return(c(local, population, covariance, residual_variance,
             variance_scale))
  probabilities <- pmax(parameters$rho, .Machine$double.xmin)
  allocation_log_odds <- if (ncol(probabilities) == 1L) numeric() else
    as.numeric(log(probabilities[, -ncol(probabilities), drop = FALSE]) -
                 log(probabilities[, ncol(probabilities)]))
  c(allocation_log_odds, parameters$alpha - 1, local, population,
    covariance, residual_variance, variance_scale)
}

# Only scheduled ELBO evaluations belong here; a final off-grid evaluation
# is reported but must not count as an additional convergence checkpoint.
update_vi_monitor <- function(state, elbo, iteration, control,
                              parameters = NULL) {
  if (length(elbo) != 1L || !is.finite(elbo)) stop("Non-finite ELBO at iteration ", iteration)
  current_parameters <- if (is.null(parameters)) NULL else
    .flatten_vi_parameters(parameters)
  if (!is.null(current_parameters) && !is.null(state$previous_parameters)) {
    if (length(current_parameters) != length(state$previous_parameters))
      stop("Variational parameter dimensions changed during fitting")
    state$parameter_change <- sqrt(sum(
      (current_parameters - state$previous_parameters)^2)) /
      (1 + sqrt(sum(state$previous_parameters^2)))
  }
  if (!is.null(current_parameters))
    state$previous_parameters <- current_parameters
  state$history <- tail(c(state$history, elbo), control$window)
  if (length(state$history) >= control$window) {
    current <- mean(state$history)
    state$smoothed_elbo <- current
    if (!is.null(state$previous)) {
      state$relative_change <- abs(current - state$previous) / (1 + abs(state$previous))
      objective_stable <- !is.null(control$epsilon) &&
        state$relative_change < control$epsilon
      parameters_stable <- is.null(control$parameter_tolerance) ||
        (!is.na(state$parameter_change) &&
         state$parameter_change < control$parameter_tolerance)
      eligible <- iteration >= control$min_iter && objective_stable &&
        parameters_stable
      state$stable <- if (eligible) state$stable + 1L else 0L
      state$converged <- state$stable >= control$patience
    }
    state$previous <- current
  }
  state$records$iteration <- c(state$records$iteration, as.integer(iteration))
  state$records$objective <- c(state$records$objective, elbo)
  state$records$smoothed_objective <- c(
    state$records$smoothed_objective, state$smoothed_elbo)
  state$records$objective_relative_change <- c(
    state$records$objective_relative_change, state$relative_change)
  state$records$parameter_relative_change <- c(
    state$records$parameter_relative_change, state$parameter_change)
  state$records$stable_checks <- c(
    state$records$stable_checks, state$stable)
  state
}

vi_monitor_history <- function(monitor) {
  as.data.frame(monitor$records, stringsAsFactors = FALSE)
}

vi_convergence_reason <- function(control, converged) {
  if (!converged) return("max_iterations")
  if (is.null(control$parameter_tolerance)) "objective_tolerance" else
    "objective_and_parameter_tolerance"
}
