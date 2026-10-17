# One evaluation path for CHLRM MCMC, CAVI and SVI. No fitting or file I/O.

check_chlrm_evaluation_budget <- function(m, K, p, S, max_output_mb) {
  S <- validate_integer_scalar(S, "S")
  if (S < 2) stop("Evaluation requires at least two draws")
  if (!is.numeric(max_output_mb) || length(max_output_mb) != 1L ||
      !is.finite(max_output_mb) || max_output_mb <= 0)
    stop("max_output_mb must be a positive finite scalar")
  # Lower bound for the main double arrays, not a process-RAM guarantee.
  if (8 * as.double(S) * (m + K * (p + 2) + 6) / 1024^2 > max_output_mb)
    stop("Evaluation arrays exceed max_output_mb; reduce S or increase the explicit budget")
  invisible(TRUE)
}

prepare_chlrm_mcmc_evaluation <- function(chains, stats, n_draws,
                                           configured_draws = NULL,
                                           seed = NULL,
                                           max_output_mb = 256) {
  available <- dim(chains$Betak)[1L]
  n_draws <- prediction_draw_count(n_draws, configured_draws, available)
  check_chlrm_evaluation_budget(stats$m, stats$K, stats$p, n_draws,
                                max_output_mb)
  scope <- inference_rng_scope(seed); on.exit(scope$restore(), add = TRUE)
  indices <- prediction_draw_indices(available, n_draws)
  list(
    draws = list(
      Betak = chains$Betak[indices, , , drop = FALSE],
      sigma2k = chains$sigma2k[indices, , drop = FALSE],
      omega = chains$omega[indices, , drop = FALSE]),
    z = chains$z[, indices, drop = FALSE],
    indices = indices,
    n_draws = n_draws)
}

prepare_chlrm_vi_evaluation <- function(object, results, n_draws,
                                         configured_draws = NULL,
                                         seed = NULL,
                                         max_output_mb = 256) {
  n_draws <- prediction_draw_count(n_draws, configured_draws)
  stats <- object$stats
  check_chlrm_evaluation_budget(stats$m, stats$K, stats$p, n_draws,
                                max_output_mb)
  scope <- inference_rng_scope(seed); on.exit(scope$restore(), add = TRUE)
  list(
    draws = vi_sample_draws_chlrm(
      object, results, n_draws, include_allocations = FALSE),
    rho = results$rho,
    n_draws = n_draws)
}

# Compact draws plus exactly one of z [m,S] or rho [m,K]. These define
# prediction for observed groups, distinct from marginal new-group scoring.
validate_allocations_chlrm <- function(z, rho, m, S, K) {
  if (is.null(z) == is.null(rho)) stop("Supply exactly one of z or rho")
  if (!is.null(z) && (!is.matrix(z) || !is.numeric(z) ||
      !identical(dim(z), as.integer(c(m, S))) || any(!is.finite(z)) ||
      any(z != floor(z) | z < 1 | z > K))) stop("z must be an integer-valued [m,S] matrix")
  if (!is.null(rho) && (!is.matrix(rho) || !is.numeric(rho) ||
      !identical(dim(rho), as.integer(c(m, K))) || any(!is.finite(rho)) ||
      any(rho < 0) || any(abs(rowSums(rho) - 1) > 1e-10)))
    stop("rho must be a nonnegative [m,K] matrix with unit row sums")
  invisible(TRUE)
}

fitted_mean_chlrm <- function(X_list, draws = NULL, z = NULL, rho = NULL,
                              component_mean = NULL) {
  m <- length(X_list); p <- ncol(X_list[[1L]])
  if (!is.null(rho)) {
    if (is.null(component_mean) && is.null(draws))
      stop("rho requires component_mean or component draws")
    if (!is.null(draws)) {
      dims <- validate_component_draws_chlrm(draws, p)
      K <- dims$K
    } else {
      component_mean <- as.matrix(component_mean)
      K <- nrow(component_mean)
      dims <- list(S = 1L, K = K, p = p)
    }
  } else {
    dims <- validate_component_draws_chlrm(draws, p)
    K <- dims$K
  }
  S <- dims$S
  validate_allocations_chlrm(z, rho, m, S, K)
  if (!is.null(component_mean) && (!is.matrix(component_mean) ||
      !is.numeric(component_mean) ||
      !identical(dim(component_mean), as.integer(c(K, p))) ||
      any(!is.finite(component_mean))))
    stop("component_mean must be a finite [K,p] matrix")
  mean_beta <- component_mean
  if (!is.null(rho) && is.null(mean_beta))
    mean_beta <- matrix(apply(draws$Betak, c(2, 3), mean), K, p)
  lapply(seq_len(m), function(j) {
    if (!is.null(rho)) b <- drop(rho[j, , drop = FALSE] %*% mean_beta)
    else {
      b <- numeric(p)
      for (s in seq_len(S)) b <- b + draws$Betak[s, z[j, s], ] / S
    }
    drop(X_list[[j]] %*% b)
  })
}

# Descriptive PPC for observed groups. A single dataset is held at a time.
# The binomial MCSE is nominal IID only; MCMC autocorrelation is not ignored
# by reporting it as an effective-sample-size corrected uncertainty estimate.
pp_check_chlrm <- function(y_list, X_list, draws, z = NULL, rho = NULL,
                            seed = NULL, return_yrep = FALSE,
                            max_output_mb = 256,
                            discrepancies = ppc_discrepancies()) {
  info <- validate_grouped_data(X_list, y_list)
  dims <- validate_component_draws_chlrm(draws, info$p)
  validate_allocations_chlrm(z, rho, info$m, dims$S, dims$K)
  simulate_draw <- function(s) lapply(seq_len(info$m), function(j) {
    k <- if (is.null(z)) sample.int(dims$K, 1, prob = rho[j, ]) else z[j, s]
    gaussian_response_draw(X_list[[j]], draws$Betak[s, k, ], draws$sigma2k[s, k])
  })
  pp_check_stream(y_list, dims$S, simulate_draw, "observed_groups_conditional_allocations",
                  seed, return_yrep, max_output_mb, X_list, discrepancies)
}

evaluate_in_sample_fit_chlrm <- function(y_list, X_list, draws, method,
                                         z = NULL, rho = NULL,
                                         component_mean = NULL,
                                         return_details = FALSE) {
  info <- validate_grouped_data(X_list, y_list)
  method <- match.arg(method, c("MCMC", "CAVI", "SVI"))
  fitted <- fitted_mean_chlrm(X_list, draws, z, rho, component_mean)
  residuals <- Map(`-`, y_list, fitted)
  new_in_sample_fit_result(
    "CHLRM", method,
    mean_squared_error(unlist(y_list), unlist(fitted)),
    coefficient_of_determination(unlist(y_list), unlist(fitted)),
    metadata = list(groups = info$m, observations = sum(info$nj),
      estimate = "posterior mean conditional on observed-group allocations",
      interpretation = "descriptive in-sample fit; not new-group prediction"),
    fitted_values = if (return_details) fitted else NULL,
    residuals = if (return_details) residuals else NULL)
}

predictive_criteria_chlrm <- function(y_list, X_list, draws, method,
                                      return_log_likelihood = FALSE,
                                      max_output_mb = 256) {
  info <- validate_grouped_data(X_list, y_list)
  dims <- validate_component_draws_chlrm(draws, info$p)
  method <- match.arg(method, c("MCMC", "CAVI", "SVI"))
  if (!is.logical(return_log_likelihood) ||
      length(return_log_likelihood) != 1L || is.na(return_log_likelihood))
    stop("return_log_likelihood must be TRUE or FALSE")
  check_chlrm_evaluation_budget(info$m, dims$K, info$p, dims$S,
                                max_output_mb)
  ll <- log_likelihood(
    "CHLRM", y_list, X_list, draws, unit = "group",
    conditioning = "group_marginal", max_output_mb = max_output_mb)
  criteria <- list(waic_group_marginal = waic_from_loglik(ll))
  new_predictive_criteria_result("CHLRM", method, "new_groups", criteria, list(
         draw_source = if (method == "MCMC") "posterior_mcmc" else "variational_q",
         likelihood = "group_marginal", unit = "group", draws = dims$S,
         log_likelihood_estimand = log_likelihood_estimand(ll),
         groups = info$m, observations = sum(info$nj),
         interpretation = "new-group allocation-marginal criterion; LOGO refitting is still required for out-of-sample assessment",
         approximation = if (method == "MCMC") "Monte Carlo; convergence must be checked separately" else
           "WAIC under q, not an exact-posterior criterion"),
    log_likelihood = if (return_log_likelihood) list(group = ll) else NULL)
}

posterior_predictive_check_chlrm <- function(y_list, X_list, draws, method,
                                             z = NULL, rho = NULL,
                                             seed = NULL,
                                             return_yrep = FALSE,
                                             max_output_mb = 256,
                                             discrepancies = ppc_discrepancies()) {
  method <- match.arg(method, c("MCMC", "CAVI", "SVI"))
  check <- pp_check_chlrm(y_list, X_list, draws, z, rho, seed,
                          return_yrep, max_output_mb, discrepancies)
  new_ppc_result(
    "CHLRM", method, "observed_groups", check,
    metadata = list(conditioning = "observed_group_allocations",
      interpretation = "posterior model adequacy for observed groups; not held-out predictive performance"))
}
