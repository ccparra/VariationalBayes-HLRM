# LRM is one Gaussian regression block; HLRM has observed-group blocks.
# These criteria condition on local coefficients/variances, NOT new-group effects.
validate_gaussian_draws <- function(draws, m, p) {
  d <- dim(draws$Betaj)
  if (!is.numeric(draws$Betaj) || length(d) != 3L || d[1] < 2L ||
      d[2] != m || d[3] != p || any(!is.finite(draws$Betaj)))
    stop("Betaj must be a finite [S,m,p] array with S >= 2")
  if (!is.matrix(draws$sigma2j) || !is.numeric(draws$sigma2j) ||
      !identical(dim(draws$sigma2j), d[1:2]) ||
      any(!is.finite(draws$sigma2j)) || any(draws$sigma2j <= 0))
    stop("sigma2j must be a positive finite [S,m] matrix")
  invisible(d[1])
}

check_gaussian_evaluation_budget <- function(nj, p, S, max_output_mb) {
  # Draws, group likelihood, one group's observation likelihood, PPC summaries
  # and pointwise summaries. Other R temporaries/process memory are not bounded.
  doubles <- as.double(S) * (length(nj) * (p + 2) + max(nj) + 6) + 3 * sum(nj)
  validate_evaluation_size(S, max_output_mb, doubles)
}

prepare_gaussian_mcmc_evaluation <- function(model, chains, stats, n_draws,
                                             configured_draws = NULL,
                                             seed = NULL,
                                             max_output_mb = 256) {
  model <- match.arg(model, c("LRM", "HLRM"))
  available <- if (model == "LRM") nrow(chains$Beta) else dim(chains$Betaj)[1L]
  n_draws <- prediction_draw_count(n_draws, configured_draws, available)
  group_sizes <- if (model == "LRM") length(stats$y) else lengths(stats$y)
  check_gaussian_evaluation_budget(group_sizes, stats$p, n_draws,
                                   max_output_mb)
  scope <- inference_rng_scope(seed); on.exit(scope$restore(), add = TRUE)
  indices <- prediction_draw_indices(available, n_draws)
  draws <- if (model == "LRM") list(
    Betaj = array(chains$Beta[indices, , drop = FALSE],
                  c(n_draws, 1L, stats$p)),
    sigma2j = matrix(chains$sigma2[indices], n_draws, 1L)) else list(
    Betaj = chains$Betaj[indices, , , drop = FALSE],
    sigma2j = chains$sigma2j[indices, , drop = FALSE])
  list(draws = draws, indices = indices, n_draws = n_draws)
}

prepare_gaussian_vi_evaluation <- function(model, object, results,
                                           n_draws,
                                           configured_draws = NULL,
                                           seed = NULL,
                                           max_output_mb = 256) {
  model <- match.arg(model, c("LRM", "HLRM"))
  n_draws <- prediction_draw_count(n_draws, configured_draws)
  stats <- object$stats
  group_sizes <- if (model == "LRM") length(stats$y) else lengths(stats$y)
  check_gaussian_evaluation_budget(group_sizes, stats$p, n_draws,
                                   max_output_mb)
  scope <- inference_rng_scope(seed); on.exit(scope$restore(), add = TRUE)
  sampled <- if (model == "LRM") vi_sample_draws_lrm(object, n_draws) else
    vi_sample_draws_hlrm(object, results, n_draws)
  draws <- if (model == "LRM") list(
    Betaj = array(sampled$Beta, c(n_draws, 1L, stats$p)),
    sigma2j = matrix(sampled$sigma2, n_draws, 1L)) else sampled
  list(draws = draws, n_draws = n_draws)
}

pp_check_gaussian <- function(y_list, X_list, draws, seed = NULL,
                                return_yrep = FALSE, max_output_mb = 256,
                                target = "observed_groups_conditional_parameters",
                                discrepancies = ppc_discrepancies()) {
  info <- validate_grouped_data(X_list, y_list)
  S <- validate_gaussian_draws(draws, info$m, info$p)
  simulate_draw <- function(s) lapply(seq_len(info$m), function(j)
    gaussian_response_draw(X_list[[j]], draws$Betaj[s, j, ], draws$sigma2j[s, j]))
  pp_check_stream(y_list, S, simulate_draw, target, seed, return_yrep,
                  max_output_mb, X_list, discrepancies)
}

evaluate_in_sample_fit_gaussian <- function(y_list, X_list, point_beta,
                                            model, method,
                                            return_details = FALSE) {
  model <- match.arg(model, c("LRM", "HLRM"))
  method <- match.arg(method, c("MCMC", "CAVI", "SVI"))
  info <- validate_grouped_data(X_list, y_list)
  if (model == "LRM" && info$m != 1L) stop("LRM requires one regression block")
  if (!is.matrix(point_beta) || !is.numeric(point_beta) ||
      !identical(dim(point_beta), as.integer(c(info$m, info$p))) || any(!is.finite(point_beta)))
    stop("point_beta must be a finite [m,p] matrix")
  fitted <- lapply(seq_len(info$m), function(j)
    drop(X_list[[j]] %*% point_beta[j, ]))
  residuals <- Map(`-`, y_list, fitted)
  new_in_sample_fit_result(
    model, method, mse(unlist(y_list), unlist(fitted)),
    r2(unlist(y_list), unlist(fitted)),
    metadata = list(groups = info$m, observations = sum(info$nj),
                    estimate = "posterior mean of the conditional response mean",
                    interpretation = "descriptive in-sample fit; not out-of-sample performance"),
    fitted_values = if (return_details) fitted else NULL,
    residuals = if (return_details) residuals else NULL)
}

predictive_criteria_gaussian <- function(y_list, X_list, draws, model, method,
                                         predictive_target = NULL,
                                         predictive_unit = NULL,
                                         return_log_likelihood = FALSE,
                                         max_output_mb = 256,
                                         groups = NULL) {
  model <- match.arg(model, c("LRM", "HLRM"))
  method <- match.arg(method, c("MCMC", "CAVI", "SVI"))
  info <- validate_grouped_data(X_list, y_list)
  if (model == "LRM" && info$m != 1L) stop("LRM requires one regression block")
  S <- validate_gaussian_draws(draws, info$m, info$p)
  check_gaussian_evaluation_budget(info$nj, info$p, S, max_output_mb)
  if (!is.logical(return_log_likelihood) ||
      length(return_log_likelihood) != 1L || is.na(return_log_likelihood))
    stop("return_log_likelihood must be TRUE or FALSE")
  if (is.null(predictive_target)) predictive_target <-
    if (model == "LRM") "observed_data" else "observed_groups"
  predictive_target <- match.arg(predictive_target,
    c("observed_data", "observed_groups", "new_groups"))
  if (model == "HLRM" && predictive_target != "observed_groups")
    stop("Use predictive_criteria_hlrm_marginal() for new HLRM groups")
  if (model == "LRM" && predictive_target != "observed_data")
    stop("LRM predictive_target must be observed_data")
  if (model == "LRM") {
    if (is.null(predictive_unit)) predictive_unit <- if (is.null(groups))
      "observation" else "group"
    predictive_unit <- match.arg(predictive_unit, c("observation", "group"))
    if (predictive_unit == "group" && is.null(groups))
      stop("groups are required when predictive_unit = 'group'")
    if (predictive_unit == "observation" && !is.null(groups))
      stop("groups must be NULL when predictive_unit = 'observation'")
  }
  if (model == "LRM" && !is.null(groups) &&
      (!is.atomic(groups) || length(groups) != sum(info$nj) || anyNA(groups) ||
       (is.numeric(groups) && any(!is.finite(groups)))))
    stop("groups must give one nonmissing label per LRM observation")
  likelihood_draws <- if (model == "LRM")
    list(Beta = matrix(draws$Betaj[, 1L, ], S, info$p),
         sigma2 = draws$sigma2j[, 1L]) else draws
  likelihood_y <- if (model == "LRM") y_list[[1L]] else y_list
  likelihood_X <- if (model == "LRM") X_list[[1L]] else X_list
  observation_ll <- log_likelihood(
    model, likelihood_y, likelihood_X, likelihood_draws,
    unit = "observation", conditioning = "parameter_conditional",
    max_output_mb = max_output_mb)
  if (model == "HLRM") {
    group_ll <- log_likelihood(
      model, likelihood_y, likelihood_X, likelihood_draws,
      unit = "group", conditioning = "parameter_conditional",
      max_output_mb = max_output_mb)
    criteria <- list(
      waic_observation_conditional = waic_from_loglik(observation_ll),
      waic_group_conditional = waic_from_loglik(group_ll))
    estimands <- list(observation = log_likelihood_estimand(observation_ll),
                      group = log_likelihood_estimand(group_ll))
    retained_log_likelihood <- list(observation = observation_ll,
                                    group = group_ll)
  } else if (predictive_unit == "observation") {
    criteria <- list(waic_observation = waic_from_loglik(observation_ll))
    estimands <- list(observation = log_likelihood_estimand(observation_ll))
    retained_log_likelihood <- list(observation = observation_ll)
  } else {
    supplied <- log_likelihood(
      "LRM", likelihood_y, likelihood_X, likelihood_draws,
      unit = "group", conditioning = "parameter_conditional",
      groups = groups, max_output_mb = max_output_mb)
    criteria <- list(waic_group = waic_from_loglik(supplied))
    estimands <- list(group = log_likelihood_estimand(supplied))
    retained_log_likelihood <- list(group = supplied)
  }
  metadata <- list(draw_source = if (method == "MCMC")
      "posterior_mcmc" else "variational_q", draws = S,
      groups = info$m, observations = sum(info$nj),
      log_likelihood_estimands = estimands,
      interpretation = "predictive criterion under the declared estimand; not a convergence diagnostic")
  if (model == "LRM") metadata$predictive_unit <- predictive_unit
  new_predictive_criteria_result(
    model, method, predictive_target, criteria, metadata = metadata,
    log_likelihood = if (return_log_likelihood)
      retained_log_likelihood else NULL)
}

posterior_predictive_check_gaussian <- function(y_list, X_list, draws,
                                                model, method, seed = NULL,
                                                return_yrep = FALSE,
                                                max_output_mb = 256,
                                                discrepancies = ppc_discrepancies()) {
  model <- match.arg(model, c("LRM", "HLRM"))
  method <- match.arg(method, c("MCMC", "CAVI", "SVI"))
  target <- if (model == "HLRM")
    "observed_groups_conditional_parameters" else "observed_design_regression"
  check <- pp_check_gaussian(y_list, X_list, draws, seed, return_yrep,
                             max_output_mb, target, discrepancies)
  new_ppc_result(model, method,
    if (model == "HLRM") "observed_groups" else "observed_data", check,
    metadata = list(interpretation =
      "posterior model adequacy for observed data; not held-out predictive performance"))
}
