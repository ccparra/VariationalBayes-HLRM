# PSIS-LOO delegates importance-weight smoothing and Pareto diagnostics to the
# reference loo package. It is intentionally restricted to MCMC posterior
# draws; applying the posterior identity directly to draws from q would not be
# exact Bayesian leave-one-out.

psis_loo <- function(
    fit, predictive_target = NULL, predictive_unit = NULL, groups = NULL,
    n_draws = NULL, seed = NULL, require_converged = TRUE,
    cores = 1L, max_output_mb = 256, marginal_control = list()) {
  if (!inherits(fit, "vb_model_fit"))
    stop("fit must be produced by fit_model()")
  if (fit$method != "MCMC")
    stop("PSIS-LOO requires MCMC posterior draws; use refitted CV for CAVI/SVI")
  if (!is.logical(require_converged) || length(require_converged) != 1L ||
      is.na(require_converged))
    stop("require_converged must be TRUE or FALSE")
  if (require_converged && fit$diagnostics$status != "converged")
    stop("PSIS-LOO requires a converged multichain fit; set require_converged = FALSE only for diagnostics")
  cores <- validate_integer_scalar(cores, "cores", 1L)
  if (!requireNamespace("loo", quietly = TRUE))
    stop("Install the suggested R package 'loo' to compute PSIS-LOO")

  if (fit$model == "LRM") {
    if (is.null(predictive_target)) predictive_target <- "observed_data"
    if (is.null(predictive_unit)) predictive_unit <- "observation"
    criteria <- predictive_criteria(
      fit, n_draws = n_draws, seed = seed,
      predictive_target = predictive_target,
      predictive_unit = predictive_unit,
      return_log_likelihood = TRUE, max_output_mb = max_output_mb,
      groups = groups)
    likelihood_name <- predictive_unit
  } else if (fit$model == "HLRM") {
    if (is.null(predictive_target)) predictive_target <- "new_groups"
    predictive_target <- match.arg(
      predictive_target, c("observed_groups", "new_groups"))
    if (is.null(predictive_unit)) predictive_unit <- "group"
    predictive_unit <- match.arg(predictive_unit, c("observation", "group"))
    if (predictive_target == "new_groups" && predictive_unit != "group")
      stop("new_groups requires predictive_unit = 'group'")
    criteria <- predictive_criteria(
      fit, n_draws = n_draws, seed = seed,
      predictive_target = predictive_target,
      return_log_likelihood = TRUE, max_output_mb = max_output_mb,
      marginal_control = marginal_control)
    likelihood_name <- predictive_unit
  } else {
    if (is.null(predictive_target)) predictive_target <- "new_groups"
    if (!identical(predictive_target, "new_groups"))
      stop("CHLRM PSIS-LOO currently targets new_groups with group-marginal likelihood")
    if (!is.null(predictive_unit) && !identical(predictive_unit, "group"))
      stop("CHLRM PSIS-LOO requires predictive_unit = 'group'")
    predictive_unit <- "group"
    criteria <- predictive_criteria(
      fit, n_draws = n_draws, seed = seed,
      predictive_target = predictive_target,
      return_log_likelihood = TRUE, max_output_mb = max_output_mb)
    likelihood_name <- "group"
  }
  likelihood <- criteria$log_likelihood[[likelihood_name]]
  if (is.null(likelihood))
    stop("The requested predictive unit is unavailable for this target")
  log_likelihood_draw_unit <- as_log_lik_matrix(likelihood)
  if (nrow(log_likelihood_draw_unit) < 4L)
    stop("PSIS-LOO requires at least four posterior draws")

  origin <- attr(.fit_backend(fit)$chains, "draw_origin", exact = TRUE)
  indices <- criteria$metadata$draw_indices
  chain_id <- if (is.null(origin)) rep.int(1L, nrow(log_likelihood_draw_unit)) else
    origin$chain[indices]
  chain_rows <- split(seq_along(chain_id), chain_id)
  balanced_count <- min(lengths(chain_rows))
  if (balanced_count < 4L)
    stop("PSIS relative-efficiency estimation requires at least four selected draws from every chain")
  efficiency_rows <- unlist(lapply(
    chain_rows, function(rows) rows[seq_len(balanced_count)]),
    use.names = FALSE)
  efficiency_likelihood <- log_likelihood_draw_unit[efficiency_rows, , drop = FALSE]
  stabilized_likelihood <- sweep(
    efficiency_likelihood, 2L,
    apply(efficiency_likelihood, 2L, max), `-`)
  relative_efficiency <- loo::relative_eff(
    exp(stabilized_likelihood), chain_id = chain_id[efficiency_rows],
    cores = cores)
  if (any(!is.finite(relative_efficiency)))
    stop("PSIS relative-efficiency estimation was not finite; use more post-warmup draws")
  loo_object <- loo::loo(
    log_likelihood_draw_unit, r_eff = relative_efficiency,
    save_psis = TRUE, cores = cores)
  pointwise <- cbind(
    likelihood$predictive_units,
    as.data.frame(loo_object$pointwise, check.names = FALSE))
  estimates <- as.data.frame(loo_object$estimates, check.names = FALSE)
  estimates$quantity <- rownames(estimates)
  rownames(estimates) <- NULL
  estimates <- estimates[, c("quantity", setdiff(names(estimates), "quantity"))]
  pareto_k <- loo::pareto_k_values(loo_object)
  structure(list(
    schema_version = "psis-loo-result-1",
    model = fit$model,
    method = fit$method,
    predictive_target = predictive_target,
    unit = likelihood$unit,
    conditioning = likelihood$conditioning,
    estimates = estimates,
    pointwise = pointwise,
    diagnostics = list(
      pareto_k = pareto_k,
      pareto_k_table = loo::pareto_k_table(loo_object),
      mcse_elpd_loo = loo::mcse_loo(loo_object),
      reliable = all(is.finite(pareto_k) & pareto_k <= 0.7)),
    metadata = list(
      draws = nrow(log_likelihood_draw_unit),
      chains = length(unique(chain_id)),
      draws_per_chain_for_relative_efficiency = balanced_count,
      convergence_status = fit$diagnostics$status,
      density_scale = "input",
      estimand = log_likelihood_estimand(likelihood),
      implementation = paste0("loo ", as.character(utils::packageVersion("loo")))),
    loo_object = loo_object),
    class = c("vb_psis_loo", "list"))
}
