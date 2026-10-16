# Generic multichain orchestration and modern convergence diagnostics.
#
# This module deliberately does not run chains in parallel.  Sequential
# execution has predictable memory use and avoids oversubscribing a laptop.
# The paper-level scheduler may parallelize independent datasets later.

.validate_multichain_data <- function(model, y, X, K) {
  if (model == "LRM") {
    X <- as.matrix(X); y <- as.numeric(y)
    if (!is.numeric(X) || !nrow(X) || !ncol(X) || length(y) != nrow(X) ||
        any(!is.finite(X)) || any(!is.finite(y)))
      stop("LRM requires finite y and X with compatible nonzero dimensions")
    return(list(m = 1L, p = ncol(X), n = length(y), X = X, y = y))
  }
  dims <- validate_grouped_data(X, y, if (model == "CHLRM") K else NULL)
  if (model == "CHLRM" && is.null(K)) stop("K is required for CHLRM")
  c(dims, list(X = X, y = y))
}

.estimate_multichain_mb <- function(model, dims, samples, n_chains, K = NULL,
                                    diagnostic_scope = "global") {
  per_draw <- .mcmc_values_per_draw(model, dims, K)
  diagnostic <- .mcmc_diagnostic_values(
    model, dims, K, diagnostic_scope)
  # Stored chains + diagnostic array + conservative list/object overhead.
  1.35 * samples * n_chains * 8 * (per_draw + diagnostic) / 1024^2
}

.new_mcmc_chain <- function(model, y, X, K, hyperparams, run_settings) {
  switch(model,
    LRM = mcmc_lrm$new(y, X, hyperparams, run_settings),
    HLRM = mcmc_hlrm$new(y, X, hyperparams, run_settings),
    CHLRM = mcmc_chlrm$new(y, X, K, hyperparams, run_settings))
}

.run_one_mcmc_chain <- function(model, y, X, K, hyperparams, run_settings,
                                seed, verbose, initialize = NULL,
                                chain_index = 1L) {
  # LRM/HLRM constructors currently draw their initial state.  Scope both
  # construction and fitting so each complete chain is reproducible and does
  # not alter the caller's ambient RNG state.
  scope <- inference_rng_scope(seed)
  on.exit(scope$restore(), add = TRUE)
  fit <- .new_mcmc_chain(model, y, X, K, hyperparams, run_settings)
  if (!is.null(initialize)) {
    fit <- initialize(fit, chain_index)
    if (!inherits(fit, "R6"))
      stop("initialize must return the fitted R6 backend")
  }
  result <- fit$fit(verbose = verbose, seed = NULL)
  list(model = fit, result = result)
}

.named_matrix_columns <- function(x, prefix) {
  x <- as.matrix(x)
  colnames(x) <- paste0(prefix, "[", seq_len(ncol(x)), "]")
  x
}

.named_array_columns <- function(x, prefix) {
  d <- dim(x)
  out <- matrix(x, nrow = d[1L], ncol = prod(d[-1L]))
  indices <- expand.grid(lapply(d[-1L], seq_len), KEEP.OUT.ATTRS = FALSE)
  colnames(out) <- apply(indices, 1L, function(idx)
    paste0(prefix, "[", paste(idx, collapse = ","), "]"))
  out
}

.sigma_columns <- function(Sigma, prefix = "Sigma") {
  p <- dim(Sigma)[2L]
  indices <- which(lower.tri(matrix(0, p, p), diag = TRUE), arr.ind = TRUE)
  out <- vapply(seq_len(nrow(indices)), function(q)
    Sigma[, indices[q, 1L], indices[q, 2L]], numeric(dim(Sigma)[1L]))
  if (is.null(dim(out))) out <- matrix(out, ncol = 1L)
  colnames(out) <- apply(indices, 1L, function(idx)
    paste0(prefix, "[", idx[1L], ",", idx[2L], "]"))
  out
}

.diagnostic_matrix_mcmc <- function(draws, model, scope = c("global", "all")) {
  scope <- match.arg(scope)
  global <- switch(model,
    LRM = cbind(.named_matrix_columns(draws$Beta, "Beta"),
                `log_sigma2` = log(draws$sigma2)),
    HLRM = cbind(.named_matrix_columns(draws$Beta, "Beta"),
                 .sigma_columns(draws$Sigma), `log_xi2` = log(draws$xi2)),
    CHLRM = {
      weights <- if (ncol(draws$omega) > 1L) {
        w <- log(draws$omega[, -ncol(draws$omega), drop = FALSE] /
                 draws$omega[, ncol(draws$omega)])
        colnames(w) <- paste0("omega_logratio[", seq_len(ncol(w)), "]")
        w
      } else NULL
      cbind(.named_array_columns(draws$Betak, "Betak"),
            .named_matrix_columns(log(draws$sigma2k), "log_sigma2k"),
            weights, .named_matrix_columns(draws$Beta, "Beta"),
            .sigma_columns(draws$Sigma), `log_xi2` = log(draws$xi2))
    })
  if (scope == "global" || model == "LRM") return(global)
  local <- cbind(.named_array_columns(draws$Betaj, "Betaj"),
                 .named_matrix_columns(log(draws$sigma2j), "log_sigma2j"))
  cbind(global, local)
}

.diagnose_mcmc_chains <- function(chains, model, elapsed_seconds,
                                  scope, thresholds) {
  if (!requireNamespace("posterior", quietly = TRUE))
    stop("Package 'posterior' is required for multichain diagnostics")
  matrices <- lapply(chains, .diagnostic_matrix_mcmc, model = model, scope = scope)
  dimensions <- vapply(matrices, nrow, integer(1L))
  if (length(unique(dimensions)) != 1L ||
      !all(vapply(matrices, function(x) identical(colnames(x), colnames(matrices[[1L]])),
                  logical(1L))))
    stop("All chains must expose identical iterations and diagnostic variables")
  draws <- array(NA_real_, c(dimensions[1L], length(matrices), ncol(matrices[[1L]])),
                 dimnames = list(iteration = NULL, chain = NULL,
                                 variable = colnames(matrices[[1L]])))
  for (chain in seq_along(matrices)) draws[, chain, ] <- matrices[[chain]]
  summary <- as.data.frame(posterior::summarise_draws(
    posterior::as_draws_array(draws), "mean", "sd", "rhat", "ess_bulk",
    "ess_tail", "mcse_mean", "mcse_sd"))
  total_seconds <- sum(elapsed_seconds)
  summary$ess_per_second <- summary$ess_bulk / total_seconds
  summary$relative_mcse_mean <- summary$mcse_mean / summary$sd
  summary$relative_mcse_sd <- summary$mcse_sd / summary$sd
  finite <- is.finite(summary$rhat) & is.finite(summary$ess_bulk) &
    is.finite(summary$ess_tail) & is.finite(summary$relative_mcse_mean) &
    is.finite(summary$relative_mcse_sd)
  summary$passed <- finite & summary$rhat <= thresholds$rhat &
    summary$ess_bulk >= thresholds$ess_bulk &
    summary$ess_tail >= thresholds$ess_tail &
    summary$relative_mcse_mean <= thresholds$relative_mcse_mean &
    summary$relative_mcse_sd <= thresholds$relative_mcse_sd
  list(
    variables = summary,
    passed = nrow(summary) > 0L && all(summary$passed),
    status = if (!all(finite)) "insufficient_draws" else
      if (all(summary$passed)) "passed" else "failed",
    thresholds = thresholds,
    scope = scope,
    implementation = paste0("posterior ",
      as.character(utils::packageVersion("posterior"))),
    total_chain_seconds = total_seconds
  )
}

.diagnose_mcmc_chain_set <- function(chains, model, elapsed_seconds,
                                     scope, thresholds) {
  reference <- NULL
  diagnostic_chains <- chains
  if (model == "CHLRM") {
    combined_raw <- .combine_chain_draws(chains, model)
    reference <- fit_label_reference_chlrm(combined_raw)
    diagnostic_chains <- lapply(
      chains, align_to_reference_chlrm, reference = reference)
  }
  list(
    reference = reference,
    chains = diagnostic_chains,
    diagnostics = .diagnose_mcmc_chains(
      diagnostic_chains, model, elapsed_seconds, scope, thresholds))
}

# Run independent MCMC chains, align CHLRM labels to one shared reference, and
# diagnose before any draws can be combined.
run_mcmc_chains <- function(model = c("LRM", "HLRM", "CHLRM"), y, X,
                            hyperparams, run_settings, K = NULL,
                            n_chains = 4L, seed = NULL, seeds = NULL,
                            diagnostic_scope = c("global", "all"),
                            rhat_max = 1.01, ess_bulk_min = 1000,
                            ess_tail_min = 1000,
                            relative_mcse_mean_max = 1,
                            relative_mcse_sd_max = 1,
                            max_memory_mb = 1024,
                            block_size = NULL, max_samples = NULL,
                            min_samples = 1L, initialize = NULL,
                            verbose = FALSE) {
  model <- match.arg(model)
  diagnostic_scope <- match.arg(diagnostic_scope)
  n_chains <- validate_integer_scalar(n_chains, "n_chains", 2L)
  validate_mcmc_settings(run_settings)
  thin <- if (is.null(run_settings$thin)) 1L else run_settings$thin
  if (!identical(as.integer(thin), 1L))
    stop("Multichain inference requires thin = 1; diagnose autocorrelation instead of discarding draws")
  if (!is.list(hyperparams)) stop("hyperparams must be a list")
  if (!is.null(initialize) && !is.function(initialize))
    stop("initialize must be NULL or a function")
  min_samples <- validate_integer_scalar(min_samples, "min_samples", 1L)
  dims <- .validate_multichain_data(model, y, X, K)
  thresholds <- list(rhat = as.numeric(rhat_max),
                     ess_bulk = as.numeric(ess_bulk_min),
                     ess_tail = as.numeric(ess_tail_min),
                     relative_mcse_mean = as.numeric(relative_mcse_mean_max),
                     relative_mcse_sd = as.numeric(relative_mcse_sd_max))
  if (any(!is.finite(unlist(thresholds))) || rhat_max < 1 ||
      ess_bulk_min <= 0 || ess_tail_min <= 0 ||
      relative_mcse_mean_max <= 0 || relative_mcse_sd_max <= 0)
    stop("Diagnostic thresholds require rhat_max >= 1 and positive ESS/relative MCSE limits")
  if (!is.numeric(max_memory_mb) || length(max_memory_mb) != 1L ||
      !is.finite(max_memory_mb) || max_memory_mb <= 0)
    stop("max_memory_mb must be a positive finite scalar")
  incremental <- !is.null(block_size) || !is.null(max_samples)
  if (incremental && (is.null(block_size) || is.null(max_samples)))
    stop("Provide both block_size and max_samples for incremental execution")
  if (incremental) {
    block_size <- validate_integer_scalar(block_size, "block_size", 1L)
    max_samples <- validate_integer_scalar(max_samples, "max_samples", 1L)
    if (block_size > max_samples) stop("block_size cannot exceed max_samples")
    if (min_samples > max_samples) stop("min_samples cannot exceed max_samples")
  } else {
    block_size <- run_settings$samples
    max_samples <- run_settings$samples
  }
  estimated_mb <- .estimate_multichain_mb(model, dims, max_samples,
    n_chains, K, diagnostic_scope)
  if (estimated_mb > max_memory_mb)
    stop(sprintf("Estimated multichain storage %.1f MB exceeds max_memory_mb %.1f",
                 estimated_mb, max_memory_mb))
  chain_seeds <- derive_independent_seeds(n_chains, seed, seeds)

  chains <- vector("list", n_chains)
  models <- vector("list", n_chains)
  elapsed <- numeric(n_chains)
  for (chain in seq_len(n_chains)) {
    settings <- run_settings
    settings$samples <- min(block_size, max_samples)
    settings$thin <- 1L
    settings$seed <- chain_seeds[chain]
    run <- .run_one_mcmc_chain(model, dims$y, dims$X, K, hyperparams,
                               settings, chain_seeds[chain], verbose,
                               initialize, chain)
    models[[chain]] <- run$model
    chains[[chain]] <- run$model$chains
    elapsed[chain] <- run$result$timing$total_seconds
  }

  finite_max <- function(x) if (any(is.finite(x)))
    max(x[is.finite(x)]) else NA_real_
  finite_min <- function(x) if (any(is.finite(x)))
    min(x[is.finite(x)]) else NA_real_
  diagnosed <- .diagnose_mcmc_chain_set(
    chains, model, elapsed, diagnostic_scope, thresholds)
  current_samples <- min(block_size, max_samples)
  passed <- diagnosed$diagnostics$passed && current_samples >= min_samples
  history <- data.frame(
    block = 1L, stored_samples = current_samples,
    status = if (current_samples < min_samples) "minimum_not_reached" else
      diagnosed$diagnostics$status,
    max_rhat = finite_max(diagnosed$diagnostics$variables$rhat),
    min_ess_bulk = finite_min(diagnosed$diagnostics$variables$ess_bulk),
    min_ess_tail = finite_min(diagnosed$diagnostics$variables$ess_tail),
    elapsed_seconds = sum(elapsed))

  while (incremental && !passed &&
         current_samples < max_samples) {
    additional <- min(block_size, max_samples - current_samples)
    for (chain in seq_len(n_chains)) {
      .resume_mcmc_backend(models[[chain]], additional, verbose = verbose)
      chains[[chain]] <- models[[chain]]$chains
      elapsed[chain] <- models[[chain]]$fit_metadata$total_seconds
    }
    current_samples <- current_samples + additional
    diagnosed <- .diagnose_mcmc_chain_set(
      chains, model, elapsed, diagnostic_scope, thresholds)
    passed <- diagnosed$diagnostics$passed && current_samples >= min_samples
    history <- rbind(history, data.frame(
      block = nrow(history) + 1L, stored_samples = current_samples,
      status = if (current_samples < min_samples) "minimum_not_reached" else
        diagnosed$diagnostics$status,
      max_rhat = finite_max(diagnosed$diagnostics$variables$rhat),
      min_ess_bulk = finite_min(diagnosed$diagnostics$variables$ess_bulk),
      min_ess_tail = finite_min(diagnosed$diagnostics$variables$ess_tail),
      elapsed_seconds = sum(elapsed)))
  }

  label_reference <- diagnosed$reference
  diagnostic_chains <- diagnosed$chains
  diagnostics <- diagnosed$diagnostics
  diagnostics$passed <- passed
  diagnostics$status <- if (passed) "passed" else
    if (current_samples < min_samples) "minimum_not_reached" else
      diagnostics$status
  out <- structure(list(
    schema_version = "mcmc-multichain-1",
    model = model,
    status = if (passed) "diagnosed" else "not_converged",
    chains = chains,
    aligned_chains = if (model == "CHLRM") diagnostic_chains else NULL,
    label_reference = label_reference,
    diagnostics = diagnostics,
    block_history = history,
    chain_metadata = data.frame(chain = seq_len(n_chains), seed = chain_seeds,
                                elapsed_seconds = elapsed),
    run_settings = modifyList(run_settings,
                              list(samples = current_samples, thin = 1L)),
    execution = list(incremental = incremental, block_size = block_size,
                     max_samples = max_samples,
                     stop_reason = if (passed) "diagnostic_thresholds" else
                       if (incremental) "max_samples" else "fixed_samples"),
    estimated_storage_mb = estimated_mb
  ), class = "mcmc_multichain_fit")
  attr(out, "models") <- models
  out
}

# Combining is a separate, guarded operation so downstream summaries cannot
# silently consume chains that failed the declared convergence rule.
combine_mcmc_chains <- function(fit, require_converged = TRUE) {
  if (!inherits(fit, "mcmc_multichain_fit"))
    stop("fit must be an mcmc_multichain_fit")
  if (!is.logical(require_converged) || length(require_converged) != 1L ||
      is.na(require_converged)) stop("require_converged must be TRUE or FALSE")
  if (require_converged && !isTRUE(fit$diagnostics$passed))
    stop("Chains did not satisfy the declared R-hat and ESS thresholds")
  source_chains <- if (fit$model == "CHLRM") fit$aligned_chains else fit$chains
  combined <- .combine_chain_draws(source_chains, fit$model)
  iterations <- vapply(source_chains, function(x)
    if (fit$model == "CHLRM") dim(x$Betak)[1L] else nrow(x$Beta), integer(1L))
  attr(combined, "draw_origin") <- data.frame(
    chain = rep(seq_along(iterations), times = iterations),
    iteration = sequence(iterations),
    seed = rep(fit$chain_metadata$seed, times = iterations))
  attr(combined, "diagnostics") <- fit$diagnostics
  attr(combined, "label_reference") <- fit$label_reference
  class(combined) <- c("mcmc_combined_draws", "list")
  combined
}

print.mcmc_multichain_fit <- function(x, ...) {
  cat(sprintf("%s multichain MCMC: %d chains x %d stored draws\n",
              x$model, nrow(x$chain_metadata), x$run_settings$samples))
  cat(sprintf("Diagnostics: %s; max R-hat %.4f; min bulk/tail ESS %.1f / %.1f\n",
              x$diagnostics$status, max(x$diagnostics$variables$rhat, na.rm = TRUE),
              min(x$diagnostics$variables$ess_bulk, na.rm = TRUE),
              min(x$diagnostics$variables$ess_tail, na.rm = TRUE)))
  invisible(x)
}
