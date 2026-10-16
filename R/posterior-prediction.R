# Unified posterior prediction for LRM, HLRM and CHLRM.
# Draw matrices are stored as observations x posterior draws.  For grouped
# models, one list element is one predictive group, so group-level latent
# quantities are shared by every row in that element.

gaussian_mean_draw <- function(X, beta) drop(as.matrix(X) %*% beta)

gaussian_response_draw <- function(X, beta, sigma2) {
  if (!is.numeric(sigma2) || length(sigma2) != 1L || !is.finite(sigma2) || sigma2 <= 0)
    stop("sigma2 must be a positive finite scalar")
  stats::rnorm(nrow(as.matrix(X)), gaussian_mean_draw(X, beta), sqrt(sigma2))
}

prediction_draw_count <- function(S, configured = NULL, available = NULL) {
  if (is.null(S)) S <- configured
  if (is.null(S)) S <- available
  if (is.null(S)) stop("Specify S or run_settings$samples_gof")
  if (!is.numeric(S) || length(S) != 1L || !is.finite(S) ||
      S < 2L || S != floor(S)) stop("S must be an integer >= 2")
  S <- as.integer(S)
  if (!is.null(available) && S > available)
    stop("S exceeds the number of stored MCMC draws")
  S
}

validate_prediction_level <- function(level) {
  if (!is.numeric(level) || length(level) != 1L || !is.finite(level) ||
      level <= 0 || level >= 1) stop("level must be a scalar strictly between 0 and 1")
  as.numeric(level)
}

validate_prediction_X <- function(X_new, p, grouped = TRUE) {
  if (!is.numeric(p) || length(p) != 1L || p < 1L) stop("p must be positive")
  if (!grouped) X_new <- list(as.matrix(X_new))
  if (!is.list(X_new) || !length(X_new))
    stop("X_new must be a non-empty list of design matrices")
  out <- lapply(seq_along(X_new), function(j) {
    Xj <- as.matrix(X_new[[j]])
    if (!is.numeric(Xj) || nrow(Xj) < 1L || ncol(Xj) != p || any(!is.finite(Xj)))
      stop(sprintf("X_new[[%d]] must be a finite matrix with %d columns", j, p))
    storage.mode(Xj) <- "double"
    Xj
  })
  out
}

resolve_prediction_groups <- function(X_new, m_train, group_ids = NULL) {
  g <- length(X_new)
  if (is.null(group_ids)) {
    if (g != m_train) {
      stop("group_ids is required unless X_new has one element per fitted group; use NA for each new group")
    }
    group_ids <- seq_len(m_train)
  }
  if (is.logical(group_ids) && length(group_ids) == g && all(is.na(group_ids)))
    group_ids <- rep(NA_integer_, g)
  if (!is.numeric(group_ids) || length(group_ids) != g ||
      any(!is.na(group_ids) & (!is.finite(group_ids) | group_ids != floor(group_ids) |
                                group_ids < 1L | group_ids > m_train))) {
    stop("group_ids must contain one fitted group index or NA per element of X_new")
  }
  as.integer(group_ids)
}

prediction_draw_indices <- function(available, S) {
  if (S == available) seq_len(S) else sort(sample.int(available, S))
}

check_prediction_budget <- function(X_new, S, max_output_mb) {
  # One returned draw matrix plus one temporary matrix used while simulating or
  # summarising. This is a local bound, not a total-process RAM guarantee.
  validate_evaluation_size(S, max_output_mb,
                           2 * as.double(sum(vapply(X_new, nrow, integer(1)))) * S)
}

prediction_result <- function(draws, model, method, target, level, group_ids,
                              return_draws, draw_source, allocation = NULL) {
  lower_prob <- (1 - level) / 2
  upper_prob <- 1 - lower_prob
  rows <- lapply(seq_along(draws), function(j) {
    value <- draws[[j]]
    qs <- t(apply(value, 1L, stats::quantile,
                  probs = c(lower_prob, 0.5, upper_prob), names = FALSE))
    data.frame(
      group = j,
      group_id = group_ids[j],
      group_status = if (model == "LRM") "not_applicable" else
        if (is.na(group_ids[j])) "new" else "observed",
      row = seq_len(nrow(value)),
      mean = rowMeans(value),
      sd = apply(value, 1L, stats::sd),
      median = qs[, 2L],
      lower = qs[, 1L],
      upper = qs[, 3L],
      stringsAsFactors = FALSE
    )
  })
  summary <- do.call(rbind, rows)
  rownames(summary) <- NULL
  new_prediction_result(
    model, method, target, level,
    metadata = list(
      draws = ncol(draws[[1L]]),
      draw_source = draw_source,
      grouping = if (model == "LRM") "not_applicable" else
        "one shared latent group state per X_new list element and posterior draw"
    ),
    summary = summary,
    allocation = allocation,
    draws = if (isTRUE(return_draws)) draws else NULL
  )
}

posterior_prediction_lrm <- function(X_new, draws, method, include_noise = FALSE,
                                     level = 0.95, return_draws = FALSE,
                                     max_output_mb = 256) {
  level <- validate_prediction_level(level)
  if (!is.matrix(draws$Beta) || !is.numeric(draws$Beta) ||
      !is.numeric(draws$sigma2) || nrow(draws$Beta) != length(draws$sigma2) ||
      any(!is.finite(draws$Beta)) || any(!is.finite(draws$sigma2)) || any(draws$sigma2 <= 0))
    stop("Invalid LRM posterior draws")
  X <- validate_prediction_X(X_new, ncol(draws$Beta), grouped = FALSE)
  S <- nrow(draws$Beta)
  check_prediction_budget(X, S, max_output_mb)
  values <- X[[1L]] %*% t(draws$Beta)
  if (include_noise) for (s in seq_len(S))
    values[, s] <- gaussian_response_draw(X[[1L]], draws$Beta[s, ], draws$sigma2[s])
  prediction_result(list(values), "LRM", method,
                    if (include_noise) "response" else "mean", level,
                    NA_integer_, return_draws,
                    if (method == "MCMC") "posterior_mcmc" else "variational_q")
}

posterior_prediction_hlrm <- function(X_new, local_draws, global_draws, group_ids,
                                      nu0, method, include_noise = FALSE,
                                      level = 0.95, return_draws = FALSE,
                                      max_output_mb = 256) {
  level <- validate_prediction_level(level)
  dims <- dim(local_draws$Betaj)
  if (length(dims) != 3L || !is.numeric(local_draws$sigma2j) ||
      !identical(dim(local_draws$sigma2j), as.integer(dims[1:2])) ||
      any(!is.finite(local_draws$Betaj)) || any(!is.finite(local_draws$sigma2j)) ||
      any(local_draws$sigma2j <= 0)) stop("Invalid HLRM local posterior draws")
  S <- dims[1L]; m_train <- dims[2L]; p <- dims[3L]
  X <- validate_prediction_X(X_new, p)
  group_ids <- resolve_prediction_groups(X, m_train, group_ids)
  check_prediction_budget(X, S, max_output_mb)
  needs_global <- anyNA(group_ids)
  if (needs_global) {
    if (is.null(global_draws) || !is.matrix(global_draws$Beta) ||
        !identical(dim(global_draws$Beta), as.integer(c(S, p))) ||
        !identical(dim(global_draws$Sigma), as.integer(c(S, p, p))) ||
        length(global_draws$xi2) != S || any(!is.finite(global_draws$Beta)) ||
        any(!is.finite(global_draws$Sigma)) || any(!is.finite(global_draws$xi2)) ||
        any(global_draws$xi2 <= 0)) stop("Invalid HLRM global posterior draws")
    if (!is.numeric(nu0) || length(nu0) != 1L || !is.finite(nu0) || nu0 <= 0)
      stop("nu0 must be positive for new-group prediction")
  }
  values <- vector("list", length(X))
  for (j in seq_along(X)) {
    id <- group_ids[j]
    if (is.na(id)) {
      beta <- matrix(vapply(seq_len(S), function(s)
        drop(mvtnorm::rmvnorm(1L, global_draws$Beta[s, ],
                              matrix(global_draws$Sigma[s, , ], p, p), method = "chol")),
        numeric(p)), nrow = S, ncol = p, byrow = TRUE)
      sigma2 <- 1 / stats::rgamma(S, shape = nu0 / 2,
                                  rate = nu0 * global_draws$xi2 / 2)
    } else {
      beta <- matrix(local_draws$Betaj[, id, ], S, p)
      sigma2 <- local_draws$sigma2j[, id]
    }
    values[[j]] <- X[[j]] %*% t(beta)
    if (include_noise) for (s in seq_len(S))
      values[[j]][, s] <- gaussian_response_draw(X[[j]], beta[s, ], sigma2[s])
  }
  prediction_result(values, "HLRM", method,
                    if (include_noise) "response" else "mean", level,
                    group_ids, return_draws,
                    if (method == "MCMC") "posterior_mcmc" else "variational_q")
}

posterior_prediction_chlrm <- function(X_new, draws, group_ids, method,
                                       z = NULL, rho = NULL,
                                       include_noise = FALSE, level = 0.95,
                                       return_draws = FALSE, max_output_mb = 256) {
  level <- validate_prediction_level(level)
  shape <- validate_component_draws_chlrm(draws)
  S <- shape$S; K <- shape$K; p <- shape$p
  m_train <- if (!is.null(z)) nrow(z) else if (!is.null(rho)) nrow(rho) else 0L
  if (m_train < 1L) stop("Provide posterior allocations z or rho for fitted groups")
  X <- validate_prediction_X(X_new, p)
  group_ids <- resolve_prediction_groups(X, m_train, group_ids)
  check_prediction_budget(X, S, max_output_mb)
  if (!is.null(z)) {
    if (!is.matrix(z) || !identical(dim(z), as.integer(c(m_train, S))) ||
        any(!is.finite(z)) || any(z < 1L | z > K) || any(z != floor(z)))
      stop("z must be an integer [fitted groups, S] matrix")
  } else {
    if (!is.matrix(rho) || ncol(rho) != K || any(!is.finite(rho)) ||
        any(rho < 0) || any(abs(rowSums(rho) - 1) > 1e-8))
      stop("rho must be a row-stochastic [fitted groups, K] matrix")
  }
  allocation <- matrix(NA_integer_, length(X), S)
  values <- vector("list", length(X))
  for (j in seq_along(X)) {
    id <- group_ids[j]
    if (is.na(id)) {
      allocation[j, ] <- vapply(seq_len(S), function(s)
        sample.int(K, 1L, prob = draws$omega[s, ]), integer(1))
    } else if (!is.null(z)) {
      allocation[j, ] <- z[id, ]
    } else {
      allocation[j, ] <- vapply(seq_len(S), function(s)
        sample.int(K, 1L, prob = rho[id, ]), integer(1))
    }
    beta <- matrix(vapply(seq_len(S), function(s)
      draws$Betak[s, allocation[j, s], ], numeric(p)),
      nrow = S, ncol = p, byrow = TRUE)
    sigma2 <- draws$sigma2k[cbind(seq_len(S), allocation[j, ])]
    values[[j]] <- X[[j]] %*% t(beta)
    if (include_noise) for (s in seq_len(S))
      values[[j]][, s] <- gaussian_response_draw(X[[j]], beta[s, ], sigma2[s])
  }
  alloc_summary <- data.frame(
    group = rep(seq_along(X), each = K),
    group_id = rep(group_ids, each = K),
    group_status = rep(ifelse(is.na(group_ids), "new", "observed"), each = K),
    component = rep(seq_len(K), times = length(X)),
    probability = as.vector(vapply(seq_along(X), function(j)
      tabulate(allocation[j, ], nbins = K) / S, numeric(K))),
    stringsAsFactors = FALSE
  )
  prediction_result(values, "CHLRM", method,
                    if (include_noise) "response" else "mean", level,
                    group_ids, return_draws,
                    if (method == "MCMC") "posterior_mcmc" else "variational_q",
                    allocation = alloc_summary)
}

run_lrm_prediction_mcmc <- function(chains, X_new, S, configured, seed, level,
                                    return_draws, max_output_mb, include_noise) {
  if (!length(chains) || is.null(chains$Beta)) stop("Call fit() successfully before prediction")
  available <- nrow(chains$Beta)
  S <- prediction_draw_count(S, configured, available)
  scope <- inference_rng_scope(seed); on.exit(scope$restore(), add = TRUE)
  idx <- prediction_draw_indices(available, S)
  ans <- posterior_prediction_lrm(X_new,
    list(Beta = chains$Beta[idx, , drop = FALSE], sigma2 = chains$sigma2[idx]),
    "MCMC", include_noise, level, return_draws, max_output_mb)
  ans$metadata$draw_indices <- idx
  ans
}

run_lrm_prediction_vi <- function(object, results, X_new, S, configured, seed,
                                  level, return_draws, max_output_mb,
                                  include_noise, stochastic = FALSE) {
  if (!length(results)) stop("Call fit() successfully before prediction")
  S <- prediction_draw_count(S, configured)
  scope <- inference_rng_scope(seed); on.exit(scope$restore(), add = TRUE)
  sampled <- vi_sample_draws_lrm(object, S)
  posterior_prediction_lrm(X_new, sampled,
                           if (isTRUE(stochastic)) "SVI" else "CAVI",
                           include_noise, level,
                           return_draws, max_output_mb)
}

run_hlrm_prediction_mcmc <- function(chains, X_new, group_ids, nu0, S,
                                     configured, seed, level, return_draws,
                                     max_output_mb, include_noise) {
  if (!length(chains) || is.null(chains$Betaj)) stop("Call fit() successfully before prediction")
  available <- dim(chains$Betaj)[1L]
  S <- prediction_draw_count(S, configured, available)
  scope <- inference_rng_scope(seed); on.exit(scope$restore(), add = TRUE)
  idx <- prediction_draw_indices(available, S)
  local <- list(Betaj = chains$Betaj[idx, , , drop = FALSE],
                sigma2j = chains$sigma2j[idx, , drop = FALSE])
  global <- list(Beta = chains$Beta[idx, , drop = FALSE],
                 Sigma = chains$Sigma[idx, , , drop = FALSE], xi2 = chains$xi2[idx])
  ans <- posterior_prediction_hlrm(X_new, local, global, group_ids, nu0,
    "MCMC", include_noise, level, return_draws, max_output_mb)
  ans$metadata$draw_indices <- idx
  ans
}

run_hlrm_prediction_vi <- function(object, results, X_new, group_ids, nu0, S,
                                   configured, seed, level, return_draws,
                                   max_output_mb, include_noise, stochastic) {
  if (!length(results)) stop("Call fit() successfully before prediction")
  S <- prediction_draw_count(S, configured)
  scope <- inference_rng_scope(seed); on.exit(scope$restore(), add = TRUE)
  local <- vi_sample_draws_hlrm(object, results, S)
  X_checked <- validate_prediction_X(X_new, object$stats$p)
  ids <- resolve_prediction_groups(X_checked, object$stats$m, group_ids)
  global <- if (anyNA(ids)) vi_global_draws_hlrm(results, S) else NULL
  posterior_prediction_hlrm(X_checked, local, global, ids, nu0,
    if (isTRUE(stochastic)) "SVI" else "CAVI", include_noise, level,
    return_draws, max_output_mb)
}

run_chlrm_prediction_mcmc <- function(chains, X_new, group_ids, S, configured,
                                      seed, level, return_draws, max_output_mb,
                                      include_noise) {
  if (!length(chains) || is.null(chains$Betak)) stop("Call fit() successfully before prediction")
  available <- dim(chains$Betak)[1L]
  S <- prediction_draw_count(S, configured, available)
  scope <- inference_rng_scope(seed); on.exit(scope$restore(), add = TRUE)
  idx <- prediction_draw_indices(available, S)
  draws <- list(Betak = chains$Betak[idx, , , drop = FALSE],
                sigma2k = chains$sigma2k[idx, , drop = FALSE],
                omega = chains$omega[idx, , drop = FALSE])
  ans <- posterior_prediction_chlrm(X_new, draws, group_ids, "MCMC",
    z = chains$z[, idx, drop = FALSE], include_noise = include_noise,
    level = level, return_draws = return_draws, max_output_mb = max_output_mb)
  ans$metadata$draw_indices <- idx
  ans
}

run_chlrm_prediction_vi <- function(object, results, X_new, group_ids, S,
                                    configured, seed, level, return_draws,
                                    max_output_mb, include_noise, stochastic) {
  if (!length(results)) stop("Call fit() successfully before prediction")
  S <- prediction_draw_count(S, configured)
  scope <- inference_rng_scope(seed); on.exit(scope$restore(), add = TRUE)
  draws <- vi_sample_draws_chlrm(object, results, S, include_allocations = FALSE)
  posterior_prediction_chlrm(X_new, draws, group_ids,
    if (isTRUE(stochastic)) "SVI" else "CAVI", rho = results$rho,
    include_noise = include_noise, level = level, return_draws = return_draws,
    max_output_mb = max_output_mb)
}
