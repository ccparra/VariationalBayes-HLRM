# Simulation-based calibration helpers owned by the evaluation suite. The
# package supplies prior-predictive simulation and posterior draws; ranks and
# truth-based performance summaries are properties of a replicated study.

.paper_rng_scope <- function(seed) {
  if (is.null(seed)) return(function() invisible(NULL))
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
      seed < 0 || seed != floor(seed))
    stop("seed must be NULL or a nonnegative integer")
  existed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  previous <- if (existed) get(".Random.seed", envir = .GlobalEnv) else NULL
  set.seed(as.integer(seed))
  function() {
    if (existed) assign(".Random.seed", previous, envir = .GlobalEnv) else if (
      exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv)
    invisible(NULL)
  }
}

sbc_ranks <- function(truth, draws, seed = NULL) {
  truth <- as.numeric(truth)
  if (!length(truth) || any(!is.finite(truth)))
    stop("truth must contain finite scalar estimands")
  draws <- as.matrix(draws)
  if (!is.numeric(draws) || nrow(draws) < 2L || ncol(draws) != length(truth) ||
      any(!is.finite(draws)))
    stop("draws must be a finite matrix with draws in rows and estimands in columns")
  quantity <- names(truth)
  if (is.null(quantity) || any(!nzchar(quantity)) || anyDuplicated(quantity)) {
    quantity <- colnames(draws)
  }
  if (is.null(quantity) || length(quantity) != length(truth) ||
      any(!nzchar(quantity)) || anyDuplicated(quantity))
    stop("estimands require unique nonempty names")
  if (!is.null(colnames(draws)) && !identical(colnames(draws), quantity))
    stop("draw column names must match truth names and order")

  restore <- .paper_rng_scope(seed)
  on.exit(restore(), add = TRUE)
  less <- vapply(seq_along(truth), function(index)
    sum(draws[, index] < truth[index]), integer(1L))
  ties <- vapply(seq_along(truth), function(index)
    sum(draws[, index] == truth[index]), integer(1L))
  randomized_ties <- vapply(ties, function(count)
    if (count) sample.int(count + 1L, 1L) - 1L else 0L, integer(1L))
  rank <- less + randomized_ties
  data.frame(
    quantity = quantity,
    truth = truth,
    posterior_draws = nrow(draws),
    less_than_truth = less,
    ties = ties,
    rank = rank,
    normalized_rank = (rank + .5) / (nrow(draws) + 1),
    stringsAsFactors = FALSE)
}

.sbc_vector_block <- function(prefix, truth, draws) {
  truth <- as.numeric(truth)
  draws <- as.matrix(draws)
  if (ncol(draws) != length(truth)) stop("Incompatible SBC block dimensions")
  names(truth) <- paste0(prefix, "[", seq_along(truth), "]")
  colnames(draws) <- names(truth)
  list(truth = truth, draws = draws)
}

chlrm_invariant_sbc_inputs <- function(simulation, posterior) {
  if (!inherits(simulation, "vb_simulation_result") ||
      !identical(simulation$model, "CHLRM"))
    stop("simulation must be a CHLRM simulation result")
  if (!inherits(posterior, "vb_posterior_draws") ||
      !identical(posterior$model, "CHLRM"))
    stop("posterior must be CHLRM posterior_draws() output")
  truth <- simulation$truth
  parameters <- posterior$parameters
  coefficient_draws <- parameters$component_coefficients
  shape <- dim(coefficient_draws)
  if (length(shape) != 3L) stop("component coefficient draws must be draw x component x coefficient")
  S <- shape[1L]; K <- shape[2L]
  if (!identical(dim(truth$component_coefficients), as.integer(shape[2:3])))
    stop("truth and posterior component coefficient dimensions differ")

  blocks <- list()
  add_block <- function(name, truth_value, draw_value) {
    blocks[[length(blocks) + 1L]] <<- .sbc_vector_block(
      name, truth_value, draw_value)
  }
  if (!is.null(truth$population_coefficients)) add_block(
    "population_coefficients", truth$population_coefficients,
    parameters$population_coefficients)
  if (!is.null(truth$variance_scale)) add_block(
    "variance_scale", truth$variance_scale,
    matrix(parameters$variance_scale, ncol = 1L))

  add_block(
    "sorted_mixture_weight", sort(truth$mixture_weights),
    t(apply(parameters$mixture_weights, 1L, sort)))
  add_block(
    "sorted_component_variance", sort(truth$component_variances),
    t(apply(parameters$component_variances, 1L, sort)))

  if (K > 1L) {
    truth_distance <- sort(as.numeric(stats::dist(
      truth$component_coefficients)))
    draw_distance <- t(vapply(seq_len(S), function(draw)
      sort(as.numeric(stats::dist(coefficient_draws[draw, , , drop = FALSE][1L, , ]))),
      numeric(length(truth_distance))))
    add_block("component_distance", truth_distance, draw_distance)
  }

  truth_vector <- unlist(lapply(blocks, `[[`, "truth"), use.names = TRUE)
  draw_matrix <- do.call(cbind, lapply(blocks, `[[`, "draws"))
  colnames(draw_matrix) <- names(truth_vector)
  list(
    truth = truth_vector,
    draws = draw_matrix,
    metadata = list(
      estimand = "label_invariant",
      excluded = c("raw component labels", "group allocations"),
      posterior_draws = S))
}
