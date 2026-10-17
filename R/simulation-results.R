new_simulation_result <- function(model, train, truth, configuration,
                                  conditioning, parameter_source,
                                  design_source, test = NULL) {
  model <- match.arg(model, c("LRM", "HLRM", "CHLRM"))
  validate_simulation_split(train, model, "train")
  if (!is.null(test)) {
    validate_simulation_split(test, model, "test")
    train_p <- if (model == "LRM") ncol(train$X) else ncol(train$X[[1L]])
    test_p <- if (model == "LRM") ncol(test$X) else ncol(test$X[[1L]])
    if (train_p != test_p)
      stop("train and test must have the same design column count")
  }
  if (!is.list(truth) || !length(truth) || is.null(names(truth)) ||
      any(!nzchar(names(truth))) || anyDuplicated(names(truth)))
    stop("truth must be a nonempty uniquely named list")
  blocks <- list(truth = truth, configuration = configuration,
                 conditioning = conditioning,
                 parameter_source = parameter_source,
                 design_source = design_source)
  if (any(!vapply(blocks, is.list, logical(1L))) ||
      any(vapply(blocks, function(value)
        is.null(names(value)) || any(!nzchar(names(value))), logical(1L))))
    stop("simulation result metadata blocks must be lists")
  out <- list(schema_version = "simulation-result-2", model = model,
              train = train)
  if (!is.null(test)) out$test <- test
  out <- c(out, blocks)
  structure(out, class = c("vb_simulation_result", "list"))
}

.truth_source <- function(simulation, quantity) {
  source <- simulation$parameter_source
  if (grepl("conditional_mean|residuals", quantity)) return("generated")
  if (grepl("allocations|mixture_proportions", quantity))
    return(if (is.null(source$allocations)) "generated" else source$allocations)
  if (grepl("^(test_)?group_", quantity))
    return(if (is.null(source$local)) "generated" else source$local)
  if (grepl("^component_", quantity))
    return(if (is.null(source$components)) "generated" else source$components)
  source$global
}

.flatten_truth_value <- function(quantity, value, source) {
  if (is.list(value) && !is.data.frame(value)) {
    rows <- lapply(seq_along(value), function(group) {
      nested <- .flatten_truth_value(quantity, value[[group]], source)
      nested$index_3 <- nested$index_2
      nested$index_2 <- nested$index_1
      nested$index_1 <- group
      nested
    })
    return(do.call(rbind, rows))
  }
  if (!is.numeric(value) && !is.integer(value))
    stop("truth values must be numeric, integer or lists of numeric values")
  dimensions <- dim(value)
  if (is.null(dimensions)) {
    indices <- cbind(seq_along(value), NA_integer_, NA_integer_)
  } else {
    if (length(dimensions) > 3L)
      stop("tidy_truth() supports truth arrays with at most three dimensions")
    indices <- arrayInd(seq_along(value), dimensions)
    indices <- cbind(indices, matrix(
      NA_integer_, nrow(indices), 3L - ncol(indices)))
  }
  data.frame(
    quantity = quantity,
    index_1 = indices[, 1L], index_2 = indices[, 2L],
    index_3 = indices[, 3L], value = as.numeric(value), source = source,
    stringsAsFactors = FALSE)
}

tidy_truth <- function(simulation, include_observation_level = FALSE) {
  if (!inherits(simulation, "vb_simulation_result") ||
      !identical(simulation$schema_version, "simulation-result-2"))
    stop("simulation must be a simulation-result-2 object")
  include_observation_level <- .logical_scalar(
    include_observation_level, "include_observation_level")
  quantities <- names(simulation$truth)
  if (!include_observation_level) quantities <- quantities[
    !grepl("conditional_mean|residuals", quantities)]
  rows <- lapply(quantities, function(quantity) .flatten_truth_value(
    quantity, simulation$truth[[quantity]],
    .truth_source(simulation, quantity)))
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  structure(out, class = c("vb_tidy_truth", "data.frame"))
}

