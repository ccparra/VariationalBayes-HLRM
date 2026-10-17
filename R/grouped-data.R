# Build the matrix/list representation shared by HLRM and CHLRM from a data
# frame. Formula support is deliberately separate: this constructor performs
# no implicit contrast expansion and keeps the design fully auditable.

.resolve_frame_columns <- function(data, columns, name, multiple = TRUE) {
  if (is.character(columns)) {
    if (!length(columns) || anyNA(columns) || any(!nzchar(columns)) ||
        anyDuplicated(columns) || any(!columns %in% names(data)))
      stop(name, " must contain unique existing column names")
    indices <- match(columns, names(data))
  } else if (is.numeric(columns)) {
    if (!length(columns) || anyNA(columns) || any(!is.finite(columns)) ||
        any(columns != floor(columns)) || any(columns < 1L) ||
        any(columns > ncol(data)) || anyDuplicated(columns))
      stop(name, " must contain unique valid column indices")
    indices <- as.integer(columns)
  } else stop(name, " must contain column names or indices")
  if (!multiple && length(indices) != 1L)
    stop(name, " must identify exactly one column")
  indices
}

#' @rdname grouped-regression-data
#' @export
grouped_regression_data <- function(data, response, predictors, group,
                                    intercept = TRUE) {
  if (!is.data.frame(data) || !nrow(data))
    stop("data must be a nonempty data.frame")
  response_index <- .resolve_frame_columns(
    data, response, "response", multiple = FALSE)
  predictor_indices <- .resolve_frame_columns(
    data, predictors, "predictors", multiple = TRUE)
  group_index <- .resolve_frame_columns(
    data, group, "group", multiple = FALSE)
  if (response_index %in% predictor_indices || group_index %in% predictor_indices ||
      response_index == group_index)
    stop("response, predictors and group must identify distinct columns")
  if (!is.logical(intercept) || length(intercept) != 1L || is.na(intercept))
    stop("intercept must be TRUE or FALSE")

  response_values <- data[[response_index]]
  predictor_frame <- data[predictor_indices]
  if (!is.numeric(response_values) || any(!is.finite(response_values)))
    stop("response must be finite numeric")
  if (!all(vapply(predictor_frame, is.numeric, logical(1L))) ||
      any(!is.finite(as.matrix(predictor_frame))))
    stop("predictors must be finite numeric columns")
  if (anyNA(data[[group_index]])) stop("group cannot contain missing values")

  labels <- as.character(data[[group_index]])
  levels <- unique(labels)
  indices <- lapply(levels, function(label) which(labels == label))
  design <- as.matrix(predictor_frame)
  storage.mode(design) <- "double"
  if (intercept)
    design <- cbind(`(Intercept)` = 1, design)
  X <- lapply(indices, function(index)
    design[index, , drop = FALSE])
  y <- lapply(indices, function(index)
    unname(as.numeric(response_values[index])))
  names(X) <- names(y) <- NULL
  validate_grouped_data(X, y)
  structure(list(
    X = X,
    y = y,
    mapping = data.frame(
      group_id = seq_along(levels), group_label = levels,
      stringsAsFactors = FALSE),
    columns = list(
      response = names(data)[response_index],
      predictors = names(data)[predictor_indices],
      group = names(data)[group_index],
      intercept = intercept)),
    class = c("vb_grouped_regression_data", "list"))
}
validate_grouped_data <- function(X, y, K = NULL) {
  if (!is.list(X) || !is.list(y) || !length(X) || length(X) != length(y)) {
    stop("X and y must be non-empty lists with the same number of groups")
  }
  coefficient_count <- ncol(as.matrix(X[[1L]]))
  if (is.null(coefficient_count) || coefficient_count < 1L) {
    stop("Each X[[j]] must be a matrix with at least one column")
  }
  for (group_index in seq_along(X)) {
    group_design <- as.matrix(X[[group_index]])
    group_response <- as.numeric(y[[group_index]])
    if (ncol(group_design) != coefficient_count ||
        nrow(group_design) != length(group_response) ||
        !length(group_response)) {
      stop(sprintf("Incompatible dimensions in group %d", group_index))
    }
    if (any(!is.finite(group_design)) || any(!is.finite(group_response))) {
      stop(sprintf("Non-finite data in group %d", group_index))
    }
  }
  if (!is.null(K) &&
      (length(K) != 1L || !is.finite(K) || K < 1L || K != as.integer(K))) {
    stop("K must be a positive integer")
  }
  invisible(list(
    m = length(X),
    p = coefficient_count,
    nj = lengths(y)
  ))
}
