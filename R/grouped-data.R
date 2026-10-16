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
