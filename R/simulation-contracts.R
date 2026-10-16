# Public contract for simulated datasets. Model-specific truth is allowed, but
# the container and the data split vocabulary are shared by all models.

validate_simulation_split <- function(data, model, name) {
  if (!is.list(data) || !all(c("X", "y") %in% names(data)))
    stop(name, " must contain X and y")
  if (model == "LRM") {
    X <- as.matrix(data$X)
    y <- as.numeric(data$y)
    if (!is.numeric(X) || !nrow(X) || nrow(X) != length(y) ||
        any(!is.finite(X)) || any(!is.finite(y)))
      stop(name, " must contain a finite design matrix and matching response")
  } else {
    validate_grouped_data(data$X, data$y)
    if (!is.null(data$group)) {
      expected <- lapply(seq_along(data$y), function(j)
        rep.int(j, length(data$y[[j]])))
      if (!identical(data$group, expected))
        stop(name, "$group must encode consecutive group identifiers")
    }
  }
  invisible(data)
}
