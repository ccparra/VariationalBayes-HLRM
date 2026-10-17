.expand_design_parameter <- function(value, q, name, positive = FALSE) {
  if (q == 0L) return(numeric())
  if (length(value) == 1L) value <- rep(value, q)
  if (!is.numeric(value) || length(value) != q || any(!is.finite(value)) ||
      (positive && any(value <= 0)))
    stop(name, " must have one ", if (positive) "positive " else "",
         "finite value per non-binary covariate")
  as.numeric(value)
}

.design_matrix <- function(design, n, coefficient_count, split, group = NULL) {
  if (design$distribution == "fixed") {
    source <- if (split == "train") design$fixed_train else design$fixed_test
    if (is.null(source)) stop("A fixed test design requires fixed_test")
    value <- if (is.list(source) && !is.data.frame(source)) source[[group]] else source
    value <- as.matrix(value)
    if (!is.numeric(value) || nrow(value) != n ||
        ncol(value) != coefficient_count || any(!is.finite(value)))
      stop("Fixed design dimensions must match the sample and coefficient count")
    if (!is.null(design$column_names)) {
      if (length(design$column_names) != ncol(value))
        stop("column_names length must equal the final design column count")
      colnames(value) <- design$column_names
    }
    return(value)
  }
  q <- coefficient_count - as.integer(design$include_intercept)
  if (q < 0L) stop("coefficient_count is incompatible with include_intercept")
  if (design$distribution == "multivariate_normal") {
    if (nrow(design$covariance) != q)
      stop("covariance dimension must match the model covariate count")
    location <- .expand_design_parameter(design$mean, q, "mean")
    covariates <- if (q) mvtnorm::rmvnorm(
      n, location, design$covariance, method = "chol") else matrix(numeric(), n, 0L)
  } else if (design$distribution == "independent_normal") {
    location <- .expand_design_parameter(design$mean, q, "mean")
    scale <- .expand_design_parameter(design$sd, q, "sd", TRUE)
    covariates <- if (q) matrix(stats::rnorm(n * q), n, q) else
      matrix(numeric(), n, 0L)
    if (q) covariates <- sweep(sweep(covariates, 2L, scale, `*`),
                               2L, location, `+`)
  } else if (design$distribution == "mixed") {
    binary_count <- length(design$binary_probabilities)
    numeric_count <- q - binary_count
    if (numeric_count < 0L)
      stop("binary_probabilities exceed the model covariate count")
    binary <- if (binary_count) vapply(
      design$binary_probabilities,
      function(probability) stats::rbinom(n, 1L, probability), numeric(n)) else
        matrix(numeric(), n, 0L)
    numeric_mean <- .expand_design_parameter(design$mean, numeric_count, "mean")
    if (!is.null(design$covariance)) {
      if (nrow(design$covariance) != numeric_count)
        stop("covariance dimension must match the continuous covariate count")
      continuous <- if (numeric_count) mvtnorm::rmvnorm(
        n, numeric_mean, design$covariance, method = "chol") else
          matrix(numeric(), n, 0L)
    } else {
      numeric_sd <- .expand_design_parameter(design$sd, numeric_count, "sd", TRUE)
      continuous <- if (numeric_count) matrix(
        stats::rnorm(n * numeric_count), n, numeric_count) else
          matrix(numeric(), n, 0L)
      if (numeric_count) continuous <- sweep(
        sweep(continuous, 2L, numeric_sd, `*`), 2L, numeric_mean, `+`)
    }
    covariates <- if (design$binary_position == "first")
      cbind(binary, continuous) else cbind(continuous, binary)
  } else {
    if (nrow(design$bounds) != q)
      stop("bounds must have one row per model covariate")
    if (!q) covariates <- matrix(numeric(), n, 0L) else {
      base <- floor(n^(1 / q)); sizes <- rep(max(1L, base), q)
      while (prod(sizes) < n) for (index in seq_len(q)) {
        if (prod(sizes) >= n) break
        sizes[index] <- sizes[index] + 1L
      }
      axes <- lapply(seq_len(q), function(index) seq(
        design$bounds[index, 1L], design$bounds[index, 2L],
        length.out = sizes[index]))
      covariates <- as.matrix(do.call(expand.grid, axes))[seq_len(n), , drop = FALSE]
      jitter <- .expand_design_parameter(
        design$jitter_sd, q, "jitter_sd")
      for (index in seq_len(q)) covariates[, index] <- pmin(pmax(
        covariates[, index] + stats::rnorm(n, sd = jitter[index]),
        design$bounds[index, 1L]), design$bounds[index, 2L])
    }
  }
  X <- if (design$include_intercept) cbind(intercept = 1, covariates) else covariates
  if (!is.null(design$column_names)) {
    if (length(design$column_names) != ncol(X))
      stop("column_names length must equal the final design column count")
    colnames(X) <- design$column_names
  } else {
    dimnames(X) <- NULL
  }
  X
}

.validate_probability_vector <- function(value, count, name) {
  value <- as.numeric(value)
  if (length(value) != count || any(!is.finite(value)) || any(value < 0) ||
      abs(sum(value) - 1) > 1e-10)
    stop(name, " must be a nonnegative vector summing to one")
  value
}

