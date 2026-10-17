.backend_state_field <- function(model, method) {
  paste0(tolower(model), if (method == "MCMC") "_mcmc_object" else
    "_vi_object")
}

.initialization_name_map <- function(model, method) {
  variational <- method != "MCMC"
  switch(model,
    LRM = if (!variational) c(
      coefficients = "Beta",
      residual_variance = "sigma2") else c(
      coefficient_mean = "MuBeta",
      coefficient_covariance = "SigmaBeta",
      residual_variance_shape = "asigma2",
      residual_variance_rate = "bsigma2"),
    HLRM = if (!variational) c(
      group_coefficients = "Betaj",
      variance_scale = "xi2",
      group_residual_variances = "sigma2j",
      population_coefficients = "Beta",
      population_covariance = "Sigma") else c(
      group_coefficient_means = "MuBetaj",
      group_coefficient_covariances = "SigmaBetaj",
      population_coefficient_mean = "MuBeta",
      population_coefficient_covariance = "SigmaBeta",
      population_covariance_df = "nuS",
      population_covariance_scale = "SS",
      group_residual_variance_shapes = "aj",
      group_residual_variance_rates = "bj",
      variance_scale_shape = "axi2",
      variance_scale_rate = "bxi2"),
    CHLRM = if (!variational) c(
      mixture_weights = "omega",
      group_allocations = "z",
      population_coefficients = "Beta",
      population_covariance = "Sigma",
      variance_scale = "xi2",
      component_coefficients = "Betak",
      component_residual_variances = "sigma2k") else c(
      allocation_probabilities = "rho",
      mixture_concentration = "alpha",
      component_coefficient_means = "MuBetak",
      component_coefficient_covariances = "SigmaBetak",
      population_coefficient_mean = "MuBeta",
      population_coefficient_covariance = "SigmaBeta",
      population_covariance_df = "nuS",
      population_covariance_scale = "SS",
      component_residual_variance_shapes = "ak",
      component_residual_variance_rates = "bk",
      variance_scale_shape = "axi2",
      variance_scale_rate = "bxi2"),
    stop("Unknown model family"))
}

.validate_initialization_state <- function(values, model, method) {
  positive <- function(x, name) {
    if (any(x <= 0)) stop(name, " must contain positive values")
  }
  spd_list <- function(x, name) {
    for (index in seq_along(x))
      chol_spd(x[[index]], paste(name, index))
  }
  if (method == "MCMC") {
    if (model == "LRM") {
      positive(values$sigma2, "residual_variance")
    } else {
      chol_spd(values$Sigma, "population_covariance")
      positive(values$xi2, "variance_scale")
      if (model == "HLRM") {
        positive(values$sigma2j, "group_residual_variances")
      } else {
        positive(values$sigma2k, "component_residual_variances")
        positive(values$omega, "mixture_weights")
        if (abs(sum(values$omega) - 1) > 1e-10)
          stop("mixture_weights must sum to one")
        if (any(values$z != floor(values$z)) ||
            any(values$z < 1 | values$z > length(values$omega)))
          stop("group_allocations must contain valid component integers")
      }
    }
    return(invisible(values))
  }
  chol_spd(values$SigmaBeta, "population_coefficient_covariance")
  if (model == "LRM") {
    positive(values$asigma2, "residual_variance_shape")
    positive(values$bsigma2, "residual_variance_rate")
    return(invisible(values))
  }
  chol_spd(values$SS, "population_covariance_scale")
  if (values$nuS <= nrow(values$SS) - 1L)
    stop("population_covariance_df must exceed coefficient_count - 1")
  positive(values$axi2, "variance_scale_shape")
  positive(values$bxi2, "variance_scale_rate")
  if (model == "HLRM") {
    spd_list(values$SigmaBetaj, "group_coefficient_covariance")
    positive(values$aj, "group_residual_variance_shapes")
    positive(values$bj, "group_residual_variance_rates")
  } else {
    spd_list(values$SigmaBetak, "component_coefficient_covariance")
    positive(values$ak, "component_residual_variance_shapes")
    positive(values$bk, "component_residual_variance_rates")
    positive(values$alpha, "mixture_concentration")
    if (any(values$rho < 0) ||
        any(abs(rowSums(values$rho) - 1) > 1e-10))
      stop("allocation_probabilities must be nonnegative with rows summing to one")
  }
  invisible(values)
}

.apply_values_initialization <- function(backend, model, method, values) {
  field <- .backend_state_field(model, method)
  state <- backend[[field]]
  template <- state$params
  name_map <- .initialization_name_map(model, method)
  if (!setequal(names(values), names(name_map)))
    stop("initialization values must contain exactly these public fields: ",
         paste(names(name_map), collapse = ", "))
  values <- values[names(name_map)]
  names(values) <- unname(name_map)
  values <- values[names(template)]
  if (method != "MCMC") {
    .validate_vi_state(values, template)
  } else {
    candidate <- unlist(values, recursive = TRUE, use.names = FALSE)
    if (!is.numeric(candidate) || any(!is.finite(candidate)))
      stop("MCMC initialization values must contain only finite numeric parameters")
    if (!identical(.vi_state_signature(values), .vi_state_signature(template)))
      stop("MCMC initialization values have incompatible dimensions")
  }
  .validate_initialization_state(values, model, method)
  state$params <- unserialize(serialize(values, NULL))
  backend[[field]] <- state
  invisible(backend)
}

.apply_ols_initialization <- function(backend, model, method, ols) {
  field <- .backend_state_field(model, method)
  state <- backend[[field]]
  if (model == "LRM") {
    if (method == "MCMC") {
      state$params$Beta <- ols$coefficients
      state$params$sigma2 <- ols$residual_variance
    } else {
      state$params$MuBeta <- ols$coefficients
      shape <- max(state$params$asigma2, 2)
      state$params$asigma2 <- shape
      state$params$bsigma2 <- (shape - 1) * ols$residual_variance
    }
  } else if (model == "HLRM") {
    if (method == "MCMC") {
      state$params$Betaj <- ols$coefficients
      state$params$sigma2j <- ols$variances
      state$params$Beta <- colMeans(ols$coefficients)
      state$params$xi2 <- stats::median(ols$variances)
    } else {
      state$params$MuBetaj <- ols$coefficients
      state$params$MuBeta <- colMeans(ols$coefficients)
    }
  } else stop("initialization_ols() is not complete for CHLRM; use initialization_kmeans()")
  backend[[field]] <- state
  invisible(backend)
}

.apply_kmeans_initialization <- function(backend, model, method, ols,
                                         component_count) {
  if (model != "CHLRM")
    stop("initialization_kmeans() applies only to CHLRM")
  # The variational constructor consumes Betaj_ols directly. MCMC needs the equivalent
  # state written explicitly because its constructor otherwise draws from the
  # prior.
  if (method != "MCMC") return(invisible(backend))
  clustering <- tryCatch(
    stats::kmeans(ols$coefficients, centers = component_count),
    error = function(error) stop(
      "k-means initialization failed: ", conditionMessage(error)))
  state <- backend$chlrm_mcmc_object
  z <- as.integer(clustering$cluster)
  centers <- unname(clustering$centers)
  component_variances <- vapply(seq_len(component_count), function(k) {
    stats::median(ols$variances[z == k])
  }, numeric(1L))
  state$params$z <- z
  state$params$omega <- tabulate(z, nbins = component_count) / length(z)
  state$params$Betak <- centers
  state$params$sigma2k <- component_variances
  state$params$Beta <- colMeans(centers)
  if (component_count > 1L) {
    dispersion <- stats::cov(centers)
    if (is.null(dim(dispersion))) dispersion <- matrix(dispersion, 1L, 1L)
    state$params$Sigma <- symmetrize_numeric(
      dispersion + 1e-8 * diag(ncol(centers)),
      "k-means component dispersion")
  }
  state$params$xi2 <- stats::median(component_variances)
  backend$chlrm_mcmc_object <- state
  invisible(backend)
}

.canonical_public_method <- function(method) {
  switch(method, mcmc = "MCMC", cavi = "CAVI", svi = "SVI")
}

.apply_initialization <- function(backend, model, method, initialization,
                                  lrm_ols = NULL, grouped_ols = NULL,
                                  component_count = NULL) {
  if (initialization$strategy == "values")
    .apply_values_initialization(
      backend, model, method, initialization$values)
  if (initialization$strategy == "ols")
    .apply_ols_initialization(
      backend, model, method, if (model == "LRM") lrm_ols else grouped_ols)
  if (initialization$strategy == "kmeans")
    .apply_kmeans_initialization(
      backend, model, method, grouped_ols, component_count)
  invisible(backend)
}

