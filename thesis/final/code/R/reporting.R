repo_root <- function() {
  root <- Sys.getenv("VB_HLRM_ROOT", unset = "")
  if (nzchar(root)) return(normalizePath(root))
  normalizePath(".")
}

ensure_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  normalizePath(path)
}

thesis_hierarchical_prior <- function(model = c("HLRM", "CHLRM"),
                                      coefficient_mean,
                                      coefficient_covariance,
                                      coefficient_dispersion_scale,
                                      component_count = NULL) {
  model <- match.arg(model)
  p <- length(coefficient_mean)
  common <- list(
    coefficient_mean = coefficient_mean,
    coefficient_covariance = coefficient_covariance,
    coefficient_dispersion_df = p + 2,
    coefficient_dispersion_scale = coefficient_dispersion_scale,
    residual_degrees_of_freedom = 6,
    residual_variance_scale_shape = 2,
    residual_variance_scale_rate = .1)
  if (model == "HLRM") return(do.call(hlrm_prior, common))
  if (is.null(component_count)) stop("component_count is required for CHLRM")
  do.call(chlrm_prior, c(common, list(
    mixture_concentration = rep(1 / component_count, component_count))))
}

thesis_execution_profile <- function() {
  profile <- Sys.getenv("VB_HLRM_THESIS_PROFILE", unset = "standard")
  if (!profile %in% c("smoke", "standard"))
    stop("VB_HLRM_THESIS_PROFILE must be 'smoke' or 'standard'")
  profile
}

thesis_smoke_controls <- function(group_count) {
  list(
    evaluation_draws = 4L,
    mcmc = mcmc_control(
      warmup = 2, post_warmup_draws = 8, chains = 2,
      evaluation_draws = 4,
      stopping = mcmc_stopping(
        min_post_warmup_draws = 4, check_every = 4,
        rhat_max = 100, ess_bulk_min = .1, ess_tail_min = .1,
        relative_mcse_mean_max = 100, relative_mcse_sd_max = 100),
      max_memory_mb = 256),
    cavi = cavi_control(
      max_iterations = 4, evaluation_draws = 4,
      stopping = cavi_stopping(
        objective_tolerance = 1e20, parameter_tolerance = 1e20,
        min_iterations = 2, patience = 1)),
    svi = svi_control(
      max_iterations = 4, minibatch_size = min(4L, group_count),
      learning_rate = power_schedule(19, .7), evaluation_draws = 4,
      stopping = svi_stopping(
        objective_tolerance = 1e20, parameter_tolerance = 1e20,
        min_iterations = 2, patience = 1, check_every = 1,
        smoothing_window = 1)))
}

rmse <- function(x, y) sqrt(mean((as.numeric(x) - as.numeric(y))^2))

ari <- function(x, y) {
  tab <- table(x, y); c2 <- function(v) v * (v - 1) / 2
  a <- sum(c2(tab)); b <- sum(c2(rowSums(tab))); cc <- sum(c2(colSums(tab))); n2 <- c2(sum(tab))
  expected <- b * cc / n2
  (a - expected) / (0.5 * (b + cc) - expected)
}

metric_row <- function(dataset, method, fit, in_sample, predictive) {
  if (!inherits(in_sample, "vb_in_sample_fit") ||
      !identical(in_sample$schema_version, "in-sample-fit-1"))
    stop("in_sample must use the current in-sample-fit-1 contract")
  if (!inherits(predictive, "vb_predictive_criteria") ||
      !identical(predictive$schema_version, "predictive-criteria-2"))
    stop("predictive must use the current predictive-criteria-2 contract")
  if (!identical(in_sample$model, predictive$model) ||
      !identical(in_sample$method, predictive$method))
    stop("Evaluation results must describe the same model and method")
  hlrm <- identical(predictive$model, "HLRM")
  criteria <- predictive$criteria
  obs <- if (hlrm) criteria$waic_observation_conditional else
    criteria$waic_observation
  grouped <- criteria$waic_group
  conditional_group <- criteria$waic_group_conditional
  marginal_waic <- criteria$waic_group_marginal
  data.frame(
    dataset = dataset, method = method,
    in_sample_schema = in_sample$schema_version,
    predictive_schema = predictive$schema_version,
    predictive_target = predictive$predictive_target,
    WAIC_observation = if (is.null(obs)) NA_real_ else obs$WAIC,
    WAIC_observation_se = if (is.null(obs)) NA_real_ else obs$se,
    WAIC_observation_likelihood = if (is.null(obs)) NA_character_ else
      if (hlrm) "local_parameter_conditional" else "gaussian_regression",
    WAIC_group = if (is.null(grouped)) NA_real_ else grouped$WAIC,
    WAIC_group_se = if (is.null(grouped)) NA_real_ else grouped$se,
    WAIC_group_conditional = if (is.null(conditional_group)) NA_real_ else
      conditional_group$WAIC,
    WAIC_group_conditional_se = if (is.null(conditional_group)) NA_real_ else
      conditional_group$se,
    WAIC_group_marginal = if (!is.null(marginal_waic)) marginal_waic$WAIC else NA_real_,
    WAIC_group_marginal_se = if (!is.null(marginal_waic)) marginal_waic$se else NA_real_,
    criterion_source = predictive$metadata$draw_source,
    MSE = in_sample$statistics$mse, R2 = in_sample$statistics$r2,
    time_seconds = fit$timing$total_seconds,
    stringsAsFactors = FALSE
  )
}

write_parameter_summary <- function(hat, method, file) {
  out <- as.data.frame(tidy_parameters(hat))
  out$method <- method
  write.csv(out, file, row.names = FALSE)
  invisible(out)
}

plot_elbo_comparison <- function(fits, file, title) {
  png(file, width = 1500, height = 900, res = 160)
  on.exit(dev.off(), add = TRUE)
  cols <- c(CAVI = "#0072B2", SVI = "#D55E00")
  xr <- range(unlist(lapply(fits, function(x) x$objective$evaluation_iterations)))
  yr <- range(unlist(lapply(fits, function(x) x$objective$elbo)), finite = TRUE)
  plot(NA, xlim = xr, ylim = yr, xlab = "Iteración", ylab = "ELBO", main = title)
  for (method in names(fits)) lines(fits[[method]]$objective$evaluation_iterations,
                                    fits[[method]]$objective$elbo,
                                    col = cols[method], lwd = 2,
                                    lty = if (method == "CAVI") 1 else 2)
  legend("bottomright", names(fits), col = cols[names(fits)], lty = c(1, 2), lwd = 2, bty = "n")
}

plot_predictive_density <- function(y, yrep, file, title) {
  observed <- unlist(y, use.names = FALSE)
  predictive <- as.numeric(yrep)
  dobs <- density(observed); dpred <- density(predictive)
  png(file, width = 1500, height = 900, res = 160)
  on.exit(dev.off(), add = TRUE)
  plot(dobs, col = "#333333", lwd = 3,
       xlim = range(dobs$x, dpred$x), ylim = c(0, 1.08 * max(dobs$y, dpred$y)),
       xlab = "Respuesta", ylab = "Densidad", main = title)
  lines(dpred, col = "#0072B2", lwd = 2.5, lty = 2)
  legend("topright", c("Observado", "Posterior predictiva"),
         col = c("#333333", "#0072B2"), lty = c(1, 2), lwd = c(3, 2.5), bty = "n")
}

plot_psm <- function(psm, file, title) {
  ord <- stats::hclust(stats::as.dist(1 - psm))$order
  png(file, width = 1100, height = 1000, res = 160)
  on.exit(dev.off(), add = TRUE)
  image(seq_len(nrow(psm)), seq_len(ncol(psm)), psm[ord, ord],
        col = rev(hcl.colors(40, "Blues 3")), xlab = "Grupo ordenado",
        ylab = "Grupo ordenado", main = title, useRaster = TRUE)
}

plot_regression_lines <- function(data, hats, file, title) {
  groups <- sort(unique(data$farm))
  png(file, width = 1900, height = 1400, res = 160)
  on.exit(dev.off(), add = TRUE)
  par(mfrow = c(4, 6), mar = c(2.2, 2.2, 2, .5), oma = c(2, 2, 2, 0))
  cols <- c(MCMC = "#333333", CAVI = "#0072B2", SVI = "#D55E00")
  for (j in seq_along(groups)) {
    d <- data[data$farm == groups[j], ]
    plot(d$N, d$size, pch = 16, cex = .5, col = adjustcolor("#555555", .45),
         xlab = "", ylab = "", main = paste("Finca", groups[j]))
    xg <- range(d$N)
    for (method in names(hats)) {
      b <- hats[[method]]$parameters$group_coefficients$mean[j, ]
      lines(xg, b[1] + b[2] * xg, col = cols[method], lwd = 2,
            lty = match(method, names(hats)))
    }
  }
  mtext("Nitrógeno", side = 1, outer = TRUE); mtext("Tamaño", side = 2, outer = TRUE)
  mtext(title, side = 3, outer = TRUE)
  legend("top", legend = names(hats), col = cols[names(hats)], lty = seq_along(hats),
         lwd = 2, bty = "n", xpd = NA)
}
