# Static integration figures built only from completed CSV artifacts. These
# plots validate the reporting pipeline; smoke values are not scientific
# estimates and every title states that explicitly.

.paper_read_csv <- function(path, required) {
  if (!file.exists(path)) stop("Missing reporting input: ", path)
  value <- utils::read.csv(path, check.names = FALSE)
  if (!all(required %in% names(value)))
    stop("Reporting input has an incompatible schema: ", path)
  value
}

.paper_png <- function(path, width = 1500, height = 1000) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  grDevices::png(path, width = width, height = height, res = 160)
}

.paper_metric_bars <- function(values, standard_errors, labels, title,
                               ylab, color) {
  positions <- graphics::barplot(
    values, names.arg = labels, col = color, border = NA,
    main = title, ylab = ylab, las = 1)
  finite <- is.finite(standard_errors) & standard_errors > 0
  if (any(finite)) graphics::arrows(
    positions[finite], values[finite] - standard_errors[finite],
    positions[finite], values[finite] + standard_errors[finite],
    angle = 90, code = 3, length = .05)
  invisible(positions)
}

plot_e1_smoke_overview <- function(input_dir, figure_dir) {
  data <- .paper_read_csv(
    file.path(input_dir, "metrics-summary.csv"),
    c("method", "ari_mean", "ari_mcse", "cpe_mean", "cpe_mcse",
      "mlpd_mean", "mlpd_mcse", "predictive_coverage_mean",
      "predictive_coverage_mcse"))
  order <- match(c("MCMC", "CAVI", "SVI"), data$method)
  data <- data[order[is.finite(order)], , drop = FALSE]
  colors <- c(MCMC = "#4D4D4D", CAVI = "#0072B2", SVI = "#D55E00")
  path <- file.path(figure_dir, "E1-overview-smoke.png")
  .paper_png(path)
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::par(mfrow = c(2, 2), mar = c(4, 4, 3, 1), oma = c(0, 0, 4, 0))
  .paper_metric_bars(data$ari_mean, data$ari_mcse, data$method,
                     "ARI", "Índice", colors[data$method])
  .paper_metric_bars(data$cpe_mean, data$cpe_mcse, data$method,
                     "Error de coclustering", "CPE", colors[data$method])
  .paper_metric_bars(data$mlpd_mean, data$mlpd_mcse, data$method,
                     "Densidad predictiva retenida", "MLPD",
                     colors[data$method])
  .paper_metric_bars(
    data$predictive_coverage_mean, data$predictive_coverage_mcse,
    data$method, "Cobertura predictiva", "Proporción",
    colors[data$method])
  graphics::mtext("E1 smoke: integración, no resultado científico",
                  outer = TRUE, line = 1, font = 2)
  path
}

plot_e1_parameter_recovery_smoke <- function(input_dir, figure_dir) {
  data <- .paper_read_csv(
    file.path(input_dir, "parameter-recovery-summary.csv"),
    c("method", "parameter", "rmse", "coverage", "mean_interval_width"))
  aggregate_metric <- function(name) stats::aggregate(
    data[[name]], list(method = data$method, parameter = data$parameter),
    function(value) mean(value[is.finite(value)]))
  rmse <- aggregate_metric("rmse")
  coverage <- aggregate_metric("coverage")
  labels <- c(component_coefficients = "Coeficientes",
              component_variances = "Varianzas", mixture_weights = "Pesos")
  methods <- c("MCMC", "CAVI", "SVI")
  matrix_for <- function(frame) {
    output <- matrix(NA_real_, length(labels), length(methods),
                     dimnames = list(names(labels), methods))
    for (i in seq_len(nrow(frame)))
      output[frame$parameter[i], frame$method[i]] <- frame$x[i]
    output
  }
  path <- file.path(figure_dir, "E1-parameter-recovery-smoke.png")
  .paper_png(path, width = 1700, height = 850)
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::par(mfrow = c(1, 2), mar = c(6, 4, 3, 1), oma = c(0, 0, 4, 0))
  colors <- c("#4D4D4D", "#0072B2", "#D55E00")
  rmse_matrix <- matrix_for(rmse)
  graphics::matplot(
    seq_along(labels), rmse_matrix, type = "b", log = "y", xaxt = "n",
    col = colors, pch = c(16, 17, 15), lty = 1, lwd = 2,
    xlab = "", ylab = "RMSE promedio por coordenada (log)",
    main = "Recuperación")
  graphics::axis(1, at = seq_along(labels), labels = unname(labels), las = 2)
  graphics::legend("topleft", methods, col = colors, pch = c(16, 17, 15),
                   lty = 1, lwd = 2, bty = "n")
  graphics::barplot(
    t(matrix_for(coverage)), beside = TRUE, col = colors, border = NA,
    names.arg = unname(labels), las = 2, ylim = c(0, 1),
    ylab = "Cobertura", main = "Intervalos marginales")
  graphics::abline(h = .95, lty = 2, col = "#666666")
  graphics::mtext("E1 smoke: sólo verifica el pipeline",
                  outer = TRUE, line = 1, font = 2)
  path
}

plot_e2_selection_smoke <- function(input_dir, figure_dir) {
  data <- .paper_read_csv(
    file.path(input_dir, "E2", "selection-results.csv"),
    c("minimum_distance", "method", "criterion", "selected_components"))
  criterion <- ifelse(data$criterion == "occupied_components_mode",
                      "ocupados", data$criterion)
  labels <- paste0("d=", data$minimum_distance, "\n", data$method,
                   "-", criterion)
  path <- file.path(figure_dir, "E2-K-selection-smoke.png")
  .paper_png(path, width = 1900, height = 1000)
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::par(mar = c(7, 4, 4, 1))
  graphics::barplot(
    data$selected_components, names.arg = labels, las = 2,
    col = ifelse(data$method == "MCMC", "#4D4D4D",
                 ifelse(data$method == "CAVI", "#0072B2", "#D55E00")),
    border = NA, ylab = "K seleccionado",
    main = "E2 smoke: selección de K (no confirmatoria)")
  graphics::abline(h = 3, lty = 2, lwd = 2)
  path
}

plot_e3_comparison_smoke <- function(input_dir, figure_dir) {
  data <- .paper_read_csv(
    file.path(input_dir, "E3", "svi-comparison.csv"),
    c("delta_elbo", "delta_mlpd", "delta_ari"))
  values <- unlist(data[1L, c("delta_elbo", "delta_mlpd", "delta_ari")])
  path <- file.path(figure_dir, "E3-SVI-comparison-smoke.png")
  .paper_png(path, width = 1700, height = 650)
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::par(mfrow = c(1, 3), mar = c(4, 4, 4, 1), oma = c(0, 0, 4, 0))
  titles <- c("ELBO: CAVI - SVI", "MLPD: SVI - CAVI", "ARI: SVI - CAVI")
  for (i in seq_along(values)) {
    graphics::barplot(values[i], col = "#D55E00", border = NA,
                      main = titles[i], ylab = "Diferencia")
    graphics::abline(h = 0, lty = 2)
  }
  graphics::mtext("E3 smoke: una celda, no tuning científico",
                  outer = TRUE, line = 1, font = 2)
  path
}

plot_e4_scaling_smoke <- function(input_dir, figure_dir) {
  data <- .paper_read_csv(
    file.path(input_dir, "E4", "fit-results.csv"),
    c("observations", "method", "runtime"))
  path <- file.path(figure_dir, "E4-scaling-smoke.png")
  .paper_png(path)
  on.exit(grDevices::dev.off(), add = TRUE)
  methods <- c("MCMC", "CAVI", "SVI")
  colors <- c(MCMC = "#4D4D4D", CAVI = "#0072B2", SVI = "#D55E00")
  symbols <- c(MCMC = 16, CAVI = 17, SVI = 15)
  graphics::plot(
    data$observations, data$runtime, log = "xy", type = "n",
    xlab = "Observaciones", ylab = "Segundos",
    main = "E4 smoke: escalamiento observado (no científico)")
  for (method in methods) {
    values <- data[data$method == method, , drop = FALSE]
    values <- values[order(values$observations), ]
    graphics::lines(values$observations, values$runtime, type = "b",
                    pch = symbols[method], lwd = 2, col = colors[method])
  }
  graphics::legend("topright", methods, col = colors, lwd = 2,
                   pch = symbols, bty = "n")
  path
}

plot_e5_k_smoke <- function(input_dir, figure_dir) {
  data <- .paper_read_csv(
    file.path(input_dir, "E5", "K-comparison.csv"),
    c("component_count", "method", "waic_group_marginal", "cv_mlpd"))
  path <- file.path(figure_dir, "E5-K-comparison-smoke.png")
  .paper_png(path, width = 1700, height = 800)
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::par(mfrow = c(1, 2), mar = c(4, 4, 4, 1), oma = c(0, 0, 4, 0))
  methods <- c("MCMC", "CAVI", "SVI")
  colors <- c(MCMC = "#4D4D4D", CAVI = "#0072B2", SVI = "#D55E00")
  for (metric in c("waic_group_marginal", "cv_mlpd")) {
    graphics::plot(
      range(data$component_count), range(data[[metric]], finite = TRUE),
      type = "n", xlab = "K", ylab = if (metric == "cv_mlpd") "CV MLPD" else
        "WAIC marginal por finca",
      main = if (metric == "cv_mlpd") "Predicción retenida" else "WAIC")
    for (method in methods) {
      values <- data[data$method == method, , drop = FALSE]
      values <- values[order(values$component_count), ]
      graphics::lines(values$component_count, values[[metric]], type = "b",
                      pch = 16, lwd = 2, col = colors[method])
    }
    graphics::legend("topright", methods, col = colors, lwd = 2, pch = 16,
                     bty = "n")
  }
  graphics::mtext("E5 smoke: K=1,2; no usar como conclusión",
                  outer = TRUE, line = 1, font = 2)
  path
}

generate_smoke_figures <- function(input_dir, figure_dir = file.path(input_dir, "figures")) {
  builders <- list(
    E1_overview = plot_e1_smoke_overview,
    E1_parameter_recovery = plot_e1_parameter_recovery_smoke,
    E2_K_selection = plot_e2_selection_smoke,
    E3_SVI_comparison = plot_e3_comparison_smoke,
    E4_scaling = plot_e4_scaling_smoke,
    E5_K_comparison = plot_e5_k_smoke)
  paths <- vapply(builders, function(builder) builder(input_dir, figure_dir),
                  character(1L))
  inputs <- c(
    file.path(input_dir, "metrics-summary.csv"),
    file.path(input_dir, "parameter-recovery-summary.csv"),
    file.path(input_dir, "E2", "selection-results.csv"),
    file.path(input_dir, "E3", "svi-comparison.csv"),
    file.path(input_dir, "E4", "fit-results.csv"),
    file.path(input_dir, "E5", "K-comparison.csv"))
  input_fingerprint <- .paper_fingerprint(as.list(stats::setNames(
    unname(tools::md5sum(inputs)), substring(inputs, nchar(input_dir) + 2L))))
  manifest <- data.frame(
    figure = names(paths), path = file.path("figures", basename(paths)),
    input_fingerprint = input_fingerprint,
    generated_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    profile = "smoke", scientific_result = FALSE,
    stringsAsFactors = FALSE)
  utils::write.csv(manifest, file.path(figure_dir, "figure-manifest.csv"),
                   row.names = FALSE)
  manifest
}

.e3_configuration_label <- function(data) {
  ifelse(data$method == "CAVI", "CAVI",
         sprintf("b=%d\nk=%.1f, t=%d", data$minibatch_size,
                 data$exponent, as.integer(data$offset)))
}

plot_e3_screen_quality <- function(input_dir, figure_dir) {
  results <- .paper_read_csv(
    file.path(input_dir, "fit-results.csv"),
    c("task_id", "method", "minibatch_size", "exponent", "offset",
      "final_elbo", "mlpd", "ari", "cpe"))
  comparison <- .paper_read_csv(
    file.path(input_dir, "svi-comparison.csv"),
    c("task_id", "delta_elbo", "delta_mlpd", "ari", "cpe", "mlpd"))
  comparison <- comparison[order(-comparison$mlpd), , drop = FALSE]
  source_rows <- match(comparison$task_id, results$task_id)
  labels <- .e3_configuration_label(results[source_rows, , drop = FALSE])
  reference <- results[results$method == "CAVI", , drop = FALSE]
  colors <- grDevices::hcl.colors(nrow(comparison), "Oranges", rev = TRUE)
  path <- file.path(figure_dir, "E3-screen-quality.png")
  .paper_png(path, width = 2100, height = 1300)
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::par(mfrow = c(2, 2), mar = c(7, 4.5, 3, 1),
                oma = c(0, 0, 4, 0))
  draw <- function(values, title, ylab, reference_value = 0,
                   lower_is_better = FALSE) {
    positions <- graphics::barplot(
      values, names.arg = labels, las = 2, col = colors, border = NA,
      main = title, ylab = ylab)
    graphics::abline(h = reference_value, lty = 2, lwd = 2,
                     col = "#0072B2")
    direction <- if (lower_is_better) "menor es mejor" else "mayor es mejor"
    graphics::mtext(direction, side = 3, line = .2, cex = .75,
                    col = "#555555")
    invisible(positions)
  }
  draw(comparison$delta_elbo, "Brecha de ELBO",
       "CAVI - SVI", 0, lower_is_better = TRUE)
  draw(comparison$delta_mlpd, "Densidad predictiva retenida",
       "SVI - CAVI", 0)
  draw(comparison$ari, "Recuperación de partición", "ARI",
       reference$ari[[1L]])
  draw(comparison$cpe, "Error de coclustering", "CPE",
       reference$cpe[[1L]], lower_is_better = TRUE)
  graphics::mtext(
    "E3 screen: sensibilidad preliminar; línea azul = referencia CAVI",
    outer = TRUE, line = 1, font = 2)
  path
}

plot_e3_screen_resources <- function(input_dir, figure_dir) {
  data <- .paper_read_csv(
    file.path(input_dir, "fit-results.csv"),
    c("method", "minibatch_size", "exponent", "offset", "runtime",
      "worker_wall_seconds", "measured_peak_rss_mb", "result_file_mb",
      "mlpd"))
  data <- data[order(data$method != "CAVI", -data$mlpd), , drop = FALSE]
  labels <- .e3_configuration_label(data)
  colors <- ifelse(data$method == "CAVI", "#0072B2", "#D55E00")
  path <- file.path(figure_dir, "E3-screen-resources.png")
  .paper_png(path, width = 2100, height = 1300)
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::par(mfrow = c(2, 2), mar = c(7, 4.5, 3, 1),
                oma = c(0, 0, 4, 0))
  panels <- list(
    list(data$worker_wall_seconds, "Tiempo total por proceso", "Segundos"),
    list(data$runtime, "Tiempo del ajuste", "Segundos"),
    list(data$measured_peak_rss_mb, "Memoria máxima observada", "MiB RSS"),
    list(data$result_file_mb, "Tamaño del artefacto", "MiB"))
  for (panel in panels) graphics::barplot(
    panel[[1L]], names.arg = labels, las = 2, col = colors, border = NA,
    main = panel[[2L]], ylab = panel[[3L]])
  graphics::mtext(
    "E3 screen: perfil secuencial, un proceso y un ajuste a la vez",
    outer = TRUE, line = 1, font = 2)
  path
}

plot_e3_screen_convergence <- function(input_dir, figure_dir) {
  results <- .paper_read_csv(
    file.path(input_dir, "fit-results.csv"),
    c("task_id", "method", "minibatch_size", "exponent", "offset",
      "final_elbo", "path"))
  reference <- results$final_elbo[results$method == "CAVI"]
  svi <- results[results$method == "SVI", , drop = FALSE]
  labels <- .e3_configuration_label(svi)
  histories <- lapply(seq_len(nrow(svi)), function(index) {
    artifact <- readRDS(file.path(input_dir, svi$path[[index]]))
    artifact$fit$diagnostics$history
  })
  colors <- grDevices::hcl.colors(nrow(svi), "Dark 3")
  gaps <- unlist(lapply(histories, function(history)
    pmax(reference[[1L]] - history$smoothed_objective, .01)))
  iterations <- unlist(lapply(histories, `[[`, "iteration"))
  finite_gaps <- gaps[is.finite(gaps)]
  if (!length(finite_gaps))
    stop("E3 screen histories contain no finite smoothed objectives")
  path <- file.path(figure_dir, "E3-screen-convergence.png")
  .paper_png(path, width = 1900, height = 1100)
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::par(mar = c(5, 5, 4, 2))
  graphics::plot(
    range(iterations), range(finite_gaps), type = "n", log = "y",
    xlab = "Iteración SVI", ylab = "CAVI ELBO - ELBO suavizado SVI (log)",
    main = "E3 screen: trayectorias respecto a la referencia CAVI")
  for (index in seq_along(histories)) {
    history <- histories[[index]]
    gap <- pmax(reference[[1L]] - history$smoothed_objective, .01)
    graphics::lines(history$iteration, gap, col = colors[[index]], lwd = 2)
    graphics::points(tail(history$iteration, 1L), tail(gap, 1L),
                     col = colors[[index]], pch = 16)
  }
  graphics::abline(v = 1200, lty = 3, col = "#555555")
  graphics::legend("bottomleft", inset = .015, legend = labels,
                   col = colors, lwd = 2, bty = "o", bg = "white",
                   ncol = 3, cex = .62)
  graphics::mtext(
    "Línea vertical: primera parada posible; estabilidad no garantiza buen óptimo",
    side = 3, line = .3, cex = .8)
  path
}

generate_e3_screen_figures <- function(
    input_dir, figure_dir = file.path(input_dir, "figures")) {
  builders <- list(
    E3_screen_quality = plot_e3_screen_quality,
    E3_screen_resources = plot_e3_screen_resources,
    E3_screen_convergence = plot_e3_screen_convergence)
  paths <- vapply(builders, function(builder) builder(input_dir, figure_dir),
                  character(1L))
  inputs <- c(
    file.path(input_dir, "fit-results.csv"),
    file.path(input_dir, "svi-comparison.csv"),
    file.path(input_dir, "screen-review.csv"),
    list.files(file.path(input_dir, "tasks"), pattern = "\\.rds$",
               full.names = TRUE))
  inputs <- inputs[file.exists(inputs) & !grepl("unmeasured", inputs, fixed = TRUE)]
  input_fingerprint <- .paper_fingerprint(as.list(stats::setNames(
    unname(tools::md5sum(inputs)), substring(inputs, nchar(input_dir) + 2L))))
  manifest <- data.frame(
    figure = names(paths), path = file.path("figures", basename(paths)),
    input_fingerprint = input_fingerprint,
    generated_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    profile = "pilot-screen", scientific_result = FALSE,
    stringsAsFactors = FALSE)
  utils::write.csv(manifest, file.path(figure_dir, "figure-manifest.csv"),
                   row.names = FALSE)
  manifest
}
