library(VariationalBayesHLRM)
source("thesis/final/code/R/reporting.R")

out <- ensure_dir(file.path(
  repo_root(), "thesis", "final", "code", "runs", "04-svi-step-tuning", "output"))
farms <- read.delim(file.path(repo_root(), "data", "farm_trees.txt"))
d <- grouped_regression_data(farms, "size", "N", "farm")
K <- 4L
model_spec <- chlrm_model(K, 2)
prior <- thesis_hierarchical_prior(
  "CHLRM", coefficient_mean = c(100, 0),
  coefficient_covariance = diag(c(50^2, 5^2)),
  coefficient_dispersion_scale = diag(c(20^2, 1^2)),
  component_count = K)
scaling <- data_scaling()
m <- length(d$y)
profile <- thesis_execution_profile()
learning_rate_at <- function(iteration, offset, exponent)
  (iteration + offset)^(-exponent)

# Screening economico: el batch queda fijo para aislar el efecto de la
# secuencia de Robbins--Monro. Cada corrida procesa igual numero de epocas.
grid <- if (profile == "smoke") data.frame(
  rho_initial = .05, kappa = .7, seed = 4101L) else expand.grid(
    rho_initial = c(0.01, 0.05, 0.10),
    kappa = c(0.55, 0.70, 0.85),
    seed = c(4101L, 4102L, 4103L),
    KEEP.OUT.ATTRS = FALSE)
batch <- max(K, ceiling(0.50 * m))
epochs <- if (profile == "smoke") 2L else 300L
iterations <- ceiling(epochs * m / batch)
evaluation_draws <- if (profile == "smoke") 4L else 50L

rows <- vector("list", nrow(grid))
traces <- vector("list", nrow(grid))
evaluations <- vector("list", nrow(grid))

for (i in seq_len(nrow(grid))) {
  g <- grid[i, ]
  tau <- g$rho_initial^(-1 / g$kappa) - 1
  fit <- fit_model(
    model_spec, d, prior, "svi", initialization_kmeans(),
    svi_control(
      max_iterations = iterations, minibatch_size = batch,
      learning_rate = power_schedule(tau, g$kappa),
      evaluation_draws = evaluation_draws,
      stopping = svi_stopping(
        objective_tolerance = 0, min_iterations = iterations,
        patience = 1, check_every = ceiling(m / batch),
        smoothing_window = 1)),
    scaling = scaling, seed = g$seed, verbose = FALSE)
  hat <- summarize_posterior(fit)
  fit_eval <- evaluate_in_sample_fit(fit)
  criteria <- predictive_criteria(
    fit, n_draws = evaluation_draws, seed = as.integer(g$seed + 10000L),
    return_log_likelihood = TRUE)
  final_elbo <- tail(fit$objective$elbo[is.finite(fit$objective$elbo)], 1L)
  effective_sizes <- colSums(responsibilities(fit))

  rows[[i]] <- data.frame(
    rho_initial = g$rho_initial, kappa = g$kappa, tau = tau,
    profile = profile,
    seed = g$seed, batch = batch, iterations = iterations,
    effective_epochs = iterations * batch / m,
    final_step = learning_rate_at(iterations, tau, g$kappa),
    final_elbo = final_elbo,
    WAIC_group_marginal = criteria$criteria$waic_group_marginal$WAIC,
    WAIC_group_marginal_se = criteria$criteria$waic_group_marginal$se,
    MSE = fit_eval$statistics$mse, R2 = fit_eval$statistics$r2,
    mean_entropy = mean(classification_uncertainty(fit)$entropy),
    min_effective_component_size = min(effective_sizes),
    occupied_MAP_components = length(unique(hat$allocations$labels)),
    finite = all(is.finite(c(final_elbo,
                             criteria$criteria$waic_group_marginal$WAIC,
                             fit_eval$statistics$mse))),
    time_seconds = fit$timing$total_seconds
  )
  traces[[i]] <- data.frame(
    rho_initial = g$rho_initial, kappa = g$kappa, seed = g$seed,
    iteration = fit$objective$evaluation_iterations, elbo = fit$objective$elbo
  )
  evaluations[[i]] <- list(
    configuration = g, in_sample = fit_eval, predictive = criteria)
  message(sprintf("[%02d/%02d] rho1=%.2f kappa=%.2f seed=%d WAIC=%.2f",
                  i, nrow(grid), g$rho_initial, g$kappa, g$seed,
                  criteria$criteria$waic_group_marginal$WAIC))
}

raw <- do.call(rbind, rows)
trace <- do.call(rbind, traces)
write.csv(raw, file.path(out, "screening-runs.csv"), row.names = FALSE)
write.csv(trace, file.path(out, "elbo-traces.csv"), row.names = FALSE)
saveRDS(evaluations,
        file.path(out, "evaluation-with-log-likelihood.rds"),
        compress = "xz")

key <- interaction(raw$rho_initial, raw$kappa, drop = TRUE)
aggregate_one <- function(idx) {
  z <- raw[idx, ]
  data.frame(
    rho_initial = z$rho_initial[1], kappa = z$kappa[1], tau = z$tau[1],
    runs = nrow(z), valid_runs = sum(z$finite),
    profile = profile,
    WAIC_mean = mean(z$WAIC_group_marginal),
    WAIC_sd = if (nrow(z) > 1L) sd(z$WAIC_group_marginal) else 0,
    ELBO_mean = mean(z$final_elbo),
    ELBO_sd = if (nrow(z) > 1L) sd(z$final_elbo) else 0,
    MSE_mean = mean(z$MSE),
    MSE_sd = if (nrow(z) > 1L) sd(z$MSE) else 0,
    R2_mean = mean(z$R2),
    entropy_mean = mean(z$mean_entropy),
    min_component_mean = mean(z$min_effective_component_size),
    occupied_components_min = min(z$occupied_MAP_components),
    time_mean = mean(z$time_seconds)
  )
}
summary <- do.call(rbind, lapply(split(seq_len(nrow(raw)), key), aggregate_one))
summary$WAIC_rank <- rank(summary$WAIC_mean, ties.method = "min")
summary$ELBO_rank <- rank(-summary$ELBO_mean, ties.method = "min")
summary$stability_rank <- rank(summary$WAIC_sd, ties.method = "min")
summary$combined_rank <- summary$WAIC_rank + summary$ELBO_rank + summary$stability_rank
summary <- summary[order(summary$combined_rank, summary$WAIC_mean), ]
rownames(summary) <- NULL
write.csv(summary, file.path(out, "screening-summary.csv"), row.names = FALSE)
write.csv(summary[1, ], file.path(out, "recommended-setting.csv"), row.names = FALSE)

cols <- c("#0072B2", "#009E73", "#D55E00")
png(file.path(out, "waic-by-step-size.png"), width = 1500, height = 900, res = 160)
boxplot(WAIC_group_marginal ~ interaction(rho_initial, kappa), data = raw,
        col = rep(cols, each = 3), las = 2,
        xlab = expression(paste(rho[1], "  /  ", kappa)), ylab = "WAIC marginal por finca",
        main = "CHLRM SVI: sensibilidad al step size")
dev.off()

png(file.path(out, "elbo-by-step-size.png"), width = 1500, height = 900, res = 160)
plot(NA, xlim = range(trace$iteration), ylim = range(trace$elbo, finite = TRUE),
     xlab = "Iteracion", ylab = "ELBO", main = "Trayectorias ELBO del screening SVI")
configs <- unique(trace[c("rho_initial", "kappa")])
palette <- hcl.colors(nrow(configs), "Dark 3")
for (i in seq_len(nrow(configs))) {
  take <- trace$rho_initial == configs$rho_initial[i] & trace$kappa == configs$kappa[i]
  for (s in unique(trace$seed[take])) {
    q <- trace[take & trace$seed == s, ]
    lines(q$iteration, q$elbo, col = adjustcolor(palette[i], 0.55), lwd = 1)
  }
}
legend("bottomright", sprintf("rho1=%.2f, k=%.2f", configs$rho_initial, configs$kappa),
       col = palette, lty = 1, lwd = 2, cex = .75, bty = "n")
dev.off()

late <- trace[trace$iteration >= 25L, ]
if (!nrow(late)) late <- trace
png(file.path(out, "elbo-final-zoom.png"), width = 1500, height = 900, res = 160)
plot(NA, xlim = range(late$iteration), ylim = range(late$elbo, finite = TRUE),
     xlab = "Iteracion", ylab = "ELBO",
     main = "CHLRM SVI: separacion entre optimos locales (desde iteracion 25)")
for (i in seq_len(nrow(configs))) {
  take <- late$rho_initial == configs$rho_initial[i] & late$kappa == configs$kappa[i]
  for (s in unique(late$seed[take])) {
    q <- late[take & late$seed == s, ]
    lines(q$iteration, q$elbo, col = adjustcolor(palette[i], 0.60), lwd = 1.5)
  }
}
dev.off()

png(file.path(out, "elbo-vs-waic.png"), width = 1300, height = 900, res = 160)
plot(raw$final_elbo, raw$WAIC_group_marginal, pch = 21, cex = 1.5,
     bg = cols[match(raw$rho_initial, sort(unique(raw$rho_initial)))],
     xlab = "ELBO final (mayor es mejor)", ylab = "WAIC marginal por finca (menor es mejor)",
     main = sprintf("ELBO identifica las mejores soluciones (correlacion = %.3f)",
                    cor(raw$final_elbo, raw$WAIC_group_marginal)))
text(raw$final_elbo, raw$WAIC_group_marginal, labels = raw$seed,
     pos = 3, cex = .65, col = "#444444")
legend("topright", legend = paste("rho1 =", sort(unique(raw$rho_initial))),
       pt.bg = cols, pch = 21, bty = "n")
dev.off()

schedule <- expand.grid(iteration = seq_len(iterations),
                        rho_initial = sort(unique(grid$rho_initial)),
                        kappa = sort(unique(grid$kappa)))
schedule$tau <- schedule$rho_initial^(-1 / schedule$kappa) - 1
schedule$step <- mapply(
  learning_rate_at, schedule$iteration, schedule$tau, schedule$kappa)
png(file.path(out, "step-size-schedules.png"), width = 1500, height = 900, res = 160)
plot(NA, xlim = c(1, iterations), ylim = range(schedule$step),
     xlab = "Iteracion", ylab = expression(rho[t]),
     main = "Secuencias de aprendizaje evaluadas")
for (i in seq_len(nrow(configs))) {
  q <- schedule[schedule$rho_initial == configs$rho_initial[i] &
                schedule$kappa == configs$kappa[i], ]
  lines(q$iteration, q$step, col = palette[i], lwd = 2)
}
legend("topright", sprintf("rho1=%.2f, k=%.2f", configs$rho_initial, configs$kappa),
       col = palette, lty = 1, lwd = 2, cex = .75, bty = "n")
dev.off()

print(summary)
