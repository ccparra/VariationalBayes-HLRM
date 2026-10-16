source('R/HLRM/hlrm.R')
source("thesis/v0/code/analysis/utils.R")

farms <- read.delim(project_path("data", "farm_trees.txt"))
runs <- run_vi_hlrm_grid(
  y, X, hypers,
  p1_vals    = c(0.05,0.1,0.15,0.2),
  kappa_vals = c(0.6,0.7,0.8,0.9),
  batch_vals = c(1, 6, 12, 18, 24),
  iterations = 20000,
  samples_gof = 1000,
  seed = 1,
  label_fun = function(p1, kappa, batch, stochastic)
    if (!stochastic) "batch" else sprintf("minibatch: %s", batch)
)

# únicos disponibles (tal cual vienen)
(unique_m      <- sort(unique(vapply(runs, function(r) r$settings$minibatch_m, numeric(1)))))
(unique_kappa  <- sort(unique(vapply(runs, function(r) r$settings$kappa,       numeric(1)))))
(unique_tau    <- sort(unique(vapply(runs, function(r) r$settings$tau,         numeric(1)))))

runs$`minibatch: 1`$settings$tau


# 1) Comparar m con κ=0.7 y p1=0.05, con límites Y y leyenda arriba en 2 columnas
plot_elbo_by_minibatch(
  runs,
  kappa = 0.7, p1 = 0.2,
  minibatch_vals = unique_m,
  k = 300, step = 300, start_at = 0,
  ylim = c(-20000, 500),                # <- tu límite manual
  legend_side = "top", legend_outside = TRUE, legend_inset = -0.08,
  legend_cex = 0.85, legend_ncol = 2, legend_nrow = 12,
  legend_labels = paste("minibatch", unique_m),
  oma = c(1.5,1.5,2.6,1)
)

# 2) Comparar (κ, p1) con m fijo = 24, leyenda a la derecha fuera, en 3 filas
plot_elbo_by_kappa_tau(
  runs,
  minibatch_m = 24,
  kappas  = c(0.6, 0.7, 0.8, 0.9),
  p1_vals = c(0.05, 0.05, 0.05, 0.05),
  ylim = c(-1000, 0), 
  k = 300, step = 300, start_at = 300,
  legend_side = "right", legend_outside = TRUE, legend_inset = -0.12,
  legend_cex = 0.9, legend_ncol = 2, legend_nrow = 4,
  legend_labels = paste0("kappa=", c(0.6,0.7,0.8,0.9)),
  oma = c(1.5,1.5,2.6,3)
)
