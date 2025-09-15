source('R/chlrm/chlrm.R')
source('Master/utils.R')
source('R/CHLRM/utils.R')

library(paletteer)
library(knitr)
library(kableExtra)
# Simulated Data ------------------------------------
K <- 3; p <- 3; m <- 15; nj <- 20
Betak   <- rbind(c(-5,  8,  3),
                 c(10, -1, -2),
                 c(35, -8, -2))
sigma2k <- c(4^2, 3^2, 2^2)

data <- sim_chlrm(m, nj, Betak, sigma2k, p,omega = c(0.4,0.3,0.3), K = K, seed = 123)
X <- data$X
y <- data$y
g <- data$group
plot(density(unlist(y)))

hypers <- hyperparams_chlrm(list(X = X, y = y),K,"OLS")
plot(density(unlist(y)))


#mcmc
run_settings_mcmc <- list(
  samples = 50000
  ,burn_in = 10000
  ,thin = 10
  ,samples_gof = 1000
)
mcmc <- mcmc_chlrm$new(y, X,K,FALSE,hypers, run_settings_mcmc)
results_mcmc <- mcmc$fit(verbose= TRUE,seed = 1)
results_mcmc$timing$total
hat_values_mcmc <- mcmc$inferences()
hat_values_mcmc
y_hat_mcmc <- mcmc$predict(X)
goodness_mcmc <- mcmc$goodness_of_fit(seed = 1)
goodness_mcmc$waic$WAIC
goodness_mcmc$dic$DIC
goodness_mcmc$mse
goodness_mcmc$r2
goodness_mcmc$ppval$ppp
goodness_mcmc$ppval$repdata

#vi
run_settings_vi <- list(
  iterations = 50000
  ,epsilon = 1e-10
  ,samples_gof = 1000
  #,stochastic = TRUE
  #,minibatch_m = 6
  #,tau = 10 
  #,kappa= 1
)
vi <- vi_chlrm$new(y, X, K,standarize=FALSE, hypers, run_settings_vi)
results_vi <- vi$fit(verbose=TRUE)
results_vi$params
results_vi$timing$total
hat_values_vi <- vi$inferences()
hat_values_vi
y_hat_vi <- vi$predict(X)
goodness_vi <- vi$goodness_of_fit(seed=1)
goodness_vi$waic$WAIC
goodness_vi$dic$DIC
goodness_vi$mse
goodness_vi$r2
goodness_vi$ppval$ppp


#svi
run_settings_svi <- list(
  iterations = 50000
  ,epsilon = 1e-10
  ,samples_gof = 1000
  ,stochastic = TRUE
  ,minibatch_m = 12
  ,tau = (0.1)^(-1/0.7)-1
  ,kappa= 0.7
)
svi <- vi_chlrm$new(y, X,K, standarize=FALSE, hypers, run_settings_svi)
results_svi <- svi$fit(verbose=TRUE)
results_svi$params
results_svi$timing$total
hat_values_svi <- svi$inferences()
hat_values_svi
y_hat_svi <- svi$predict(X)
goodness_svi <- svi$goodness_of_fit(seed=1)
goodness_svi$waic$WAIC
goodness_svi$dic$DIC
goodness_svi$mse
goodness_svi$r2
goodness_svi$ppval$ppp


plot_y_predictive_densities(
  y          = y,
  y_hat_vi   = y_hat_vi,
  y_hat_svi  = y_hat_svi,
  y_hat_mcmc = y_hat_mcmc,
  col_obs = "blue",   lty_obs = 1, lwd_obs = 2,
  col_vi  = "black", lty_vi  = 5, lwd_vi  = 4,
  col_svi = "darkgreen",   lty_svi = 3, lwd_svi = 4,
  col_mcmc= "red", lty_mcmc= 1, lwd_mcmc=1,
  S = 12,
  adjust = 0.5,
  legend_pos = "topleft",
  legend_cex = 0.85,
  main = "",
  add_rug = FALSE,
)


# Supón que ya tienes X (lista de matrices) y y (lista o vector) como en tu ejemplo
# 3 grupos => personalizamos colores y formas
colores <- paletteer_c("grDevices::Viridis", 30)[0:15]  # azul, verde, rojo
formas  <- 0:15                       # círculo, triángulo, cuadrado

# Leyenda en 2 columnas, llenado por filas, fuera a la derecha
plot_3d_by_group(
  X, y,
  col_groups = colores,
  pch_groups = formas,
  legend_title = "Groups",
  legend_pos   = "bottomright",     # "right" + inset negativo la empuja fuera
  legend_inset = c(0.05,-0.09),
  legend_ncol  = 3,           # 2 columnas
  legend_byrow = TRUE,        # llenar por filas (1,2 en la primera fila; 3 en la segunda)
  main = "",
  xlab = expression(x[1]), ylab = expression(x[2]), zlab = "y"
)


# Mapeos: map[j] = índice actual que corresponde al grupo canónico j
map_mcmc <- c(2, 3, 1)   # MCMC ya está en el orden canónico
map_vi   <- c(2, 1, 3)   # en VI, el grupo canónico 2 está en la fila 3; el 3 en la fila 2
map_svi  <- c(3, 2, 1)   # igual que VI

labs <- c("Cluster 1","Cluster 2","Cluster 3")

# β1 por grupos 1:3, SIN alinear Y
plot_betak_ci_maps(
  zero_line = FALSE,main = "",
  hat_values_mcmc, map_mcmc,
  hat_values_vi,   map_vi,
  hat_values_svi,  map_svi,
  beta_index = 1,
  groups = 1:3, layout = c(1,3),
  group_labels = labs, xlab_common = "",
  sync_ylim = FALSE,
  legend_side = "right", legend_inset = -0.25,
  oma =  c(1.8, 1.5, 2, 12)
)

# β2 por grupos 1:3, SIN alinear Y
plot_betak_ci_maps(main = "",
  hat_values_mcmc, map_mcmc,
  hat_values_vi,   map_vi,
  hat_values_svi,  map_svi,
  beta_index = 2,
  groups = 1:3, layout = c(1,3),
  group_labels = labs, xlab_common = "",
  sync_ylim = FALSE,
  legend_side = "right", legend_inset = -0.25,
  oma =  c(1.8, 1.5, 2, 12)
)

# σ_k^2 por grupos 1:3, SIN alinear Y
plot_sigma2k_ci_maps(main = "",
  hat_values_mcmc, map_mcmc,
  hat_values_vi,   map_vi,
  hat_values_svi,  map_svi,
  groups = 1:3, layout = c(1,3),
  group_labels = labs, xlab_common = "",
  sync_ylim = FALSE,
  legend_side = "right", legend_inset = -0.25,
  oma =  c(1.8, 1.5, 2, 12)
)

# Tus objetos:
# goodness_mcmc, goodness_vi, goodness_svi
# Tiempos (en segundos, por ejemplo):
times <- c(MCMC = results_mcmc$timing$total, VI = results_vi$timing$total, SVI = results_svi$timing$total)

tabs <- make_latex_goodness_tables(
  goodness_mcmc,
  goodness_vi,
  goodness_svi,
  times = times,
  digits_fit = c(WAIC = 3, DIC = 3, MSE = 3, R2 = 3,Tiempo=3))
tabs$latex_fit   # código LaTeX tabla 1
tabs$latex_ppp

# Mostrar LaTeX (en consola) o guardados si diste save_tex_prefix
tabs$latex_fit
tabs$latex_ppp



plot(hat_values_mcmc$z$k_vals,hat_values_mcmc$z$prob_nc,type='h')
A <- hat_values_mcmc$z$incidence_matrix
clust_A <- mclust::Mclust(A)
indices <- order(clust_A$classification)
A_ord<- A[indices,indices]
par(mar = c(2.75,2.75,0.5,0.5), mgp = c(1.7,0.7,0))
corrplot::corrplot(corr = A, is.corr = FALSE, addgrid.col = NA, method = "color", tl.pos = "both", tl.cex = 0.3, tl.col = "black")
corrplot::corrplot(corr = A_ord, is.corr = FALSE, addgrid.col = NA, method = "color", tl.pos = "both", tl.cex = 0.3, tl.col = "black")

A_vi <- hat_values_vi$z$incidence_matrix
clust_A_vi <- mclust::Mclust(A_vi)
indices_vi <- order(clust_A_vi$classification)
A_vi_ord<- A_vi[indices_vi,indices_vi]
par(mar = c(2.75,2.75,0.5,0.5), mgp = c(1.7,0.7,0))
corrplot::corrplot(corr = A_vi, is.corr = FALSE, addgrid.col = NA, method = "color", tl.pos = "both", tl.cex = 0.3, tl.col = "black")
corrplot::corrplot(corr = A_vi_ord, is.corr = FALSE, addgrid.col = NA, method = "color", tl.pos = "both", tl.cex = 0.3, tl.col = "black")

A_svi <- hat_values_svi$z$incidence_matrix
clust_A_svi <- mclust::Mclust(A_svi)
indices_svi <- order(clust_A_svi$classification)
A_svi_ord<- A_svi[indices_svi,indices_svi]
par(mar = c(2.75,2.75,0.5,0.5), mgp = c(1.7,0.7,0))
corrplot::corrplot(corr = A_svi, is.corr = FALSE, addgrid.col = NA, method = "color", tl.pos = "both", tl.cex = 0.3, tl.col = "black")
corrplot::corrplot(corr = A_svi_ord, is.corr = FALSE, addgrid.col = NA, method = "color", tl.pos = "both", tl.cex = 0.3, tl.col = "black")

hypers_byK <- setNames(lapply(1:14, function(K) hyperparams_chlrm(list(X = X, y = y),K,"OLS")), paste0("K=", 1:14))

res <- run_vi_over_K(
  y = y, X = X,
  Ks = 1:14,
  hypers_byK = hypers_byK,   # << lista por K
  iterations = 20000,
  samples_gof = 1000,
  epsilon = 1e-10,
  standarize = FALSE,
  verbose = TRUE
)


# results_grid es la lista que te devolvió run_vi_over_K(...)
# (con elementos tipo "K=14", cada uno con $params$elbo)

df_elbo <- plot_elbo_last_by_K(
  results_list = res,
  # Ks = 1:15,                # opcional; por defecto usa todos los presentes
  type = "b", pch = 16, col = "#1f78b4", lwd = 2,
  main = "",
  xlab = "K (clusters)", ylab = "ELBO",
  grid = TRUE,
  mark_best = TRUE,
  best_point_col = "firebrick",
  best_label = TRUE,
  best_label_fmt = ""
)

# df_elbo tiene las parejas (K, ELBO_last) usadas en la gráfica
df_elbo


#mcmc
run_settings_mcmc <- list(
  samples = 50000
  ,burn_in = 10000
  ,thin = 10
  ,samples_gof = 1000
)
hypers_k <- hyperparams_chlrm(list(X = X, y = y),14,"OLS")
mcmc_k <- mcmc_chlrm$new(y, X,14,FALSE,hypers_k, run_settings_mcmc)
results_mcmc_k <- mcmc_k$fit(verbose= TRUE,seed = 1)
results_mcmc_k$timing$total
hat_values_mcmc_k <- mcmc_k$inferences()
hat_values_mcmc_k

plot(hat_values_mcmc_k$z$k_vals,
     hat_values_mcmc_k$z$prob_nc,
     type='h',
     lwd=5,
     xlab='K* (clusters)',
     ylab='Density')

plot(results_vi$timing$up_to_iter,results_vi$params$elbo,type='l',lwd=2)
lines(results_svi$timing$up_to_iter,results_svi$params$elbo, col='red',lwd=3,lty=13)

elbo_vi  <- results_vi$params$elbo
elbo_svi <- results_svi$params$elbo

plot(seq_along(elbo_vi), elbo_vi, type = "l", lwd = 2,
     xlab = "Iteration", ylab = "ELBO")
lines(seq_along(elbo_svi), elbo_svi, col = "red", lwd = 2, lty = 2)
points(seq_along(elbo_svi), elbo_svi, col = "red", pch = 16)
legend("bottomright", c("VI","SVI"), col = c("black","red"), lwd = 2, lty = c(1,2), pch = c(NA,16))


plot(results_mcmc$chains$logjoint, type = "l", col = "lightgray",
     xlab = "Iteration", ylab = "Log joint")

m <- mean(results_mcmc$chains$logjoint, na.rm = TRUE)
abline(h = m, col = "red", lty = 2, lwd = 2)

legend("bottomright",
       legend = c("log joint", sprintf("mean = %.3f", m)),
       col = c("lightgray", "red"),
       lty = c(1, 2), lwd = c(1, 2), bty = "n")

