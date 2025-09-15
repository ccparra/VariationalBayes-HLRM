source('R/HLRM/hlrm.R')
source('Master/utils.R')

library(paletteer)

farms <- read.delim("Data/farm_trees.txt")
farms_l <- df_to_hlrm_lists(
  data  = farms,
  y     = "size",
  x     = c("N"),
  group = "farm"
)
X <- farms_l$X
y <- farms_l$y
hypers <- hyperparams_hlrm(farms_l,type='OLS')

#mcmc
run_settings_mcmc <- list(
  samples = 20000
  ,burn_in = 5000
  ,thin = 10
  ,samples_gof = 1000
)
mcmc <- mcmc_hlrm$new(y, X,hypers, run_settings_mcmc)
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
vi <- vi_hlrm$new(y, X, standarize=FALSE, hypers, run_settings_vi)
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
  ,minibatch_m = 20
  ,tau = (0.05)^(-1/0.7)-1
  ,kappa= 0.7
)
svi <- vi_hlrm$new(y, X, standarize=FALSE, hypers, run_settings_svi)
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
  col_obs = "orange",   lty_obs = 1, lwd_obs = 2,
  col_vi  = "blue", lty_vi  = 5, lwd_vi  = 2,
  col_svi = "black",   lty_svi = 2, lwd_svi = 2,
  col_mcmc= "red", lty_mcmc= 1, lwd_mcmc=2,
  S = 18,
  adjust = 0.5,
  legend_pos = "topright",
  legend_cex = 0.85,
  main = "",
  add_rug = FALSE,
)


# Ejemplo con tus objetos
cat(
  latex_ic_counts_table(hat_mcmc = hat_values_mcmc,
                        hat_vi   = hat_values_vi,
                        hat_svi  = hat_values_svi,
                        caption = "Conteo de IC que contienen 0 y la media, por método",
                        label   = "tab:ic_metodos",
                        show_totals = TRUE)  # usa FALSE si quieres solo el número
)

times <- c(MCMC = results_mcmc$timing$total, VI = results_vi$timing$total, SVI = results_svi$timing$total)

tabs <- make_latex_goodness_tables(
  goodness_mcmc,
  goodness_vi,
  goodness_svi,
  times = times,
  digits_fit = c(WAIC = 3, DIC = 3, MSE = 3, R2 = 3,Tiempo=3))
tabs$latex_fit   # código LaTeX tabla 1
tabs$latex_ppp

plot(results_vi$timing$up_to_iter,results_vi$params$elbo)
lines(results_svi$timing$up_to_iter,results_svi$params$elbo,col='red')


m <- nrow(hat_values_mcmc$Betaj$mean)

plot_hlrm_betas_by_group(
  results_mcmc    = results_mcmc,
  hat_values_mcmc = hat_values_mcmc,
  hat_values_vi   = hat_values_vi,
  hat_values_svi  = hat_values_svi,
  groups = 1:24, layout = c(4,6), S = 4000,
  main = "",
  group_label = "Farm",
  mcmc_col = "red", vi_col = "blue", svi_col = "black",
  pch_mcmc = 20, pch_vi = 20, pch_svi = 25, cex_mean = 1.9,
  lty_mcmc_ci = "13", lty_vi_ci = "42", lty_svi_ci = "44",
  ci_nudge_frac = 0.01, legend_cex = 0.7,
  legend_outside = TRUE, legend_side = "top", legend_inset = c(0, -0.045),
  oma = c(1.5, 1.5, 2.8, 0.6)
)

plot_hlrm_sigma2j_by_group(
  hat_values_mcmc = hat_values_mcmc,
  hat_values_vi   = hat_values_vi,
  hat_values_svi  = hat_values_svi,
  groups = 1:24, layout = c(4,6),
  group_label = "Farm",
  main = "",
  col_mcmc = "red", col_vi = "blue", col_svi = "black",
  pch_mcmc = 20, pch_vi = 20, pch_svi = 25, cex_mean = 0.8,
  legend_cex = 0.7,
  legend_outside = TRUE, legend_side = "top", legend_inset = c(0, -0.045),
  oma = c(1.5, 1.5, 2.8, 0.4)
)



plot_hlrm_lines_by_group(
  data = farms, x = N, y = size, group = farm,
  hat_values_mcmc = hat_values_mcmc,
  hat_values_vi   = hat_values_vi,
  hat_values_svi  = hat_values_svi,
  layout = c(4, 6),
  main = "",
  xlab = "Nitrogen (N)",
  ylab = "Size",
  group_label = "Farm",
  legend_outside = TRUE, legend_side = "top", legend_cex = 0.85,legend_inset=c(0, -0.08),
  lty_mcmc = "13", col_mcmc = "red",   lwd_mcmc = 1.5,
  lty_vi   = "solid", col_vi   = "blue", lwd_vi   = 2,
  lty_svi  = "1343", col_svi  = "black", lwd_svi = 3,
  sync_ylim = TRUE, y_pad_frac = 0,
)


# tu scatter:
plot(farms[, c(1, 2)],
     pch = c(0:24)[farms$farm],
     col = paletteer_d("colorBlindness::SteppedSequential5Steps")[farms$farm],
     xlab = "Nitrogen", ylab = "Size")

# recta gris global (OLS):
abline(lm(size ~ N, data = farms), col = "gray", lwd = 1.2, lty = 1)

# etiquetas por grupo (en centro mediano, con colocación inteligente si hay plotrix):
add_group_labels(
  data = farms, x = N, y = size, group = farm,
  fun = median,                 # prueba mean si prefieres
  cex = 1, offset = 0.6,
  col_map = paletteer_d("colorBlindness::SteppedSequential5Steps"),
  use_plotrix = TRUE,           # si no está plotrix, hace fallback a text()
  draw_boxes = FALSE            # pon TRUE si quieres cajitas blancas
)




