source('R/LRM/lrm.R')
source('Master/utils.R')

iris <- iris
plot(iris)
X <- as.matrix(cbind(1,iris$Petal.Length))
y <- as.matrix(iris$Sepal.Length)
hypers <- hyperparams_lrm(list(X = X, y = y),type = 'NON_INFORMATIVE')
summary(lm(y~X))
#vi
run_settings_vi <- list(
  iterations = 20000,
  #epsilon = 0.00001,
  samples_gof = 1000
)
vi <- vi_lrm$new(y, X, hypers,run_settings_vi)
results_vi <- vi$fit(verbose=TRUE)
results_vi$params
results_vi$timing$total
hat_values_vi <- vi$inferences()
hat_values_vi
y_hat_vi <-  vi$predict(X)
goodness_vi <- vi$goodness_of_fit(seed=1)
goodness_vi$waic$WAIC
goodness_vi$dic$DIC
goodness_vi$mse
goodness_vi$r2
goodness_vi$ppval$ppp

#mcmc
run_settings_mcmc <- list(
  samples = 20000,
  burn_in = 5000,
  thin = 10
)
mcmc <- mcmc_lrm$new(y, X, hypers,run_settings_mcmc)
results_mcmc <- mcmc$fit(verbose=TRUE,seed = 1)
results_mcmc$timing$total
hat_values_mcmc <- mcmc$inferences()
hat_values_mcmc
y_hat_mcmc <-  mcmc$predict(X)
goodness_mcmc <- mcmc$goodness_of_fit()
goodness_mcmc$waic$WAIC
goodness_mcmc$dic$DIC
goodness_mcmc$mse
goodness_mcmc$r2
goodness_mcmc$ppval$ppp

plot(X[,2], y,
     ylab = "Sepal Length",
     xlab = "Petal Length",
     pch  = 4, col = "darkgray")

# recta MCMC (roja, sólida)
abline(coef = hat_values_mcmc$Beta$mean, col = "red",   lty = 1, lwd = 2)

# recta VI (negra, discontinua)
abline(coef = hat_values_vi$Beta$mean,   col = "black", lty = 2, lwd = 2)

# leyenda
legend("topleft",
       legend = c("MCMC", "VI"),
       col    = c("red", "black"),
       lty    = c(1, 2),
       lwd    = c(2, 2),
       bty    = "n",      # sin recuadro
       inset  = 0.02,
       seg.len = 2)

# Asignas las covarianzas completas p×p
hat_values_vi$Beta$cov    <- results_vi$params$SigmaBeta
hat_values_mcmc$Beta$cov  <- cov(results_mcmc$chains$Beta)

# (opcional) si tienes las covarianzas completas p×p:
# hat_values_mcmc$Beta$cov <- cov(results_mcmc$chains$Beta)   # empírica MCMC
# hat_values_vi$Beta$cov    <- results_vi$params$SigmaBeta    # o la que estimes en VI

# Elipses 50% y 95%; leyenda compacta
plot_betas_joint(
  results_mcmc, hat_values_mcmc, hat_values_vi,
  j1 = 1, j2 = 2,
  probs_ellipses = c(0.5, 0.95),
  ci_label_level = 95,
  legend_pos = "topright",
  legend_cex = 0.4, legend_pt_cex = 0.75,
  legend_seg_len = 1.1, legend_inset = 0.02,
  legend_ncol = 1, legend_box_lwd = 0,
)

# Más pequeña y en dos columnas
plot_betas_joint(
  results_mcmc, hat_values_mcmc, hat_values_vi,
  j1 = 1, j2 = 2,
  probs_ellipses = c(0.5, 0.9),
  legend_pos = "topright",
  legend_cex = 0.55, legend_pt_cex = 0.6,
  legend_ncol = 2, legend_bg = "white", legend_box_lwd = 0.4, legend_box_col = "grey70"
)

plot_sigma2_density_ig(
  results_mcmc, hat_values_mcmc, hat_values_vi,
  adjust = 2,
  main_title = "",
  legend_pos = "topright",
  legend_cex = 0.55, legend_pt_cex = 0.75
)

scatterplot3d::scatterplot3d(results_mcmc$chains$Beta[,1],results_mcmc$chains$Beta[,2],results_mcmc$chains$sigma2)

plot_vi_vs_mcmc_time(
  results_vi, results_mcmc,
  main_title = "",
  xlab='Time (sec)',
  ylab = expression(E[q(bold(Theta))] * "[" * log * phantom(.) * p(y, X, bold(Theta))-log * phantom(.) * q(bold(Theta)) * "]"),
  legend_pos = "right",
  legend_cex = 0.85,
  lty_vi = 1,
  lwd_vi = 1,
  pch_vi = NA
)


cat( latex_fit_table(goodness_mcmc, goodness_vi,results_mcmc,results_vi,
                     digits = 3,
                     caption = "Linear Regression: fit comparison",
                     label   = "tab:fit_lr",
                     use_booktabs = TRUE) )

cat(latex_param_table_both(hat_values_mcmc, hat_values_vi,
                           digits = 3,
                           caption = "Posterior summaries (MCMC vs VI)",
                           label   = "tab:post_mcmc_vi",
                           use_booktabs = TRUE))



