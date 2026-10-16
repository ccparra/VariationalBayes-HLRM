source('R/CHLRM/chlrm.R')
source("thesis/v0/code/analysis/utils.R")
source('R/CHLRM/utils.R')


library(mclust)
# Simulated Data ------------------------------------
farms <- read.delim(project_path("data", "farm_trees.txt"))
farms_l <- df_to_chlrm_lists(
  data  = farms,
  y     = "size",
  x     = c("N"),
  group = "farm"
)
X <- farms_l$X
y <- farms_l$y
K <- 23
hypers <- hyperparams_chlrm(list(X = X, y = y),K,"OLS")

#mcmc
run_settings_mcmc <- list(
  samples = 20000
  ,burn_in = 10000
  ,thin = 10
  ,samples_gof = 100
)
K <- 23
hypers_k23 <- hyperparams_chlrm(list(X = X, y = y),K,"OLS")
mcmc_k23 <- mcmc_chlrm$new(y, X, K, hypers_k23, run_settings_mcmc)
results_mcmc_k23 <- mcmc_k23$fit(verbose= TRUE,seed = 1)
hat_values_mcmc_k23 <- mcmc_k23$inferences()
plot(hat_values_mcmc_k23$z$k_vals,hat_values_mcmc_k23$z$prob_nc,type='h',xlab='K* (clusters)',ylab='Density',lwd=5)

K <- 5
#mcmc
run_settings_mcmc <- list(
  samples = 120000
  ,burn_in = 30000
  ,thin = 20
  ,samples_gof = 1000
)
hypers <- hyperparams_chlrm(list(X = X, y = y),K,"OLS")
mcmc <- mcmc_chlrm$new(y, X, K, hypers, run_settings_mcmc)
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


hypers_byK <- setNames(lapply(1:23, function(K) hyperparams_chlrm(list(X = X, y = y),K,"OLS")), paste0("K=", 1:23))

res <- run_vi_over_K(
  y = y, X = X,
  Ks = 1:23,
  hypers_byK = hypers_byK,   # << lista por K
  iterations = 20000,
  samples_gof = 1000,
  epsilon = 1e-10,
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
  best_label_fmt = " "
)

K <- 5
  
#vi
run_settings_vi <- list(
    iterations = 20000
    ,epsilon = 1e-10
    ,samples_gof = 1000
    #,stochastic = TRUE
    #,minibatch_m = 6
    #,tau = 10 
    #,kappa= 1
  )
hypers <- hyperparams_chlrm(list(X = X, y = y),K,"OLS")
vi <- vi_chlrm$new(y, X, K, hypers, run_settings_vi)
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
  iterations = 20000
  ,epsilon = 1e-4
  ,samples_gof = 1000
  ,stochastic = TRUE
  ,minibatch_m = 18
  ,tau = (0.05)^(-1/0.7)-1
  ,kappa= 0.7
)
svi <- vi_chlrm$new(y, X, K, hypers, run_settings_svi)
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


x_vi  <- results_vi$timing$up_to_iter
x_svi <- results_svi$timing$up_to_iter
elbo_vi  <- results_vi$params$elbo
elbo_svi <- results_svi$params$elbo

plot(x_vi, elbo_vi, type = "l", lwd = 2, log = "x",
     xlab = "Time (s, log scale)", ylab = "ELBO",
     xlim = range(c(x_vi, x_svi)))
lines(x_svi, elbo_svi, col = "red", lwd = 2, lty = 2)
points(x_svi, elbo_svi, col = "red", pch = 16)
legend("bottomright", c("VI","SVI"), col = c("black","red"), lwd = 2, lty = c(1,2), pch = c(NA,16))


elbo_vi  <- results_vi$params$elbo
elbo_svi <- results_svi$params$elbo

plot(seq_along(elbo_svi)[1:50], elbo_svi[1:50], type = "l", lwd = 2,
     xlab = "Iteration", ylab = "ELBO")
lines(seq_along(elbo_vi), elbo_vi, col = "red", lwd = 2, lty = 2)
points(seq_along(elbo_svi), elbo_svi, col = "red", pch = 16)
legend("bottomright", c("VI","SVI"), col = c("black","red"), lwd = 2, lty = c(1,2), pch = c(NA,16))


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
  adjust = 0.8,
  legend_pos = "topright",
  legend_cex = 0.85,
  main = "",
  add_rug = FALSE,
)

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
tabs$latex_fit  
tabs$latex_ppp
