source('R/LRM/lrm.R')

# Simulated Data --------------------------------------------
beta   <- c(2, 0.5)
beta   <- c(25, 10,-30)
sigma2 <- 100^2
x_range <-matrix(rbind(c(0,100),c(-5,5)),nrow = 2, ncol =2, byrow = TRUE)
data <- sim_lrm(n = 200, beta = beta, sigma2 = sigma2,
                x_range = x_range,seed=123,jitter_frac = 0.0001)
X <- data$X
y <- data$y
scatterplot3d::scatterplot3d(X[,2],X[,3],y)
plot(data$X[,2], data$y, pch = 16, cex = 0.6)
plot(data$X[,3], data$y, pch = 16, cex = 0.6)
hypers <- hyperparams_lrm(data, type = "NON_INFORMATIVE")

#vi
run_settings_vi <- list(
  iterations = 200000
  ,epsilon = 0.00001
  ,samples_gof = 1000
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
  samples = 2000,
  burn_in = 10000,
  thin = 10
)
mcmc <- mcmc_lrm$new(y, X, hypers,run_settings_mcmc)
results_mcmc <- mcmc$fit(verbose=TRUE,seed = 1)
hat_values_mcmc <- mcmc$inferences()
hat_values_mcmc
y_hat_mcmc <-  mcmc$predict(X)
goodness_mcmc <- mcmc$goodness_of_fit()
goodness_mcmc$waic$WAIC
goodness_mcmc$dic$DIC
goodness_mcmc$mse
goodness_mcmc$r2
goodness_mcmc$ppval$ppp


plot(y~X[,2], col= 'pink')
abline((lm(y~X[,2])),col='red')
points(y_hat_vi~X[,2],col='darkgreen')
abline(reg = vi$hat_values$Beta$mean, col='green')
points(y_hat_mcmc~X[,2],col='blue')
abline(reg = mcmc$hat_values$Beta$mean,col='lightblue')

# Farms Data --------------------------------------------------------------
farms <- read.delim(project_path("data", "farm_trees.txt"))
X <- as.matrix(cbind(1,farms$N))
y <- as.matrix(farms$size)
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
