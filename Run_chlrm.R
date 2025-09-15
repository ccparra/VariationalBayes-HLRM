source('R/CHLRM/chlrm.R')
source('R/CHLRM/utils.R')

# Simulated Data ------------------------------------
K <- 3; p <- 3; m <- 300; nj <- 200
Betak   <- rbind(c(5,  4,  0.5),
                 c(-3, 1.5, 30),
                 c(12, -2, -1))
sigma2k <- c(12^2, 3^2, 20^2)

data <- sim_chlrm(m, nj, Betak, sigma2k, p,omega = c(0.3,0.1,0.6), K = K, seed = 123)
X <- data$X
y <- data$y
g <- data$group
hypers <- hyperparams_chlrm(list(X = X, y = y),K,"NON_INFORMATIVE")


#vi
run_settings_vi <- list(
  iterations = 20000
  ,epsilon = 10^-5
  ,samples_gof = 1000
  #,stochastic = TRUE
  #,minibatch_m = 20
  #,tau = 10
  #,kappa = 1
)
vi <- vi_chlrm$new(y,X,K,standardize= FALSE,hypers,run_settings_vi)
results_vi <- vi$fit()
results_vi$params
results_vi$timing$total
hat_values_vi <- vi$inferences()
hat_values_vi
y_hat_vi <- vi$predict(X,as_vector=TRUE)
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
  thin = 10,
  samples_gof = 1000
)
mcmc <- mcmc_chlrm$new(y,X,K,standardize= FALSE,hypers,run_settings_mcmc)
results_mcmc <- mcmc$fit(verbose=TRUE,seed = 123)
results_mcmc$timing$total
hat_values_mcmc <- mcmc$inferences()
hat_values_mcmc
y_hat_mcmc <- mcmc$predict(X,as_vector = TRUE)
goodness_mcmc <- mcmc$goodness_of_fit()
goodness_mcmc$waic$WAIC
goodness_mcmc$dic$DIC
goodness_mcmc$mse
goodness_mcmc$r2
goodness_mcmc$ppval$ppp




# Farms Data --------------------------------------------------------------
farms <- read.delim("Data/farm_trees.txt")
farms_l <- df_to_chlrm_lists(
  data  = farms,
  y     = "size",
  x     = c("N"),
  group = "farm"
)
X <- farms_l$X
y <- farms_l$y
K <- 6
hypers <- hyperparams_chlrm(list(X = X, y = y),K,"NON_INFORMATIVE")
#vi
run_settings_vi <- list(
  iterations = 10000
  ,epsilon = 0.00001
  ,samples_gof = 1000
  #,stochastic = TRUE
  #,minibatch_m = 8
  #,tau = 10
  #,kappa = 1
)
vi <- vi_chlrm$new(y,X,K,standardize= FALSE,hypers,run_settings_vi)
results_vi <- vi$fit(verbose = TRUE)
results_vi$params
results_vi$timing$total
hat_values_vi <- vi$inferences()
hat_values_vi
y_hat_vi <- vi$predict(X,as_vector=TRUE)
goodness_vi <- vi$goodness_of_fit(seed=1)
goodness_vi$waic$WAIC
goodness_vi$dic$DIC
goodness_vi$mse
goodness_vi$r2
goodness_vi$ppval$ppp

#mcmc
hypers <- hyperparams_chlrm(list(X = X, y = y),K,"NON_INFORMATIVE")
run_settings_mcmc <- list(
  samples = 10000,
  burn_in = 10000,
  thin = 10,
  samples_gof = 1000
)
mcmc <- mcmc_chlrm$new(y,X,K,standardize= FALSE,hypers,run_settings_mcmc)
results_mcmc <- mcmc$fit(verbose=TRUE,seed = 123)
results_mcmc$timing$total
hat_values_mcmc <- mcmc$inferences()
hat_values_mcmc
y_hat_mcmc <- mcmc$predict(X,as_vector = TRUE)
goodness_mcmc <- mcmc$goodness_of_fit()
goodness_mcmc$waic$WAIC
goodness_mcmc$dic$DIC
goodness_mcmc$mse
goodness_mcmc$r2
goodness_mcmc$ppval$ppp



