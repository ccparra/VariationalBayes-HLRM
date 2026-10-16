source('R/HLRM/hlrm.R')

# Simulated Data ------------------------------------
# ----- Parámetros del ejemplo -----
m <- 4
Betaj <- rbind(
  c(-7, 9,  35),
  c(  10, -20,  -4),
  c(  15,  3,  12),
  c( 1,  8,  -10)
)
nj <- c(150,100,300,100)
sigma2j <- rep(3^2, m)
p <- ncol(Betaj)


data <- sim_hlrm(
  m        = m,
  nj       = nj,
  Beta_j   = Betaj,
  p        = p ,
  dummies  = 0,
  sigma2_j = sigma2j,
  seed     = 123
)

X <- data$X
y <- data$y
g <- data$group

# ----- Chequeo rápido -----
sapply(X, dim)        # dimensiones de X_j
sapply(y, length)     # tamaños de y_j

# ----- Densidades de y por grupo (visualiza la separación) -----
yy <- unlist(y)
gg <- as.integer(unlist(g))
K  <- length(unique(gg))
cols <- c("#1b9e77","#d95f02","#7570b3","#e7298a")[seq_len(K)]

plot(density(yy[gg==1]), col = cols[1], lwd = 2,
     xlab = "y", main = "Densidades de y por grupo (HLRM)",
     xlim = range(yy))
for (k in 2:K) lines(density(yy[gg==k]), col = cols[k], lwd = 2)
legend("topright", legend = paste("Grupo", 1:K), col = cols, lwd = 2, bty = "n")

hypers <- hyperparams_hlrm(list(X = X, y = y),type='OLS')

#vi
run_settings_vi <- list(
  iterations = 20000
  #,epsilon =1e-10
  ,samples_gof = 1000
  #,stochastic = TRUE
  #,minibatch_m = 2
  #,tau = 10
  #,kappa= 1
)
vi <- vi_hlrm$new(y, X, hypers, run_settings_vi)
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
goodness_vi$dic$DIC
goodness_vi$r2
goodness_vi$ppval$ppp

#mcmc
hypers <- hyperparams_hlrm(list(X = X, y = y),type='OLS')
run_settings_mcmc <- list(
  samples = 1500,
  burn_in = 5000,
  thin = 10
)
mcmc <- mcmc_hlrm$new(y,X,hypers,run_settings_mcmc)
results_mcmc <- mcmc$fit(verbose=TRUE,seed = 1)
results_mcmc$timing$total
hat_values_mcmc <- mcmc$inferences()
hat_values_mcmc
y_hat_mcmc <- mcmc$predict(X)
goodness_mcmc <- mcmc$goodness_of_fit()
goodness_mcmc$waic$WAIC
goodness_mcmc$dic$DIC
goodness_mcmc$mse
goodness_mcmc$r2
goodness_mcmc$ppval$ppp


# Farms Data --------------------------------------------------------------
farms <- read.delim(project_path("data", "farm_trees.txt"))
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
  samples = 15000,
  burn_in = 5000,
  thin = 10
)
mcmc <- mcmc_hlrm$new(y, X,hypers, run_settings_mcmc)
results_mcmc <- mcmc$fit(verbose= TRUE,seed = 1)
results_mcmc$timing$total
hat_values_mcmc <- mcmc$inferences()
hat_values_mcmc
y_hat_vi <- mcmc$predict(X)
goodness_mcmc <- mcmc$goodness_of_fit()
goodness_mcmc$waic$WAIC
goodness_mcmc$dic$DIC
goodness_mcmc$mse
goodness_mcmc$r2
goodness_mcmc$ppval$ppp


#vi
run_settings_vi <- list(
  iterations = 200,
  epsilon = 0.00001,
  samples_gof = 1000,
  stochastic = TRUE,
  minibatch_m = 6,
  tau = 10 ,
  kappa= 1
)
vi <- vi_hlrm$new(y, X, hypers, run_settings_vi)
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

