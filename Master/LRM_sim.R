source('R/LRM/lrm.R')
source('Master/utils.R')

sim_specs <- list(
  p3_n50   = list(beta=c(25,10,-30), sigma2=100^2, n=50, x_range=rbind(c(0,100),c(-5,5))),
  p3_n100   = list(beta=c(25,10,-30), sigma2=100^2, n=100, x_range=rbind(c(0,100),c(-5,5))),
  p3_n1000   = list(beta=c(25,10,-30), sigma2=100^2, n=1000, x_range=rbind(c(0,100),c(-5,5))),
  p2_n50 = list(beta=c( 2, -12), sigma2= 4^2, n=50, x_range=matrix(c(0,100),1,2)),
  p2_n100 = list(beta=c( 2, -12), sigma2= 4^2, n=100, x_range=matrix(c(0,100),1,2)),
  p2_n1000 = list(beta=c( 2, -12), sigma2= 4^2, n=1000, x_range=matrix(c(0,100),1,2))
)

# Ajustes VI y MCMC (puedes dejarlos por default si quieres)
run_settings_vi <- list(iterations = 2000, epsilon = 1e-10, samples_gof = 1000)
run_settings_mcmc <- list(samples = 20000, burn_in = 10000, thin = 10)

# Ejecutar:
all_runs_lrm <- run_simulations_lrm(
  sim_specs = sim_specs,
  vi_settings = run_settings_vi,
  mcmc_settings = run_settings_mcmc,
  hyper_type = "OLS",
  verbose = TRUE,
  save_rds_path = "Master/sim_lrm_runs.rds"
)

tex <- latex_methods_overview(
  all_runs_lrm,
  caption = "Coverage, zero-in-CI and runtime by method and scenario",
  label   = "tab:method_quality"
)
cat(tex)
