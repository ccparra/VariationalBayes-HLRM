# VariationalBayesHLRM 0.1.0

- Unified LRM, HLRM, and CHLRM under a typed public API.
- Added MCMC, CAVI, and SVI engines for all three model families.
- Added posterior draws, prediction, posterior predictive checks, predictive
  scoring, WAIC, PSIS-LOO for MCMC, and cross-validation.
- Added explicit prior, initialization, convergence, scaling, simulation, and
  resource-control contracts.
- Added multichain MCMC, multistart CAVI/SVI, checkpointing, and CHLRM label
  alignment.
- Separated the reusable package from thesis and paper-specific research code.
