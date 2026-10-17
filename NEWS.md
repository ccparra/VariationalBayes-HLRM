# VariationalBayesHLRM 0.2.0

- Established `docs/API.md` and `docs/ROADMAP.md` as the canonical English
  contract and execution plan.
- Reorganized internal modules by responsibility without changing statistical
  estimands or public result schemas.
- Froze the pre-1.0 public symbol manifest with automated tests.
- Added English workflow vignettes, continuous integration, repository
  metadata, and a software citation.
- Standardized source comments, messages, labels, and generated package output
  on English.

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
