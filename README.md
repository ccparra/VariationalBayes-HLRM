# VariationalBayesHLRM

`VariationalBayesHLRM` is an R package for Bayesian Gaussian regression with
three model families:

- LRM: linear regression for one population;
- HLRM: hierarchical regression with group-specific coefficients and residual
  variances;
- CHLRM: a finite hierarchical mixture that clusters groups by shared
  regression profiles.

Each model supports Markov chain Monte Carlo (MCMC), coordinate-ascent
variational inference (CAVI), and stochastic variational inference (SVI). SVI
uses observation minibatches for LRM and complete-group minibatches for HLRM
and CHLRM.

This repository contains only reusable package code. Thesis sources, applied
data, paper experiments, and E1--E5 outputs live in the private companion
repository `ccparra/VariationalBayes-HLRM-research`.

## Documentation

- [Public API contract](docs/API.md): normative statistical and software
  behavior.
- [Package roadmap](docs/ROADMAP.md): current status, implementation order,
  release gates, and deferred work.

The package reference in `man/` documents individual functions. Vignettes show
complete fitting, simulation, prediction, assessment, and CHLRM relabeling
workflows.

## Installation

Install released dependencies and then the local source package:

```sh
Rscript -e 'install.packages(c("R6", "mvtnorm", "LaplacesDemon", "gtools", "posterior", "loo", "testthat", "clue"), repos="https://cloud.r-project.org")'
R CMD INSTALL .
```

## Minimal example

```r
library(VariationalBayesHLRM)

set.seed(10)
X <- cbind(intercept = 1, predictor = seq(-2, 2, length.out = 40))
data <- list(X = X, y = drop(X %*% c(1, 2)) + rnorm(40, sd = 0.5))

prior <- lrm_prior(
  coefficient_mean = c(0, 0),
  coefficient_covariance = 10 * diag(2),
  residual_degrees_of_freedom = 4,
  residual_variance_scale = 1
)

fit <- fit_model(
  model = lrm_model(2),
  data = data,
  prior = prior,
  method = "cavi",
  initialization = initialization_ols(),
  control = cavi_control(max_iterations = 2000),
  seed = 10,
  verbose = FALSE
)

summary(fit)
coef(fit)
predict(fit, X[1:3, , drop = FALSE], n_draws = 500, seed = 11)
posterior_predict(fit, X[1:3, , drop = FALSE], n_draws = 500, seed = 12)
```

## Public workflow

The main reusable operations are:

- `simulate_data()` for conditional, hierarchical, and prior-predictive DGPs;
- `fit_model()` for MCMC, CAVI, and SVI;
- `summarize_posterior()` and `posterior_draws()` for posterior output;
- `predict()` and `posterior_predict()` for mean and future-response
  prediction;
- `evaluate_in_sample_fit()` for descriptive training fit;
- `posterior_predictive_check()` for generative adequacy;
- `predictive_criteria()`, `predictive_score()`, `predictive_performance()`,
  `cross_validate_model()`, and `psis_loo()` for predictive assessment;
- CHLRM-specific accessors for components, responsibilities,
  classification uncertainty, and co-clustering.

Common result contracts are shared where the mathematics agrees. Quantities
that do not exist for a model are omitted rather than filled with placeholders.
MCMC retains posterior draws; CAVI/SVI retain variational factors. Their
log-likelihood, log-joint, and ELBO objectives remain explicitly distinct.

## Development validation

Run from the repository root:

```sh
Rscript dev/run-validation.R
R CMD build .
R CMD check --no-manual VariationalBayesHLRM_0.1.0.tar.gz
```

`tests/` contains fast deterministic regression tests for algebraic identities,
CAVI/SVI updates, MCMC orchestration, scaling, simulation, prediction,
likelihood estimands, and CHLRM relabeling. Generated validation output, source
archives, and check directories are ignored by Git.

## Statistical safeguards

- Covariance, precision, and scale matrices must be symmetric positive
  definite; invalid matrices are not silently repaired.
- Scaling is explicit, learned from training data only, and reversible for
  parameters and predictions.
- WAIC, PSIS-LOO, and log scores are compared only under the same dataset,
  predictive unit, target, conditioning convention, and likelihood.
- CHLRM component summaries use relabeling or label-invariant quantities.
- Resource estimates are conservative preflight guards, not operating-system
  memory limits.
