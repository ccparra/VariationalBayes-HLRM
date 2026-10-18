# Public API contract

## Authority and purpose

This document is the normative public contract for `VariationalBayesHLRM`.
Code, tests, generated reference documentation, and external consumers must
agree with it. [ROADMAP.md](ROADMAP.md) is the sole source for incomplete and
deferred work.

The package provides reusable simulation, inference, posterior analysis,
prediction, and assessment for:

- LRM: Gaussian linear regression for one population;
- HLRM: hierarchical Gaussian regression with group-specific coefficients and
  residual variances;
- CHLRM: a finite hierarchical mixture in which groups share component-level
  regression profiles.

The exported namespace is the public surface. R6 backends, Gibbs steps,
variational updates, sufficient statistics, and native parameter names are
implementation details.

## Contract principles

1. The API is organized around statistical questions rather than algorithms
   or source files.
2. LRM, HLRM, and CHLRM use the same public name when an operation has the
   same meaning.
3. Model-specific quantities are omitted when they do not exist; they are not
   represented by `NULL`, zero, or artificial estimates.
4. Every predictive result records its target, predictive unit, and
   conditioning or marginalization convention.
5. Evaluation operations do not mutate fitted objects.
6. Potentially large draws and pointwise matrices are retained only when
   requested and are guarded by explicit memory budgets.
7. A supplied seed makes an operation reproducible and restores the caller's
   RNG state.
8. The package does not know paper scenarios, experimental grids, paths,
   figures, or manuscript-specific sample sizes.
9. English is used for public names, messages, labels, metadata, and output.

## Models and inference engines

| Model | MCMC | CAVI | SVI | SVI sampling unit |
|---|---:|---:|---:|---|
| LRM | yes | yes | yes | observation |
| HLRM | yes | yes | yes | complete group |
| CHLRM | yes | yes | yes | complete group |

`method` accepts `"mcmc"`, `"cavi"`, or `"svi"`. CAVI is the canonical name
of deterministic coordinate-ascent variational inference. VI denotes the
statistical family containing CAVI and SVI and is not an accepted ambiguous
method value.

SVI samples complete groups in HLRM and CHLRM because the local state of a
group must be updated jointly. LRM samples observations. Full-batch SVI with a
unit learning rate reproduces the corresponding CAVI update.

## Typed construction and fitting

The public construction layer separates six responsibilities:

| Specification | Responsibility |
|---|---|
| model | likelihood and model dimensions |
| prior | probability distributions for model parameters |
| initialization | construction of an initial algorithmic state |
| inference control | execution limits, diagnostics, and stopping |
| data scaling | explicit preprocessing learned from training data |
| simulation | DGP, covariate design, sample size, and split |

Algorithm controls such as learning rate, minibatch size, iteration limit,
chain count, and tolerances are not model hyperparameters. The word
hyperparameter is reserved for the prior or probabilistic hierarchy.

The main fitting interface is:

```r
fit_model(
  model,
  data,
  prior,
  method = c("mcmc", "cavi", "svi"),
  initialization,
  control,
  scaling = NULL,
  seed = NULL,
  verbose = TRUE
)
```

The fitted object is an S3 `vb_model_fit`. It retains validated
specifications, the training-only scaling map, seeds, provenance, diagnostics,
timing, and the internal result needed by public methods. A backend represents
one fit and is never silently refitted. MCMC continuation uses
`resume_mcmc()`.

### Model and prior constructors

- `lrm_model()`, `hlrm_model()`, and `chlrm_model()` define a model and its
  likelihood.
- `lrm_prior()`, `hlrm_prior()`, and `chlrm_prior()` define explicit priors.
- `prior_moments()` reports whether requested moments exist and their values.

Public prior names describe distribution and scale:

| Prior block | LRM | HLRM | CHLRM |
|---|---:|---:|---:|
| `coefficient_mean`, `coefficient_covariance` | yes | yes | yes |
| `coefficient_dispersion_df`, `coefficient_dispersion_scale` | — | yes | yes |
| `residual_degrees_of_freedom` | yes | yes | yes |
| `residual_variance_scale` | yes | — | — |
| `residual_variance_scale_shape`, `residual_variance_scale_rate` | — | yes | yes |
| `mixture_concentration` | — | — | yes |

The corresponding hierarchy is

$$
\begin{aligned}
\boldsymbol\beta &\sim N(m_0,V_0),\\
\boldsymbol\Sigma &\sim IW(n_0,S_0),\\
\sigma_j^2\mid\xi^2 &\sim IG(\nu/2,\nu\xi^2/2),\\
\xi^2 &\sim Gamma(a,b),\\
\boldsymbol\omega &\sim Dirichlet(\boldsymbol\alpha).
\end{aligned}
$$

Inverse-gamma uses shape--scale, Gamma uses shape--rate, and inverse-Wishart
uses a scale matrix. Covariance, precision, and scale are not interchangeable.
An inverse-Wishart prior is proper when its degrees of freedom exceed
`p - 1`; existence of its mean requires the stronger distribution-specific
condition.

The package does not provide a prior called non-informative and does not use
an empirical-Bayes prior by default. OLS is an initialization strategy, not a
prior. Supplying truth to an inference prior is leakage except in an explicitly
defined algebraic test. In simulation-based calibration, the DGP and fit share
the prior distribution but the fit never receives the realized truth.

### Initialization

`initialization_from_prior()`, `initialization_ols()`,
`initialization_kmeans()`, and `initialization_values()` create validated
initialization specifications. K-means initialization is CHLRM-specific.
Incompatible combinations fail explicitly. Initialization never changes the
prior. Multiple chains or starts receive independent reproducible states.

`variational_initialization_set()` materializes a versioned, serializable set
of CAVI/SVI states from a model, dataset, prior, base initialization strategy,
optional scaling specification, and required seed. Passing the same set as the
`initialization` argument of separate CAVI and SVI calls guarantees identical
starting variational factors. The `seed` supplied to each `fit_model()` call
then controls optimization independently. The number of states must equal
`control$starts`; model dimensions and the learned scaling map must match.
Native state field names are implementation details and are not part of this
public contract.

### Controls and convergence

- `mcmc_control()` uses `mcmc_stopping()` and monitors rank-normalized R-hat,
  bulk ESS, tail ESS, and MCSE for multiple chains.
- `cavi_control()` uses `cavi_stopping()` and monitors deterministic full-data
  ELBO and variational-parameter change.
- `svi_control()` uses `svi_stopping()` and monitors periodically comparable
  full-data ELBO evaluations plus natural-parameter change.
- `power_schedule()` defines the SVI learning-rate schedule using descriptive
  `offset` and `exponent` arguments.

SVI never compares objectives computed from different minibatches as a
stopping rule. The learning-rate schedule controls step size; the stopper
assesses stabilization. Every engine also has a hard work limit.

The common convergence contract is `convergence-1`:

```r
list(
  status = "converged",
  reason = "objective_and_parameters_stable",
  work_completed = 4200L,
  checks_completed = 64L,
  metrics = list(...),
  thresholds = list(...),
  history = data.frame(...)
)
```

Only applicable diagnostics are present. CAVI/SVI do not receive fictitious
R-hat or ESS, and MCMC does not receive a variational convergence label. A
single MCMC chain reaches its declared budget but cannot establish
multi-chain convergence. With `starts > 1`, CAVI/SVI retain every attempt and
select the largest finite terminal ELBO among attempts satisfying the declared
convergence rule. If none qualifies, the status is `no_valid_start`.

## Result schemas

| Schema | Class | Producer |
|---|---|---|
| `fit-result-2` | `vb_fit_result` | `fit_model()` |
| `posterior-summary-1` | `vb_posterior_summary` | `summarize_posterior()` |
| `posterior-draws-1` | `vb_posterior_draws` | `posterior_draws()` |
| `in-sample-fit-1` | `vb_in_sample_fit` | `evaluate_in_sample_fit()` |
| `posterior-predictive-check-2` | `vb_posterior_predictive_check` | `posterior_predictive_check()` |
| `predictive-criteria-2` | `vb_predictive_criteria` | `predictive_criteria()` |
| `posterior-prediction-2` | `vb_prediction_result` | `predict()`, `posterior_predict()` |
| `log-likelihood-1` | `vb_log_likelihood` | `log_likelihood()` |
| `simulation-result-2` | `vb_simulation_result` | `simulate_data()` |
| `fit-resource-estimate-1` | `vb_fit_resource_estimate` | `estimate_fit_resources()` |

Changing the structure or meaning of a schema requires a version increment;
renaming only the schema label is not a migration.

`fit-result-2` contains `parameters`, `objective`, `diagnostics`, `timing`, and
`metadata`. CHLRM additionally contains `latent_variables`. MCMC stores draws;
CAVI/SVI store variational-factor parameters. MCMC log likelihood/log joint and
variational ELBO are separate algorithmic objectives and are not comparable as
one metric.

## Posterior summaries and draws

`summarize_posterior(fit)` returns `posterior-summary-1`. Applicable parameter
summaries contain `mean`, `sd`, `lower`, `upper`, `level`, `distribution`,
`mean_exists`, and `contains_zero`. CHLRM also reports classification and
relabeling information where requested by its accessors.

`posterior_draws(fit, n_draws, seed, include_latent, max_output_mb)` returns
joint draws with parameter names matching the fit contract. The first array
dimension is always draw. For MCMC, chain and source iteration are preserved.
For CAVI/SVI, draws are explicitly labeled as samples from the variational
approximation and never receive MCMC diagnostics. Latent group allocations
appear only for CHLRM and only when requested.

`posterior::as_draws()` converts a fit or posterior-draw result to
`draws_df`. `tidy_parameters()` returns one row per scalar parameter element.

## Prediction and model assessment

The package separates five questions.

### Descriptive in-sample fit

`evaluate_in_sample_fit(fit, return_details = FALSE)` asks how close the fitted
mean is to the training response. It returns MSE and R-squared, and optionally
fitted values and residuals. It is neither posterior calibration nor
out-of-sample generalization.

### Posterior predictive checks

`posterior_predictive_check()` asks whether replicated responses generated by
the fitted model reproduce selected characteristics of the observed data. The
standard discrepancies cover global mean, global standard deviation, MSE
against a fixed observed-data OLS reference, and applicable between/within
group variation. `ppc_discrepancies()` selects a subset.

The result includes observed and replicated discrepancies, posterior
predictive probabilities, and Monte Carlo uncertainty. For MCMC, indicator
ESS respects chain order; for CAVI/SVI, it describes independent draws from
the approximation. `y_rep` is returned only when requested. An inapplicable
discrepancy fails explicitly rather than producing a placeholder.

`pp_check_from_prediction()` reuses a response-prediction result. Posterior
prediction remains a separate responsibility; PPC consumes it.

### Posterior mean prediction

`predict(fit, newdata, group_ids, n_draws, level, seed, return_draws,
max_output_mb)` returns `posterior-prediction-2` with `target = "mean"`. Its
interval is a pointwise credible interval for the conditional mean.

### Future-response prediction

`posterior_predict()` uses the same data arguments and returns
`posterior-prediction-2` with `target = "response"`. It includes observational
noise and therefore produces pointwise predictive intervals. With
`return_draws = TRUE`, it retains response draws in `y_rep`.

HLRM and CHLRM distinguish observed from new groups. New-group HLRM prediction
draws new local parameters from the population hierarchy. New-group CHLRM
prediction integrates or samples a component assignment under the mixture.

### Predictive criteria and held-out performance

`predictive_criteria()` returns WAIC and applicable diagnostics for an explicit
predictive target. With `return_log_likelihood = TRUE`, it also returns the
pointwise `log-likelihood-1` object used in the calculation.

`predictive_log_likelihood()` and `predictive_score()` evaluate held-out
responses using ELPD/MLPD. `predictive_performance()` reports applicable RMSE,
MAE, R-squared, interval coverage/width, and logarithmic scores in a contract
separate from in-sample fit and PPC.

`cross_validation_plan()`, `cross_validation_split()`, and
`cross_validate_model()` support:

- observation K-fold for new LRM observations;
- within-group K-fold for new observations in known HLRM/CHLRM groups;
- group K-fold for new groups;
- leave-one-group-out for one new group per refit.

Cross-validation is sequential and requires an explicit maximum fit count.
Massive grids, parallel scheduling, and experiment checkpoints belong to the
research project.

`psis_loo()` delegates Pareto smoothing to `loo` and accepts only converged
MCMC fits. The leave-one-out importance identity requires posterior draws;
CAVI/SVI use refitted cross-validation rather than treating draws from `q` as
exact posterior draws.

DIC is not part of the API. Its dependence on a posterior point estimate and
parametrization is particularly problematic for hierarchical and mixture
models. WAIC, PSIS-LOO when valid, refitted CV, and held-out log scores retain
an explicit predictive estimand.

## Likelihood ownership and predictive estimands

The likelihood belongs to the model specification, not to an inference engine
or GOF function. The same validated probability model drives fitting,
`log_likelihood()`, posterior prediction, and conditional simulation. Changing
from MCMC to CAVI or SVI does not change the observational model.

The high-level backend selects conditioning from the requested predictive
target:

| Model | Target | Unit | Likelihood convention |
|---|---|---|---|
| LRM | observed data | observation | conditional on coefficients and variance |
| LRM | supplied group | group | same conditional likelihood, aggregated |
| HLRM | observed groups | observation/group | conditional on local parameters |
| HLRM | new groups | group | local coefficients and variance marginalized |
| CHLRM | new groups | group | component allocation marginalized |

The advanced `log_likelihood()` interface allows explicit `unit` and
`conditioning`. Its matrix stores predictive units in rows and draws in
columns:

$$L_{us}=\log p(y_u\mid\theta^{(s)}).$$

`as_log_lik_matrix()` converts to the draws-by-unit convention used by external
tools while preserving metadata. `waic_from_loglik()` is the transparent WAIC
primitive; `predictive_criteria()` is the recommended user interface.

WAIC, PSIS-LOO, and log scores are comparable only when dataset, predictive
unit, predictive target, conditioning/marginalization, and likelihood
convention match. Observation-level and group-level WAIC are not directly
comparable. Conditional HLRM and marginal new-group HLRM scores are not the
same estimand. WAIC based on draws from `q` evaluates the fitted variational
approximation; it does not prove that CAVI/SVI reproduce MCMC.

## Data scaling

`data_scaling()` declares explicitly selected predictors, response scaling,
and the intercept column. `fit_scaling()` learns a `data_scaling_map` using
training data only. `transform_data()` and `transform_prior()` apply the map;
`scaling_map(fit)` exposes it for audit.

`predictors = "all"` selects every non-intercept column. Numeric columns are
never guessed to be continuous: categorical indicators remain unchanged unless
the user explicitly selects them. Names, order, and intercept placement are
validated for new data.

With `fit_model(..., scaling = data_scaling(...))`, scaling is automatic but
not implicit: the fitted object retains the specification and learned map.
Priors are declared in the input scale. Public parameters, draws, summaries,
predictions, PPC results, and log likelihoods are restored to the input scale.
Every cross-validation fold learns its map from fold training data only.

If coefficients transform as
$\boldsymbol\beta=A+T\boldsymbol\beta^*$, then prior parameters transform as
$m^*=T^{-1}(m-A)$ and $V^*=T^{-1}VT^{-\mathsf T}$. Hierarchical covariance
uses the same linear map. Residual variances and $\xi^2$ use the squared
response scale. With a Gamma shape--rate prior and
$\xi^{2*}=\xi^2/s_y^2$, the transformed rate is $b^*=b s_y^2$. With an
inverse-gamma shape--scale prior, the transformed scale is divided by
$s_y^2$.

The public log likelihood is on the input scale. Returning from a transformed
response applies the appropriate Jacobian per observation or group.

## Simulation contract

Simulation separates the DGP, covariate design, and sample structure:

```r
simulate_data(
  generator,
  design,
  sample,
  seed = NULL,
  max_output_mb = 256
)
```

Three DGPs answer different questions:

| DGP | Fixed | Regenerated | Primary use |
|---|---|---|---|
| `conditional_dgp()` | likelihood parameters | remaining latent variables and response | recovery/stress at fixed truth |
| `hierarchical_dgp()` | population parameters | local/component parameters, latent state, response | integrated performance for new units |
| `prior_predictive_dgp()` | design and prior | all parameters, latent state, response | prior prediction and SBC |

LRM has no intermediate local hierarchy, so `hierarchical_dgp()` does not
apply. HLRM and CHLRM support all three modes. In CHLRM, conditional truth
contains component coefficients, component variances, and mixture weights;
group allocations are sampled by default but may be supplied for paired
designs.

`simulation_design()` declares fixed or random covariates, distribution,
covariance, intercept handling, and split behavior. `simulation_sample()`
declares observation/group sizes and optional test sizes. The simulator
generates both `X` and `y`, validates support and SPD matrices, isolates RNG,
and returns `simulation-result-2` with `train`, optional `test`, `truth`,
`configuration`, `conditioning`, `parameter_source`, and `design_source`.

`truth` is output, not a simulation strategy. Observation-level truth is
retained only when requested to avoid large redundant objects. `tidy_truth()`
converts retained truth to a long representation.

Conditional repeated simulation estimates frequentist performance at a fixed
scenario. Hierarchical simulation evaluates population-level recovery and new
groups. SBC requires exactly the same prior and likelihood in generation and
fitting; changing either creates a robustness study, not SBC. Raw CHLRM labels
are not identifiable, so recovery uses relabeling or label-invariant
quantities.

The package does not define paper coefficient geometries, scenario names,
replication counts, or master seeds. A research project constructs those
values and supplies them to the typed DGP.

## Model-specific accessors and label switching

All fits implement `print()`, `summary()`, `coef()`, `fitted()`, `residuals()`,
and `nobs()` where their conventional meaning is defined.

HLRM and CHLRM implement `population_parameters()`, `group_parameters()`, and
`group_deviations()`. Group deviations summarize the joint quantity
$\beta_j-\beta$ rather than manufacturing intervals from unrelated marginal
summaries.

CHLRM additionally implements `components()`, `responsibilities()`,
`classify_groups()`, `classification_uncertainty()`, and `co_clustering()`.
Component-specific posterior summaries use an explicit relabeling reference.
Co-clustering probabilities and other permutation-invariant summaries remain
valid without assigning semantic meaning to raw labels. LRM and HLRM do not
expose mixture accessors.

`relabel_fit(fit, reference)` applies a compatible reference obtained from
`fit_orchestration()` and returns a new CHLRM fit. It consistently transforms
canonical component parameters, latent allocations, posterior draws,
summaries, and mixture accessors without mutating the source fit. References
record their coordinate scale. A standardized-backend reference is accepted
only when the fitted scaling maps are identical.

## Built-in plots

The package provides dependency-light base-R `plot()` methods for:

- fitted-model convergence histories and coefficient intervals;
- scalar posterior densities from `posterior_draws()`;
- CHLRM allocation-probability heatmaps;
- posterior predictions and their pointwise intervals;
- posterior predictive discrepancy distributions;
- pointwise log-likelihood summaries;
- chain/start elapsed times from `fit_orchestration()`;
- simulated response-versus-predictor views.

Plot methods perform no inference and invisibly return the data used in the
figure. Publication-specific layouts and multi-panel paper figures belong to
the research repository.

## Robust execution and resource control

`fit_orchestration()` returns a read-only audit of chains or starts without
exposing mutable R6 backends. MCMC runs multiple chains sequentially to avoid
automatic CPU/RAM oversubscription. CHLRM chains are aligned to a common label
reference before joint diagnostics.

Timing fields have deliberately narrow meanings. For a single fit,
`fit$timing$total_seconds` is backend fitting time. For variational multistart
it is the selected attempt's time; all attempt times, including `NA` for
failures without a valid duration, are in
`fit_orchestration(fit)$summary$elapsed_seconds`. For multichain MCMC it is the
sum of per-chain compute times, also available in `chain_metadata`; because
chains run sequentially, this is not declared to be process-level wall time.
`update_seconds` and `objective_seconds` partition timed variational fitting
work when available. Downstream prediction and assessment time, complete
`fit_model()` call wall time, and operating-system process time are not stored
and must be measured externally when an experiment requires them.

`resume_mcmc()` returns a new continued fit and recalculates applicable
diagnostics and CHLRM alignment. `save_mcmc_checkpoint()` and
`read_mcmc_checkpoint()` persist validated flat state rather than serializing
R6 environments. Saving refuses accidental overwrite.

`estimate_fit_resources()` uses the same preflight calculation as
`fit_model()`. It distinguishes estimated fitting peak, retained state, and
potential evaluation storage. `max_memory_mb` is available in every inference
control. The estimate includes safety overhead but is not an operating-system
hard limit or a measured RSS value.

## Interoperability and I/O boundary

The package adopts useful conventions without depending on `lme4`, `brms`,
`bayesplot`, `mclust`, or `mixtools`:

- `posterior` is the canonical draws representation;
- `as_yrep()` provides draws-by-observation response matrices;
- `as_log_lik_matrix()` provides draws-by-unit log-likelihood matrices;
- `loo` performs PSIS-LOO;
- standard S3 fit methods follow base R conventions;
- CHLRM accessors follow established mixture-model questions.

The typed matrix/list API is authoritative. A `formula + data` interface is not
part of the current contract.

The package returns serializable objects and writes files only through an
explicit operation such as `save_mcmc_checkpoint()`. Consumers own paths,
filenames, compression, manifests, cache retention, plots, tables, and
experiment checkpoints.

## Non-public implementation

The following remain internal and may change without breaking the API when
numerical and statistical contracts are preserved:

- `sample_*`, `update_q_*`, and `initialize_*` kernels;
- total log-density and specialized integration helpers;
- Cholesky and validation helpers;
- native draw preparation and `new_*_result` constructors;
- relabeling internals;
- backend R6 classes and native `hyperparams`/`run_settings` lists.
