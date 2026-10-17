# Package roadmap

## Objective

Deliver an installable, tested, and documented R package for simulation,
inference, posterior analysis, prediction, and model assessment with LRM, HLRM,
and CHLRM. The package exposes reusable statistical mechanisms; datasets,
experimental grids, replication counts, tables, and manuscripts belong to the
separate research repository.

The normative public surface is described in [API.md](API.md). A feature is
complete only when its implementation, documentation, and tests agree.

## Current state

| Area | Status | Remaining work |
|---|---|---|
| MCMC, CAVI, and SVI for LRM/HLRM/CHLRM | Complete | Profile only where research experiments identify a bottleneck |
| Typed model, prior, initialization, and control specifications | Complete | Preserve backward-incompatible changes for a later major version |
| Simulation contracts and conditional/hierarchical/prior-predictive DGPs | Complete | Add examples to a vignette |
| Posterior draws and summaries | Complete | Expand interoperability examples |
| Prediction and posterior predictive checks | Complete | Document grouped prediction targets in a vignette |
| WAIC, MLPD, CV/LOGO, and MCMC PSIS-LOO | Complete | Add a model-assessment workflow vignette |
| CHLRM relabeling and invariant clustering summaries | Complete | Document interpretation and limitations |
| Reversible data scaling | Complete | Add end-to-end examples |
| Resource estimation and execution budgets | Implemented | Validate estimates against confirmatory benchmarks |
| Package documentation and release automation | In progress | Roxygen consolidation, vignettes, CI, and release checks |

## Package boundary

The package may contain:

- general model, prior, initialization, simulation, and inference code;
- posterior draws, summaries, predictions, log-likelihoods, and diagnostics;
- generic cross-validation and resource controls;
- deterministic unit and integration tests.

It must not contain:

- thesis or manuscript sources;
- review-specific scenarios, sample sizes, or seeds;
- Farms or other analysis-specific datasets unless deliberately exposed as
  documented package data;
- experiment schedulers, result caches, figures, or reviewer-response tables.

The companion research project is maintained separately at
`ccparra/VariationalBayes-HLRM-research` and consumes only the installed public
API.

## Next increments

1. **Documentation:** consolidate manual Rd generation with roxygen2 and add
   concise vignettes for fitting, simulation, prediction, and assessment.
2. **Idiomatic R interface:** decide whether a `formula + data` layer provides
   enough value to justify its maintenance; keep the typed API authoritative.
3. **Continuous integration:** run package tests and `R CMD check` on supported
   R versions and operating systems.
4. **Performance:** use measured E1--E5 profiles to optimize only demonstrated
   bottlenecks; preserve numerical equivalence tests around each change.
5. **Release readiness:** complete citation metadata, lifecycle policy,
   exported-symbol audit, reverse-dependency checks, and a clean source build.

## Definition of done

An increment closes only when:

- there is one canonical public name and one schema version;
- LRM, HLRM, and CHLRM behave analogously where their mathematics permits;
- non-applicable model-specific quantities are omitted rather than filled with
  placeholders;
- theory, implementation, tests, and documentation use the same estimand;
- package checks pass without relying on research files or local paths.
