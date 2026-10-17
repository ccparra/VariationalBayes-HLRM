# Package roadmap

## Authority and objective

This is the only execution roadmap for `VariationalBayesHLRM`. The current
public behavior is defined by [API.md](API.md); this document records work that
is not yet complete. A milestone is complete only when implementation,
generated reference documentation, tests, and the API contract agree.

The objective is an installable, tested, documented, and reproducible R
package for simulation, Bayesian inference, posterior analysis, prediction,
and model assessment with LRM, HLRM, and CHLRM. The package provides reusable
statistical mechanisms. Datasets, experimental grids, replication counts,
figures, tables, and manuscripts belong to the companion research repository.

## Verified baseline

Baseline verified on 2026-10-17:

- package version `0.1.0`;
- 65 exported functions and generics;
- 29 registered S3 methods;
- MCMC, CAVI, and SVI implemented for LRM, HLRM, and CHLRM;
- 150 test blocks with zero failures, errors, warnings, or skips;
- source package builds successfully;
- `R CMD check --no-manual` reports `Status: OK`;
- package and research repositories have separate histories and boundaries.

## Package boundary

The package contains general model, prior, initialization, simulation,
inference, posterior, prediction, diagnostic, cross-validation, and resource
control mechanisms. It may contain deterministic unit and integration tests
and small synthetic examples.

The package does not contain:

- thesis or manuscript sources;
- reviewer-specific scenarios, sample sizes, seeds, or coefficient grids;
- Farms or other analysis-specific data unless deliberately released as
  documented package data;
- experiment schedulers, result caches, plots, or reviewer-response tables;
- parallel research runners or research checkpoint policies.

The companion project `ccparra/VariationalBayes-HLRM-research` consumes the
installed public API and owns E1--E5, the thesis, the paper, and their outputs.

## Language policy

English is the canonical language for source code, comments, public names,
errors, warnings, documentation, tables, plots, labels, captions, and other
package-generated output. Mathematical symbols may remain inside kernels when
they map directly to the theory, but public metadata and labels are descriptive
English. Research artifacts follow the same output-language rule.

## Milestones

### M0 — Canonical documentation

**Status:** complete.

- Keep `API.md` as the present-tense normative contract.
- Keep this file as the only plan and status document.
- Maintain both documents in English.
- Remove transitional, historical, and unresolved wording from the API.
- Keep README concise and link to these two sources of truth.

**Gate:** the API describes only implemented behavior; future work appears
only here; all public estimands and result schemas are explicit.

### M1 — Structural hygiene

**Status:** complete.

- Preserve the idiomatic flat `R/` layout.
- Remove the mixed-purpose `R/utils.R` by moving each helper to its owning
  subsystem.
- Split only modules that currently combine distinct responsibilities:
  simulation, fit orchestration, result contracts, specifications, MCMC
  continuation/persistence, and CHLRM variational updates/objective/draws.
- Use English comments and remove historical wording, dead branches, obsolete
  aliases, and unused helpers.
- Use descriptive names at public and orchestration boundaries. Mathematical
  notation may remain local to kernels when it maps directly to the theory.
- Make the documented source-tree validation command declare all development
  dependencies.

**Gate:** no public behavior or result schema changes; fixed-seed numerical
equivalence tests pass for all nine model--method combinations; package check
remains clean.

### M2 — API freeze and generated documentation

**Status:** in progress.

- Audit and freeze the public symbol manifest.
- Preserve the current 65 exports and 29 S3 methods unless a separately
  documented pre-1.0 correction is required.
- Move public reference documentation to roxygen2 and generate `NAMESPACE`
  and `man/` from source.
- Add fast executable examples for construction, fitting, simulation,
  prediction, assessment, scaling, and CHLRM accessors.
- Add vignettes for fitting and convergence, simulation and recovery,
  prediction and assessment, and CHLRM label switching.

**Gate:** no undocumented exports; no manually divergent Rd signatures;
examples and vignettes build in a clean session; the exported-symbol contract
test passes.

### M3 — Automation and release preparation

**Status:** pending M2.

- Translate the README to English.
- Add continuous integration for supported R versions on Linux, macOS, and
  Windows.
- Add repository/support metadata, vignette metadata, and a software citation.
- Audit Imports and Suggests.
- Run source installation, clean-session examples, vignette builds, and
  `R CMD check --as-cran`.
- Release the work as a pre-1.0 package version.

**Gate:** CI is green on the declared platform matrix and the source tarball
passes release checks without local paths or research files.

### M4 — Evidence-driven performance

**Status:** waiting for confirmatory E1--E5 evidence.

- Optimize only bottlenecks demonstrated by research profiles.
- Add numerical-equivalence tests around every optimization.
- Compare resource estimates with observed confirmatory runtime and memory.
- Document the conservatism and limitations of `estimate_fit_resources()`.

**Gate:** optimizations preserve the same estimand, numerical contract, and
public result schema; resource-estimation error is reported from confirmatory
measurements.

### M5 — Stable release

**Status:** deferred until M0--M4 and confirmatory evaluation are complete.

- Freeze the v1 public API and schema versions.
- Complete citation and lifecycle policy.
- Perform exported-symbol, reverse-dependency, and clean-source audits.
- Tag the stable release only after the research consumer passes against the
  release candidate.

## Compatibility policy before v1

The M1 hygiene work does not intentionally break the public API. A controlled
breaking correction is permitted before v1 only when the existing contract is
demonstrably inconsistent. It must update `API.md`, schema/version decisions,
implementation, tests, `NEWS.md`, and the research consumer in the same
increment.

Internal R6 backends and kernel notation are not public. Public function names,
S3 methods, result schemas, predictive targets, units, and conditioning rules
are public.

## Deferred work

The following items are intentionally outside the first stable release unless
publication makes them necessary:

- a `formula + data` interface;
- a pkgdown site;
- optional bayesplot-specific convenience methods;
- a general empirical-Bayes prior constructor;
- parallel experiment execution.

DIC, implicit standardization, a prior called “non-informative,” and
paper-specific scenario constructors are not planned public features.

## Validation after every milestone

1. Run the complete package test suite with no failures, errors, warnings, or
   skips.
2. Run fixed-seed smoke/equivalence cases for LRM, HLRM, and CHLRM under MCMC,
   CAVI, and SVI.
3. Verify result classes, schema versions, parameter names, model-specific
   fields, RNG restoration, and numerical equivalence.
4. Confirm that non-applicable quantities are omitted rather than represented
   by placeholders.
5. Build and check the source package in a clean session.
6. Install to a temporary library and run the companion research API tests and
   smoke consumers without starting confirmatory experiments.
7. Confirm that Git tracks no check directories, archives, RDS files,
   experimental outputs, thesis, paper, or machine-specific paths.

## Definition of done

The cleanup is complete when `API.md`, this roadmap, implementation, generated
reference documentation, tests, namespace, and the research consumer describe
and use the same statistical contract. Stable v1 additionally requires the M4
confirmatory evidence and a successful release-candidate audit.
