# Package roadmap

## Purpose

This file contains only unfinished work, execution order, and acceptance
criteria for `VariationalBayesHLRM`. The implemented statistical and software
contract is defined in [API.md](API.md).

The package owns reusable simulation, inference, posterior analysis,
prediction, model assessment, and resource controls. The companion research
repository owns the thesis, paper, reviewer scenarios, experimental grids,
Farms data analysis, caches, figures, and tables.

English is required for source code, public names, messages, documentation,
and generated output.

## Status

| Area | Status | Remaining gate |
|---|---|---|
| Canonical documentation | complete | Keep code, tests, and API contract synchronized |
| Structural cleanup | complete | Preserve numerical equivalence |
| Generated reference documentation and vignettes | complete | Keep examples executable |
| Pre-1.0 public workflows | complete | Preserve installed-API tests |
| Release automation | locally complete | Verify remote CI and tag the release |
| Performance optimization | waiting | Requires confirmatory E1--E5 evidence |
| Stable v1 | deferred | Requires every preceding gate |

## Release preparation

The local release candidate has completed the dependency audit, clean
installation, executable examples, vignette builds, source build, all package
tests, installed-package research validation, and the Farms E5 integration
smoke. `R CMD check --as-cran` reports zero errors, zero warnings, and only the
expected new-submission note.

Remaining release actions:

1. Commit and push the release candidate.
2. Verify GitHub Actions on supported R versions for Linux, macOS, and Windows.
3. Tag the pre-1.0 release after remote CI is green.

The source tree must not track `.Rcheck` directories, archives, RDS files,
experimental outputs, thesis or paper sources, or machine-specific paths.

## Evidence-driven performance work

Performance changes follow confirmatory E1--E5 measurements from the research
repository. Every optimization must have a numerical-equivalence test and
must preserve the same estimand and result schema. Validate
`estimate_fit_resources()` against observed runtime and memory, and report its
error and conservatism.

Parallel experiment execution, research caching, scenario grids, and paper
figures remain outside the package.

## Compatibility before v1

The cleanup does not intentionally break the public API. A controlled
pre-1.0 correction is allowed only when the current contract is inconsistent.
The same increment must update `API.md`, schema/version decisions,
implementation, tests, and the research consumer. Git history and tagged
releases record completed changes; this pre-release package does not maintain
a separate changelog.

Internal R6 backends and kernel notation are private. Public function names,
S3 methods, result schemas, predictive targets, units, and conditioning rules
are public.

## Deferred or excluded

Deferred until after v1 unless publication requires them:

- a `formula + data` interface;
- a pkgdown site;
- optional bayesplot convenience methods;
- a general empirical-Bayes prior constructor;
- parallel experiment execution.

DIC, implicit standardization, a prior called "non-informative", paper-specific
scenario constructors, and manuscript workflows are not planned package
features.

## Validation gate

For each completed increment:

1. Run the complete package test suite without failures, errors, warnings, or
   skips.
2. Run deterministic public-simulation fixtures for all nine combinations of
   LRM, HLRM, and CHLRM with MCMC, CAVI, and SVI.
3. Exercise fitting, posterior draws, in-sample fit, posterior predictive
   checks, and predictive criteria.
4. Verify result classes, schema versions, parameter names, model-specific
   fields, RNG restoration, and numerical equivalence.
5. Confirm that non-applicable quantities are omitted rather than represented
   by placeholders.
6. Build, check, and install the source package in a clean temporary library.
7. Run the companion research API tests and Farms smoke consumer against that
   installed package.

## Definition of done

The package is ready for a stable v1 when `API.md`, implementation, generated
reference documentation, tests, namespace, vignettes, and the installed
research consumer use the same statistical contract; CI and release checks are
green; and confirmatory evidence supports the published resource claims.
