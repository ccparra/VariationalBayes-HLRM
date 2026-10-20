# VariationalBayesHLRM 0.2.2

- `posterior_draws()` and `posterior::as_draws()` now accept exact ordered
  `draw_indices` and selected `parameters` for stored MCMC fits. Extraction
  copies only requested blocks and retains chain/iteration provenance, with
  support for latent-only CHLRM allocations and a selected-output budget.
- Parameter subsets identify schema `posterior-draws-2`; default calls retain
  the complete schema, output values, and random-stream behavior. The selectors
  do not apply to CAVI/SVI or provide disk streaming/checkpoint compaction.
- Standard draws conversion adds `source_chain` beside the existing source
  iteration/stored-draw columns, preserving original chain identifiers when
  `posterior` normalizes a selected subset's reserved indices.

# VariationalBayesHLRM 0.2.1

- `fit_model()` now starts single-start variational optimization from the
  explicit fitting seed when given a reusable initialization set, matching
  the existing multistart seed contract. Discarded constructor draws no
  longer change the optimization stream across initialization strategies.
  This changes single-start SVI trajectories with supplied states; retain
  the original package version and commit to reproduce historical runs.
- Added regression coverage for identical-state/seed equivalence across
  constructor strategies, agreement with the native seeded start, caller
  RNG preservation, and prior immutability.
