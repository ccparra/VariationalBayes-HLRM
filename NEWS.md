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
