# VariationalBayesHLRM 0.2.6

- CHLRM label validation scans every retained allocation in bounded chunks.
  Posterior allocation probabilities and occupied-component counts accumulate
  exact counts in draw blocks rather than allocating full-history logical
  matrices or a full transpose. Labels, component matching, posterior summaries,
  chain counts, RNG state and sampling transitions are unchanged. Complete
  histories remain in RAM; this does not provide disk-streamed sampling.

# VariationalBayesHLRM 0.2.5

- Checkpoint-store decoding no longer retains encoded numeric chunks alongside
  reconstructed arrays, and releases temporary assignment arguments promptly.
  Object checksums, validation, formats and decoded values remain unchanged.
- Typed checkpoint restoration binds validated plain state to fresh backends
  using R copy-on-modify semantics instead of serializing histories again.
  Legacy snapshots with mutable metadata retain the previous cloning path.
  Transactional sampler cloning and checkpoint writing are unchanged. Restored
  histories remain in RAM; this is not disk-streamed sampling.

# VariationalBayesHLRM 0.2.4

- Object-store checkpoint writing now reads a validated plain-state view of
  backend histories without first serializing and deep-cloning all chains.
  It retains R copy-on-modify protection and the existing manifest/object
  formats, checksums, writer lock and publication order. The default single-file
  checkpoint and transactional sampler/restoration copies are unchanged.
- Saving still retains the full fit and copies individual encoded chunks;
  restoration still reconstructs complete state in RAM. This does not provide
  disk-streamed sampling or guarantee a workload's RSS budget.

# VariationalBayesHLRM 0.2.3

- `save_mcmc_checkpoint()` accepts an optional immutable object `store` and
  `chunk_size`. Identical compressed blocks and metadata nodes are shared
  across complete checkpoint manifests; `read_mcmc_checkpoint()` restores both
  the original single-file format and the new store format without recomputing
  alignment, diagnostics or sampler state. Default serialization is unchanged.
- Stores validate object checksums, serialize writers with a lock and publish
  manifests last. They preserve all raw/aligned histories and RNG states but
  still build and restore full in-memory checkpoints; no RSS reduction or
  large-workload resource safety is claimed.
- Checkpoint reconstruction now restores the original complete `run_settings`,
  including its seed provenance, after using seed-free constructor settings.
  Continuation still uses the saved terminal RNG state.

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
