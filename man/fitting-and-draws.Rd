\name{fitting-and-draws}
\alias{fit_model}
\alias{summarize_posterior}
\alias{posterior_draws}
\alias{as_draws.vb_model_fit}
\alias{as_draws.vb_posterior_draws}
\alias{fit_orchestration}
\alias{data_scaling}
\alias{fit_scaling}
\alias{transform_data}
\alias{transform_prior}
\alias{scaling_map}
\title{Fit Models and Extract Posterior Information}
\description{
Fit a typed model with a selected inference engine, summarize only parameters
that belong to that model, extract consistently oriented draws, and inspect
multichain or multistart provenance. Optional scaling is learned from training
data and all public outputs are restored to the input scale.
}
\details{
\code{posterior::as_draws()} converts either a fitted model or a
\code{posterior-draws-1} object to a standard \code{draws_df}. Tensor indices
are preserved in bracketed variable names. MCMC chain provenance is retained
in the reserved column and original gapped/subsampled indices are retained as
\code{source_iteration} and \code{source_stored_draw}; the standard
\code{.iteration} column is necessarily reindexed by \pkg{posterior}.

\code{fit_orchestration()} returns a read-only
\code{fit-orchestration-1} audit object. Multichain fits expose chain metadata,
diagnostics, and block history. Multistart fits expose the start summary and
selected start; CHLRM also exposes a long table of hard classifications for
every completed start. Mutable backends, native parameter lists, and partial
failed states are not exposed.
}
\seealso{\code{lrm_model()}, \code{predictive_criteria()}}
