\name{orchestration}
\alias{save_mcmc_checkpoint}
\alias{read_mcmc_checkpoint}
\alias{resume_mcmc}
\title{Checkpoint and Resume a Typed MCMC Fit}
\usage{
save_mcmc_checkpoint(fit, path)

read_mcmc_checkpoint(path)

resume_mcmc(fit, additional_samples, seed = NULL, verbose = FALSE)
}
\arguments{
  \item{fit}{An MCMC \code{vb_model_fit} produced by \code{fit_model()} or
  restored by \code{read_mcmc_checkpoint()}.}
  \item{path}{One nonempty path. Saving refuses to overwrite an existing
  checkpoint; use a versioned path.}
  \item{additional_samples}{Positive number of post-warmup draws to append to
  each chain.}
  \item{seed}{Optional seed for an explicit new branch. If omitted, every
  chain resumes from its saved RNG state.}
  \item{verbose}{Whether to print sampler progress.}
}
\description{
Save a typed MCMC fit as plain serializable sampler state, reconstruct it, or
continue every chain without repeating warm-up. Multichain continuation
preserves each terminal chain/RNG state and recomputes joint relabeling and
diagnostics. \code{resume_mcmc()} returns a new fit and does not mutate its
input.
}
\value{
\code{save_mcmc_checkpoint()} invisibly returns the normalized path.
\code{read_mcmc_checkpoint()} and \code{resume_mcmc()} return an MCMC
\code{vb_model_fit}.
}
\seealso{\code{fit_model()}, \code{mcmc_control()}}
