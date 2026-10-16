\name{inference-controls}
\alias{mcmc_control}
\alias{mcmc_stopping}
\alias{cavi_control}
\alias{cavi_stopping}
\alias{svi_control}
\alias{svi_stopping}
\alias{power_schedule}
\title{Inference and Stopping Controls}
\usage{
mcmc_control(warmup = 1000L, post_warmup_draws = 1000L,
  thinning = 1L, chains = 1L, evaluation_draws = NULL,
  stopping = NULL, diagnostic_scope = NULL, max_memory_mb = NULL)

cavi_control(max_iterations = 10000L, evaluation_draws = 500L,
  stopping = cavi_stopping(), starts = 1L,
  require_converged = NULL, max_memory_mb = NULL)

svi_control(max_iterations = 50000L, minibatch_size,
  learning_rate = power_schedule(), evaluation_draws = 500L,
  stopping = svi_stopping(), starts = 1L,
  require_converged = NULL, max_memory_mb = NULL)
}
\arguments{
  \item{warmup}{Number of MCMC warm-up iterations per chain.}
  \item{post_warmup_draws}{Maximum retained MCMC draws per chain.}
  \item{thinning}{MCMC thinning interval. Multiple chains require one so that
  autocorrelation is diagnosed instead of discarded.}
  \item{chains}{Number of independent MCMC chains.}
  \item{evaluation_draws}{Number of posterior draws made available by default
  to downstream evaluation operations.}
  \item{stopping}{A method-specific stopping specification from
  \code{mcmc_stopping()}, \code{cavi_stopping()} or \code{svi_stopping()}.}
  \item{diagnostic_scope}{For multichain MCMC, either \code{"global"} or
  \code{"all"}. The latter also diagnoses all local parameters.}
  \item{max_memory_mb}{Positive preflight memory budget in megabytes. If
  \code{NULL}, 1024 MB is used.}
  \item{max_iterations}{Maximum number of CAVI or SVI iterations.}
  \item{starts}{Number of independent variational starts.}
  \item{require_converged}{For multiple starts, whether only starts satisfying
  their stopping rule are eligible for selection.}
  \item{minibatch_size}{Number of observations for LRM or complete groups for
  HLRM/CHLRM sampled at each SVI iteration.}
  \item{learning_rate}{An SVI schedule created by a learning-rate
  constructor such as \code{power_schedule()}.}
}
\description{
Declare computational budgets, convergence rules, evaluation draw counts,
multichain or multistart settings, and the power learning-rate schedule. In
SVI, minibatches contain observations for LRM and complete groups for
HLRM/CHLRM.
}
\details{
\code{max_memory_mb} is a positive preflight budget for every fit, including
a single chain or variational start. It defaults to 1024 MB. It is checked
against \code{estimate_fit_resources()} before backend allocation. The
\code{require_converged} selection rule applies only to multiple variational
starts, and \code{diagnostic_scope} applies only to multiple MCMC chains.
}
\seealso{\code{fit_model()}}
