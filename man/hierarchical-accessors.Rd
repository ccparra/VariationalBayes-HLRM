\name{hierarchical-accessors}
\alias{population_parameters}
\alias{population_parameters.vb_model_fit}
\alias{group_parameters}
\alias{group_parameters.vb_model_fit}
\alias{group_deviations}
\alias{group_deviations.vb_model_fit}
\title{Inspect Hierarchical Parameter Blocks}
\description{
Extract population-level and observed-group parameter summaries using names
that match the package's hierarchical model rather than importing fixed/random
effect terminology from a different parameterization.
}
\usage{
population_parameters(fit, ...)

\method{population_parameters}{vb_model_fit}(fit, ...)

group_parameters(fit, ...)

\method{group_parameters}{vb_model_fit}(fit, ...)

group_deviations(fit, ...)

\method{group_deviations}{vb_model_fit}(fit, ...)
}
\arguments{
  \item{fit}{An HLRM or CHLRM object returned by \code{fit_model()}.}
  \item{...}{Arguments forwarded to \code{summarize_posterior()}.}
}
\details{
Population parameters apply to HLRM and CHLRM. Group parameters and group
deviations apply only to HLRM because CHLRM models groups through component
allocations and exposes those quantities through the CHLRM accessors.

\code{group_deviations()} returns
\eqn{E(\beta_j-\beta\mid y)=E(\beta_j\mid y)-E(\beta\mid y)}. It deliberately
does not attach intervals: intervals for the difference require its joint
posterior draws, not subtraction of marginal interval endpoints.
}
\value{
\code{population_parameters()} and \code{group_parameters()} return named
lists of parameter summaries. \code{group_deviations()} returns a
group-by-coefficient matrix of posterior means.
}
\seealso{\code{components()}, \code{coef.vb_model_fit()},
\code{summarize_posterior()}}
