\name{tidy-parameters}
\alias{tidy_parameters}
\alias{tidy_parameters.vb_model_fit}
\alias{tidy_parameters.vb_posterior_summary}
\title{Represent Posterior Parameter Summaries in Long Form}
\description{
Convert all and only the applicable parameter summaries of a fitted model to a
dependency-free long data frame suitable for tables, plots, and replicated
simulation studies.
}
\usage{
tidy_parameters(x, ...)

\method{tidy_parameters}{vb_model_fit}(x, ...)

\method{tidy_parameters}{vb_posterior_summary}(x, ...)
}
\arguments{
  \item{x}{A fitted model or \code{posterior-summary-1} object.}
  \item{...}{For a fitted model, arguments forwarded to
  \code{summarize_posterior()}.}
}
\details{
Vector elements use \code{index_1}; matrix-valued parameters use
\code{index_1} and \code{index_2}. The output retains the distribution label,
interval level, existence of the posterior mean, and whether the pointwise
interval contains zero. Nonapplicable model quantities are absent rather than
encoded as null rows.
}
\value{A data frame with class \code{vb_tidy_parameters}.}
\seealso{\code{summarize_posterior()}, \code{tidy_truth()}}
