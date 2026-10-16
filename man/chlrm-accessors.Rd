\name{chlrm-accessors}
\alias{components}
\alias{components.vb_model_fit}
\alias{responsibilities}
\alias{responsibilities.vb_model_fit}
\alias{classify_groups}
\alias{classify_groups.vb_model_fit}
\alias{classification_uncertainty}
\alias{co_clustering}
\alias{co_clustering.vb_model_fit}
\title{Inspect Components and Group Classification in CHLRM Fits}
\description{
Extract model-specific mixture quantities from a fitted CHLRM without
navigating internal result lists. These operations reject LRM and HLRM fits
instead of returning nonapplicable null fields.
}
\usage{
components(fit, ...)

\method{components}{vb_model_fit}(fit, ...)

responsibilities(fit, ...)

\method{responsibilities}{vb_model_fit}(fit, ...)

classify_groups(fit, ...)

\method{classify_groups}{vb_model_fit}(fit, ...)

classification_uncertainty(fit, ...)

co_clustering(fit, ...)

\method{co_clustering}{vb_model_fit}(fit, ...)
}
\arguments{
  \item{fit}{A CHLRM object returned by \code{fit_model()}.}
  \item{...}{Arguments forwarded to \code{summarize_posterior()}, such as a
  common label reference.}
}
\details{
Responsibilities are posterior allocation probabilities for MCMC and the
variational probabilities for CAVI/SVI. Classification uses their maximum;
uncertainty is one minus that maximum, while entropy is reported separately.
Co-clustering is the posterior similarity matrix for MCMC and
\eqn{\rho\rho^T} under the mean-field variational approximation.
}
\value{
\code{components()} returns the three component-parameter summaries;
\code{responsibilities()} returns a group-by-component matrix;
\code{classify_groups()} and \code{classification_uncertainty()} return a data
frame; and \code{co_clustering()} returns a square group similarity matrix.
}
\seealso{\code{summarize_posterior()}, \code{coef.vb_model_fit()}}
