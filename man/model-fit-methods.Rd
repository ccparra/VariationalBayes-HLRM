\name{model-fit-methods}
\alias{print.vb_model_fit}
\alias{summary.vb_model_fit}
\alias{coef.vb_model_fit}
\alias{fitted.vb_model_fit}
\alias{residuals.vb_model_fit}
\alias{nobs.vb_model_fit}
\title{Standard Methods for Fitted Bayesian Regression Models}
\description{
Inspect a fitted model through conventional R generics without exposing its
mutable inference backend.
}
\usage{
\method{print}{vb_model_fit}(x, ...)

\method{summary}{vb_model_fit}(object, ...)

\method{coef}{vb_model_fit}(object, ...)

\method{fitted}{vb_model_fit}(object, ...)

\method{residuals}{vb_model_fit}(object, ...)

\method{nobs}{vb_model_fit}(object, ...)
}
\arguments{
  \item{x, object}{An object returned by \code{fit_model()}.}
  \item{...}{Additional arguments passed to the underlying versioned operation
  where applicable.}
}
\details{
\code{summary()} returns the same \code{posterior-summary-1} object as
\code{summarize_posterior()}. \code{coef()} returns posterior or variational
means: a named vector for LRM; population and group coefficient blocks for
HLRM; and population and component coefficient blocks for CHLRM. It does not
invent mixture or latent fields for models where they do not apply.

\code{fitted()} and \code{residuals()} return numeric vectors in the original
row order. Grouped inputs are concatenated in group order and retain names that
identify group and within-group observation. These are descriptive in-sample
quantities, not out-of-sample predictions.
}
\value{
\code{print()} returns \code{x} invisibly; \code{summary()} returns a posterior
summary; \code{coef()} returns a vector or model-specific named list;
\code{fitted()} and \code{residuals()} return numeric vectors; and
\code{nobs()} returns the total observation count.
}
\seealso{\code{fit_model()}, \code{summarize_posterior()}, \code{predict()}}
