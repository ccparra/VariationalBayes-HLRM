\name{interoperability}
\alias{as_yrep}
\alias{as_yrep.vb_prediction_result}
\alias{as_yrep.vb_posterior_predictive_check}
\alias{as_log_lik_matrix}
\alias{as_log_lik_matrix.vb_log_likelihood}
\title{Extract Matrices for Bayesian R Ecosystem Tools}
\description{
Convert versioned prediction and log-likelihood results to the draw-by-unit
matrix orientation used by packages such as \pkg{bayesplot} and \pkg{loo}.
}
\usage{
as_yrep(x, ...)

\method{as_yrep}{vb_prediction_result}(x, ...)

\method{as_yrep}{vb_posterior_predictive_check}(x, ...)

as_log_lik_matrix(x, ...)

\method{as_log_lik_matrix}{vb_log_likelihood}(x, ...)
}
\arguments{
  \item{x}{A response prediction, posterior predictive check, or pointwise
  log-likelihood object.}
  \item{...}{Reserved for methods.}
}
\details{
Package-native prediction and likelihood objects store predictive units in
rows and draws in columns. These adapters return draws in rows and observations
or groups in columns. \code{as_yrep()} requires retained response draws;
latent-mean predictions are rejected. \code{as_log_lik_matrix()} preserves the
predictive-unit table and full estimand as matrix attributes.
}
\value{A numeric matrix with posterior or variational draws in rows.}
\seealso{\code{posterior_predict()}, \code{posterior_predictive_check()},
\code{predictive_criteria()}}
