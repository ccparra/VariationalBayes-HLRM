\name{estimate_fit_resources}
\alias{estimate_fit_resources}
\alias{print.vb_fit_resource_estimate}
\title{Estimate Resources Required by a Model Fit}
\usage{
estimate_fit_resources(model, data, control, scaling = NULL)

\method{print}{vb_fit_resource_estimate}(x, ...)
}
\arguments{
  \item{model}{A specification from \code{lrm_model()},
  \code{hlrm_model()} or \code{chlrm_model()}.}
  \item{data}{A list containing \code{X} and \code{y} in the representation
  required by the selected model.}
  \item{control}{An inference control from \code{mcmc_control()},
  \code{cavi_control()} or \code{svi_control()}.}
  \item{scaling}{Optional preprocessing specification from
  \code{data_scaling()}.}
  \item{x}{A resource estimate.}
  \item{...}{Unused.}
}
\value{
An object of class \code{vb_fit_resource_estimate}. It reports conservative
peak memory for fitting, retained fit storage, potential downstream evaluation
storage, the declared memory budget and whether the request is within budget.
}
\description{
Performs a dimension-based preflight estimate before fitting. The estimate
includes headroom for R object overhead and distinguishes memory used by the
fit from storage that a later retained posterior evaluation could request.
It is a safety guard, not a measured or guaranteed operating-system memory
bound.
}
\seealso{\code{fit_model()}, \code{mcmc_control()}, \code{cavi_control()},
\code{svi_control()}}
