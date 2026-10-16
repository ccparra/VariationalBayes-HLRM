\name{model-specifications}
\alias{lrm_model}
\alias{hlrm_model}
\alias{chlrm_model}
\alias{lrm_prior}
\alias{hlrm_prior}
\alias{chlrm_prior}
\alias{prior_moments}
\alias{initialization_from_prior}
\alias{initialization_ols}
\alias{initialization_kmeans}
\alias{initialization_values}
\title{Model, Prior, and Initialization Specifications}
\description{
Construct validated, data-independent model and prior specifications and an
explicit initialization strategy. Priors and initial values are separate
objects so an OLS start never silently changes the Bayesian model.
}
\seealso{\code{fit_model()}, \code{simulate_data()}}
