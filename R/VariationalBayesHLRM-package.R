#' VariationalBayesHLRM: Bayesian Gaussian regression
#'
#' Fit LRM, HLRM, and CHLRM models with MCMC, CAVI, or SVI through one typed
#' interface. The package separates inference, simulation, posterior
#' prediction, model adequacy checks, and out-of-sample predictive evaluation.
#'
#' The primary entry point is [fit_model()]. Model, prior, initialization,
#' algorithm control, and optional scaling are declared explicitly before
#' fitting. See `docs/API.md` in the source repository for the complete public
#' contract.
#'
#' @importFrom R6 R6Class
#' @importFrom posterior as_draws
#' @importFrom stats coef fitted nobs predict residuals rgamma setNames var
#' @importFrom utils head modifyList object.size setTxtProgressBar tail txtProgressBar
#' @keywords internal
"_PACKAGE"
