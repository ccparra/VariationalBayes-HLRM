\name{predictive-evaluation}
\alias{predict.vb_model_fit}
\alias{posterior_predict}
\alias{evaluate_in_sample_fit}
\alias{posterior_predictive_check}
\alias{ppc_discrepancies}
\alias{pp_check_from_prediction}
\alias{log_likelihood}
\alias{predictive_criteria}
\alias{waic_from_loglik}
\alias{predictive_log_likelihood}
\alias{predictive_score}
\alias{predictive_performance}
\alias{cross_validation_plan}
\alias{cross_validation_split}
\alias{cross_validate_model}
\alias{psis_loo}
\title{Prediction and Bayesian Model Evaluation}
\description{
Separate conditional-mean prediction, response prediction, descriptive
in-sample fit, posterior predictive checking, pointwise log likelihood,
predictive criteria, held-out scoring, refitted cross-validation, and
PSIS-LOO. Every result records its predictive target, unit, and conditioning.
}
\details{
PSIS-LOO is restricted to diagnosed MCMC draws. CAVI and SVI should be compared
with refitted cross-validation. Posterior predictive checks are generative
adequacy diagnostics and are not interchangeable with held-out performance.
}
