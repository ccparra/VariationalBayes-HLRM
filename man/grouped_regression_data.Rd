\name{grouped_regression_data}
\alias{grouped_regression_data}
\title{Build Grouped Regression Data from a Data Frame}
\usage{
grouped_regression_data(data, response, predictors, group, intercept = TRUE)
}
\arguments{
  \item{data}{A nonempty data frame.}
  \item{response}{One response column name or index.}
  \item{predictors}{One or more numeric predictor column names or indices.}
  \item{group}{One grouping column name or index.}
  \item{intercept}{Whether to prepend an explicit intercept column.}
}
\value{
A \code{vb_grouped_regression_data} list containing grouped \code{X} and
\code{y}, a stable mapping from integer group ids to input labels, and column
provenance. The same object can be passed to HLRM or CHLRM operations.
}
\description{
Creates the common list-of-matrices representation used by hierarchical and
clustered hierarchical regressions. It performs no implicit factor contrast
expansion; predictors must already be finite numeric columns.
}
\seealso{\code{hlrm_model()}, \code{chlrm_model()}, \code{fit_model()}}
