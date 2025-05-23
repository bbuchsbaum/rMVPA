% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/allgeneric.R
\name{performance}
\alias{performance}
\title{Compute Performance Metrics}
\usage{
performance(x, ...)
}
\arguments{
\item{x}{The classification/regression result object to evaluate.}

\item{...}{Additional arguments passed to method-specific performance functions.}
}
\value{
Named numeric vector of performance metrics.
}
\description{
Generic function to compute performance metrics from result objects.
}
\examples{
cres <- binary_classification_result(
  observed  = factor(c("a", "b")),
  predicted = factor(c("a", "b")),
  probs     = matrix(c(0.8, 0.2, 0.3, 0.7), ncol = 2,
                     dimnames = list(NULL, c("a", "b")))
)
performance(cres)
}
