% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/allgeneric.R
\name{compute_performance}
\alias{compute_performance}
\title{Compute Performance for an Object}
\usage{
compute_performance(obj, result)
}
\arguments{
\item{obj}{The input object.}

\item{result}{The classification/regression result object to evaluate.}
}
\value{
A list of performance metrics.
}
\description{
Delegates calculation of performance metrics to the appropriate method.
}
