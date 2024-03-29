% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/allgeneric.R
\name{merge_predictions}
\alias{merge_predictions}
\title{Merge predictions from multiple models}
\usage{
merge_predictions(obj1, rest, ...)
}
\arguments{
\item{obj1}{The first object containing predictions.}

\item{rest}{The rest of the objects containing predictions.}

\item{...}{Additional arguments to be passed to the function.}
}
\value{
A combined object containing merged predictions from multiple models.
}
\description{
Combine predictions from several models applied to the same test set.
}
