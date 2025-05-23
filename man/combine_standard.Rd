% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/searchlight.R
\name{combine_standard}
\alias{combine_standard}
\title{Combine standard classifier results}
\usage{
combine_standard(model_spec, good_results, bad_results)
}
\arguments{
\item{model_spec}{A list containing the model specification}

\item{good_results}{A data frame containing the successful classifier results}

\item{bad_results}{A data frame containing the unsuccessful classifier results}
}
\value{
A list containing the combined performance matrix and other information
}
\description{
This function combines the standard classifier results from a good results data frame
by binding the performance rows together and optionally computes the observed probabilities.
}
\keyword{internal}
