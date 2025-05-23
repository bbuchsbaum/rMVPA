% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/custom.R
\name{combine_custom_randomized}
\alias{combine_custom_randomized}
\title{Combine Custom Randomized Searchlight Results}
\usage{
combine_custom_randomized(dataset, iteration_results)
}
\arguments{
\item{dataset}{The original mvpa_dataset object.}

\item{iteration_results}{The raw tibble output from *all* iterations of `mvpa_iterate`.}
}
\value{
A `searchlight_result` object.
}
\description{
Internal function to combine results from a randomized custom searchlight run.
Averages results for each metric across overlapping spheres.
}
\keyword{internal}
