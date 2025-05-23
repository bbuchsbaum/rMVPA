% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/searchlight.R
\name{do_randomized}
\alias{do_randomized}
\title{Perform randomized searchlight analysis}
\usage{
do_randomized(
  model_spec,
  radius,
  niter,
  mvpa_fun = mvpa_iterate,
  combiner = pool_randomized,
  ...
)
}
\arguments{
\item{model_spec}{An object specifying the model to be used in the searchlight analysis.}

\item{radius}{The radius of the searchlight sphere.}

\item{niter}{The number of iterations for randomized searchlight.}

\item{mvpa_fun}{The MVPA function to be used in the searchlight analysis (default is \code{mvpa_iterate}).}

\item{combiner}{The function to be used to combine results (default is \code{pool_randomized}).}

\item{...}{Additional arguments to be passed to the MVPA function.}
}
\description{
This function performs randomized searchlight analysis using a specified model, radius, and number of iterations.
It can be customized with different MVPA functions, combiners, and permutation options.
}
\keyword{internal}
