% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/searchlight.R
\name{run_searchlight.vector_rsa}
\alias{run_searchlight.vector_rsa}
\title{run_searchlight method for vector_rsa_model}
\usage{
\method{run_searchlight}{vector_rsa}(
  model_spec,
  radius = 8,
  method = c("randomized", "standard"),
  niter = 4,
  ...
)
}
\arguments{
\item{model_spec}{A \code{vector_rsa_model} object.}

\item{radius}{Numeric searchlight radius (1 to 100).}

\item{method}{Character: "standard" or "randomized".}

\item{niter}{Number of iterations if \code{method="randomized"}.}

\item{...}{Additional arguments passed on to \code{do_standard} or \code{do_randomized}.}
}
\description{
This sets a custom \code{mvpa_fun} (e.g., \code{vector_rsa_iterate}) or 
different combiners for standard vs. randomized, etc.
}
