% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/regional.R
\name{run_regional.default}
\alias{run_regional.default}
\title{Default method for run_regional}
\usage{
\method{run_regional}{default}(model_spec, region_mask, ...)
}
\arguments{
\item{model_spec}{A generic model specification object.}

\item{region_mask}{A mask representing different brain regions.}

\item{...}{Additional arguments passed to \code{run_regional_base}.}
}
\value{
A \code{regional_mvpa_result} object.
}
\description{
By default, if an object does not implement a specific \code{run_regional.*} method,
this fallback will attempt to run the \code{run_regional_base} pipeline.
}
