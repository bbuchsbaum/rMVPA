% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/regional.R
\name{run_regional.mvpa_model}
\alias{run_regional.mvpa_model}
\title{Run regional MVPA analysis on a specified MVPA model}
\usage{
\method{run_regional}{mvpa_model}(
  model_spec,
  region_mask,
  coalesce_design_vars = FALSE,
  processor = NULL,
  verbose = FALSE,
  ...
)
}
\arguments{
\item{model_spec}{An object of type \code{mvpa_model} specifying the MVPA model to be used.}

\item{region_mask}{A mask representing different regions in the brain image.}

\item{processor}{}

\item{...}{Additional arguments to be passed to the function.}
}
\value{
A \code{list} of type \code{regional_mvpa_result} containing a named list of \code{NeuroVol} objects,
where each element contains a performance metric and is labeled according to the metric used (e.g. Accuracy, AUC).
}
\description{
This function runs a regional MVPA analysis using a specified MVPA model and
region mask. The analysis can be customized to return model fits, predictions,
and performance measures.
}
