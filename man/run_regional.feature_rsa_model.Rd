% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/regional.R, R/feature_rsa_model.R
\name{run_regional.feature_rsa_model}
\alias{run_regional.feature_rsa_model}
\title{Specialized run_regional method for feature_rsa_model}
\usage{
\method{run_regional}{feature_rsa_model}(
  model_spec,
  region_mask,
  coalesce_design_vars = FALSE,
  processor = NULL,
  verbose = FALSE,
  ...
)

\method{run_regional}{feature_rsa_model}(
  model_spec,
  region_mask,
  coalesce_design_vars = FALSE,
  processor = NULL,
  verbose = FALSE,
  ...
)
}
\arguments{
\item{model_spec}{A \code{feature_rsa_model} object.}

\item{region_mask}{A mask representing different brain regions.}

\item{coalesce_design_vars}{If TRUE, merges design variables into prediction table.}

\item{processor}{A custom processor function for ROIs. If NULL, uses defaults.}

\item{verbose}{Print progress messages.}

\item{...}{Additional arguments}
}
\description{
This function runs a regional analysis using a feature RSA model and region mask.
}
\details{
This integrates `feature_rsa_model` with the MVPA framework, similar to `run_regional.mvpa_model`.
}
