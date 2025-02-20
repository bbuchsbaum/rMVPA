% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/feature_rsa_model.R
\name{merge_results.feature_rsa_model}
\alias{merge_results.feature_rsa_model}
\title{Merge Multiple Results for Feature RSA Model}
\usage{
\method{merge_results}{feature_rsa_model}(obj, result_set, indices, id, ...)
}
\arguments{
\item{obj}{A \code{feature_rsa_model} object}

\item{result_set}{A data frame of results from cross-validation folds}

\item{indices}{The voxel indices used (may not be relevant for feature_rsa_model)}

\item{id}{An identifier for the merged result (e.g., ROI id)}

\item{...}{Additional arguments}
}
\value{
A tibble with merged results
}
\description{
Merge Multiple Results for Feature RSA Model
}
