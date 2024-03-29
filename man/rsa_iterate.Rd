% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mvpa_iterate.R
\name{rsa_iterate}
\alias{rsa_iterate}
\title{rsa_iterate}
\usage{
rsa_iterate(
  mod_spec,
  vox_list,
  ids = 1:length(vox_list),
  permute = FALSE,
  regtype = c("pearson", "spearman", "lm", "rfit"),
  distmethod = c("spearman", "pearson")
)
}
\arguments{
\item{mod_spec}{An object of class \code{rsa_model} specifying the RSA model.}

\item{vox_list}{A \code{list} of voxel indices or coordinates for each voxel set.}

\item{ids}{A \code{vector} of IDs for each voxel set (defaults to 1:length(vox_list)).}

\item{permute}{Logical, whether to permute the labels (defaults to FALSE).}

\item{regtype}{A character string specifying the analysis method. One of: \code{"pearson"}, \code{"spearman"}, \code{"lm"}, or \code{"rfit"} (defaults to "pearson").}

\item{distmethod}{A character string specifying the method used to compute distances between observations. One of: \code{"pearson"} or \code{"spearman"} (defaults to "spearman").}
}
\description{
Runs representational similarity analysis (RSA) for each voxel set in a list.
}
