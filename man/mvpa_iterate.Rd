% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mvpa_iterate.R
\name{mvpa_iterate}
\alias{mvpa_iterate}
\title{MVPA Iteration for Voxel Sets with Parallelization}
\usage{
mvpa_iterate(
  mod_spec,
  vox_list,
  ids = 1:length(vox_list),
  compute_performance = TRUE,
  return_predictions = TRUE,
  return_fits = FALSE,
  batch_size = as.integer(0.1 * length(ids)),
  permute = FALSE,
  verbose = TRUE
)
}
\arguments{
\item{mod_spec}{An object of class \code{mvpa_model} specifying the model.}

\item{vox_list}{A \code{list} of voxel indices or coordinates.}

\item{ids}{A \code{vector} of IDs for each voxel set (defaults to 1:length(vox_list)).}

\item{compute_performance}{A \code{logical} indicating whether to compute and store performance measures for each voxel set (defaults to TRUE).}

\item{return_predictions}{A \code{logical} indicating whether to return row-wise predictions for each voxel set (defaults to TRUE).}

\item{return_fits}{A \code{logical} indicating whether to return the model fit for each voxel set (defaults to FALSE).}

\item{batch_size}{An \code{integer} specifying the number of voxel sets to process in each batch (defaults to 10% of the total voxel sets).}

\item{permute}{A \code{logical} indicating whether to permute the labels (defaults to FALSE).}

\item{verbose}{A \code{logical} indicating whether to print progress messages (defaults to TRUE).}
}
\value{
A \code{data.frame} containing the results for each voxel set, including performance measures, predictions, and model fits, as specified by the input parameters.
}
\description{
This function fits a classification or regression model for each voxel set in a list using parallelization.
It can compute and store performance measures, return row-wise predictions, and return the model fit for each voxel set.
}
\details{
This function utilizes parallel processing to speed up the process of fitting the specified model for each voxel set in a list.
The parallelization is achieved using the `furrr` package, which provides a parallel backend for the `purrr` package.
By default, it divides the voxel sets into batches and processes them in parallel according to the specified batch size.
The function provides options to control the return of performance measures, predictions, and model fits for each voxel set.
}
