% Please edit documentation in R/spatial_nmf.R
\name{spatial_nmf}
\alias{spatial_nmf}
\title{Spatial Non-negative Matrix Factorization}
\usage{
spatial_nmf(
  x,
  group_B = NULL,
  groups = NULL,
  k,
  mask = NULL,
  dims = NULL,
  transform = c("none", "shift", "auc", "auc_raw", "zscore", "relu", "abs"),
  min_val = 0,
  floor = -0.5,
  lambda = 0,
  fast = FALSE,
  graph = NULL,
  neighbors = 6,
  na_action = c("zero", "error"),
  return_maps = .is_map_input(x),
  return_data = FALSE,
  component_test = NULL,
  global_test = NULL,
  stability = NULL,
  voxelwise_stats = c("none", "stability_zp"),
  parallel = NULL,
  progress = FALSE,
  ...
)
}
\arguments{
\item{x}{A numeric subject-by-voxel matrix/data frame, a single map, or a
list of NeuroVol/NeuroSurface maps for group A.}

\item{group_B}{Optional list of NeuroVol/NeuroSurface maps for group B. Use
this only when \code{x} is map input.}

\item{groups}{Optional group labels for matrix input. Required for
\code{component_test = TRUE} or \code{global_test = TRUE} with matrix input.}

\item{k}{Number of components.}

\item{mask}{Mask object for volumetric/surface maps, or optional spatial
metadata for matrix input when \code{return_maps = TRUE} or \code{lambda > 0}.}

\item{dims}{Optional spatial dimensions for volumetric masks.}

\item{transform}{Input transformation before NMF: "none" requires
non-negative input; other values use the same semantics as
\code{\link{nmf_preprocess_maps}}.}

\item{min_val}{Minimum value after transformation.}

\item{floor}{Lower floor used by the "auc" and "auc_raw" transformations.}

\item{lambda}{Spatial regularization strength (0 = none).}

\item{fast}{Logical; use faster, lower-iteration defaults in the NMF fit.}

\item{graph}{Optional graph Laplacian list from \code{build_graph_laplacian}.}

\item{neighbors}{Neighborhood size for volumetric adjacency (6/18/26).}

\item{na_action}{How to handle NA values: "zero" or "error".}

\item{return_maps}{Logical; return component maps when spatial metadata is
available.}

\item{return_data}{Logical; include the subject-by-voxel data matrix.}

\item{component_test}{NULL to skip, TRUE for defaults, or a list of arguments
passed to \code{\link{spatial_nmf_component_test}}.}

\item{global_test}{NULL to skip, TRUE for defaults, or a list of arguments
passed to \code{\link{spatial_nmf_global_test}}.}

\item{stability}{NULL to skip, TRUE for defaults, or a list of arguments
passed to \code{\link{spatial_nmf_stability}}.}

\item{voxelwise_stats}{Optional voxelwise statistics. Use "stability_zp" to
return bootstrap z- and p-value maps derived from stability.}

\item{parallel}{Logical; enable parallel processing for inference functions.}

\item{progress}{Logical; report progress via progressr package.}

\item{...}{Additional arguments passed to the NMF optimizer.}
}
\value{
A \code{spatial_nmf_result} object with \code{fit}, optional \code{components},
optional \code{data}, and requested inference results.
}
\description{
A single high-level interface for spatial NMF. Use \code{spatial_nmf()} with a
non-negative subject-by-voxel matrix, or with a list of NeuroVol/NeuroSurface
maps. Optional preprocessing can make signed maps or matrices non-negative
before factorization.
}
\examples{
set.seed(1)
W <- matrix(runif(12 * 2), nrow = 12)
H <- matrix(runif(2 * 20), nrow = 2)
X <- W \%*\% H
result <- spatial_nmf(X, k = 2, max_iter = 20)

\dontrun{
result <- spatial_nmf(group_A_vols, group_B = group_B_vols, mask = mask, k = 5)
}
}
