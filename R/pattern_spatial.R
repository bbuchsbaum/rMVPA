# Feature-aligned spatial graphs
#
# A spatial_graph is the anatomical adjacency structure used by spatially
# regularized estimators. Its single invariant is that graph vertex j refers
# to column j of the feature matrix produced by get_feature_matrix() for the
# dataset it was built from. Everything else (Laplacian, degree, edge lists)
# is derived from the adjacency matrix.

#' @keywords internal
#' @noRd
new_spatial_graph <- function(A, feature_ids, domain_type, geometry_id,
                              basis = NULL, weighted = FALSE,
                              parent_index = NULL) {
  if (!inherits(A, "sparseMatrix")) {
    A <- Matrix::Matrix(as.matrix(A), sparse = TRUE)
  }
  # Normalize to a numeric general column-sparse matrix (pattern matrices such
  # as the edgeless ngCMatrix returned for isolated voxels have no @x slot).
  A <- methods::as(methods::as(methods::as(A, "generalMatrix"), "CsparseMatrix"), "dMatrix")
  if (nrow(A) != ncol(A)) {
    stop("spatial_graph: adjacency matrix must be square.", call. = FALSE)
  }
  n <- nrow(A)
  feature_ids <- as.integer(feature_ids)
  if (length(feature_ids) != n) {
    stop(sprintf("spatial_graph: %d feature ids supplied for a %d-vertex graph.",
                 length(feature_ids), n), call. = FALSE)
  }
  if (!is.null(basis) && length(basis) != n) {
    stop("spatial_graph: 'basis' must have one entry per vertex.", call. = FALSE)
  }

  # Validate raw weights before any symmetrization can mask them.
  A <- Matrix::drop0(A)
  if (length(A@x) && any(A@x < 0)) {
    stop("spatial_graph: edge weights must be non-negative.", call. = FALSE)
  }
  # Symmetrize (an edge present in either direction is kept, with the larger
  # weight), drop self-loops, and binarize unless weights are meaningful.
  A <- Matrix::drop0(pmax(A, Matrix::t(A)))
  Matrix::diag(A) <- 0
  A <- Matrix::drop0(A)
  if (!isTRUE(weighted) && length(A@x)) A@x[] <- 1
  A <- methods::as(methods::as(A, "generalMatrix"), "CsparseMatrix")

  degree <- as.numeric(Matrix::rowSums(A))
  L <- Matrix::Diagonal(n = n, x = degree) - A

  structure(
    list(
      A = A,
      degree = degree,
      L = L,
      feature_ids = feature_ids,
      n_features = n,
      domain_type = domain_type,
      geometry_id = geometry_id,
      basis = basis,
      weighted = isTRUE(weighted),
      n_edges = as.integer(Matrix::nnzero(A) / 2),
      parent_index = parent_index
    ),
    class = c("spatial_graph", "list")
  )
}

#' @export
print.spatial_graph <- function(x, ...) {
  cat(sprintf("spatial_graph [%s]: %d features, %d edges%s\n",
              x$domain_type, x$n_features, x$n_edges,
              if (x$weighted) " (weighted)" else ""))
  if (!is.null(x$basis)) {
    cat(sprintf("  basis channels: %d\n", length(unique(x$basis))))
  }
  if (!is.null(x$parent_index)) {
    cat("  restricted from a parent graph\n")
  }
  cat("  geometry:", x$geometry_id, "\n")
  invisible(x)
}

#' @keywords internal
#' @noRd
.mask_active_vec <- function(mask) {
  vals <- if (inherits(mask, c("NeuroVol", "NeuroVec"))) as.numeric(neuroim2::values(mask)) else as.numeric(mask)
  vals[is.na(vals)] <- 0
  vals > 0
}

#' @keywords internal
#' @noRd
.volume_geometry_id <- function(mask, dims, neighbors) {
  spacing <- tryCatch(
    paste(signif(neuroim2::spacing(neuroim2::space(mask)), 6), collapse = "x"),
    error = function(e) "NA"
  )
  sprintf("volume:%s:spacing=%s:nbr=%d", paste(dims, collapse = "x"), spacing, as.integer(neighbors))
}

#' @keywords internal
#' @noRd
.volume_adjacency_from_mask <- function(mask, neighbors) {
  dims <- .infer_spatial_dims(mask, NULL)
  active <- .mask_active_vec(mask)
  if (length(active) != prod(dims)) {
    stop("spatial_graph: mask length does not match its spatial dimensions.", call. = FALSE)
  }
  A <- build_voxel_adjacency(active, dims = dims, neighbors = neighbors)
  list(A = A, feature_ids = which(active), dims = dims)
}

#' @rdname spatial_graph
#' @param neighbors Voxel neighbourhood for volumetric grids: 6 (faces),
#'   18 (faces and edges), or 26 (faces, edges, and corners).
#' @export
spatial_graph.mvpa_image_dataset <- function(x, neighbors = 6, ...) {
  vol <- .volume_adjacency_from_mask(x$mask, neighbors)
  new_spatial_graph(
    vol$A,
    feature_ids = vol$feature_ids,
    domain_type = "volume",
    geometry_id = .volume_geometry_id(x$mask, vol$dims, neighbors)
  )
}

#' @rdname spatial_graph
#' @param connect_basis Logical; for multibasis datasets, also connect each
#'   voxel to the same voxel in every other basis channel (default
#'   \code{FALSE}: channels form disconnected components, so smoothing never
#'   crosses channels).
#' @export
spatial_graph.mvpa_multibasis_image_dataset <- function(x, neighbors = 6,
                                                        connect_basis = FALSE, ...) {
  vol <- .volume_adjacency_from_mask(x$mask, neighbors)
  k <- as.integer(x$basis_count)
  n_vox <- length(vol$feature_ids)
  A <- Matrix::bdiag(replicate(k, vol$A, simplify = FALSE))
  if (isTRUE(connect_basis) && k > 1L) {
    pairs <- utils::combn(k, 2)
    ii <- unlist(lapply(seq_len(ncol(pairs)), function(p) (pairs[1, p] - 1L) * n_vox + seq_len(n_vox)))
    jj <- unlist(lapply(seq_len(ncol(pairs)), function(p) (pairs[2, p] - 1L) * n_vox + seq_len(n_vox)))
    cross <- Matrix::sparseMatrix(i = ii, j = jj, x = 1, dims = c(k * n_vox, k * n_vox))
    A <- A + cross + Matrix::t(cross)
  }
  new_spatial_graph(
    A,
    feature_ids = rep(vol$feature_ids, times = k),
    domain_type = "multibasis_volume",
    geometry_id = paste0(.volume_geometry_id(x$mask, vol$dims, neighbors), ":k=", k),
    basis = rep(seq_len(k), each = n_vox)
  )
}

#' @rdname spatial_graph
#' @export
spatial_graph.mvpa_surface_dataset <- function(x, ...) {
  require_package("neurosurf", "for surface adjacency graphs")
  geom <- x$train_data@geometry
  A_full <- neurosurf::adjacency(geom)
  idx <- which(as.numeric(x$mask) > 0)
  if (max(idx) > nrow(A_full)) {
    stop("spatial_graph: surface mask indexes nodes beyond the mesh.", call. = FALSE)
  }
  new_spatial_graph(
    A_full[idx, idx, drop = FALSE],
    feature_ids = idx,
    domain_type = "surface",
    geometry_id = sprintf("surface:nodes=%d:edges=%d", nrow(A_full),
                          as.integer(Matrix::nnzero(A_full) / 2))
  )
}

#' @rdname spatial_graph
#' @details
#' For clustered datasets the graph vertices are the parcels, in the column
#' order of the cluster time series, and \code{feature_ids} holds the actual
#' cluster labels (which need not be \code{1..K}). Column positions, not
#' labels, index the feature matrix.
#' @export
spatial_graph.mvpa_clustered_dataset <- function(x, neighbors = 6, ...) {
  cvol <- x$train_data@cvol
  vol <- .volume_adjacency_from_mask(cvol@mask, neighbors)
  clusters <- as.integer(cvol@clusters)
  if (length(clusters) != length(vol$feature_ids)) {
    stop("spatial_graph: cluster labels do not align with the clustered mask.", call. = FALSE)
  }
  cl_ids <- sort(unique(clusters))
  K <- ncol(x$train_data@ts)
  if (length(cl_ids) != K) {
    stop(sprintf("spatial_graph: %d cluster labels but %d cluster time series.",
                 length(cl_ids), K), call. = FALSE)
  }
  M <- Matrix::sparseMatrix(
    i = seq_along(clusters),
    j = match(clusters, cl_ids),
    x = 1,
    dims = c(length(clusters), K)
  )
  A_c <- Matrix::crossprod(M, vol$A %*% M)
  new_spatial_graph(
    A_c,
    feature_ids = cl_ids,
    domain_type = "cluster",
    geometry_id = paste0(.volume_geometry_id(cvol@mask, vol$dims, neighbors), ":clusters=", K)
  )
}

#' @rdname spatial_graph
#' @param feature_ids Integer identifiers mapping graph vertices to dataset
#'   locations when a raw adjacency is supplied. Defaults to \code{seq_len(n)}.
#' @param weighted Logical; keep edge weights of a raw adjacency instead of
#'   binarizing them.
#' @param domain_type Label for a raw adjacency's domain (default
#'   \code{"custom"}).
#' @export
spatial_graph.default <- function(x, feature_ids = NULL, weighted = FALSE,
                                  domain_type = "custom", ...) {
  A <- if (is.list(x) && !is.data.frame(x) && !is.null(x$A)) {
    # Fields carried by the list are defaults; explicit arguments win.
    if (is.null(feature_ids) && !is.null(x$feature_ids)) feature_ids <- x$feature_ids
    if (missing(weighted) && !is.null(x$weighted)) weighted <- isTRUE(x$weighted)
    x$A
  } else if (is.matrix(x) || inherits(x, "Matrix")) {
    x
  } else {
    stop(sprintf("spatial_graph: no method for class '%s'.", paste(class(x), collapse = "/")),
         call. = FALSE)
  }
  n <- nrow(A)
  if (is.null(feature_ids)) feature_ids <- seq_len(n)
  new_spatial_graph(
    A,
    feature_ids = feature_ids,
    domain_type = domain_type,
    geometry_id = sprintf("%s:n=%d", domain_type, n),
    weighted = weighted
  )
}

#' Restrict a Spatial Graph to a Subset of Features
#'
#' Returns the induced subgraph on the selected feature columns, keeping the
#' feature-to-column alignment. Use it after dropping invalid or constant
#' columns from the feature matrix so the graph continues to describe exactly
#' the columns the estimator sees.
#'
#' @param graph A \code{spatial_graph}.
#' @param keep Either a logical vector with one entry per feature column or an
#'   integer vector of column positions to keep (in the order they should
#'   appear in the restricted matrix).
#' @return A \code{spatial_graph} over the kept columns. Its
#'   \code{parent_index} field records the position of each kept column in
#'   the parent graph.
#' @examples
#' ds <- gen_sample_dataset(c(4, 4, 4), 8)
#' g <- spatial_graph(ds$dataset)
#' g2 <- restrict_graph(g, seq_len(10))
#' g2$n_features
#' @seealso \code{\link{spatial_graph}}
#' @export
restrict_graph <- function(graph, keep) {
  stopifnot(inherits(graph, "spatial_graph"))
  if (is.logical(keep)) {
    if (length(keep) != graph$n_features) {
      stop("restrict_graph: logical 'keep' must have one entry per feature.", call. = FALSE)
    }
    keep <- which(keep)
  }
  keep <- as.integer(keep)
  if (any(keep < 1L | keep > graph$n_features) || anyDuplicated(keep)) {
    stop("restrict_graph: 'keep' must be unique column positions within the graph.", call. = FALSE)
  }
  parent <- if (is.null(graph$parent_index)) keep else graph$parent_index[keep]
  new_spatial_graph(
    graph$A[keep, keep, drop = FALSE],
    feature_ids = graph$feature_ids[keep],
    domain_type = graph$domain_type,
    geometry_id = graph$geometry_id,
    basis = if (!is.null(graph$basis)) graph$basis[keep] else NULL,
    weighted = graph$weighted,
    parent_index = parent
  )
}

#' Edge List of a Spatial Graph
#'
#' @param graph A \code{spatial_graph}.
#' @return A data frame with integer columns \code{from} and \code{to}
#'   (column positions, \code{from < to}) and a numeric \code{weight}.
#' @examples
#' ds <- gen_sample_dataset(c(4, 4, 4), 8)
#' head(graph_edges(spatial_graph(ds$dataset)))
#' @seealso \code{\link{spatial_graph}}
#' @export
graph_edges <- function(graph) {
  stopifnot(inherits(graph, "spatial_graph"))
  s <- Matrix::summary(graph$A)
  s <- s[s$i < s$j, , drop = FALSE]
  data.frame(from = as.integer(s$i), to = as.integer(s$j), weight = as.numeric(s$x))
}

#' @keywords internal
#' @noRd
.assert_graph_aligned <- function(graph, n_columns, what = "feature matrix") {
  if (!inherits(graph, "spatial_graph")) {
    stop("A 'spatial_graph' is required.", call. = FALSE)
  }
  if (!identical(as.integer(graph$n_features), as.integer(n_columns))) {
    stop(sprintf("spatial_graph has %d vertices but the %s has %d columns.",
                 graph$n_features, what, as.integer(n_columns)), call. = FALSE)
  }
  invisible(TRUE)
}
