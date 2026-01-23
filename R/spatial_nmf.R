#' @keywords internal
#' @noRd
spatial_nmf_fit <- function(X,
                            k,
                            graph = NULL,
                            lambda = 0,
                            init = c("nndsvd", "random"),
                            W_init = NULL,
                            H_init = NULL,
                            max_iter = 200,
                            min_iter = 5,
                            tol = 1e-4,
                            check_every = 1,
                            eps = 1e-8,
                            normalize = c("none", "H"),
                            seed = NULL,
                            verbose = FALSE) {
  X <- as.matrix(X)
  if (any(!is.finite(X))) stop("X contains non-finite values.")
  if (any(X < 0)) stop("X must be non-negative.")
  if (!is.numeric(k) || length(k) != 1L || k < 1) stop("k must be a positive integer.")
  k <- as.integer(k)
  if (!is.numeric(lambda) || length(lambda) != 1L || !is.finite(lambda) || lambda < 0) {
    stop("lambda must be a finite non-negative scalar.")
  }
  if (!is.numeric(max_iter) || length(max_iter) != 1L || max_iter < 1) {
    stop("max_iter must be a positive integer.")
  }
  if (!is.numeric(min_iter) || length(min_iter) != 1L || min_iter < 0) {
    stop("min_iter must be a non-negative integer.")
  }
  if (!is.numeric(check_every) || length(check_every) != 1L || check_every < 1) {
    stop("check_every must be a positive integer.")
  }
  if (!is.numeric(tol) || length(tol) != 1L || !is.finite(tol) || tol <= 0) {
    stop("tol must be a positive scalar.")
  }
  if (!is.numeric(eps) || length(eps) != 1L || !is.finite(eps) || eps <= 0) {
    stop("eps must be a positive scalar.")
  }
  max_iter <- as.integer(max_iter)
  min_iter <- as.integer(min_iter)
  check_every <- as.integer(check_every)
  if (min_iter > max_iter) stop("min_iter must be <= max_iter.")

  init <- match.arg(init)
  normalize <- match.arg(normalize)

  if (!is.null(seed)) set.seed(seed)

  n <- nrow(X)
  p <- ncol(X)
  if (k > min(n, p)) stop("k must be <= min(nrow(X), ncol(X)).")

  graph_info <- .prepare_graph(graph, lambda, p)
  use_graph <- graph_info$use_graph

  init_res <- .nmf_initialize(X, k, init, W_init, H_init, eps)
  W <- init_res$W
  H <- init_res$H

  obj_hist <- numeric(0)
  iter_hist <- integer(0)
  converged <- FALSE

  for (iter in seq_len(max_iter)) {
    XHt <- tcrossprod(X, H)
    HHT <- tcrossprod(H)
    WHHT <- W %*% HHT
    W <- W * (XHt / pmax(WHHT, eps))

    WtX <- crossprod(W, X)
    WtW <- crossprod(W)
    denom <- WtW %*% H
    numer <- WtX

    if (use_graph) {
      HA <- as.matrix(H %*% graph_info$A)
      numer <- numer + lambda * HA
      denom <- denom + lambda * (H * graph_info$degree)
    }

    H <- H * (numer / pmax(denom, eps))

    if (normalize == "H") {
      scale <- pmax(rowSums(H), eps)
      H <- H / scale
      W <- W * scale
    }

    if (iter %% check_every == 0 || iter == 1 || iter == max_iter) {
      obj <- .nmf_objective(X, W, H, graph_info, lambda, eps)
      obj_hist <- c(obj_hist, obj)
      iter_hist <- c(iter_hist, iter)
      if (length(obj_hist) > 1L) {
        rel_change <- abs(obj_hist[length(obj_hist)] - obj_hist[length(obj_hist) - 1L]) /
          (obj_hist[length(obj_hist) - 1L] + eps)
        if (verbose) {
          message(sprintf("iter %d | obj %.6g | rel_change %.3g", iter, obj, rel_change))
        }
        if (iter >= min_iter && rel_change < tol) {
          converged <- TRUE
          break
        }
      } else if (verbose) {
        message(sprintf("iter %d | obj %.6g", iter, obj))
      }
    }
  }

  structure(
    list(
      W = W,
      H = H,
      k = k,
      n = n,
      p = p,
      lambda = lambda,
      graph = graph_info$meta,
      objective = obj_hist,
      objective_iter = iter_hist,
      iter = iter,
      converged = converged
    ),
    class = "spatial_nmf_fit"
  )
}

#' @keywords internal
#' @noRd
spatial_nmf_project <- function(X,
                                H,
                                W_init = NULL,
                                max_iter = 200,
                                min_iter = 5,
                                tol = 1e-4,
                                check_every = 1,
                                eps = 1e-8,
                                seed = NULL,
                                verbose = FALSE) {
  X <- as.matrix(X)
  if (any(!is.finite(X))) stop("X contains non-finite values.")
  if (any(X < 0)) stop("X must be non-negative.")
  H <- as.matrix(H)
  if (any(!is.finite(H))) stop("H contains non-finite values.")
  if (any(H < 0)) stop("H must be non-negative.")
  if (ncol(X) != ncol(H)) stop("ncol(X) must match ncol(H).")
  if (!is.numeric(max_iter) || length(max_iter) != 1L || max_iter < 1) {
    stop("max_iter must be a positive integer.")
  }
  if (!is.numeric(min_iter) || length(min_iter) != 1L || min_iter < 0) {
    stop("min_iter must be a non-negative integer.")
  }
  if (!is.numeric(check_every) || length(check_every) != 1L || check_every < 1) {
    stop("check_every must be a positive integer.")
  }
  if (!is.numeric(tol) || length(tol) != 1L || !is.finite(tol) || tol <= 0) {
    stop("tol must be a positive scalar.")
  }
  if (!is.numeric(eps) || length(eps) != 1L || !is.finite(eps) || eps <= 0) {
    stop("eps must be a positive scalar.")
  }
  max_iter <- as.integer(max_iter)
  min_iter <- as.integer(min_iter)
  check_every <- as.integer(check_every)
  if (min_iter > max_iter) stop("min_iter must be <= max_iter.")

  n <- nrow(X)
  k <- nrow(H)

  if (!is.null(seed)) set.seed(seed)

  if (is.null(W_init)) {
    W <- matrix(runif(n * k), nrow = n, ncol = k)
  } else {
    W <- as.matrix(W_init)
    if (!all(dim(W) == c(n, k))) stop("W_init has incompatible dimensions.")
  }

  obj_hist <- numeric(0)
  iter_hist <- integer(0)
  converged <- FALSE

  HHT <- tcrossprod(H)
  XHt <- tcrossprod(X, H)

  for (iter in seq_len(max_iter)) {
    WHHT <- W %*% HHT
    W <- W * (XHt / pmax(WHHT, eps))

    if (iter %% check_every == 0 || iter == 1 || iter == max_iter) {
      obj <- sum((X - W %*% H)^2)
      obj_hist <- c(obj_hist, obj)
      iter_hist <- c(iter_hist, iter)
      if (length(obj_hist) > 1L) {
        rel_change <- abs(obj_hist[length(obj_hist)] - obj_hist[length(obj_hist) - 1L]) /
          (obj_hist[length(obj_hist) - 1L] + eps)
        if (verbose) {
          message(sprintf("iter %d | obj %.6g | rel_change %.3g", iter, obj, rel_change))
        }
        if (iter >= min_iter && rel_change < tol) {
          converged <- TRUE
          break
        }
      } else if (verbose) {
        message(sprintf("iter %d | obj %.6g", iter, obj))
      }
    }
  }

  structure(
    list(
      W = W,
      H = H,
      n = n,
      k = k,
      p = ncol(H),
      objective = obj_hist,
      objective_iter = iter_hist,
      iter = iter,
      converged = converged
    ),
    class = "spatial_nmf_projection"
  )
}

.nmf_initialize <- function(X, k, init, W_init, H_init, eps) {
  if (!is.null(W_init) && !is.null(H_init)) {
    W <- as.matrix(W_init)
    H <- as.matrix(H_init)
    if (nrow(W) != nrow(X) || ncol(W) != k) stop("W_init has incompatible dimensions.")
    if (nrow(H) != k || ncol(H) != ncol(X)) stop("H_init has incompatible dimensions.")
    if (any(W < 0) || any(H < 0)) stop("W_init and H_init must be non-negative.")
    return(list(W = W, H = H))
  }

  if (init == "nndsvd") {
    return(.nndsvd_initialize(X, k, eps))
  }

  n <- nrow(X)
  p <- ncol(X)
  scale <- max(mean(X), eps)
  W <- matrix(runif(n * k), nrow = n, ncol = k) * scale
  H <- matrix(runif(k * p), nrow = k, ncol = p) * scale
  list(W = W, H = H)
}

.nndsvd_initialize <- function(X, k, eps) {
  n <- nrow(X)
  p <- ncol(X)
  min_dim <- min(n, p)

  svd_res <- NULL
  if (k >= min_dim || min_dim <= 5L || (k / min_dim) > 0.5) {
    svd_res <- base::svd(X)
  } else {
    svd_res <- tryCatch(
      irlba::irlba(X, nv = k, nu = k),
      error = function(e) base::svd(X)
    )
  }

  U <- as.matrix(svd_res$u[, seq_len(k), drop = FALSE])
  V <- as.matrix(svd_res$v[, seq_len(k), drop = FALSE])
  d <- svd_res$d[seq_len(k)]

  W <- matrix(0, nrow = n, ncol = k)
  H <- matrix(0, nrow = k, ncol = p)

  W[, 1] <- sqrt(d[1]) * pmax(U[, 1], 0)
  H[1, ] <- sqrt(d[1]) * pmax(V[, 1], 0)

  if (k > 1L) {
    for (i in 2:k) {
      u <- U[, i]
      v <- V[, i]
      u_pos <- pmax(u, 0)
      u_neg <- pmax(-u, 0)
      v_pos <- pmax(v, 0)
      v_neg <- pmax(-v, 0)

      norm_pos <- sqrt(sum(u_pos^2)) * sqrt(sum(v_pos^2))
      norm_neg <- sqrt(sum(u_neg^2)) * sqrt(sum(v_neg^2))

      if (norm_pos >= norm_neg && norm_pos > 0) {
        W[, i] <- sqrt(d[i] * norm_pos) * (u_pos / (sqrt(sum(u_pos^2)) + eps))
        H[i, ] <- sqrt(d[i] * norm_pos) * (v_pos / (sqrt(sum(v_pos^2)) + eps))
      } else if (norm_neg > 0) {
        W[, i] <- sqrt(d[i] * norm_neg) * (u_neg / (sqrt(sum(u_neg^2)) + eps))
        H[i, ] <- sqrt(d[i] * norm_neg) * (v_neg / (sqrt(sum(v_neg^2)) + eps))
      } else {
        W[, i] <- runif(n) * eps
        H[i, ] <- runif(p) * eps
      }
    }
  }

  list(W = W, H = H)
}

.nmf_objective <- function(X, W, H, graph_info, lambda, eps) {
  resid <- X - W %*% H
  obj <- sum(resid * resid)
  if (graph_info$use_graph) {
    HA <- as.matrix(H %*% graph_info$A)
    smooth <- sum(H * (H * graph_info$degree - HA))
    obj <- obj + lambda * smooth
  }
  obj + eps
}

.prepare_graph <- function(graph, lambda, p) {
  if (is.null(graph) || lambda <= 0) {
    return(list(use_graph = FALSE, A = NULL, degree = NULL, meta = NULL))
  }

  if (!is.list(graph) || is.null(graph$A)) {
    stop("graph must be a list with at least an adjacency matrix 'A'.")
  }

  A <- graph$A
  if (!inherits(A, "sparseMatrix")) {
    A <- Matrix::Matrix(A, sparse = TRUE)
  }
  if (nrow(A) != p || ncol(A) != p) {
    stop("graph$A must be a square sparse matrix matching ncol(X).")
  }
  A <- (A + Matrix::t(A)) * 0.5
  Matrix::diag(A) <- 0
  weighted <- isTRUE(graph$weighted)
  if (!weighted) {
    if (!is.null(A@x) && any(A@x != 0 & abs(A@x - 1) > 1e-12)) {
      warning("graph$A has non-binary weights; binarizing. Set graph$weighted=TRUE to preserve weights.")
    }
    if (!is.null(A@x)) A@x[A@x != 0] <- 1
  }

  degree <- graph$degree
  if (is.null(degree) || length(degree) != p) {
    if (!is.null(degree) && length(degree) != p) {
      warning("graph$degree length does not match ncol(X); recomputing degree from A.")
    }
    degree <- Matrix::rowSums(A)
  }

  meta <- list(n_edges = Matrix::nnzero(A) / 2L, weighted = weighted)

  list(use_graph = TRUE, A = A, degree = degree, meta = meta)
}

#' @keywords internal
#' @noRd
build_voxel_adjacency <- function(mask, dims = NULL, neighbors = 6) {
  dims <- .infer_spatial_dims(mask, dims)
  if (length(dims) < 2L || length(dims) > 3L) {
    stop("dims must have length 2 or 3.")
  }

  mask_vec <- .resolve_mask_vec(mask, prod(dims))
  mask_idx <- which(mask_vec)
  if (length(mask_idx) == 0L) {
    return(Matrix::sparseMatrix(i = integer(0), j = integer(0), dims = c(0, 0)))
  }

  coords <- arrayInd(mask_idx, dims)
  n <- nrow(coords)

  index_map <- integer(prod(dims))
  index_map[mask_idx] <- seq_len(n)

  active_dims <- which(dims > 1L)
  if (length(active_dims) == 2L) {
    if (neighbors %in% c(6, 18, 26)) {
      neighbors <- if (neighbors == 6) 4 else 8
    }
    offsets2 <- .neighbor_offsets(neighbors, 2L)
    offsets <- matrix(0L, nrow = nrow(offsets2), ncol = length(dims))
    offsets[, active_dims] <- offsets2
  } else {
    offsets <- .neighbor_offsets(neighbors, length(dims))
  }
  edges_i <- vector("list", nrow(offsets))
  edges_j <- vector("list", nrow(offsets))

  for (i in seq_len(nrow(offsets))) {
    off <- offsets[i, ]
    neigh_coords <- sweep(coords, 2, off, "+")

    in_bounds <- rep(TRUE, n)
    for (d in seq_along(dims)) {
      in_bounds <- in_bounds & neigh_coords[, d] >= 1L & neigh_coords[, d] <= dims[d]
    }

    if (!any(in_bounds)) next

    base_rows <- which(in_bounds)
    neigh_lin <- .linear_index(neigh_coords[base_rows, , drop = FALSE], dims)
    neigh_row <- index_map[neigh_lin]
    ok <- neigh_row > 0L
    if (!any(ok)) next

    edges_i[[i]] <- base_rows[ok]
    edges_j[[i]] <- neigh_row[ok]
  }

  ii <- unlist(edges_i, use.names = FALSE)
  jj <- unlist(edges_j, use.names = FALSE)
  if (!length(ii)) {
    return(Matrix::sparseMatrix(i = integer(0), j = integer(0), dims = c(n, n)))
  }

  A <- Matrix::sparseMatrix(i = ii, j = jj, x = 1, dims = c(n, n))
  A <- A + Matrix::t(A)
  Matrix::diag(A) <- 0
  if (!is.null(A@x)) A@x[A@x != 0] <- 1
  A
}

#' @keywords internal
#' @noRd
build_graph_laplacian <- function(A, normalized = FALSE) {
  if (!inherits(A, "sparseMatrix")) {
    A <- Matrix::Matrix(A, sparse = TRUE)
  }
  A <- (A + Matrix::t(A)) * 0.5
  Matrix::diag(A) <- 0
  if (!is.null(A@x)) A@x[A@x != 0] <- 1

  degree <- Matrix::rowSums(A)
  if (normalized) {
    deg_inv_sqrt <- 1 / sqrt(pmax(degree, .Machine$double.eps))
    D_inv_sqrt <- Matrix::Diagonal(x = deg_inv_sqrt)
    L <- Matrix::Diagonal(n = nrow(A)) - (D_inv_sqrt %*% A %*% D_inv_sqrt)
  } else {
    L <- Matrix::Diagonal(x = degree) - A
  }

  list(A = A, degree = degree, L = L, normalized = normalized)
}

.infer_spatial_dims <- function(mask, dims) {
  if (!is.null(dims)) {
    return(as.integer(dims))
  }

  if (inherits(mask, c("NeuroVol", "NeuroVec"))) {
    space_obj <- neuroim2::space(mask)
    dims <- dim(space_obj)
    if (length(dims) > 3L) dims <- dims[seq_len(3L)]
    return(as.integer(dims))
  }

  stop("dims must be provided when mask is not a NeuroVol/NeuroVec.")
}

.resolve_mask_vec <- function(mask, spatial_length) {
  if (inherits(mask, c("NeuroVol", "NeuroVec"))) {
    vals <- neuroim2::values(mask)
    if (is.matrix(vals)) {
      vals <- vals[, 1, drop = TRUE]
    }
    return(as.logical(as.numeric(vals)))
  }
  if (is.numeric(mask) || is.logical(mask)) {
    mask_vec <- as.logical(mask)
    if (!is.null(spatial_length) && length(mask_vec) != spatial_length) {
      stop("mask length does not match spatial dimensions.")
    }
    return(mask_vec)
  }

  stop("mask must be a NeuroVol/NeuroVec or a logical/numeric vector.")
}

.neighbor_offsets <- function(neighbors, ndim) {
  if (!is.numeric(neighbors) || length(neighbors) != 1L) {
    stop("neighbors must be a single numeric value.")
  }

  if (ndim == 2L) {
    if (!neighbors %in% c(4, 8)) stop("neighbors must be 4 or 8 for 2D.")
  } else if (ndim == 3L) {
    if (!neighbors %in% c(6, 18, 26)) stop("neighbors must be 6, 18, or 26 for 3D.")
  } else {
    stop("ndim must be 2 or 3.")
  }

  offsets <- as.matrix(expand.grid(rep(list(-1:1), ndim)))
  offsets <- offsets[rowSums(abs(offsets)) > 0, , drop = FALSE]

  if (ndim == 2L) {
    if (neighbors == 4) {
      offsets <- offsets[rowSums(abs(offsets)) == 1, , drop = FALSE]
    }
  } else {
    if (neighbors == 6) {
      offsets <- offsets[rowSums(abs(offsets)) == 1, , drop = FALSE]
    } else if (neighbors == 18) {
      offsets <- offsets[rowSums(abs(offsets)) <= 2, , drop = FALSE]
    }
  }

  offsets
}

.linear_index <- function(coords, dims) {
  if (length(dims) == 2L) {
    coords[, 1] + (coords[, 2] - 1L) * dims[1]
  } else {
    coords[, 1] +
      (coords[, 2] - 1L) * dims[1] +
      (coords[, 3] - 1L) * dims[1] * dims[2]
  }
}

#' Spatial NMF on Map Lists
#'
#' \strong{Raison d'etre.}
#' Spatial NMF factorizes a non-negative subject-by-voxel matrix \eqn{X} into
#' \eqn{W} (subject loadings) and \eqn{H} (spatial components) such that
#' \eqn{X \approx W H}. \eqn{W} captures how strongly each subject expresses each
#' component, while \eqn{H} encodes spatially interpretable, additive patterns.
#' An optional graph Laplacian penalty encourages smooth, neuroanatomically
#' plausible components. The resulting representation is compact and supports
#' downstream inference (component tests, global CV tests, and stability).
#'
#' \strong{Why not just voxelwise t-tests?}
#' Spatial NMF targets multivariate structure by learning \emph{patterns of
#' covarying voxels} rather than testing each voxel independently. This yields a
#' low-dimensional representation that is often more interpretable, statistically
#' powerful (fewer multiple comparisons), and robust to noise. With optional
#' spatial regularization, components are anatomically smoother than voxelwise
#' maps, and the framework provides component-level inference and stability
#' analyses that reflect network-level effects.
#'
#' \strong{Interpreting inference outputs.}
#' \itemize{
#'   \item \emph{Component tests} (`spatial_nmf_component_test`): use `p_fwer` to
#'   identify components whose loadings differ reliably between groups; `p_unc`
#'   is uncorrected. `p_global` tests whether any component differs.
#'   \item \emph{Global test} (`spatial_nmf_global_test`): `stat` is cross-validated
#'   performance (AUC or accuracy) on component loadings; a significant `p_value`
#'   indicates a multivariate group difference without pinpointing a component.
#'   \item \emph{Stability} (`spatial_nmf_stability`): interpret `selection` as
#'   voxelwise stability (frequency among top-loadings), and `cv` as variability;
#'   higher `component_similarity` indicates stable component identity.
#' }
#'
#' Convenience wrapper that converts lists of NeuroVol/NeuroSurface maps into
#' a subject-by-voxel matrix and fits spatially regularized NMF.
#'
#' @param group_A List of NeuroVol/NeuroSurface maps for group A. All maps must share
#'   the same spatial grid/geometry; if maps carry sparse indices, those indices must match.
#' @param group_B Optional list of NeuroVol/NeuroSurface maps for group B (same requirements as group_A).
#' @param mask Mask object for volumetric data (NeuroVol or logical/numeric vector).
#'   Required for volumetric inputs; optional for surface inputs.
#' @param dims Optional spatial dimensions for volumetric masks given as vectors.
#' @param k Number of components.
#' @param lambda Spatial regularization strength (0 = none).
#' @param graph Optional graph Laplacian list (from build_graph_laplacian). If `graph$A`
#'   is weighted, set `graph$weighted=TRUE` to preserve weights; otherwise edges are binarized.
#' @param neighbors Neighborhood size for volumetric adjacency (6/18/26).
#' @param na_action How to handle NA values in maps: "zero" or "error".
#' @param return_maps Logical; return component maps as NeuroVol/NeuroSurface.
#' @param return_data Logical; include the data matrix in the result.
#' @param component_test NULL to skip, TRUE for defaults, or a list of arguments
#'   passed to spatial_nmf_component_test.
#' @param global_test NULL to skip, TRUE for defaults, or a list of arguments
#'   passed to spatial_nmf_global_test.
#' @param stability NULL to skip, TRUE for defaults, or a list of arguments
#'   passed to spatial_nmf_stability.
#' @param ... Additional arguments passed to spatial_nmf_fit.
#'
#' @return A list with fields:
#'   \itemize{
#'     \item fit: spatial_nmf_fit object.
#'     \item components: list of spatial component maps (if return_maps).
#'     \item groups: factor of group labels for each subject.
#'     \item mask_indices: indices used to vectorize maps.
#'     \item map_type: "volume" or "surface".
#'     \item data: subject-by-voxel matrix (if return_data).
#'   }
#' @export
spatial_nmf_maps <- function(group_A,
                             group_B = NULL,
                             mask = NULL,
                             dims = NULL,
                             k,
                             lambda = 0,
                             graph = NULL,
                             neighbors = 6,
                             na_action = c("zero", "error"),
                             return_maps = TRUE,
                             return_data = FALSE,
                             component_test = NULL,
                             global_test = NULL,
                             stability = NULL,
                             ...) {
  na_action <- match.arg(na_action)

  group_A <- .ensure_map_list(group_A, "group_A")
  group_B <- .ensure_map_list(group_B, "group_B")

  map_type <- .infer_map_type(group_A[[1]])
  .validate_map_list(group_A, map_type, "group_A")
  if (!is.null(group_B)) {
    .validate_map_list(group_B, map_type, "group_B")
  }
  .validate_map_indices(group_A, "group_A")
  if (!is.null(group_B)) {
    .validate_map_indices(group_B, "group_B")
  }

  if (map_type == "volume") {
    if (is.null(mask)) {
      stop("mask is required for volumetric inputs.")
    }
    dims <- .infer_spatial_dims(mask, dims)
    mask_vec <- .resolve_mask_vec(mask, prod(dims))
    mask_idx <- which(mask_vec)
    if (!length(mask_idx)) stop("mask contains no active voxels.")
    .validate_volume_dims(group_A, dims)
    if (!is.null(group_B)) .validate_volume_dims(group_B, dims)
    full_length <- length(mask_vec)
  } else {
    full_length <- .infer_surface_length(group_A[[1]], mask)
    mask_vec <- .resolve_surface_mask(mask, full_length)
    mask_idx <- which(mask_vec)
    if (!length(mask_idx)) stop("mask contains no active vertices.")
    .validate_surface_length(group_A, full_length)
    if (!is.null(group_B)) .validate_surface_length(group_B, full_length)
  }

  X_A <- .maps_to_matrix(group_A, full_length, mask_idx, na_action)
  X <- X_A
  groups <- rep("A", nrow(X_A))
  if (!is.null(group_B)) {
    X_B <- .maps_to_matrix(group_B, full_length, mask_idx, na_action)
    X <- rbind(X_A, X_B)
    groups <- c(groups, rep("B", nrow(X_B)))
  }
  groups <- factor(groups)

  if (lambda > 0 && is.null(graph)) {
    if (map_type == "volume") {
      if (length(dims) == 2L && neighbors == 6) {
        neighbors <- 4
      }
      A <- build_voxel_adjacency(mask, dims = dims, neighbors = neighbors)
      graph <- build_graph_laplacian(A)
    } else {
      stop("graph is required for surface inputs when lambda > 0.")
    }
  }

  fit <- spatial_nmf_fit(
    X = X,
    k = k,
    graph = graph,
    lambda = lambda,
    ...
  )

  components <- NULL
  if (isTRUE(return_maps)) {
    components <- .components_to_maps(
      fit$H,
      mask_idx = mask_idx,
      map_type = map_type,
      mask = mask,
      dims = dims,
      full_length = full_length,
      ref_map = group_A[[1]]
    )
  }

  component_res <- NULL
  component_args <- .as_arg_list(component_test, "component_test")
  if (!is.null(component_args)) {
    if (is.null(component_args$fit) && is.null(component_args$W)) {
      component_args$fit <- fit
    }
    if (is.null(component_args$groups) && !is.null(group_B)) {
      component_args$groups <- groups
    }
    if (is.null(component_args$test)) {
      component_args$test <- if (!is.null(group_B)) "two_group" else "one_group"
    }
    if (component_args$test == "one_group" && is.null(component_args$null_W)) {
      stop("component_test requires null_W for one_group inference.")
    }
    component_res <- do.call(spatial_nmf_component_test, component_args)
  }

  global_res <- NULL
  global_args <- .as_arg_list(global_test, "global_test")
  if (!is.null(global_args)) {
    if (is.null(group_B)) stop("global_test requires group_B.")
    if (is.null(global_args$X)) global_args$X <- X
    if (is.null(global_args$groups)) global_args$groups <- groups
    if (is.null(global_args$k)) global_args$k <- k
    if (is.null(global_args$lambda)) global_args$lambda <- lambda
    if (is.null(global_args$graph)) global_args$graph <- graph
    if (map_type == "volume" && is.null(global_args$neighbors)) {
      global_args$neighbors <- neighbors
    }
    global_res <- do.call(spatial_nmf_global_test, global_args)
  }

  stability_res <- NULL
  stability_args <- .as_arg_list(stability, "stability")
  if (!is.null(stability_args)) {
    if (is.null(stability_args$fit)) stability_args$fit <- fit
    if (is.null(stability_args$X)) stability_args$X <- X
    if (is.null(stability_args$graph)) stability_args$graph <- graph
    if (is.null(stability_args$lambda)) stability_args$lambda <- lambda
    stability_res <- do.call(spatial_nmf_stability, stability_args)
  }

  res <- list(
    fit = fit,
    components = components,
    groups = groups,
    mask_indices = mask_idx,
    map_type = map_type,
    mask = mask,
    dims = dims,
    full_length = full_length,
    ref_map = group_A[[1]]
  )
  if (return_data) res$data <- X
  if (!is.null(component_res)) res$component_test <- component_res
  if (!is.null(global_res)) res$global_test <- global_res
  if (!is.null(stability_res)) res$stability <- stability_res

  structure(res, class = "spatial_nmf_maps_result")
}

.ensure_map_list <- function(x, name) {
  if (is.null(x)) return(NULL)
  if (inherits(x, c("NeuroVol", "NeuroVec", "NeuroSurface", "NeuroSurfaceVector"))) {
    return(list(x))
  }
  if (!is.list(x) || length(x) == 0L) {
    stop(sprintf("%s must be a non-empty list of map objects.", name))
  }
  x
}

.as_arg_list <- function(x, name) {
  if (is.null(x) || identical(x, FALSE)) return(NULL)
  if (isTRUE(x)) return(list())
  if (!is.list(x)) {
    stop(sprintf("%s must be NULL, TRUE/FALSE, or a list of arguments.", name))
  }
  x
}

.infer_map_type <- function(x) {
  if (inherits(x, c("NeuroVol", "NeuroVec"))) return("volume")
  if (inherits(x, c("NeuroSurface", "NeuroSurfaceVector"))) return("surface")
  stop("Maps must inherit from NeuroVol/NeuroVec or NeuroSurface/NeuroSurfaceVector.")
}

.validate_map_list <- function(lst, map_type, name) {
  if (is.null(lst)) return(invisible(NULL))
  ok <- vapply(lst, function(x) {
    if (map_type == "volume") {
      inherits(x, c("NeuroVol", "NeuroVec"))
    } else {
      inherits(x, c("NeuroSurface", "NeuroSurfaceVector"))
    }
  }, logical(1))
  if (!all(ok)) {
    stop(sprintf("%s contains maps with incompatible types.", name))
  }
  invisible(NULL)
}

.validate_map_indices <- function(lst, name) {
  if (is.null(lst)) return(invisible(NULL))
  idx_list <- lapply(lst, function(x) {
    tryCatch(neuroim2::indices(x), error = function(e) NULL)
  })
  has_idx <- vapply(idx_list, function(idx) !is.null(idx), logical(1))
  if (any(has_idx) && !all(has_idx)) {
    stop(sprintf("%s maps must either all have indices or none do.", name))
  }
  if (any(has_idx)) {
    ref <- as.integer(idx_list[[which(has_idx)[1L]]])
    for (i in seq_along(idx_list)) {
      idx <- as.integer(idx_list[[i]])
      if (length(idx) != length(ref) || any(idx != ref)) {
        stop(sprintf("%s maps have inconsistent indices; resample to a common grid.", name))
      }
    }
  }
  invisible(NULL)
}

.validate_volume_dims <- function(lst, dims) {
  for (x in lst) {
    space_obj <- tryCatch(neuroim2::space(x), error = function(e) NULL)
    if (is.null(space_obj)) next
    xdims <- dim(space_obj)
    if (length(xdims) > 3L) xdims <- xdims[seq_len(3L)]
    if (!all(as.integer(xdims) == as.integer(dims))) {
      stop("Map spatial dimensions do not match mask/dims.")
    }
  }
  invisible(NULL)
}

.infer_surface_length <- function(ref_map, mask) {
  if (!is.null(mask)) {
    if (inherits(mask, "NeuroSurface")) {
      geom <- neurosurf::geometry(mask)
      return(length(neurosurf::nodes(geom)))
    }
    if (is.numeric(mask) || is.logical(mask)) {
      return(length(mask))
    }
  }

  geom <- neurosurf::geometry(ref_map)
  length(neurosurf::nodes(geom))
}

.resolve_surface_mask <- function(mask, full_length) {
  if (is.null(mask)) {
    return(rep(TRUE, full_length))
  }
  if (inherits(mask, "NeuroSurface")) {
    vals <- as.numeric(neuroim2::values(mask))
    if (length(vals) != full_length) {
      stop("Surface mask length does not match surface geometry.")
    }
    return(as.logical(vals))
  }
  if (is.numeric(mask) || is.logical(mask)) {
    if (length(mask) != full_length) {
      stop("Surface mask length does not match surface geometry.")
    }
    return(as.logical(mask))
  }
  stop("mask must be NULL, a NeuroSurface, or a logical/numeric vector for surfaces.")
}

.validate_surface_length <- function(lst, full_length) {
  for (x in lst) {
    geom <- neurosurf::geometry(x)
    if (length(neurosurf::nodes(geom)) != full_length) {
      stop("Surface map geometry length does not match mask geometry.")
    }
  }
  invisible(NULL)
}

.maps_to_matrix <- function(lst, full_length, mask_idx, na_action) {
  mats <- lapply(lst, function(map) {
    vec <- .map_to_masked_vector(map, full_length, mask_idx)
    if (anyNA(vec)) {
      if (na_action == "error") {
        stop("Map contains NA values.")
      }
      vec[is.na(vec)] <- 0
    }
    vec
  })
  do.call(rbind, mats)
}

.map_to_masked_vector <- function(map, full_length, mask_idx) {
  vals <- neuroim2::values(map)
  if (is.matrix(vals)) {
    if (inherits(map, c("NeuroVec", "NeuroSurfaceVector")) && ncol(vals) > 1L) {
      stop("Each map must contain a single volume of data.")
    }
    if (ncol(vals) == 1L) {
      vals <- vals[, 1]
    }
    vals <- as.numeric(vals)
  } else if (is.array(vals)) {
    vals <- as.numeric(vals)
  } else {
    vals <- as.numeric(vals)
  }

  idx <- tryCatch(neuroim2::indices(map), error = function(e) NULL)
  if (!is.null(idx) && length(idx) == length(vals)) {
    full <- numeric(full_length)
    full[idx] <- vals
    return(full[mask_idx])
  }

  if (length(vals) == full_length) {
    return(vals[mask_idx])
  }

  if (length(vals) == length(mask_idx)) {
    return(vals)
  }

  stop("Map length does not match mask or indices.")
}

.components_to_maps <- function(H, mask_idx, map_type, mask, dims, full_length, ref_map) {
  k <- nrow(H)
  comps <- vector("list", k)
  names(comps) <- paste0("comp", seq_len(k))

  if (map_type == "volume") {
    if (inherits(mask, c("NeuroVol", "NeuroVec"))) {
      space_obj <- spatial_only_space(neuroim2::space(mask))
    } else {
      space_obj <- neuroim2::NeuroSpace(dims)
    }
    for (i in seq_len(k)) {
      comps[[i]] <- neuroim2::NeuroVol(
        data = as.numeric(H[i, ]),
        space = space_obj,
        indices = mask_idx
      )
    }
  } else {
    geom <- if (inherits(mask, "NeuroSurface")) {
      neurosurf::geometry(mask)
    } else {
      neurosurf::geometry(ref_map)
    }
    for (i in seq_len(k)) {
      comps[[i]] <- neurosurf::NeuroSurface(
        geometry = geom,
        indices = mask_idx,
        data = as.numeric(H[i, ])
      )
    }
  }

  comps
}
