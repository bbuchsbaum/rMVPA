#' Create an MVPA Dataset for Clustered/Parcellated Data
#'
#' Creates a dataset object for MVPA analysis on parcellated brain data stored
#' as \code{ClusteredNeuroVec} objects from \pkg{neuroim2}. A "clustered searchlight"
#' iterates over all K parcels, pooling nearby parcels as features for each analysis.
#'
#' @param train_data A \code{ClusteredNeuroVec} instance (training data).
#' @param test_data  An optional \code{ClusteredNeuroVec} instance (test data). Default: NULL.
#' @param mask       An optional mask. If NULL, derived from \code{train_data}'s cluster volume.
#'
#' @return An \code{mvpa_clustered_dataset} object (S3 class) containing:
#'   \describe{
#'     \item{train_data}{The training data as a \code{ClusteredNeuroVec}.}
#'     \item{test_data}{The test data (if provided).}
#'     \item{mask}{The mask.}
#'     \item{cvol}{The \code{ClusteredNeuroVol} defining cluster assignments.}
#'     \item{region_mask}{A \code{NeuroVol} where each voxel's value is its cluster ID
#'       (for broadcasting results via \code{map_values}).}
#'     \item{has_test_set}{Logical flag.}
#'   }
#'
#' @examples
#' \dontrun{
#'   ds <- gen_clustered_sample_dataset(K=5, nobs=20)
#'   print(ds$dataset)
#' }
#' @importFrom assertthat assert_that
#' @importFrom neuroim2 NeuroVol space num_clusters
#' @export
mvpa_clustered_dataset <- function(train_data, test_data = NULL, mask = NULL) {
  assert_that(inherits(train_data, "ClusteredNeuroVec"),
              msg = "train_data must be a ClusteredNeuroVec")

  if (!is.null(test_data)) {
    assert_that(inherits(test_data, "ClusteredNeuroVec"),
                msg = "test_data must be a ClusteredNeuroVec")
  }

  cvol <- train_data@cvol

  # Build the region_mask: a NeuroVol where each voxel stores its cluster ID.
  # This is used by build_output_map to broadcast per-cluster results to voxels.
  sp3 <- neuroim2::space(cvol@mask)
  cluster_ids <- cvol@clusters
  mask_idx <- which(cvol@mask > 0)
  region_arr <- integer(prod(dim(sp3)))
  region_arr[mask_idx] <- cluster_ids
  region_mask <- neuroim2::NeuroVol(array(region_arr, dim(sp3)), sp3)

  if (is.null(mask)) {
    mask <- region_mask
  }

  has_test <- !is.null(test_data)

  structure(
    list(
      train_data   = train_data,
      test_data    = test_data,
      mask         = mask,
      cvol         = cvol,
      region_mask  = region_mask,
      has_test_set = has_test
    ),
    class = c("mvpa_clustered_dataset", "mvpa_dataset", "list")
  )
}


#' @export
nobs.mvpa_clustered_dataset <- function(x) {
  nrow(x$train_data@ts)
}

#' @export
#' @method print mvpa_clustered_dataset
print.mvpa_clustered_dataset <- function(x, ...) {
  K <- neuroim2::num_clusters(x$cvol)
  T_obs <- nrow(x$train_data@ts)
  cat("\n  Clustered MVPA Dataset\n\n")
  cat("  - Clusters:     ", K, "\n")
  cat("  - Observations: ", T_obs, "\n")
  cat("  - Test set:     ", if (x$has_test_set) "yes" else "no", "\n\n")
}

#' @export
has_test_set.mvpa_clustered_dataset <- function(obj) {
  isTRUE(obj$has_test_set)
}


# ---- Searchlight generics ----

#' @export
get_center_ids.mvpa_clustered_dataset <- function(dataset, ...) {
  seq_len(neuroim2::num_clusters(dataset$cvol))
}

#' @export
build_output_map.mvpa_clustered_dataset <- function(dataset, metric_vector, ids, ...) {
  lookup <- cbind(as.numeric(ids), as.numeric(metric_vector))
  neuroim2::map_values(dataset$region_mask, lookup)
}

#' @export
searchlight_scope.mvpa_clustered_dataset <- function(dataset, ...) {
  "regional"
}


# ---- Searchlight iterator ----

#' Clustered ROI specification (lightweight descriptor)
#'
#' @return A \code{clustered_roi_spec} object.
#' @keywords internal
clustered_roi_spec <- function(seed, neighbors) {
  structure(
    list(seed = seed, neighbors = neighbors),
    class = "clustered_roi_spec"
  )
}

#' @export
length.clustered_roi_spec <- function(x) {
  length(x$neighbors)
}

#' @keywords internal
#' @noRd
.clustered_nn_fastpath_enabled <- function() {
  TRUE
}

#' @keywords internal
#' @noRd
.clustered_neighbor_specs_dense <- function(cents, radius = NULL, k = 6) {
  K <- nrow(cents)
  dmat <- as.matrix(stats::dist(cents))

  specs <- vector("list", K)
  for (i in seq_len(K)) {
    if (!is.null(radius) && !is.null(k)) {
      ord <- order(dmat[i, ])
      top_k <- ord[seq_len(min(k + 1L, K))]
      neighbors <- top_k[dmat[i, top_k] <= radius]
    } else if (!is.null(radius)) {
      neighbors <- which(dmat[i, ] <= radius)
    } else {
      ord <- order(dmat[i, ])
      neighbors <- ord[seq_len(min(k + 1L, K))]
    }
    specs[[i]] <- clustered_roi_spec(seed = i, neighbors = as.integer(neighbors))
  }

  specs
}

#' @keywords internal
#' @noRd
.clustered_neighbors_blockwise <- function(cents, radius = NULL, k = 6, block_size = 256L) {
  cents <- as.matrix(cents)
  storage.mode(cents) <- "double"

  K <- nrow(cents)
  if (K == 0L) {
    return(vector("list", 0L))
  }

  if (!is.null(k)) {
    k <- as.integer(k)
    if (is.na(k) || k < 0L) {
      stop("`k` must be NULL or a non-negative integer.", call. = FALSE)
    }
    k_keep <- min(K, k + 1L)
  } else {
    k_keep <- NULL
  }

  if (!is.null(radius)) {
    radius <- as.numeric(radius)
    if (!is.finite(radius) || radius < 0) {
      stop("`radius` must be NULL or a non-negative finite number.", call. = FALSE)
    }
    radius_sq <- radius * radius
    radius_tol <- max(1e-12, radius_sq * 1e-12)
  } else {
    radius_sq <- NULL
    radius_tol <- NULL
  }

  block_size <- suppressWarnings(as.integer(block_size))
  if (!is.finite(block_size) || is.na(block_size) || block_size < 16L) {
    block_size <- 256L
  }

  norms <- rowSums(cents * cents)
  neighbors <- vector("list", K)

  for (start in seq.int(1L, K, by = block_size)) {
    stop_idx <- min(K, start + block_size - 1L)
    q <- cents[start:stop_idx, , drop = FALSE]
    q_norms <- rowSums(q * q)

    cross <- tcrossprod(q, cents)
    dist_sq <- outer(q_norms, norms, "+") - 2 * cross
    dist_sq[dist_sq < 0] <- 0

    for (ii in seq_len(nrow(q))) {
      idx <- start + ii - 1L
      row_dist_sq <- dist_sq[ii, ]

      keep <- if (!is.null(radius_sq)) {
        which(row_dist_sq <= radius_sq + radius_tol)
      } else {
        ord <- order(row_dist_sq, seq_len(K))
        ord[seq_len(k_keep)]
      }

      if (!is.null(k_keep) && length(keep) > k_keep) {
        ord <- order(row_dist_sq[keep], keep)
        keep <- keep[ord[seq_len(k_keep)]]
      } else if (!is.null(k_keep) && length(keep) > 1L) {
        ord <- order(row_dist_sq[keep], keep)
        keep <- keep[ord]
      }

      if (!(idx %in% keep)) {
        keep <- c(idx, keep)
        if (!is.null(k_keep) && length(keep) > k_keep) {
          ord <- order(row_dist_sq[keep], keep)
          keep <- keep[ord[seq_len(k_keep)]]
        }
      }

      neighbors[[idx]] <- as.integer(unique(keep))
    }
  }

  neighbors
}

#' @keywords internal
#' @noRd
.clustered_neighbors_rann <- function(cents, radius = NULL, k = 6) {
  if (!requireNamespace("RANN", quietly = TRUE)) {
    return(NULL)
  }

  cents <- as.matrix(cents)
  K <- nrow(cents)
  if (K == 0L) {
    return(vector("list", 0L))
  }

  k <- as.integer(k)
  if (is.na(k) || k < 0L) {
    stop("`k` must be a non-negative integer when using RANN path.", call. = FALSE)
  }
  k_keep <- min(K, k + 1L)

  if (is.null(radius)) {
    nn <- RANN::nn2(data = cents, query = cents, k = k_keep, searchtype = "standard")
  } else {
    radius <- as.numeric(radius)
    if (!is.finite(radius) || radius < 0) {
      stop("`radius` must be a non-negative finite number.", call. = FALSE)
    }
    nn <- RANN::nn2(
      data = cents,
      query = cents,
      k = k_keep,
      searchtype = "radius",
      radius = radius
    )
  }

  idx <- nn$nn.idx
  dists <- nn$nn.dists
  neighbors <- vector("list", K)

  radius_tol <- if (is.null(radius)) NULL else max(1e-12, radius * 1e-12)

  for (i in seq_len(K)) {
    row_idx <- as.integer(idx[i, ])
    row_dist <- as.numeric(dists[i, ])

    valid <- row_idx > 0L & is.finite(row_dist)
    if (!is.null(radius)) {
      valid <- valid & (row_dist <= radius + radius_tol)
    }

    keep <- row_idx[valid]
    if (length(keep) > 0L) {
      ord <- order(row_dist[valid], keep)
      keep <- keep[ord]
    }

    if (length(keep) > k_keep) {
      keep <- keep[seq_len(k_keep)]
    }

    if (!(i %in% keep)) {
      if (length(keep) < k_keep) {
        keep <- c(i, keep)
      } else if (length(keep) > 0L) {
        keep[length(keep)] <- i
      } else {
        keep <- i
      }
    }

    if (length(keep) > 1L) {
      keep <- unique(keep)
      if (length(keep) > k_keep) {
        keep <- keep[seq_len(k_keep)]
      }
    }

    neighbors[[i]] <- as.integer(keep)
  }

  neighbors
}

#' @keywords internal
#' @noRd
.clustered_neighbor_specs_fast <- function(cents, radius = NULL, k = 6) {
  neighbors <- if (is.null(k)) {
    .clustered_neighbors_blockwise(cents, radius = radius, k = NULL)
  } else if (is.null(radius)) {
    .clustered_neighbors_rann(cents, radius = NULL, k = k)
  } else {
    .clustered_neighbors_rann(cents, radius = radius, k = k)
  }

  if (is.null(neighbors)) {
    neighbors <- .clustered_neighbors_blockwise(cents, radius = radius, k = k)
  }

  lapply(seq_len(length(neighbors)), function(i) {
    clustered_roi_spec(seed = i, neighbors = as.integer(neighbors[[i]]))
  })
}

#' @export
#' @method get_searchlight mvpa_clustered_dataset
get_searchlight.mvpa_clustered_dataset <- function(obj,
                                                    type = c("standard"),
                                                    radius = NULL,
                                                    k = 6,
                                                    ...) {
  type <- match.arg(type)

  cvol <- obj$cvol
  cents <- neuroim2::centroids(cvol)

  if (.clustered_nn_fastpath_enabled()) {
    .clustered_neighbor_specs_fast(cents = cents, radius = radius, k = k)
  } else {
    .clustered_neighbor_specs_dense(cents = cents, radius = radius, k = k)
  }
}


# ---- Data sampling for clustered datasets ----

#' Clustered data sample descriptor
#' @return A \code{clustered_data_sample} object.
#' @keywords internal
clustered_data_sample <- function(seed, neighbors) {
  structure(
    list(
      data = NULL,
      vox  = neighbors,
      seed = seed,
      neighbors = neighbors
    ),
    class = c("clustered_data_sample", "data_sample")
  )
}

#' @export
get_samples.mvpa_clustered_dataset <- function(obj, vox_list) {
  ret <- lapply(vox_list, function(spec) {
    if (inherits(spec, "clustered_roi_spec")) {
      clustered_data_sample(seed = spec$seed, neighbors = spec$neighbors)
    } else {
      # Fallback for plain integer vectors
      clustered_data_sample(seed = spec[1], neighbors = as.integer(spec))
    }
  })

  n <- length(ret)
  df <- tibble::tibble(sample = ret)
  df[[".id"]] <- sprintf(paste0("%0", max(2, nchar(n)), "d"), seq_len(n))
  df
}


#' @export
#' @method as_roi clustered_data_sample
#' @importFrom neuroim2 ROIVec NeuroSpace
as_roi.clustered_data_sample <- function(obj, data, ...) {
  neighbors <- obj$neighbors

  # Extract T x N matrix from the ClusteredNeuroVec time-series matrix
  train_ts <- data$train_data@ts[, neighbors, drop = FALSE]

  # Build coords: round centroids to integer for ROIVec compatibility
  cents <- neuroim2::centroids(data$cvol)
  roi_coords <- round(cents[neighbors, , drop = FALSE])

  # Create 3D space for the ROIVec
  sp3 <- neuroim2::space(data$cvol@mask)

  train_roi <- neuroim2::ROIVec(sp3, coords = roi_coords, data = train_ts)

  test_roi <- if (has_test_set(data)) {
    test_ts <- data$test_data@ts[, neighbors, drop = FALSE]
    neuroim2::ROIVec(sp3, coords = roi_coords, data = test_ts)
  } else {
    NULL
  }

  list(train_roi = train_roi, test_roi = test_roi)
}


# ---- Test data helper ----

#' Generate a Synthetic Clustered MVPA Dataset
#'
#' Creates a synthetic \code{ClusteredNeuroVec} and associated design for testing.
#'
#' @param D Spatial dimensions (e.g. \code{c(10, 10, 10)}).
#' @param nobs Number of observations (time points).
#' @param K Number of clusters.
#' @param nlevels Number of category levels.
#' @param blocks Number of cross-validation blocks.
#' @param external_test Whether to generate an external test set.
#'
#' @return A list with \code{dataset} (mvpa_clustered_dataset) and \code{design} (mvpa_design).
#' @examples
#' ds <- gen_clustered_sample_dataset(K=10, nobs=20, nlevels=2)
#' @export
gen_clustered_sample_dataset <- function(D = c(10, 10, 10), nobs = 20, K = 5,
                                          nlevels = 2, blocks = 3,
                                          external_test = FALSE) {
  sp3 <- neuroim2::NeuroSpace(D, c(1, 1, 1))

  # Create a mask covering the interior
  mask_arr <- array(FALSE, D)
  r <- lapply(D, function(d) max(2, floor(d * 0.2)):min(d, floor(d * 0.8)))
  mask_arr[r[[1]], r[[2]], r[[3]]] <- TRUE
  mask_vol <- neuroim2::LogicalNeuroVol(mask_arr, sp3)

  n_vox <- sum(mask_arr)
  # K-means on coordinates to make clusters
  mask_idx <- which(mask_arr)
  coords <- neuroim2::index_to_grid(sp3, mask_idx)
  set.seed(42)
  km <- stats::kmeans(coords, centers = min(K, n_vox))
  cvol <- neuroim2::ClusteredNeuroVol(mask_vol, km$cluster)

  # Create 4D data and reduce to clustered representation
  sp4 <- neuroim2::NeuroSpace(c(D, nobs), c(1, 1, 1))
  arr4d <- array(stats::rnorm(prod(D) * nobs), c(D, nobs))
  vec4d <- neuroim2::NeuroVec(arr4d, sp4)

  train_data <- neuroim2::ClusteredNeuroVec(vec4d, cvol)

  test_data <- if (external_test) {
    arr_test <- array(stats::rnorm(prod(D) * nobs), c(D, nobs))
    vec_test <- neuroim2::NeuroVec(arr_test, sp4)
    neuroim2::ClusteredNeuroVec(vec_test, cvol)
  } else {
    NULL
  }

  dset <- mvpa_clustered_dataset(train_data, test_data)

  # Design
  Y <- factor(rep(letters[seq_len(nlevels)], length.out = nobs))
  block_var <- as.integer(cut(seq_len(nobs), blocks, labels = seq_len(blocks)))

  if (external_test) {
    Ytest <- factor(rep(letters[seq_len(nlevels)], length.out = nobs))
    des <- mvpa_design(
      data.frame(Y = Y, block_var = block_var),
      test_design = data.frame(Ytest = Ytest),
      block_var = "block_var",
      y_train = ~ Y,
      y_test = ~ Ytest
    )
  } else {
    des <- mvpa_design(
      data.frame(Y = Y, block_var = block_var),
      block_var = "block_var",
      y_train = ~ Y
    )
  }

  list(dataset = dset, design = des)
}
