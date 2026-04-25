#' Model-space representational connectivity from a fitted rMVPA result
#'
#' Unified entry point for second-order representational connectivity. Accepts
#' a matrix/list/tibble of per-unit RDM vectors, or a fitted regional result
#' carrying per-ROI model-space fingerprints / RDM vectors, and forwards to
#' \code{\link{rdm_model_space_connectivity}} for the projection and similarity
#' summaries.
#'
#' This is the seamless front-end advocated in the pair-observation
#' model-space RSA design: within-unit RSA fits and across-unit
#' representational connectivity become two views of the same fitted object.
#'
#' @param result Either a numeric matrix / list / tibble of per-unit RDM
#'   vectors (forwarded directly to \code{rdm_model_space_connectivity}), or a
#'   fitted regional rMVPA object that carries fingerprints / RDM vectors per
#'   ROI (e.g. \code{regional_mvpa_result} from \code{run_regional()} on an
#'   \code{rsa_model(..., return_fingerprint = TRUE)} or a
#'   \code{feature_rsa_model(..., return_rdm_vectors = TRUE)}).
#' @param model_rdms Named list of square symmetric model RDMs / \code{dist}
#'   objects, or a numeric matrix with RDM cells in rows. Required when the
#'   per-unit summaries are RDM vectors. Optional when fingerprints already
#'   live on \code{result}: in that case the fingerprints are taken as the
#'   ROI-by-axis matrix \code{F} and the function returns
#'   \code{tcrossprod(F)} plus its decompositions, without re-projecting.
#' @param method,basis,use,tol,return_projected Forwarded to
#'   \code{rdm_model_space_connectivity()}. See that function for details.
#' @param prefer Either \code{"fingerprint"} (default) or \code{"rdm"} when
#'   both representations are available on \code{result}.
#' @param ... Additional arguments forwarded to underlying methods.
#'
#' @return An object of class \code{rdm_model_space_connectivity} (see
#'   \code{rdm_model_space_connectivity()} for details).
#'
#' @seealso \code{\link{rdm_model_space_connectivity}},
#'   \code{\link{pair_rsa_design}}, \code{\link{rsa_model}},
#'   \code{\link{feature_rsa_rdm_vectors}}
#'
#' @export
model_space_connectivity <- function(result, model_rdms = NULL, ...) {
  UseMethod("model_space_connectivity")
}


#' @rdname model_space_connectivity
#' @export
model_space_connectivity.default <- function(result,
                                             model_rdms = NULL,
                                             method = c("pearson", "spearman"),
                                             basis = c("pca", "qr"),
                                             use = c("complete.obs", "everything"),
                                             tol = 1e-8,
                                             return_projected = FALSE,
                                             ...) {
  if (is.null(model_rdms)) {
    stop("model_space_connectivity: `model_rdms` is required when `result` ",
         "does not carry fingerprints.", call. = FALSE)
  }
  method <- match.arg(method)
  basis  <- match.arg(basis)
  use    <- match.arg(use)
  rdm_model_space_connectivity(
    roi_rdms = result,
    model_rdms = model_rdms,
    method = method,
    basis = basis,
    use = use,
    tol = tol,
    return_projected = return_projected
  )
}


#' @rdname model_space_connectivity
#' @export
model_space_connectivity.regional_mvpa_result <- function(result,
                                                          model_rdms = NULL,
                                                          method = c("pearson", "spearman"),
                                                          basis = c("pca", "qr"),
                                                          use = c("complete.obs", "everything"),
                                                          tol = 1e-8,
                                                          return_projected = FALSE,
                                                          prefer = c("fingerprint", "rdm"),
                                                          ...) {
  prefer <- match.arg(prefer)
  fp <- attr(result, "fingerprints", exact = TRUE)

  if (identical(prefer, "fingerprint") && !is.null(fp) && !is.null(fp$scores)) {
    return(.model_space_connectivity_from_fingerprints(fp, tol = tol))
  }

  # Fall back to feature-RSA-style RDM vectors when present.
  if (is.null(model_rdms)) {
    if (!is.null(fp) && !is.null(fp$scores)) {
      return(.model_space_connectivity_from_fingerprints(fp, tol = tol))
    }
    stop("model_space_connectivity: `model_rdms` is required when fitting ",
         "from per-unit RDM vectors. Either supply `model_rdms` or refit with ",
         "`return_fingerprint = TRUE`.", call. = FALSE)
  }

  method <- match.arg(method)
  basis  <- match.arg(basis)
  use    <- match.arg(use)

  rdm_model_space_connectivity(
    roi_rdms = result,
    model_rdms = model_rdms,
    method = method,
    basis = basis,
    use = use,
    tol = tol,
    return_projected = return_projected
  )
}


#' @rdname model_space_connectivity
#'
#' @details For \code{searchlight_result} inputs the full \eqn{F F^\top}
#' matrix is intentionally never materialized. Instead, \code{k} anchor
#' searchlights are chosen (default: k-means clustering of the per-center
#' fingerprint matrix, with the searchlight closest to each centroid as the
#' anchor) and the connectivity is summarized as an \code{n_centers x k}
#' similarity matrix plus one brain map per anchor. Memory is
#' \eqn{\mathcal{O}(n_{centers} \cdot k)} rather than
#' \eqn{\mathcal{O}(n_{centers}^2)}.
#'
#' @param k Integer; number of anchors to select. Default \code{20}, clamped
#'   to the number of available searchlight centers.
#' @param scale One of \code{"norm"} (default; row-normalize fingerprints so
#'   similarity is cosine) or \code{"raw"} (use raw fingerprint inner
#'   products).
#' @param seeds Anchor selection strategy: \code{"kmeans"} (default; cluster
#'   fingerprints and pick the searchlight closest to each centroid) or an
#'   integer vector of explicit center IDs to use as anchors.
#' @param nstart,iter.max Forwarded to \code{stats::kmeans()}. Default
#'   \code{nstart = 10}, \code{iter.max = 30}.
#' @param random_seed Integer seed for reproducible k-means starts. Default
#'   \code{NULL} (no seeding). When supplied, the caller's global RNG state is
#'   restored after anchor selection.
#' @param build_maps Logical; if \code{TRUE} (default) build one
#'   \code{NeuroVol}/\code{NeuroSurface} per anchor via
#'   \code{\link{build_output_map}}. Skip with \code{FALSE} if you only need
#'   the numerical similarity matrix.
#'
#' @return For \code{searchlight_result} inputs, an object of class
#'   \code{model_space_anchor_connectivity}. This is an anchor summary, not a
#'   full searchlight-by-searchlight matrix. It contains an
#'   \code{n_centers x k} \code{similarity} matrix, anchor center IDs,
#'   optional anchor maps, and clustering metadata when k-means anchors are
#'   used.
#'
#' @export
model_space_connectivity.searchlight_result <- function(result,
                                                        model_rdms = NULL,
                                                        k = 20L,
                                                        scale = c("norm", "raw"),
                                                        seeds = c("kmeans"),
                                                        nstart = 10L,
                                                        iter.max = 30L,
                                                        random_seed = NULL,
                                                        build_maps = TRUE,
                                                        ...) {
  scale <- match.arg(scale)

  if (!is.null(model_rdms)) {
    stop("model_space_connectivity.searchlight_result uses stored fingerprints only; ",
         "`model_rdms` is not used for searchlight anchor summaries.",
         call. = FALSE)
  }
  if (!is.logical(build_maps) || length(build_maps) != 1L || is.na(build_maps)) {
    stop("`build_maps` must be TRUE or FALSE.", call. = FALSE)
  }

  fp <- attr(result, "fingerprints", exact = TRUE)
  if (is.null(fp) || is.null(fp$scores) || !is.matrix(fp$scores)) {
    stop("model_space_connectivity: searchlight_result has no stored fingerprints. ",
         "Refit with `rsa_model(..., return_fingerprint = TRUE)`.",
         call. = FALSE)
  }

  scores <- fp$scores
  ids    <- fp$ids %||% as.integer(rownames(scores))
  if (is.null(ids) || length(ids) != nrow(scores)) {
    ids <- seq_len(nrow(scores))
  }
  ids <- as.integer(ids)
  if (anyNA(ids) || anyDuplicated(ids)) {
    stop("Searchlight fingerprint ids must be unique, non-missing integers.",
         call. = FALSE)
  }

  finite <- stats::complete.cases(scores)
  if (!any(finite)) {
    stop("All searchlight fingerprints are non-finite.", call. = FALSE)
  }
  if (!all(finite)) {
    scores <- scores[finite, , drop = FALSE]
    ids    <- ids[finite]
  }

  F_use <- if (identical(scale, "norm")) {
    norms <- sqrt(rowSums(scores^2))
    keep  <- is.finite(norms) & norms > 0
    out <- scores
    out[keep, ]  <- out[keep, , drop = FALSE] / norms[keep]
    out[!keep, ] <- 0
    out
  } else {
    scores
  }

  if (is.numeric(seeds) && length(seeds) > 0L) {
    anchors_id <- as.integer(seeds)
    if (anyNA(anchors_id)) {
      stop("Explicit `seeds` must be non-missing searchlight center ids.",
           call. = FALSE)
    }
    if (anyDuplicated(anchors_id)) {
      stop("Explicit `seeds` must be unique.", call. = FALSE)
    }
    anchor_pos <- match(anchors_id, ids)
    if (any(is.na(anchor_pos))) {
      stop("Some explicit seeds are not present in the searchlight fingerprint ids.",
           call. = FALSE)
    }
    centroids  <- F_use[anchor_pos, , drop = FALSE]
    cluster_id <- rep(NA_integer_, nrow(F_use))
    anchor_method <- "explicit"
  } else {
    seeds <- match.arg(as.character(seeds), c("kmeans"))
    k <- .model_space_integer_scalar(k, "k", min = 1L)
    nstart <- .model_space_integer_scalar(nstart, "nstart", min = 1L)
    iter.max <- .model_space_integer_scalar(iter.max, "iter.max", min = 1L)
    n_distinct <- nrow(unique(F_use))
    k_eff <- max(1L, min(k, n_distinct))

    if (k_eff <= 1L) {
      grand <- colMeans(F_use)
      d <- rowSums((F_use - matrix(grand, nrow(F_use), ncol(F_use),
                                    byrow = TRUE))^2)
      anchor_pos <- which.min(d)
      centroids  <- matrix(grand, nrow = 1L)
      cluster_id <- rep(1L, nrow(F_use))
    } else if (k_eff >= nrow(F_use)) {
      anchor_pos <- seq_len(nrow(F_use))
      centroids  <- F_use
      cluster_id <- seq_len(nrow(F_use))
    } else {
      km <- tryCatch(
        .model_space_with_seed(
          random_seed,
          stats::kmeans(F_use, centers = k_eff, nstart = nstart,
                        iter.max = iter.max)
        ),
        error = function(e) {
          stop("k-means clustering of fingerprints failed: ", conditionMessage(e),
               call. = FALSE)
        }
      )
      centroids  <- km$centers
      cluster_id <- km$cluster

      anchor_pos <- integer(k_eff)
      for (j in seq_len(k_eff)) {
        members <- which(cluster_id == j)
        if (length(members) == 0L) {
          anchor_pos[j] <- which.min(rowSums((F_use - matrix(centroids[j, ],
                                                              nrow(F_use),
                                                              ncol(F_use),
                                                              byrow = TRUE))^2))
        } else {
          d <- rowSums((F_use[members, , drop = FALSE] -
                          matrix(centroids[j, ], length(members),
                                 ncol(F_use), byrow = TRUE))^2)
          anchor_pos[j] <- members[which.min(d)]
        }
      }
    }
    anchors_id <- ids[anchor_pos]
    anchor_method <- "kmeans"
  }

  similarity <- tcrossprod(F_use, F_use[anchor_pos, , drop = FALSE])
  rownames(similarity) <- as.character(ids)
  colnames(similarity) <- paste0("anchor", seq_along(anchors_id))

  vol_results <- NULL
  if (isTRUE(build_maps) && !is.null(fp$dataset)) {
    vol_results <- vector("list", ncol(similarity))
    for (j in seq_len(ncol(similarity))) {
      vol_results[[j]] <- tryCatch(
        build_output_map(fp$dataset, similarity[, j], ids),
        error = function(e) {
          futile.logger::flog.warn(
            "model_space_connectivity: build_output_map failed for anchor %d: %s",
            j, conditionMessage(e))
          NULL
        }
      )
    }
    names(vol_results) <- colnames(similarity)
  }

  out <- list(
    method        = paste0(anchor_method, "_anchors"),
    anchor_method = anchor_method,
    n_centers     = nrow(F_use),
    k             = length(anchors_id),
    anchors       = anchors_id,
    anchor_pos    = anchor_pos,
    cluster_id    = cluster_id,
    centroids     = centroids,
    similarity    = similarity,
    vol_results   = vol_results,
    scale         = scale,
    fingerprint_dim = ncol(F_use)
  )
  structure(out, class = c("model_space_anchor_connectivity", "list"))
}


#' @keywords internal
#' @noRd
.model_space_integer_scalar <- function(x, name, min = 1L) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || is.na(x) ||
      x < min || x != as.integer(x)) {
    if (min <= 0L) {
      stop(sprintf("`%s` must be an integer scalar >= %d.", name, min),
           call. = FALSE)
    }
    stop(sprintf("`%s` must be a positive integer scalar.", name), call. = FALSE)
  }
  as.integer(x)
}


#' @keywords internal
#' @noRd
.model_space_with_seed <- function(seed, expr) {
  if (is.null(seed)) {
    return(force(expr))
  }
  seed <- .model_space_integer_scalar(seed, "random_seed", min = 0L)
  has_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_seed <- if (has_seed) get(".Random.seed", envir = .GlobalEnv) else NULL
  on.exit({
    if (has_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(seed)
  force(expr)
}


#' @export
print.model_space_anchor_connectivity <- function(x, ...) {
  cat("Model-space anchor connectivity\n")
  cat("  method:           ", x$method, "\n", sep = "")
  cat("  n_centers:        ", x$n_centers, "\n", sep = "")
  cat("  k anchors:        ", x$k, "\n", sep = "")
  cat("  fingerprint dim:  ", x$fingerprint_dim, "\n", sep = "")
  cat("  scale:            ", x$scale, "\n", sep = "")
  if (!is.null(x$vol_results)) {
    n_maps <- sum(!vapply(x$vol_results, is.null, logical(1)))
    cat("  anchor maps:      ", n_maps, " of ", length(x$vol_results), "\n", sep = "")
  }
  invisible(x)
}


#' @keywords internal
#' @noRd
.model_space_connectivity_from_fingerprints <- function(fp, tol = 1e-8) {
  scores <- fp$scores
  if (!is.matrix(scores) || nrow(scores) < 1L || ncol(scores) < 1L) {
    stop("Stored fingerprints are empty.", call. = FALSE)
  }

  axis_names <- colnames(scores)
  if (is.null(axis_names) || any(!nzchar(axis_names))) {
    axis_names <- paste0("PC", seq_len(ncol(scores)))
    colnames(scores) <- axis_names
  }

  roi_labels <- rownames(scores)
  if (is.null(roi_labels) || any(!nzchar(roi_labels))) {
    roi_labels <- paste0("unit", seq_len(nrow(scores)))
    rownames(scores) <- roi_labels
  }

  similarity <- tcrossprod(scores)
  dimnames(similarity) <- list(roi_labels, roi_labels)

  norms <- sqrt(rowSums(scores^2))
  profile_scores <- scores
  keep <- is.finite(norms) & norms > tol
  profile_scores[keep, ] <- profile_scores[keep, , drop = FALSE] / norms[keep]
  profile_scores[!keep, ] <- NA_real_
  profile_similarity <- tcrossprod(profile_scores)
  dimnames(profile_similarity) <- list(roi_labels, roi_labels)

  component_similarity <- lapply(seq_len(ncol(scores)), function(i) {
    mat <- tcrossprod(scores[, i, drop = FALSE])
    dimnames(mat) <- list(roi_labels, roi_labels)
    mat
  })
  names(component_similarity) <- axis_names

  common_similarity <- component_similarity[[1L]]
  difference_similarity <- if (length(component_similarity) > 1L) {
    Reduce(`+`, component_similarity[-1L])
  } else {
    matrix(0, nrow = nrow(similarity), ncol = ncol(similarity),
           dimnames = dimnames(similarity))
  }

  out <- list(
    similarity = similarity,
    profile_similarity = profile_similarity,
    component_similarity = component_similarity,
    common_similarity = common_similarity,
    difference_similarity = difference_similarity,
    raw_similarity = NULL,
    residual_similarity = NULL,
    scores = scores,
    raw_model_scores = NULL,
    model_axis_cor = NULL,
    model_cor = NULL,
    eigenvalues = NULL,
    rank = ncol(scores),
    basis = "fingerprint",
    method = "fingerprint",
    n_pairs_total = NA_integer_,
    n_pairs_used = NA_integer_,
    roi_labels = roi_labels,
    model_labels = axis_names,
    decomposition_available = FALSE
  )
  structure(out, class = c("rdm_model_space_connectivity", "list"))
}
