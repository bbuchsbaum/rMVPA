#' Compute item-level global similarity / RDMs for ERA nuisance modeling
#'
#' Extracts trial patterns from the full \code{dataset$mask}, builds item-level
#' prototypes, and returns:
#' \itemize{
#'   \item \code{S_cross}: K x K cross-state similarity (target rows by source columns)
#'   \item \code{D_enc}, \code{D_ret}: K x K within-state RDMs computed with \code{distfun}
#'   \item \code{common_keys}: sorted character vector of item keys
#' }
#'
#' Used internally by \code{era_rsa_model()} and \code{era_partition_model()}
#' when \code{global_nuisance = TRUE}: the resulting matrices are added as
#' nuisance regressors so that ROI/searchlight similarity can be residualized
#' against the whole-mask similarity, removing pervasive non-ROI-specific
#' structure (e.g., shared temporal autocorrelation, arousal, motion).
#'
#' @section Caveat:
#' The full mask contains the ROI/sphere voxels themselves. For small ROIs and
#' searchlight spheres in a whole-brain mask this is negligible; for large
#' regional ROIs that occupy most of the mask, the residualization partially
#' removes the local signal too. v1 does not implement leave-one-out
#' (mask-minus-current-ROI) global computation.
#'
#' @param dataset An \code{mvpa_dataset} with \code{train_data},
#'   \code{test_data}, and \code{mask}.
#' @param design  An \code{mvpa_design} whose train/test design tables provide
#'   \code{key_var} levels.
#' @param key_var Item key column or formula (matched against
#'   \code{design$train_design} / \code{design$test_design}).
#' @param distfun A distfun used for within-state RDMs.
#' @return A list with \code{common_keys}, \code{S_cross}, \code{D_enc},
#'   \code{D_ret}, or \code{NULL} if extraction fails or fewer than two common
#'   items are present.
#' @noRd
.era_compute_global_rdms <- function(dataset, design, key_var, distfun) {
  if (is.null(dataset$mask)) {
    return(NULL)
  }
  vox_full <- try(compute_mask_indices(dataset$mask), silent = TRUE)
  if (inherits(vox_full, "try-error") || !length(vox_full)) {
    return(NULL)
  }

  pull <- function(src) {
    if (is.null(src)) return(NULL)
    roi <- try(neuroim2::series_roi(src, vox_full), silent = TRUE)
    if (inherits(roi, "try-error")) {
      return(NULL)
    }
    out <- try(as.matrix(neuroim2::values(roi)), silent = TRUE)
    if (inherits(out, "try-error")) {
      out <- try(as.matrix(roi), silent = TRUE)
    }
    if (inherits(out, "try-error")) NULL else out
  }

  Xenc_full <- pull(dataset$train_data)
  Xret_full <- pull(dataset$test_data)
  if (is.null(Xenc_full) || is.null(Xret_full)) {
    return(NULL)
  }

  key_tr <- factor(parse_variable(key_var, design$train_design))
  key_te <- factor(parse_variable(key_var, design$test_design),
                   levels = levels(key_tr))

  if (nrow(Xenc_full) != length(key_tr) || nrow(Xret_full) != length(key_te)) {
    return(NULL)
  }

  E_full <- group_means(Xenc_full, margin = 1, group = key_tr)
  R_full <- group_means(Xret_full, margin = 1, group = key_te)

  common <- sort(intersect(rownames(E_full), rownames(R_full)))
  if (length(common) < 2L) {
    return(NULL)
  }
  E <- E_full[common, , drop = FALSE]
  R <- R_full[common, , drop = FALSE]

  sd_E <- apply(E, 2, stats::sd, na.rm = TRUE)
  sd_R <- apply(R, 2, stats::sd, na.rm = TRUE)
  keep <- (sd_E > 0) | (sd_R > 0)
  if (!any(keep)) {
    return(NULL)
  }
  E <- E[, keep, drop = FALSE]
  R <- R[, keep, drop = FALSE]

  S_cross <- suppressWarnings(stats::cor(t(R), t(E), use = "pairwise.complete.obs"))
  rownames(S_cross) <- colnames(S_cross) <- common

  D_enc <- as.matrix(pairwise_dist(distfun, E))
  D_ret <- as.matrix(pairwise_dist(distfun, R))
  rownames(D_enc) <- colnames(D_enc) <- common
  rownames(D_ret) <- colnames(D_ret) <- common

  list(
    common_keys = common,
    S_cross     = S_cross,
    D_enc       = D_enc,
    D_ret       = D_ret
  )
}

#' Resolve a user-facing global_nuisance argument into a list of K x K matrices
#'
#' Accepts \code{TRUE} (auto-compute via \code{.era_compute_global_rdms()}),
#' \code{FALSE}/\code{NULL} (no-op), or a list with elements
#' \code{S_cross}/\code{first}, \code{D_enc}/\code{enc}, \code{D_ret}/\code{ret}
#' (pre-supplied matrices for non-standard dataset backends).
#'
#' @param global_nuisance User-supplied value.
#' @param dataset,design,key_var,distfun Forwarded to the auto-compute helper.
#' @return A list with \code{S_cross}, \code{D_enc}, \code{D_ret},
#'   \code{common_keys}, or \code{NULL} if disabled / extraction failed.
#' @noRd
.era_resolve_global_nuisance <- function(global_nuisance,
                                         dataset, design, key_var, distfun) {
  if (is.null(global_nuisance) || isFALSE(global_nuisance)) {
    return(NULL)
  }
  if (isTRUE(global_nuisance)) {
    return(.era_compute_global_rdms(dataset, design, key_var, distfun))
  }
  if (is.list(global_nuisance)) {
    pick <- function(...) {
      nms <- c(...)
      for (n in nms) {
        if (!is.null(global_nuisance[[n]])) return(as.matrix(global_nuisance[[n]]))
      }
      NULL
    }
    S <- pick("S_cross", "first", "global_first")
    De <- pick("D_enc", "enc", "global_enc")
    Dr <- pick("D_ret", "ret", "global_ret")
    if (is.null(S) && is.null(De) && is.null(Dr)) {
      return(NULL)
    }
    keys <- if (!is.null(De)) rownames(De) else
            if (!is.null(Dr)) rownames(Dr) else
            if (!is.null(S))  rownames(S)  else NULL
    return(list(common_keys = keys, S_cross = S, D_enc = De, D_ret = Dr))
  }
  NULL
}
