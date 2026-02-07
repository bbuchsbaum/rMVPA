# Representational Network Analysis: Shared Design Helpers
#
# This file defines design helper functions shared by the ReNA models:
# - repnet_design(): seed-based representational connectivity (ReNA-RC)
# - repmed_design(): representational mediation (ReNA-RM)
# - repmap_design(): representational mapping from seed features (ReNA-Map)


#' Representational connectivity design helper (ReNA-RC)
#'
#' Build a design object for representational connectivity, specifying item keys,
#' a seed RDM, and optional confound RDMs.
#'
#' @param design An \code{mvpa_design} object (as used elsewhere in rMVPA).
#' @param key_var Column name or formula giving item identity (e.g. \code{~ ImageID}).
#' @param seed_rdm A K x K matrix or \code{"dist"} object; rows/cols labelled by item IDs.
#' @param confound_rdms Optional named list of K x K matrices or \code{"dist"} objects,
#'   used as nuisance RDMs (e.g. block, lag, behavior).
#'
#' @return A list with fields:
#'   \itemize{
#'     \item \code{key}: factor of item IDs (length = nrow(design$train_design))
#'     \item \code{seed_rdm}: seed RDM as a matrix with row/colnames
#'     \item \code{confound_rdms}: named list of confound RDM matrices
#'   }
#' @export
repnet_design <- function(design,
                          key_var,
                          seed_rdm,
                          confound_rdms = NULL) {

  assertthat::assert_that(inherits(design, "mvpa_design"))
  d <- design$train_design

  key <- parse_variable(key_var, d)
  key_fac <- factor(key)

  # Coerce seed_rdm to a square matrix with item labels
  if (inherits(seed_rdm, "dist")) {
    lab <- attr(seed_rdm, "Labels")
    if (is.null(lab)) {
      # Fall back to key levels if available
      lab <- levels(key_fac)[seq_len(attr(seed_rdm, "Size"))]
    }
    S <- as.matrix(seed_rdm)
    rownames(S) <- colnames(S) <- lab
  } else {
    S <- as.matrix(seed_rdm)
    if (is.null(rownames(S)) || is.null(colnames(S))) {
      # Default labels from key levels if possible; otherwise simple sequence
      default_lab <- levels(key_fac)
      if (length(default_lab) < nrow(S)) {
        default_lab <- as.character(seq_len(nrow(S)))
      }
      rownames(S) <- colnames(S) <- default_lab[seq_len(nrow(S))]
    }
  }

  conf_list <- list()
  if (!is.null(confound_rdms)) {
    assertthat::assert_that(
      is.list(confound_rdms),
      !is.null(names(confound_rdms)),
      msg = "confound_rdms must be a named list."
    )

    conf_list <- lapply(confound_rdms, function(M) {
      if (inherits(M, "dist")) {
        lab <- attr(M, "Labels")
        if (is.null(lab)) {
          lab <- rownames(S)
          if (is.null(lab)) {
            lab <- levels(key_fac)[seq_len(attr(M, "Size"))]
          }
        }
        Mmat <- as.matrix(M)
        rownames(Mmat) <- colnames(Mmat) <- lab
      } else {
        Mmat <- as.matrix(M)
        if (is.null(rownames(Mmat)) || is.null(colnames(Mmat))) {
          default_lab <- rownames(S)
          if (is.null(default_lab) || length(default_lab) < nrow(Mmat)) {
            default_lab <- as.character(seq_len(nrow(Mmat)))
          }
          rownames(Mmat) <- colnames(Mmat) <- default_lab[seq_len(nrow(Mmat))]
        }
      }
      Mmat
    })
  }

  list(
    key           = key_fac,
    seed_rdm      = S,
    confound_rdms = conf_list
  )
}


#' Representational mediation design helper (ReNA-RM)
#'
#' Construct a design object for representational mediation, with predictor (X)
#' and outcome (Y) RDMs over items and optional confound RDMs.
#'
#' @param items Character vector of item IDs (keys).
#' @param X_rdm Predictor RDM (matrix or \code{"dist"}).
#' @param Y_rdm Outcome RDM (matrix or \code{"dist"}).
#' @param confound_rdms Optional named list of confound RDMs (matrix or \code{"dist"}).
#'
#' @return A list with elements:
#'   \itemize{
#'     \item \code{items}: character vector of items
#'     \item \code{X_rdm}, \code{Y_rdm}: predictor and outcome RDMs as matrices
#'     \item \code{confound_rdms}: named list of confound RDM matrices
#'   }
#' @export
repmed_design <- function(items,
                          X_rdm,
                          Y_rdm,
                          confound_rdms = NULL) {

  items <- as.character(items)

  as_mat <- function(M, default_items) {
    if (inherits(M, "dist")) {
      lab <- attr(M, "Labels")
      if (is.null(lab)) {
        lab <- default_items[seq_len(attr(M, "Size"))]
      }
      M <- as.matrix(M)
      rownames(M) <- colnames(M) <- lab
    } else {
      M <- as.matrix(M)
      if (is.null(rownames(M)) || is.null(colnames(M))) {
        rownames(M) <- colnames(M) <- default_items[seq_len(nrow(M))]
      }
    }
    M
  }

  X <- as_mat(X_rdm, items)
  Y <- as_mat(Y_rdm, items)

  conf_list <- list()
  if (!is.null(confound_rdms)) {
    assertthat::assert_that(
      is.list(confound_rdms),
      !is.null(names(confound_rdms)),
      msg = "confound_rdms must be a named list."
    )
    conf_list <- lapply(confound_rdms, as_mat, default_items = items)
  }

  list(
    items         = items,
    X_rdm         = X,
    Y_rdm         = Y,
    confound_rdms = conf_list
  )
}


#' Representational mapping design helper (ReNA-Map)
#'
#' Construct a design object for representational mapping, specifying item-wise
#' seed feature vectors.
#'
#' @param items Character vector of item IDs (keys).
#' @param seed_features K x P matrix of seed features; rows = items.
#'
#' @return A list with elements:
#'   \itemize{
#'     \item \code{items}: character vector of items
#'     \item \code{seed_features}: K x P feature matrix with rownames = items
#'   }
#' @export
repmap_design <- function(items,
                          seed_features) {

  items <- as.character(items)
  F <- as.matrix(seed_features)

  if (is.null(rownames(F))) {
    rownames(F) <- items[seq_len(nrow(F))]
  }

  list(
    items         = items,
    seed_features = F
  )
}

