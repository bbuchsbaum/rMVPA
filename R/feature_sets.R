#' Feature sets: grouped predictor matrices
#'
#' Many rMVPA analyses use continuous predictors (rows = TRs/observations, columns
#' = stimulus features) to explain continuous neural responses (rows = TRs,
#' columns = voxels/parcels). In practice, stimulus predictors often come in
#' \emph{multiple correlated blocks} (e.g. VGG low/mid/high/semantic PCs, or audio
#' vs vision features), and it is useful to:
#' \itemize{
#'   \item regularize each block differently (banded/grouped ridge), and
#'   \item compute attribution/competition measures such as leave-one-set-out \eqn{\Delta R^2}.
#' }
#'
#' The `feature_sets()` constructor wraps a predictor matrix (or a list of per-set
#' matrices) into a `feature_sets` object that carries set labels and column
#' indices so that downstream models do not need to manually track slices.
#'
#' @section Feature set specifications:
#' Use `blocks()` when your columns are already concatenated into consecutive blocks
#' (e.g. \code{[low|mid|high|sem]}). Use `by_set()` when you have a per-column set
#' label vector (non-contiguous groups).
#'
#' @section Row weights:
#' A `feature_sets` object also stores `row_weights`. These are optional observation
#' weights (length = \code{nrow(X)}). They are primarily intended for soft-alignment
#' use cases where recall predictors are computed as \code{gamma \%*\% X_enc} and
#' the remaining posterior mass can be used to down-weight uncertain recall TRs
#' (see `expected_features()`).
#'
#' @name feature_sets
#' @seealso
#'   \code{\link{blocks}}, \code{\link{by_set}}, \code{\link{feature_sets}},
#'   \code{\link{expected_features}}, \code{\link{feature_sets_design}},
#'   \code{\link{banded_ridge_da_model}}, \code{\link{grouped_ridge_da_model}},
#'   \code{\link{banded_ridge_da}}, \code{\link{grouped_ridge_da}}
#' @examples
#' X <- matrix(rnorm(20 * 8), 20, 8)
#' fs <- feature_sets(X, blocks(low = 3, sem = 5))
#' fs
NULL

#' Define consecutive column blocks as feature sets
#'
#' Convenience constructor for the common case where columns are arranged as
#' consecutive blocks (e.g., low/mid/high/sem PCs concatenated).
#'
#' @details
#' `blocks()` does not touch any data; it only declares how columns should be
#' interpreted when passed to `feature_sets()`.
#'
#' @param ... Named integers giving number of columns per set (must sum to `ncol(X)`).
#' @return An object of class `feature_set_spec_blocks`.
#' @export
#' @examples
#' spec <- blocks(low = 100, mid = 100, high = 100, sem = 100)
blocks <- function(...) {
  sizes <- c(...)
  nm <- names(sizes)
  if (length(sizes) < 1L || is.null(nm) || any(nm == "")) {
    stop("blocks(): provide one or more named integer sizes, e.g. blocks(low=100, mid=100).", call. = FALSE)
  }
  if (!is.numeric(sizes) || any(!is.finite(sizes)) || any(sizes <= 0)) {
    stop("blocks(): all sizes must be positive numbers.", call. = FALSE)
  }
  sizes <- as.integer(sizes)
  names(sizes) <- nm
  structure(
    list(sizes = sizes, order = nm),
    class = c("feature_set_spec_blocks", "list")
  )
}

#' Define feature sets by a per-column set label
#'
#' Use this when columns are not arranged as contiguous blocks.
#'
#' @details
#' This is useful when:
#' \itemize{
#'   \item you have interleaved predictors (e.g. time-lagged features),
#'   \item you want to group columns by an external annotation, or
#'   \item you already have column names encoding the set membership.
#' }
#'
#' @param set Character or factor vector of length ncol(X), naming the set for each column.
#' @param order Optional character vector giving the desired set order.
#' @return An object of class `feature_set_spec_by`.
#' @export
#' @examples
#' X <- matrix(rnorm(10 * 6), 10, 6)
#' set <- rep(c("audio", "vision"), each = 3)
#' fs <- feature_sets(X, by_set(set, order = c("vision", "audio")))
#' fs
by_set <- function(set, order = NULL) {
  if (is.factor(set)) set <- as.character(set)
  if (!is.character(set) || length(set) < 1L) {
    stop("by_set(): 'set' must be a character/factor vector naming the set for each column.", call. = FALSE)
  }
  if (is.null(order)) order <- unique(set)
  if (!is.character(order) || any(order == "")) {
    stop("by_set(): 'order' must be a character vector of set names.", call. = FALSE)
  }
  structure(
    list(set = set, order = order),
    class = c("feature_set_spec_by", "list")
  )
}

#' Construct a `feature_sets` object
#'
#' @description
#' Wrap a predictor representation into a `feature_sets` instance.
#'
#' @details
#' You can supply predictors in two equivalent ways:
#' \enumerate{
#'   \item A single matrix \code{X} (observations x features) plus a set specification
#'         via `blocks()` or `by_set()`.
#'   \item A named list of matrices, one per set, each with the same number of rows.
#'         The matrices are column-bound in `set_order`.
#' }
#'
#' The returned object includes:
#' \itemize{
#'   \item `X`: the concatenated numeric predictor matrix,
#'   \item `set`: a factor of length `ncol(X)` giving per-column set membership,
#'   \item `indices`: a named list mapping set name -> integer column indices,
#'   \item `dims`: number of columns per set,
#'   \item `row_weights`: optional observation weights (length `nrow(X)`).
#' }
#'
#' @param x Either (1) a numeric matrix (observations x features) with a set spec,
#'   or (2) a named list of numeric matrices, one per set (all with the same number
#'   of rows), which will be column-bound in `set_order`.
#' @param spec A set spec as created by `blocks()` or `by_set()` (required when `x` is a matrix).
#' @param set_order Optional character vector giving the desired set order (defaults to
#'   the order implied by `spec` or the names of the list).
#' @param row_weights Optional numeric vector of length nrow(X), used as observation weights
#'   by downstream models (default: all 1).
#'
#' @return An object of class `feature_sets`.
#' @export
#' @seealso \code{\link{blocks}}, \code{\link{by_set}}, \code{\link{expected_features}}
#' @examples
#' # 1) Matrix input + blocks()
#' X <- matrix(rnorm(20 * 8), 20, 8)
#' fs <- feature_sets(X, blocks(low = 3, sem = 5))
#'
#' # 2) List input (already split per set)
#' Xlist <- list(low = X[, 1:3], sem = X[, 4:8])
#' fs2 <- feature_sets(Xlist)
feature_sets <- function(x, spec = NULL, set_order = NULL, row_weights = NULL) {
  if (is.matrix(x)) {
    if (is.null(spec)) {
      stop("feature_sets(): 'spec' is required when 'x' is a matrix (use blocks(...) or by_set(...)).", call. = FALSE)
    }
    X <- as.matrix(x)
    if (!is.numeric(X)) stop("feature_sets(): 'x' must be a numeric matrix.", call. = FALSE)

    if (inherits(spec, "feature_set_spec_blocks")) {
      sizes <- spec$sizes
      order <- spec$order
      if (!is.null(set_order)) order <- set_order
      if (sum(sizes) != ncol(X)) {
        stop(sprintf("feature_sets(): blocks() sizes sum to %d but ncol(X) is %d.",
                     sum(sizes), ncol(X)), call. = FALSE)
      }
      indices <- list()
      set_vec <- character(ncol(X))
      start <- 1L
      for (nm in order) {
        k <- as.integer(sizes[nm])
        if (length(k) != 1L || is.na(k)) {
          stop(sprintf("feature_sets(): no size found for set '%s' in blocks() spec.", nm), call. = FALSE)
        }
        idx <- start:(start + k - 1L)
        indices[[nm]] <- idx
        set_vec[idx] <- nm
        start <- start + k
      }
      set_fac <- factor(set_vec, levels = order)
    } else if (inherits(spec, "feature_set_spec_by")) {
      set_vec <- spec$set
      if (length(set_vec) != ncol(X)) {
        stop(sprintf("feature_sets(): by_set() provided %d set labels but ncol(X) is %d.",
                     length(set_vec), ncol(X)), call. = FALSE)
      }
      order <- spec$order
      if (!is.null(set_order)) order <- set_order
      if (any(!set_vec %in% order)) {
        bad <- unique(set_vec[!set_vec %in% order])
        stop("feature_sets(): by_set() contains set labels not present in 'order': ",
             paste(bad, collapse = ", "), call. = FALSE)
      }
      indices <- split(seq_len(ncol(X)), factor(set_vec, levels = order))
      # Drop empty groups while preserving order
      indices <- indices[order[order %in% names(indices)]]
      set_fac <- factor(set_vec, levels = order)
    } else {
      stop("feature_sets(): unknown 'spec' type (use blocks(...) or by_set(...)).", call. = FALSE)
    }

    if (is.null(row_weights)) {
      row_weights <- rep(1, nrow(X))
    }
    if (!is.numeric(row_weights) || length(row_weights) != nrow(X)) {
      stop("feature_sets(): 'row_weights' must be a numeric vector of length nrow(X).", call. = FALSE)
    }
    if (any(!is.finite(row_weights)) || any(row_weights < 0)) {
      stop("feature_sets(): 'row_weights' must be finite and non-negative.", call. = FALSE)
    }

    ret <- list(
      X = X,
      set = set_fac,
      indices = indices,
      dims = vapply(indices, length, integer(1)),
      set_order = levels(set_fac),
      row_weights = as.numeric(row_weights)
    )
    class(ret) <- c("feature_sets", "list")
    return(ret)
  }

  if (is.list(x)) {
    if (length(x) < 1L) stop("feature_sets(): list input must have at least one element.", call. = FALSE)
    if (is.null(names(x)) || any(names(x) == "")) {
      stop("feature_sets(): list input must be a *named* list of matrices.", call. = FALSE)
    }
    if (is.null(set_order)) set_order <- names(x)
    if (!all(set_order %in% names(x))) {
      stop("feature_sets(): 'set_order' contains names not present in the input list.", call. = FALSE)
    }
    mats <- lapply(set_order, function(nm) {
      X <- x[[nm]]
      if (!is.matrix(X) || !is.numeric(X)) {
        stop(sprintf("feature_sets(): element '%s' must be a numeric matrix.", nm), call. = FALSE)
      }
      X
    })
    n_rows <- vapply(mats, nrow, integer(1))
    if (length(unique(n_rows)) != 1L) {
      stop("feature_sets(): all matrices must have the same number of rows (observations).", call. = FALSE)
    }
    X <- do.call(cbind, mats)
    set_vec <- rep(set_order, vapply(mats, ncol, integer(1)))
    set_fac <- factor(set_vec, levels = set_order)
    indices <- split(seq_len(ncol(X)), set_fac)

    if (is.null(row_weights)) {
      row_weights <- rep(1, nrow(X))
    }
    if (!is.numeric(row_weights) || length(row_weights) != nrow(X)) {
      stop("feature_sets(): 'row_weights' must be a numeric vector of length nrow(X).", call. = FALSE)
    }
    if (any(!is.finite(row_weights)) || any(row_weights < 0)) {
      stop("feature_sets(): 'row_weights' must be finite and non-negative.", call. = FALSE)
    }

    ret <- list(
      X = X,
      set = set_fac,
      indices = indices,
      dims = vapply(indices, length, integer(1)),
      set_order = levels(set_fac),
      row_weights = as.numeric(row_weights)
    )
    class(ret) <- c("feature_sets", "list")
    return(ret)
  }

  stop("feature_sets(): 'x' must be a matrix or a named list of matrices.", call. = FALSE)
}

#' Print method for feature_sets
#'
#' @param x feature_sets object
#' @param ... ignored
#' @return Invisibly returns the input object \code{x} (called for side effects).
#' @examples
#' \dontrun{
#'   # print method called on feature_sets object
#' }
#' @export
print.feature_sets <- function(x, ...) {
  cat("feature_sets\n")
  cat("===========\n\n")
  cat(sprintf("Observations: %d\n", nrow(x$X)))
  cat(sprintf("Features:     %d\n", ncol(x$X)))
  cat(sprintf("Sets:         %d (%s)\n", length(x$indices), paste(names(x$indices), collapse = ", ")))
  invisible(x)
}

#' Build expected-domain features from a soft alignment matrix
#'
#' Given encoding-domain predictors and a recall->encoding alignment posterior,
#' compute expected recall-domain predictors:
#' \deqn{X_{rec} = \Gamma X_{enc}.}
#'
#' This is the core "soft label" trick for bringing recall into regression when
#' recall TRs do not have a known one-to-one correspondence with encoding TRs.
#'
#' @details
#' \strong{Gamma shapes.}
#' `gamma` should be a numeric matrix where rows index recall TRs and columns
#' index encoding TRs:
#' \itemize{
#'   \item without a NULL state: \code{(T_rec x T_enc)}
#'   \item with a NULL state in the first column: \code{(T_rec x (T_enc+1))}
#' }
#'
#' When a NULL column is present and `drop_null = TRUE`, the NULL column is dropped.
#' If `renormalize = FALSE` (default), the remaining row mass is stored as
#' \code{row_weights} (so uncertain TRs with high NULL probability can be
#' down-weighted by downstream models). If `renormalize = TRUE`, rows are
#' renormalized to sum to 1 and `row_weights` is set to 1.
#'
#' @param train feature_sets object for encoding-domain predictors.
#' @param gamma Numeric matrix of shape (T_rec x T_enc) or (T_rec x (T_enc+1)) if null column present.
#' @param drop_null Logical; if TRUE and gamma has T_enc+1 columns, drop the first column.
#' @param renormalize Logical; if TRUE, renormalize rows to sum to 1 after dropping null.
#' @param eps Small constant to avoid division by zero in renormalization.
#'
#' @return A feature_sets object for recall-domain predictors with the same set layout.
#' @export
#' @seealso \code{\link{feature_sets}}, \code{\link{feature_sets_design}}
#' @examples
#' X <- matrix(rnorm(10 * 5), 10, 5)
#' fs_enc <- feature_sets(X, blocks(a = 2, b = 3))
#' gamma <- matrix(runif(6 * 10), 6, 10)
#' gamma <- gamma / rowSums(gamma)
#' fs_rec <- expected_features(fs_enc, gamma, drop_null = FALSE, renormalize = TRUE)
expected_features <- function(train,
                              gamma,
                              drop_null = TRUE,
                              renormalize = FALSE,
                              eps = 1e-12) {
  if (!inherits(train, "feature_sets")) stop("expected_features(): 'train' must be a feature_sets object.", call. = FALSE)
  G <- as.matrix(gamma)
  if (!is.numeric(G) || nrow(G) < 1L || ncol(G) < 1L) stop("expected_features(): 'gamma' must be a numeric matrix.", call. = FALSE)

  T_enc <- nrow(train$X)
  if (drop_null && ncol(G) == (T_enc + 1L)) {
    G <- G[, -1L, drop = FALSE]
  }
  if (ncol(G) != T_enc) {
    stop(sprintf("expected_features(): gamma has %d columns but expected %d (encoding TRs).",
                 ncol(G), T_enc), call. = FALSE)
  }

  row_mass <- rowSums(G)
  if (renormalize) {
    denom <- pmax(row_mass, eps)
    G <- G / denom
    row_weights <- rep(1, nrow(G))
  } else {
    row_weights <- row_mass
  }

  X_rec <- G %*% train$X

  ret <- train
  ret$X <- X_rec
  ret$row_weights <- as.numeric(row_weights)
  ret
}
