#' Construct a pair-observation RSA design
#'
#' Generalizes \code{\link{rsa_design}} to support arbitrary pair-observation
#' geometries: lower-triangle within-domain pairs (the classical RSA layout),
#' rectangular between-domain pairs (e.g. items in domain A vs. items in
#' domain B), and function-valued model RDM entries. Function entries may
#' either operate on item identifiers, or on item identifiers plus feature
#' rows when \code{features_a} / \code{features_b} are supplied.
#'
#' The returned object inherits from \code{rsa_design} so it is a drop-in
#' replacement when used with \code{\link{rsa_model}} and the regional /
#' searchlight engines. Within-domain pair designs are fully interoperable
#' with the existing \code{train_model.rsa_model} path; between-domain
#' designs are dispatched on the \code{pair_kind} field, which causes
#' \code{train_model.rsa_model} to compute a rectangular neural-pair
#' dissimilarity block instead of the lower triangle.
#'
#' @param items_a Vector of item identifiers in domain A. Length defines
#'   \code{n_a}.
#' @param items_b Optional vector of item identifiers in domain B. Required
#'   when \code{pairs = "between"}; ignored otherwise. Length defines
#'   \code{n_b}.
#' @param model Named list of model RDM specifications. Each entry may be a
#'   \code{dist} object, a numeric square matrix (within mode) or
#'   \code{n_a x n_b} matrix (between mode), a numeric vector of length
#'   \code{n_pairs}, or a function \code{function(a, b)} returning a numeric
#'   vector of pairwise values for the requested pair table. If
#'   \code{features_a} is supplied and the function has at least four
#'   formal arguments (or \code{...}), it is called as
#'   \code{function(a, b, features_a, features_b)} where the feature
#'   arguments are row-aligned pair tables.
#' @param nuisance Optional named list of nuisance pair predictors using the
#'   same accepted forms as \code{model}. Included in the RSA design matrix
#'   but excluded from model-space fingerprints returned by
#'   \code{rsa_model(..., return_fingerprint = TRUE)}.
#' @param pairs Either \code{"within"} (default) or \code{"between"}.
#' @param features_a Optional data frame or matrix of item features for
#'   \code{items_a}. Function-valued model/nuisance entries can use these
#'   rows to define feature-pair dissimilarities.
#' @param features_b Optional data frame or matrix of item features for
#'   \code{items_b}. Defaults to \code{features_a} in within mode.
#' @param block_var_a Optional vector of block labels (length \code{n_a}). If
#'   supplied and \code{keep_intra_run = FALSE}, within-block pairs are
#'   excluded.
#' @param block_var_b Optional vector of block labels for domain B, length
#'   \code{n_b}. Defaults to \code{block_var_a} in within mode.
#' @param keep_intra_run Logical; if \code{TRUE}, do not drop within-block
#'   pairs.
#' @param row_idx_a Optional integer vector of dataset row indices
#'   corresponding to \code{items_a}. In within mode, this lets a pair design
#'   address a subset/reordering of dataset rows. Required for
#'   \code{pairs = "between"} so the per-ROI engine can extract the correct
#'   neural sub-blocks.
#' @param row_idx_b Integer vector of dataset row indices for \code{items_b}.
#'   Required for \code{pairs = "between"}.
#'
#' @return A list with class \code{c("pair_rsa_design", "rsa_design", "list")}
#'   containing all fields produced by \code{rsa_design()} plus
#'   \describe{
#'     \item{pair_kind}{Either \code{"within"} or \code{"between"}.}
#'     \item{items_a, items_b}{Item identifier vectors.}
#'     \item{n_a, n_b}{Item counts.}
#'     \item{pair_index}{A data frame describing each retained pair.}
#'     \item{row_idx_a, row_idx_b}{Dataset row indices when \code{pairs = "between"}.}
#'   }
#'
#' @seealso \code{\link{rsa_design}}, \code{\link{rsa_model}},
#'   \code{\link{model_space_connectivity}}
#'
#' @export
#' @examples
#' set.seed(1)
#' items <- paste0("item", 1:8)
#' R1 <- as.matrix(dist(matrix(rnorm(8 * 4), 8, 4)))
#' R2 <- as.matrix(dist(matrix(rnorm(8 * 4), 8, 4)))
#' rownames(R1) <- colnames(R1) <- rownames(R2) <- colnames(R2) <- items
#' des <- pair_rsa_design(items, model = list(rdm1 = R1, rdm2 = R2))
#' lengths(des$model_mat)
pair_rsa_design <- function(items_a,
                            items_b = NULL,
                            model = list(),
                            nuisance = list(),
                            pairs = c("within", "between"),
                            features_a = NULL,
                            features_b = NULL,
                            block_var_a = NULL,
                            block_var_b = NULL,
                            keep_intra_run = FALSE,
                            row_idx_a = NULL,
                            row_idx_b = NULL) {
  pairs <- match.arg(pairs)

  if (length(items_a) < 2L) {
    stop("`items_a` must contain at least two items.", call. = FALSE)
  }

  if (identical(pairs, "between")) {
    if (is.null(items_b) || length(items_b) < 1L) {
      stop("`items_b` must be supplied when pairs = 'between'.", call. = FALSE)
    }
    if (is.null(row_idx_a) || is.null(row_idx_b)) {
      stop("`row_idx_a` and `row_idx_b` are required when pairs = 'between'.",
           call. = FALSE)
    }
    if (length(row_idx_a) != length(items_a)) {
      stop("`row_idx_a` must have the same length as `items_a`.", call. = FALSE)
    }
    if (length(row_idx_b) != length(items_b)) {
      stop("`row_idx_b` must have the same length as `items_b`.", call. = FALSE)
    }
  } else {
    if (is.null(items_b)) items_b <- items_a
    if (!identical(length(items_a), length(items_b)) || !all(items_a == items_b)) {
      stop("`items_b` must equal `items_a` when pairs = 'within'.", call. = FALSE)
    }
    if (is.null(block_var_b)) block_var_b <- block_var_a
    if (!is.null(row_idx_a) && length(row_idx_a) != length(items_a)) {
      stop("`row_idx_a` must have the same length as `items_a`.", call. = FALSE)
    }
    if (!is.null(row_idx_b) && !identical(row_idx_b, row_idx_a)) {
      stop("`row_idx_b` is ignored in within mode; use `row_idx_a` only.",
           call. = FALSE)
    }
    row_idx_b <- row_idx_a
  }

  n_a <- length(items_a)
  n_b <- length(items_b)

  if (!is.null(features_a)) {
    if (!is.data.frame(features_a) && !is.matrix(features_a)) {
      stop("`features_a` must be a data frame or matrix.", call. = FALSE)
    }
    if (nrow(features_a) != n_a) {
      stop("`features_a` must have one row per `items_a` entry.", call. = FALSE)
    }
  }
  if (is.null(features_b) && identical(pairs, "within")) {
    features_b <- features_a
  }
  if (!is.null(features_b)) {
    if (!is.data.frame(features_b) && !is.matrix(features_b)) {
      stop("`features_b` must be a data frame or matrix.", call. = FALSE)
    }
    if (nrow(features_b) != n_b) {
      stop("`features_b` must have one row per `items_b` entry.", call. = FALSE)
    }
  }

  if ((!is.list(model) && !is.null(model)) ||
      (!is.list(nuisance) && !is.null(nuisance))) {
    stop("`model` and `nuisance` must be (possibly empty) named lists.",
         call. = FALSE)
  }
  model <- if (is.null(model)) list() else as.list(model)
  nuisance <- if (is.null(nuisance)) list() else as.list(nuisance)
  if (length(model) == 0L && length(nuisance) == 0L) {
    stop("Provide at least one entry in `model` or `nuisance`.", call. = FALSE)
  }

  if (length(model) > 0L && (is.null(names(model)) || any(!nzchar(names(model))))) {
    stop("`model` entries must be named.", call. = FALSE)
  }
  if (length(nuisance) > 0L && (is.null(names(nuisance)) || any(!nzchar(names(nuisance))))) {
    stop("`nuisance` entries must be named.", call. = FALSE)
  }

  pair_index <- .pair_design_pair_index(items_a, items_b, pairs)
  expected_n <- nrow(pair_index)

  vec_model <- .pair_design_vectorize_list(
    model, items_a, items_b, pairs, expected_n, label = "model",
    features_a = features_a, features_b = features_b
  )
  vec_nuis <- .pair_design_vectorize_list(
    nuisance, items_a, items_b, pairs, expected_n, label = "nuisance",
    features_a = features_a, features_b = features_b
  )

  if (!is.null(block_var_a)) {
    if (length(block_var_a) != n_a) {
      stop("`block_var_a` must have length n_a.", call. = FALSE)
    }
    if (identical(pairs, "between")) {
      if (is.null(block_var_b) || length(block_var_b) != n_b) {
        stop("`block_var_b` must have length n_b when pairs = 'between' and ",
             "`block_var_a` is supplied.", call. = FALSE)
      }
    }
  }

  include <- NULL
  if (!is.null(block_var_a) && !isTRUE(keep_intra_run)) {
    include <- if (identical(pairs, "within")) {
      as.vector(stats::dist(as.numeric(as.factor(block_var_a)))) != 0
    } else {
      as.vector(outer(block_var_a, block_var_b, FUN = function(x, y) x != y))
    }
  }

  model_mat_raw <- c(vec_model, vec_nuis)
  model_names <- sanitize(names(vec_model))
  nuisance_names <- sanitize(names(vec_nuis))

  if (!is.null(include)) {
    model_mat <- lapply(model_mat_raw, function(v) v[include])
  } else {
    model_mat <- model_mat_raw
  }
  names(model_mat) <- sanitize(names(model_mat))
  if (anyDuplicated(names(model_mat))) {
    stop("Sanitized model/nuisance predictor names must be unique.",
         call. = FALSE)
  }

  if (length(model_mat) > 0L) {
    rhs <- paste(names(model_mat), collapse = " + ")
    formula <- stats::as.formula(paste("~", rhs))
  } else {
    formula <- NULL
  }

  des <- list(
    formula      = formula,
    data         = model_mat,
    split_by     = NULL,
    split_groups = NULL,
    block_var    = block_var_a,
    include      = include,
    model_mat    = model_mat,
    model_predictors = model_names,
    nuisance_predictors = nuisance_names,
    predictor_roles = stats::setNames(
      c(rep("model", length(model_names)), rep("nuisance", length(nuisance_names))),
      c(model_names, nuisance_names)
    ),
    pair_kind    = pairs,
    items_a      = items_a,
    items_b      = if (identical(pairs, "between")) items_b else items_a,
    features_a   = features_a,
    features_b   = features_b,
    n_a          = n_a,
    n_b          = n_b,
    pair_index   = pair_index,
    row_idx_a    = row_idx_a,
    row_idx_b    = row_idx_b
  )
  class(des) <- c("pair_rsa_design", "rsa_design", "list")
  des
}


#' @export
print.pair_rsa_design <- function(x, ...) {
  cat("Pair RSA design\n")
  cat("  pair_kind:    ", x$pair_kind, "\n", sep = "")
  cat("  n_a x n_b:    ", x$n_a, " x ", x$n_b, "\n", sep = "")
  cat("  raw pairs:    ", nrow(x$pair_index), "\n", sep = "")
  if (!is.null(x$include)) {
    cat("  retained:     ", sum(x$include), "\n", sep = "")
  }
  cat("  predictors:   ", paste(names(x$model_mat), collapse = ", "), "\n", sep = "")
  invisible(x)
}


#' @keywords internal
#' @noRd
.pair_design_pair_index <- function(items_a, items_b, pairs) {
  if (identical(pairs, "within")) {
    n <- length(items_a)
    if (n < 2L) {
      stop("Need at least two items for within-domain pairs.", call. = FALSE)
    }
    M <- matrix(0L, n, n)
    M[lower.tri(M)] <- seq_len(n * (n - 1L) / 2L)
    idx <- which(lower.tri(M), arr.ind = TRUE)
    data.frame(
      i      = idx[, 1L],
      j      = idx[, 2L],
      item_a = items_a[idx[, 1L]],
      item_b = items_a[idx[, 2L]],
      stringsAsFactors = FALSE
    )
  } else {
    n_a <- length(items_a); n_b <- length(items_b)
    grid <- expand.grid(i = seq_len(n_a), j = seq_len(n_b),
                        KEEP.OUT.ATTRS = FALSE)
    data.frame(
      i      = grid$i,
      j      = grid$j,
      item_a = items_a[grid$i],
      item_b = items_b[grid$j],
      stringsAsFactors = FALSE
    )
  }
}


#' @keywords internal
#' @noRd
.pair_design_vectorize_list <- function(lst, items_a, items_b, pairs, expected_n,
                                        label, features_a = NULL,
                                        features_b = NULL) {
  if (length(lst) == 0L) return(list())
  out <- vector("list", length(lst))
  names(out) <- names(lst)
  for (nm in names(lst)) {
    out[[nm]] <- .pair_design_make_block(lst[[nm]], nm, pairs,
                                          items_a, items_b, label = label,
                                          features_a = features_a,
                                          features_b = features_b)
    if (length(out[[nm]]) != expected_n) {
      stop(sprintf("%s entry '%s' produced %d pair values; expected %d.",
                   label, nm, length(out[[nm]]), expected_n),
           call. = FALSE)
    }
  }
  out
}


#' @keywords internal
#' @noRd
.pair_design_make_block <- function(value, name, pairs, items_a, items_b,
                                     label = "model", features_a = NULL,
                                     features_b = NULL) {
  n_a <- length(items_a); n_b <- length(items_b)

  if (is.function(value)) {
    pair_idx <- .pair_design_pair_index(items_a, items_b, pairs)
    fa <- .pair_design_feature_rows(features_a, pair_idx$i)
    fb <- .pair_design_feature_rows(features_b, pair_idx$j)
    use_features <- !is.null(fa) || !is.null(fb)
    fmls <- names(formals(value))
    accepts_features <- use_features &&
      ("..." %in% fmls || length(formals(value)) >= 4L)
    out <- tryCatch(
      {
        if (accepts_features) {
          as.numeric(value(pair_idx$item_a, pair_idx$item_b, fa, fb))
        } else {
          as.numeric(value(pair_idx$item_a, pair_idx$item_b))
        }
      },
      error = function(e) {
        stop(sprintf("%s function '%s' failed: %s", label, name, conditionMessage(e)),
             call. = FALSE)
      }
    )
    return(out)
  }

  if (inherits(value, "dist")) {
    sz <- attr(value, "Size")
    if (identical(pairs, "within")) {
      if (sz != n_a) {
        stop(sprintf("%s '%s': dist size %d does not match length(items_a)=%d.",
                     label, name, sz, n_a), call. = FALSE)
      }
      return(as.vector(value))
    }
    stop(sprintf("%s '%s': dist objects cannot represent rectangular between-domain pair blocks; supply a %d x %d matrix, vector, or function.",
                 label, name, n_a, n_b), call. = FALSE)
  }

  if (is.matrix(value)) {
    if (identical(pairs, "within")) {
      if (nrow(value) != n_a || ncol(value) != n_a) {
        stop(sprintf("%s '%s': matrix is %d x %d; expected %d x %d.",
                     label, name, nrow(value), ncol(value), n_a, n_a),
             call. = FALSE)
      }
      if (!isSymmetric(value)) {
        return(as.vector(stats::dist(value)))
      }
      return(value[lower.tri(value)])
    }
    if (nrow(value) != n_a || ncol(value) != n_b) {
      stop(sprintf("%s '%s': matrix is %d x %d; expected %d x %d for between pairs.",
                   label, name, nrow(value), ncol(value), n_a, n_b),
           call. = FALSE)
    }
    return(as.vector(value))
  }

  if (is.numeric(value)) {
    expected <- if (identical(pairs, "within")) n_a * (n_a - 1L) / 2L else n_a * n_b
    if (length(value) != expected) {
      stop(sprintf("%s '%s': numeric vector length %d; expected %d.",
                   label, name, length(value), expected), call. = FALSE)
    }
    return(as.numeric(value))
  }

  stop(sprintf("%s entry '%s' has unsupported type: %s.",
               label, name, paste(class(value), collapse = "/")), call. = FALSE)
}


#' @keywords internal
#' @noRd
.pair_design_feature_rows <- function(features, idx) {
  if (is.null(features)) return(NULL)
  features[idx, , drop = FALSE]
}
