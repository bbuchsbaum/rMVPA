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
#'   \code{rsa_model(..., return_fingerprint = TRUE)}. Alternatively, a
#'   right-hand-side formula evaluated on pair metadata, e.g.
#'   \code{~ factor(b.row)} for retrieval-observation intercepts. Its intercept
#'   is supplied by the regression engine, rather than duplicated here.
#' @param modulation Optional named list of right-hand-side formulas, keyed
#'   by names in \code{model}. Each formula's design columns multiply that
#'   relationship template. For example, \code{list(item = ~ b.precision *
#'   b.vividness)} expands \code{item} into an intercept and three modulated
#'   columns. Templates omitted from this list are unchanged.
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
#'     \item{block_var_b}{Block labels for domain B (equal to
#'       \code{block_var_a} in within mode); used by
#'       \code{\link{permute_labels}}.}
#'   }
#'
#' @seealso \code{\link{rsa_design}}, \code{\link{rsa_model}},
#'   \code{\link{model_space_connectivity}}
#' @details Formula metadata exposes \code{a.row}/\code{b.row} (observation
#'   positions), \code{a.item}/\code{b.item} (item IDs), and feature columns
#'   prefixed with \code{a.}/\code{b.}. Repeated item IDs retain separate
#'   observations. Within-domain formulas must be invariant to swapping the
#'   two sides; declare a symmetric rule such as
#'   \code{~ I((a.vividness + b.vividness) / 2)}. Formula evaluation retains
#'   missing rows; coefficient fitting uses complete eligible pairs.
#'   Formula transformations operate on pair metadata before masking. To
#'   standardize a trial attribute over observations rather than pairs,
#'   standardize it in \code{features_a}/\code{features_b} first. Background
#'   columns must remain identifiable on the complete eligible pairs.
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
#'
#' # Item correspondence modulated by retrieval attributes.
#' # Neural observations occupy encoding rows 1:8 and retrieval rows 9:16.
#' retrieval <- data.frame(precision = rnorm(8), vividness = rnorm(8))
#' relational <- pair_rsa_design(
#'   items_a = 1:8, items_b = 1:8, pairs = "between",
#'   row_idx_a = 1:8, row_idx_b = 9:16, features_b = retrieval,
#'   model = list(item = function(a, b) as.numeric(a == b)),
#'   modulation = list(item = ~ b.precision * b.vividness),
#'   nuisance = ~ factor(b.row)
#' )
#' relational$model_predictors
#' # Pass to rsa_model(dataset, relational, distmethod = "pearson",
#' #   measure = "similarity", regtype = "lm", statistic = "beta").
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
                            row_idx_b = NULL,
                            modulation = NULL) {
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

  nuisance_formula <- if (inherits(nuisance, "formula")) nuisance else NULL
  if ((!is.list(model) && !is.null(model)) ||
      (is.null(nuisance_formula) && !is.list(nuisance) && !is.null(nuisance))) {
    stop("`model` and `nuisance` must be (possibly empty) named lists.",
         call. = FALSE)
  }
  model <- if (is.null(model)) list() else as.list(model)
  nuisance <- if (is.null(nuisance) || !is.null(nuisance_formula)) list() else as.list(nuisance)
  if (!is.null(modulation) &&
      (!is.list(modulation) || (length(modulation) && (is.null(names(modulation)) ||
       anyNA(names(modulation)) || any(!nzchar(names(modulation))) ||
       anyDuplicated(names(modulation)) || any(!names(modulation) %in% names(model)))))) {
    stop("`modulation` must be a named list keyed by unique names in `model`.", call. = FALSE)
  }
  if (length(model) == 0L && length(nuisance) == 0L && is.null(nuisance_formula)) {
    stop("Provide at least one entry in `model` or `nuisance`.", call. = FALSE)
  }

  if (length(model) > 0L && (is.null(names(model)) || anyNA(names(model)) || any(!nzchar(names(model))))) {
    stop("`model` entries must be named.", call. = FALSE)
  }
  if (length(nuisance) > 0L && (is.null(names(nuisance)) || anyNA(names(nuisance)) || any(!nzchar(names(nuisance))))) {
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

  formula_terms <- list()
  if (length(modulation) || !is.null(nuisance_formula)) {
    metadata <- .pair_design_metadata(pair_index, features_a, features_b)
    swapped <- if (identical(pairs, "within")) {
      reverse_index <- pair_index
      reverse_index$i <- pair_index$j
      reverse_index$j <- pair_index$i
      reverse_index$item_a <- pair_index$item_b
      reverse_index$item_b <- pair_index$item_a
      .pair_design_metadata(reverse_index, features_a, features_b)
    } else NULL
    expanded <- list()
    for (nm in names(vec_model)) {
      if (!nm %in% names(modulation)) {
        expanded <- c(expanded, stats::setNames(list(vec_model[[nm]]), nm))
        next
      }
      basis <- .pair_design_formula(modulation[[nm]], metadata, swapped,
                                    paste0("modulation for '", nm, "'"))
      columns <- ifelse(colnames(basis) == "(Intercept)", nm,
                        paste(nm, colnames(basis), sep = "."))
      columns <- make.names(sanitize(columns), unique = FALSE)
      expanded <- c(expanded, stats::setNames(lapply(seq_len(ncol(basis)), function(j) {
        as.numeric(vec_model[[nm]] * basis[, j])
      }), columns))
      formula_terms[[nm]] <- list(formula = modulation[[nm]], columns = columns,
                                  assign = attr(basis, "assign"),
                                  contrasts = attr(basis, "contrasts"))
    }
    vec_model <- expanded
    if (!is.null(nuisance_formula)) {
      basis <- .pair_design_formula(nuisance_formula, metadata, swapped, "nuisance")
      assignment <- attr(basis, "assign")
      basis <- basis[, assignment != 0L, drop = FALSE]
      columns <- make.names(sanitize(paste0("background.", colnames(basis))))
      vec_nuis <- stats::setNames(lapply(seq_len(ncol(basis)), function(j) as.numeric(basis[, j])), columns)
      nuisance_terms <- list(formula = nuisance_formula, columns = columns,
                             assign = assignment[assignment != 0L])
    }
  }

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
  if (!length(model_mat_raw)) {
    stop("Provide at least one non-intercept model or nuisance predictor.", call. = FALSE)
  }
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
  if (length(modulation) || !is.null(nuisance_formula)) {
    .rsa_coefficient_query(model_mat, model_names, validate_only = TRUE)
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
    block_var_b  = block_var_b,
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
    row_idx_b    = row_idx_b,
    templates    = model,
    modulation   = modulation,
    nuisance_formula = nuisance_formula,
    formula_terms = formula_terms,
    nuisance_terms = if (!is.null(nuisance_formula)) nuisance_terms else NULL
  )
  class(des) <- c("pair_rsa_design", "rsa_design", "list")
  des
}

# Pair metadata is built only for formula specifications. Item identity and
# observation identity are deliberately separate, even when item IDs repeat.
.pair_design_metadata <- function(index, features_a, features_b) {
  out <- data.frame(a.row = index$i, b.row = index$j,
                    a.item = index$item_a, b.item = index$item_b)
  for (side in c("a", "b")) {
    features <- if (side == "a") features_a else features_b
    if (is.null(features)) next
    nms <- colnames(features)
    if (is.null(nms) || anyNA(nms) || any(!nzchar(nms)) || anyDuplicated(nms) ||
        any(nms %in% c("row", "item"))) {
      stop("Formula feature columns must have unique names other than 'row' and 'item'.",
           call. = FALSE)
    }
    rows <- if (side == "a") index$i else index$j
    values <- as.data.frame(features[rows, , drop = FALSE])
    names(values) <- paste0(side, ".", nms)
    out <- cbind(out, values)
  }
  out
}

.pair_design_formula <- function(formula, metadata, swapped = NULL, label) {
  if (!inherits(formula, "formula") || length(formula) != 2L) {
    stop(label, " must be a right-hand-side formula.", call. = FALSE)
  }
  missing <- setdiff(all.vars(formula), names(metadata))
  if (length(missing)) {
    stop(label, " references unknown pair metadata: ", paste(missing, collapse = ", "),
         ". Use a. or b. prefixes for feature columns.", call. = FALSE)
  }
  evaluate <- function(data) {
    frame <- stats::model.frame(formula, data, na.action = stats::na.pass)
    matrix <- stats::model.matrix(formula, frame)
    if (nrow(matrix) != nrow(data)) {
      stop(label, " lost pair rows during formula evaluation.", call. = FALSE)
    }
    matrix
  }
  basis <- evaluate(metadata)
  if (!is.null(swapped)) {
    reverse <- evaluate(swapped)
    if (!identical(colnames(basis), colnames(reverse)) ||
        !isTRUE(all.equal(unname(basis), unname(reverse), check.attributes = FALSE,
                          tolerance = 1e-12))) {
      stop(label, " must be symmetric under exchanging a and b for within-domain pairs.",
           call. = FALSE)
    }
  }
  basis
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
