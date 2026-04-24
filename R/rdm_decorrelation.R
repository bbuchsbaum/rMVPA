#' Decorrelation of correlated model RDMs
#'
#' Estimates layer- or model-specific RDM score vectors by removing shared
#' representational structure across a set of model RDMs. This is intended as a
#' preprocessing step before RSA when several model RDMs are strongly correlated
#' and the goal is to evaluate their more distinct representational components.
#'
#' @details
#' This function works in lower-triangular RDM-vector space. With
#' `method = "ordered_innovation"`, RDMs are processed in the order supplied:
#' each RDM vector is decomposed into the component predicted by earlier
#' innovation vectors and a residual innovation. A shrinkage parameter
#' `gamma` controls how much of the predicted/shared component is retained.
#'
#' The default objective is constrained preservation: preserve as much of each
#' original standardized RDM vector as possible while forcing the mean pairwise
#' association among adjusted RDM vectors below `epsilon` when feasible.
#'
#' With `method = "zca"`, the standardized RDM vectors are symmetrically
#' whitened. This gives exact joint decorrelation in the Pearson sense when the
#' RDM-vector covariance matrix is full rank, but each output RDM is a linear
#' mixture of the original RDMs.
#'
#' Ordinary identity shrinkage of a correlation matrix does not change
#' rank-based RSA correlations among off-diagonal RDM vectors, except at the
#' degenerate endpoint. This function instead shrinks shared RDM-vector
#' components.
#'
#' @param rdms Named list of square symmetric matrices or `dist` objects. All
#'   RDMs must have the same size and contain finite lower-triangular values.
#' @param method Character string, either `"ordered_innovation"` or `"zca"`.
#' @param similarity Character string. `"spearman"` rank-standardizes each RDM
#'   vector before decorrelation. `"pearson"` uses centered/scaled numeric RDM
#'   values directly.
#' @param objective Character string. `"constraint"` maximizes preservation
#'   subject to the requested decorrelation tolerance when feasible.
#'   `"penalty"` minimizes preservation loss plus `mu` times the decorrelation
#'   penalty. `"full_residual"` uses the fully residualized ordered innovations
#'   and ignores `epsilon`.
#' @param epsilon Non-negative tolerance for the mean pairwise adjusted-RDM
#'   association under `objective = "constraint"`.
#' @param gamma_grid Numeric vector in `[0, 1]` used for ordered innovation
#'   shrinkage. The first RDM is fixed at `gamma = 1`.
#' @param mu Non-negative penalty multiplier used when `objective = "penalty"`.
#' @param target Character string. `"mean_abs"` penalizes mean absolute
#'   pairwise correlation; `"mean_squared"` penalizes mean squared pairwise
#'   correlation.
#' @param return Character string. `"dist"` returns adjusted RDMs as `dist`
#'   objects, `"matrix"` returns symmetric matrices, and `"vector"` returns
#'   lower-triangular vectors.
#' @param project_correlation Logical. If `TRUE`, each adjusted RDM-score matrix
#'   is interpreted as `1 - correlation`, projected to the nearest positive
#'   semidefinite correlation matrix with unit diagonal using `Matrix::nearPD`,
#'   and converted back to an RDM. This is off by default because projection can
#'   reintroduce RDM correlations.
#' @param zca_ridge Numeric ridge added only when the RDM-vector covariance is
#'   rank deficient under `method = "zca"`.
#' @param max_combinations Maximum number of exhaustive `gamma_grid`
#'   combinations for ordered innovation. Larger problems use coordinate search.
#' @param tol Numerical tolerance for degeneracy and rank checks.
#'
#' @return An object of class `rdm_decorrelation_result` with adjusted RDMs,
#'   lower-triangular vectors, selected shrinkage parameters, and diagnostics.
#' @export
#'
#' @examples
#' set.seed(1)
#' X <- matrix(rnorm(12 * 4), 12, 4)
#' r1 <- dist(X)
#' r2 <- dist(X + matrix(rnorm(12 * 4, sd = 0.2), 12, 4))
#' dec <- rdm_decorrelate(list(layer1 = r1, layer2 = r2), epsilon = 0.05)
#' dec$mean_abs_cor_after
rdm_decorrelate <- function(rdms,
                            method = c("ordered_innovation", "zca"),
                            similarity = c("spearman", "pearson"),
                            objective = c("constraint", "penalty", "full_residual"),
                            epsilon = 0.05,
                            gamma_grid = seq(0, 1, by = 0.05),
                            mu = 1,
                            target = c("mean_abs", "mean_squared"),
                            return = c("dist", "matrix", "vector"),
                            project_correlation = FALSE,
                            zca_ridge = 1e-8,
                            max_combinations = 2e5,
                            tol = 1e-10) {
  method <- match.arg(method)
  similarity <- match.arg(similarity)
  objective <- match.arg(objective)
  target <- match.arg(target)
  return <- match.arg(return)

  parsed <- .rdm_decorrelation_parse_rdms(rdms)
  Y <- .rdm_decorrelation_standardize(parsed$vectors, similarity = similarity, tol = tol)

  if (!is.numeric(epsilon) || length(epsilon) != 1L || !is.finite(epsilon) || epsilon < 0) {
    stop("`epsilon` must be a single non-negative finite number.", call. = FALSE)
  }
  if (!is.numeric(mu) || length(mu) != 1L || !is.finite(mu) || mu < 0) {
    stop("`mu` must be a single non-negative finite number.", call. = FALSE)
  }

  fit <- switch(method,
    ordered_innovation = .rdm_decorrelation_ordered(
      Y = Y,
      objective = objective,
      epsilon = epsilon,
      gamma_grid = gamma_grid,
      mu = mu,
      target = target,
      max_combinations = max_combinations,
      tol = tol
    ),
    zca = .rdm_decorrelation_zca(Y = Y, zca_ridge = zca_ridge, tol = tol)
  )

  X <- fit$X
  dimnames(X) <- dimnames(Y)

  mats <- .rdm_decorrelation_vectors_to_matrices(
    X,
    n = parsed$n,
    labels = parsed$labels,
    names = parsed$names
  )
  if (isTRUE(project_correlation)) {
    mats <- lapply(mats, .rdm_decorrelation_project_correlation)
    X <- .rdm_decorrelation_matrix_list_to_vectors(mats)
  }

  out_rdms <- switch(return,
    vector = X,
    matrix = mats,
    dist = lapply(mats, stats::as.dist)
  )
  if (return != "vector") {
    names(out_rdms) <- parsed$names
  }

  before <- .rdm_decorrelation_cor_matrix(Y)
  after <- .rdm_decorrelation_cor_matrix(X)
  preservation <- .rdm_decorrelation_col_cor(X, Y)

  structure(
    list(
      rdms = out_rdms,
      vectors = X,
      original_vectors = parsed$vectors,
      standardized_vectors = Y,
      gamma = fit$gamma,
      method = method,
      similarity = similarity,
      objective = objective,
      epsilon = epsilon,
      target = target,
      preservation = preservation,
      cross_rdm_cor_before = before,
      cross_rdm_cor_after = after,
      mean_abs_cor_before = .rdm_decorrelation_mean_abs_cor(before),
      mean_abs_cor_after = .rdm_decorrelation_mean_abs_cor(after),
      mean_squared_cor_before = .rdm_decorrelation_mean_sq_cor(before),
      mean_squared_cor_after = .rdm_decorrelation_mean_sq_cor(after),
      projected = isTRUE(project_correlation),
      diagnostics = fit$diagnostics
    ),
    class = c("rdm_decorrelation_result", "list")
  )
}

#' @export
print.rdm_decorrelation_result <- function(x, ...) {
  cat("RDM Decorrelation Result\n")
  cat("  method:              ", x$method, "\n", sep = "")
  cat("  similarity:          ", x$similarity, "\n", sep = "")
  cat("  mean |cor| before:   ", format(signif(x$mean_abs_cor_before, 4)), "\n", sep = "")
  cat("  mean |cor| after:    ", format(signif(x$mean_abs_cor_after, 4)), "\n", sep = "")
  cat("  mean preservation:   ", format(signif(mean(x$preservation, na.rm = TRUE), 4)), "\n", sep = "")
  if (!is.null(x$gamma)) {
    cat("  gamma:               ", paste(format(signif(x$gamma, 4)), collapse = ", "), "\n", sep = "")
  }
  invisible(x)
}

#' ROI-by-ROI similarity through a correlated model-RDM space
#'
#' Projects ROI RDM vectors into the subspace spanned by one or more model RDMs,
#' decorrelates that model space, and compares ROIs through their model-space
#' fingerprints.
#'
#' @details
#' This helper implements a projection-based second-order RSA summary for the
#' common situation where several model RDMs are correlated. Let `Y` be the
#' lower-triangle ROI RDM matrix, with RDM cells in rows and ROIs in columns,
#' and let `X` be the corresponding matrix of model RDM vectors. After optional
#' rank transformation and column standardization, the function builds an
#' orthonormal basis `Q` for the column space of `X` and computes
#'
#' \deqn{F = Y' Q / \sqrt{p - 1}}
#'
#' where `p` is the number of retained RDM cells. Rows of `F` are decorrelated
#' ROI model-space fingerprints. The strength-sensitive ROI-by-ROI similarity is
#'
#' \deqn{S_{model} = F F'}
#'
#' which is equivalent to \eqn{Y' P_X Y / (p - 1)}. The profile-only similarity
#' is the cosine similarity between rows of `F`.
#'
#' With `basis = "pca"`, the first component is often the common model axis
#' when all model RDMs are positively correlated, and later components describe
#' differentiating model contrasts. With `basis = "qr"`, the axes are an
#' orthonormal basis generated by QR decomposition of the supplied model order.
#'
#' @param roi_rdms ROI RDM vectors. Accepted forms are a matrix with RDM cells
#'   in rows and ROIs in columns, a named list of numeric vectors, a tibble/data
#'   frame with columns `roinum` and `rdm_vec` as returned by
#'   \code{\link{feature_rsa_rdm_vectors}}, or a `regional_mvpa_result`
#'   containing stored feature-RSA RDM vectors.
#' @param model_rdms Named list of square symmetric model RDM matrices or
#'   `dist` objects, or a numeric matrix with RDM cells in rows and model RDMs
#'   in columns.
#' @param method Similarity scale used before projection. `"pearson"` centers
#'   and scales numeric RDM values. `"spearman"` rank-transforms each RDM vector
#'   before centering and scaling.
#' @param basis Basis for the model subspace. `"pca"` orders orthogonal axes by
#'   variance in the model-RDM correlation matrix. `"qr"` uses a QR basis tied
#'   to the supplied model order.
#' @param use Missing-value policy. `"complete.obs"` removes RDM cells with any
#'   non-finite model or ROI value before fitting the subspace. `"everything"`
#'   requires all model and ROI vectors to be finite.
#' @param tol Numerical tolerance for dropping near-zero model-space dimensions
#'   and detecting zero-norm ROI fingerprints.
#' @param return_projected Logical; if `TRUE`, include the projected ROI RDM
#'   vectors and residual ROI RDM vectors in the returned object.
#'
#' @return An object of class `rdm_model_space_connectivity`, a list with:
#' \describe{
#'   \item{similarity}{Strength-sensitive ROI-by-ROI similarity through the
#'     full model-RDM subspace.}
#'   \item{profile_similarity}{Cosine similarity between ROI model-space
#'     fingerprints, ignoring overall model-alignment strength.}
#'   \item{component_similarity}{Named list of rank-1 ROI-by-ROI matrices, one
#'     per orthogonal model-space axis.}
#'   \item{common_similarity}{The first component matrix, useful for the shared
#'     model axis under `basis = "pca"`.}
#'   \item{difference_similarity}{Sum of all components after the first.}
#'   \item{raw_similarity}{Second-order ROI similarity before model projection.}
#'   \item{residual_similarity}{ROI similarity outside the model-RDM subspace.}
#'   \item{scores}{ROI-by-axis fingerprint matrix `F`.}
#'   \item{raw_model_scores}{ROI-by-original-model correlation matrix.}
#'   \item{model_axis_cor}{Correlation of orthogonal model axes with original
#'     model RDMs, used to interpret the axes.}
#' }
#' @export
#'
#' @examples
#' set.seed(1)
#' n <- 8
#' p <- n * (n - 1) / 2
#' model_vecs <- matrix(rnorm(p * 4), p, 4)
#' colnames(model_vecs) <- paste0("A", 1:4)
#' roi_vecs <- model_vecs %*% matrix(rnorm(4 * 5), 4, 5) +
#'   matrix(rnorm(p * 5, sd = 0.2), p, 5)
#' colnames(roi_vecs) <- paste0("ROI", 1:5)
#'
#' conn <- rdm_model_space_connectivity(roi_vecs, model_vecs)
#' conn$similarity
rdm_model_space_connectivity <- function(roi_rdms,
                                         model_rdms,
                                         method = c("pearson", "spearman"),
                                         basis = c("pca", "qr"),
                                         use = c("complete.obs", "everything"),
                                         tol = 1e-8,
                                         return_projected = FALSE) {
  method <- match.arg(method)
  basis <- match.arg(basis)
  use <- match.arg(use)

  if (!is.numeric(tol) || length(tol) != 1L || !is.finite(tol) || tol <= 0) {
    stop("`tol` must be a positive finite scalar.", call. = FALSE)
  }
  if (!is.logical(return_projected) || length(return_projected) != 1L || is.na(return_projected)) {
    stop("`return_projected` must be TRUE or FALSE.", call. = FALSE)
  }

  roi <- .rdm_model_space_parse_roi_rdms(roi_rdms)
  model <- .rdm_model_space_parse_model_rdms(model_rdms)

  if (nrow(roi$vectors) != nrow(model$vectors)) {
    stop(
      "`roi_rdms` and `model_rdms` must have the same lower-triangle vector length.",
      call. = FALSE
    )
  }

  all_vecs <- cbind(model$vectors, roi$vectors)
  finite_rows <- stats::complete.cases(all_vecs)
  if (identical(use, "everything") && !all(finite_rows)) {
    stop("All model and ROI RDM vector entries must be finite when `use = \"everything\"`.", call. = FALSE)
  }
  if (!all(finite_rows)) {
    model$vectors <- model$vectors[finite_rows, , drop = FALSE]
    roi$vectors <- roi$vectors[finite_rows, , drop = FALSE]
  }

  p <- nrow(model$vectors)
  if (p < 2L) {
    stop("At least two finite RDM cells are required after missing-value handling.", call. = FALSE)
  }

  X <- .rdm_model_space_standardize_matrix(model$vectors, method = method, tol = tol, label = "model_rdms")
  Y <- .rdm_model_space_standardize_matrix(roi$vectors, method = method, tol = tol, label = "roi_rdms")

  raw_model_scores <- crossprod(Y, X) / (p - 1L)
  dimnames(raw_model_scores) <- list(roi$labels, model$labels)

  basis_fit <- .rdm_model_space_basis(X, basis = basis, tol = tol)
  Q <- basis_fit$Q
  axis_names <- basis_fit$axis_names
  colnames(Q) <- axis_names

  scores <- crossprod(Y, Q) / sqrt(p - 1L)
  dimnames(scores) <- list(roi$labels, axis_names)

  similarity <- tcrossprod(scores)
  dimnames(similarity) <- list(roi$labels, roi$labels)

  raw_similarity <- crossprod(Y) / (p - 1L)
  dimnames(raw_similarity) <- list(roi$labels, roi$labels)

  projected <- Q %*% crossprod(Q, Y)
  residual <- Y - projected
  residual_similarity <- crossprod(residual) / (p - 1L)
  dimnames(residual_similarity) <- list(roi$labels, roi$labels)

  component_similarity <- lapply(seq_len(ncol(scores)), function(i) {
    mat <- tcrossprod(scores[, i, drop = FALSE])
    dimnames(mat) <- list(roi$labels, roi$labels)
    mat
  })
  names(component_similarity) <- axis_names

  norms <- sqrt(rowSums(scores^2))
  profile_scores <- scores
  keep <- is.finite(norms) & norms > tol
  profile_scores[keep, ] <- profile_scores[keep, , drop = FALSE] / norms[keep]
  profile_scores[!keep, ] <- NA_real_
  profile_similarity <- tcrossprod(profile_scores)
  dimnames(profile_similarity) <- list(roi$labels, roi$labels)

  model_axis_cor <- .rdm_model_space_col_cor(Q, X)
  dimnames(model_axis_cor) <- list(axis_names, model$labels)

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
    raw_similarity = raw_similarity,
    residual_similarity = residual_similarity,
    scores = scores,
    raw_model_scores = raw_model_scores,
    model_axis_cor = model_axis_cor,
    model_cor = stats::cor(X),
    eigenvalues = basis_fit$eigenvalues,
    rank = ncol(Q),
    basis = basis,
    method = method,
    n_pairs_total = nrow(all_vecs),
    n_pairs_used = p,
    roi_labels = roi$labels,
    model_labels = model$labels
  )

  if (isTRUE(return_projected)) {
    colnames(projected) <- roi$labels
    colnames(residual) <- roi$labels
    out$projected_vectors <- projected
    out$residual_vectors <- residual
  }

  structure(out, class = c("rdm_model_space_connectivity", "list"))
}

#' @export
print.rdm_model_space_connectivity <- function(x, ...) {
  cat("RDM Model-Space Connectivity\n")
  cat("  method:       ", x$method, "\n", sep = "")
  cat("  basis:        ", x$basis, "\n", sep = "")
  cat("  ROIs:         ", length(x$roi_labels), "\n", sep = "")
  cat("  model RDMs:   ", length(x$model_labels), "\n", sep = "")
  cat("  rank:         ", x$rank, "\n", sep = "")
  cat("  RDM cells:    ", x$n_pairs_used, "/", x$n_pairs_total, "\n", sep = "")
  invisible(x)
}

.rdm_model_space_parse_roi_rdms <- function(roi_rdms) {
  if (is.data.frame(roi_rdms) && "rdm_vec" %in% names(roi_rdms)) {
    if ("observation_index" %in% names(roi_rdms)) {
      .feature_rsa_validate_obs_order(
        roi_rdms$observation_index,
        fn_label = "rdm_model_space_connectivity"
      )
    }
    lens <- vapply(roi_rdms$rdm_vec, length, integer(1))
    if (length(unique(lens)) != 1L) {
      stop("`roi_rdms$rdm_vec` entries must all have the same length.", call. = FALSE)
    }
    labels <- if ("roinum" %in% names(roi_rdms)) {
      as.character(roi_rdms$roinum)
    } else {
      as.character(seq_len(nrow(roi_rdms)))
    }
    mat <- .feature_rsa_bind_numeric_columns(lapply(roi_rdms$rdm_vec, as.numeric))
    colnames(mat) <- labels
    return(list(vectors = mat, labels = labels))
  }

  if (inherits(roi_rdms, "regional_mvpa_result")) {
    return(.rdm_model_space_parse_roi_rdms(feature_rsa_rdm_vectors(roi_rdms)))
  }

  if (is.matrix(roi_rdms) || is.data.frame(roi_rdms)) {
    mat <- as.matrix(roi_rdms)
    if (!is.numeric(mat)) {
      stop("Matrix `roi_rdms` must be numeric.", call. = FALSE)
    }
    labels <- colnames(mat)
    if (is.null(labels) || any(!nzchar(labels))) {
      labels <- paste0("ROI", seq_len(ncol(mat)))
    }
    colnames(mat) <- labels
    return(list(vectors = mat, labels = labels))
  }

  if (is.list(roi_rdms) && length(roi_rdms) > 0L) {
    lens <- vapply(roi_rdms, length, integer(1))
    if (length(unique(lens)) != 1L) {
      stop("ROI RDM vectors must all have the same length.", call. = FALSE)
    }
    mat <- .feature_rsa_bind_numeric_columns(lapply(roi_rdms, as.numeric))
    labels <- names(roi_rdms)
    if (is.null(labels) || any(!nzchar(labels))) {
      labels <- paste0("ROI", seq_along(roi_rdms))
    }
    colnames(mat) <- labels
    return(list(vectors = mat, labels = labels))
  }

  stop(
    "`roi_rdms` must be a matrix, list of numeric vectors, RDM-vector tibble, or regional_mvpa_result.",
    call. = FALSE
  )
}

.rdm_model_space_parse_model_rdms <- function(model_rdms) {
  if (is.matrix(model_rdms) || is.data.frame(model_rdms)) {
    mat <- as.matrix(model_rdms)
    if (!is.numeric(mat)) {
      stop("Matrix `model_rdms` must be numeric.", call. = FALSE)
    }
    labels <- colnames(mat)
    if (is.null(labels) || any(!nzchar(labels))) {
      labels <- paste0("model", seq_len(ncol(mat)))
    }
    colnames(mat) <- labels
    return(list(vectors = mat, labels = labels))
  }

  parsed <- .rdm_decorrelation_parse_rdms(model_rdms)
  list(vectors = parsed$vectors, labels = colnames(parsed$vectors))
}

.rdm_model_space_standardize_matrix <- function(mat, method, tol, label) {
  mat <- as.matrix(mat)
  if (identical(method, "spearman")) {
    mat <- apply(mat, 2L, rank, ties.method = "average")
  }
  out <- apply(mat, 2L, function(x) {
    x <- as.numeric(x)
    x <- x - mean(x)
    s <- stats::sd(x)
    if (!is.finite(s) || s <= tol) {
      stop(label, " contains a constant or degenerate RDM vector.", call. = FALSE)
    }
    x / s
  })
  out <- as.matrix(out)
  dimnames(out) <- dimnames(mat)
  out
}

.rdm_model_space_basis <- function(X, basis, tol) {
  if (identical(basis, "pca")) {
    G <- crossprod(X) / (nrow(X) - 1L)
    eg <- eigen(G, symmetric = TRUE)
    keep <- eg$values > tol
    if (!any(keep)) {
      stop("Model RDM space has no non-degenerate dimensions.", call. = FALSE)
    }
    values <- eg$values[keep]
    vectors <- eg$vectors[, keep, drop = FALSE]
    Q <- X %*% vectors %*% diag(1 / sqrt((nrow(X) - 1L) * values), nrow = length(values))
    axis_names <- paste0("PC", seq_along(values))
    list(Q = Q, eigenvalues = values, axis_names = axis_names)
  } else {
    qrx <- qr(X, tol = tol)
    rank_x <- qrx$rank
    if (rank_x < 1L) {
      stop("Model RDM space has no non-degenerate dimensions.", call. = FALSE)
    }
    Q <- qr.Q(qrx, complete = FALSE)[, seq_len(rank_x), drop = FALSE]
    axis_names <- paste0("QR", seq_len(rank_x))
    list(Q = Q, eigenvalues = rep(NA_real_, rank_x), axis_names = axis_names)
  }
}

.rdm_model_space_col_cor <- function(A, B) {
  denom <- sqrt(colSums(A^2)) %o% sqrt(colSums(B^2))
  out <- crossprod(A, B) / denom
  out[!is.finite(out)] <- NA_real_
  out
}

.rdm_decorrelation_parse_rdms <- function(rdms) {
  if (!is.list(rdms) || length(rdms) < 2L) {
    stop("`rdms` must be a list of at least two matrices or dist objects.", call. = FALSE)
  }
  nm <- names(rdms)
  if (is.null(nm) || any(!nzchar(nm))) {
    nm <- paste0("rdm", seq_along(rdms))
  }

  mats <- lapply(rdms, function(x) {
    if (inherits(x, "dist")) {
      as.matrix(x)
    } else if (is.matrix(x) || is.data.frame(x)) {
      as.matrix(x)
    } else {
      stop("Each RDM must be a square matrix or a dist object.", call. = FALSE)
    }
  })

  dims <- vapply(mats, nrow, integer(1))
  if (any(vapply(mats, ncol, integer(1)) != dims) || length(unique(dims)) != 1L) {
    stop("All RDMs must be square matrices with the same dimensions.", call. = FALSE)
  }
  n <- dims[[1]]
  if (n < 3L) {
    stop("RDMs must contain at least three items.", call. = FALSE)
  }

  labels <- rownames(mats[[1]])
  if (is.null(labels)) labels <- as.character(seq_len(n))

  for (i in seq_along(mats)) {
    M <- mats[[i]]
    if (!is.numeric(M)) {
      stop("All RDMs must be numeric.", call. = FALSE)
    }
    if (!isTRUE(all.equal(M, t(M), tolerance = sqrt(.Machine$double.eps), check.attributes = FALSE))) {
      stop("All matrix RDMs must be symmetric.", call. = FALSE)
    }
    vals <- M[lower.tri(M)]
    if (any(!is.finite(vals))) {
      stop("All lower-triangular RDM values must be finite.", call. = FALSE)
    }
    if (!is.null(rownames(M)) && !identical(rownames(M), labels)) {
      stop("RDM row names must match across inputs when supplied.", call. = FALSE)
    }
    if (!is.null(colnames(M)) && !identical(colnames(M), labels)) {
      stop("RDM column names must match row names and match across inputs.", call. = FALSE)
    }
  }

  V <- vapply(mats, function(M) M[lower.tri(M)], numeric(n * (n - 1L) / 2L))
  colnames(V) <- nm
  list(mats = mats, vectors = V, n = n, labels = labels, names = nm)
}

.rdm_decorrelation_standardize <- function(V, similarity, tol) {
  Y <- V
  if (identical(similarity, "spearman")) {
    Y <- apply(V, 2L, function(x) rank(x, ties.method = "average"))
  }
  Y <- apply(Y, 2L, .rdm_decorrelation_zscore, tol = tol)
  Y <- as.matrix(Y)
  colnames(Y) <- colnames(V)
  Y
}

.rdm_decorrelation_zscore <- function(x, tol = 1e-10) {
  x <- as.numeric(x)
  x <- x - mean(x)
  s <- stats::sd(x)
  if (!is.finite(s) || s <= tol) {
    stop("RDM vectors must be non-constant after the requested transformation.", call. = FALSE)
  }
  x / s
}

.rdm_decorrelation_ordered <- function(Y,
                                       objective,
                                       epsilon,
                                       gamma_grid,
                                       mu,
                                       target,
                                       max_combinations,
                                       tol) {
  L <- ncol(Y)
  dec <- .rdm_decorrelation_ordered_decompose(Y, tol = tol)
  E <- dec$E
  P <- dec$P

  if (identical(objective, "full_residual")) {
    gamma <- c(1, rep(0, L - 1L))
    X <- .rdm_decorrelation_ordered_apply(E, P, gamma, tol = tol)
    return(list(
      X = X,
      gamma = stats::setNames(gamma, colnames(Y)),
      diagnostics = c(dec$diagnostics, list(search = "full_residual"))
    ))
  }

  gamma_grid <- sort(unique(as.numeric(gamma_grid)))
  gamma_grid <- gamma_grid[is.finite(gamma_grid) & gamma_grid >= 0 & gamma_grid <= 1]
  if (!length(gamma_grid)) {
    stop("`gamma_grid` must contain at least one finite value in [0, 1].", call. = FALSE)
  }
  if (!1 %in% gamma_grid) gamma_grid <- sort(unique(c(gamma_grid, 1)))
  if (!0 %in% gamma_grid) gamma_grid <- sort(unique(c(gamma_grid, 0)))

  grids <- c(list(1), rep(list(gamma_grid), L - 1L))
  ncomb <- prod(vapply(grids, length, integer(1)))
  if (ncomb <= max_combinations) {
    opt <- .rdm_decorrelation_ordered_exhaustive(
      E = E,
      P = P,
      Y = Y,
      grids = grids,
      objective = objective,
      epsilon = epsilon,
      mu = mu,
      target = target,
      tol = tol
    )
    search <- "exhaustive"
  } else {
    opt <- .rdm_decorrelation_ordered_coordinate(
      E = E,
      P = P,
      Y = Y,
      gamma_grid = gamma_grid,
      objective = objective,
      epsilon = epsilon,
      mu = mu,
      target = target,
      tol = tol
    )
    search <- "coordinate"
  }
  opt$gamma <- stats::setNames(opt$gamma, colnames(Y))
  opt$diagnostics <- c(dec$diagnostics, opt$diagnostics, list(search = search))
  opt
}

.rdm_decorrelation_ordered_decompose <- function(Y, tol) {
  E <- matrix(0, nrow(Y), ncol(Y), dimnames = dimnames(Y))
  P <- matrix(0, nrow(Y), ncol(Y), dimnames = dimnames(Y))
  degenerate <- rep(FALSE, ncol(Y))
  E[, 1L] <- Y[, 1L]
  for (j in seq_len(ncol(Y))[-1L]) {
    prev <- E[, seq_len(j - 1L), drop = FALSE]
    qr_prev <- qr(prev, tol = tol)
    if (qr_prev$rank > 0L) {
      Q <- qr.Q(qr_prev)[, seq_len(qr_prev$rank), drop = FALSE]
      P[, j] <- Q %*% crossprod(Q, Y[, j])
    }
    E[, j] <- Y[, j] - P[, j]
    degenerate[j] <- stats::sd(E[, j]) <= tol
  }
  list(E = E, P = P, diagnostics = list(degenerate_innovation = stats::setNames(degenerate, colnames(Y))))
}

.rdm_decorrelation_ordered_apply <- function(E, P, gamma, tol) {
  X <- E
  deg <- rep(FALSE, ncol(E))
  for (j in seq_len(ncol(E))) {
    x <- E[, j] + gamma[[j]] * P[, j]
    s <- stats::sd(x)
    if (!is.finite(s) || s <= tol) {
      X[, j] <- 0
      deg[j] <- TRUE
    } else {
      X[, j] <- (x - mean(x)) / s
    }
  }
  attr(X, "degenerate_adjusted") <- deg
  X
}

.rdm_decorrelation_ordered_exhaustive <- function(E, P, Y, grids, objective, epsilon, mu, target, tol) {
  combos <- expand.grid(grids, KEEP.OUT.ATTRS = FALSE)
  best <- NULL
  best_score <- Inf
  best_feasible <- NULL
  best_feasible_pres <- -Inf
  best_feasible_dec <- Inf

  for (i in seq_len(nrow(combos))) {
    gamma <- as.numeric(combos[i, ])
    X <- .rdm_decorrelation_ordered_apply(E, P, gamma, tol = tol)
    score <- .rdm_decorrelation_score(X, Y, objective, epsilon, mu, target)
    if (identical(objective, "constraint") && isTRUE(score$feasible)) {
      if (score$preservation > best_feasible_pres ||
          (isTRUE(all.equal(score$preservation, best_feasible_pres)) && score$decorrelation < best_feasible_dec)) {
        best_feasible <- list(X = X, gamma = gamma, diagnostics = score)
        best_feasible_pres <- score$preservation
        best_feasible_dec <- score$decorrelation
      }
    }
    if (score$objective < best_score) {
      best <- list(X = X, gamma = gamma, diagnostics = score)
      best_score <- score$objective
    }
  }
  if (!is.null(best_feasible)) best_feasible else best
}

.rdm_decorrelation_ordered_coordinate <- function(E, P, Y, gamma_grid, objective, epsilon, mu, target, tol) {
  gamma <- c(1, rep(1, ncol(Y) - 1L))
  for (iter in seq_len(5L)) {
    changed <- FALSE
    for (j in seq_len(ncol(Y))[-1L]) {
      candidates <- lapply(gamma_grid, function(g) {
        gg <- gamma
        gg[[j]] <- g
        X <- .rdm_decorrelation_ordered_apply(E, P, gg, tol = tol)
        sc <- .rdm_decorrelation_score(X, Y, objective, epsilon, mu, target)
        list(gamma = gg, X = X, score = sc)
      })
      ord <- order(vapply(candidates, function(x) x$score$objective, numeric(1)))
      chosen <- candidates[[ord[[1L]]]]
      changed <- changed || !isTRUE(all.equal(chosen$gamma, gamma))
      gamma <- chosen$gamma
    }
    if (!changed) break
  }
  X <- .rdm_decorrelation_ordered_apply(E, P, gamma, tol = tol)
  list(X = X, gamma = gamma, diagnostics = .rdm_decorrelation_score(X, Y, objective, epsilon, mu, target))
}

.rdm_decorrelation_score <- function(X, Y, objective, epsilon, mu, target) {
  C <- .rdm_decorrelation_cor_matrix(X)
  dec <- if (identical(target, "mean_squared")) {
    .rdm_decorrelation_mean_sq_cor(C)
  } else {
    .rdm_decorrelation_mean_abs_cor(C)
  }
  pres <- mean(.rdm_decorrelation_col_cor(X, Y), na.rm = TRUE)
  feasible <- dec <= epsilon
  obj <- if (identical(objective, "constraint")) {
    if (feasible) -pres else dec + (1 - pres)
  } else {
    (1 - pres) + mu * dec
  }
  list(objective = obj, preservation = pres, decorrelation = dec, feasible = feasible)
}

.rdm_decorrelation_zca <- function(Y, zca_ridge, tol) {
  C <- crossprod(Y) / (nrow(Y) - 1L)
  ee <- eigen((C + t(C)) / 2, symmetric = TRUE)
  vals <- ee$values
  ridge_used <- 0
  if (any(vals <= tol)) {
    ridge_used <- zca_ridge
  }
  inv_sqrt <- ee$vectors %*% diag(1 / sqrt(pmax(vals, 0) + ridge_used), nrow = length(vals)) %*% t(ee$vectors)
  X <- Y %*% inv_sqrt
  colnames(X) <- colnames(Y)
  list(
    X = X,
    gamma = NULL,
    diagnostics = list(eigenvalues = vals, ridge_used = ridge_used)
  )
}

.rdm_decorrelation_col_cor <- function(A, B) {
  out <- rep(NA_real_, ncol(A))
  for (j in seq_len(ncol(A))) {
    if (stats::sd(A[, j]) > 0 && stats::sd(B[, j]) > 0) {
      out[[j]] <- stats::cor(A[, j], B[, j])
    }
  }
  stats::setNames(out, colnames(A))
}

.rdm_decorrelation_cor_matrix <- function(X, tol = 1e-12) {
  X <- as.matrix(X)
  L <- ncol(X)
  out <- diag(1, L)
  colnames(out) <- rownames(out) <- colnames(X)
  if (L < 2L) return(out)
  sds <- apply(X, 2L, stats::sd)
  for (a in seq_len(L - 1L)) {
    for (b in (a + 1L):L) {
      val <- if (is.finite(sds[[a]]) && is.finite(sds[[b]]) && sds[[a]] > tol && sds[[b]] > tol) {
        stats::cor(X[, a], X[, b])
      } else {
        0
      }
      out[a, b] <- out[b, a] <- val
    }
  }
  out
}

.rdm_decorrelation_mean_abs_cor <- function(C) {
  if (ncol(C) < 2L) return(0)
  mean(abs(C[lower.tri(C)]), na.rm = TRUE)
}

.rdm_decorrelation_mean_sq_cor <- function(C) {
  if (ncol(C) < 2L) return(0)
  mean(C[lower.tri(C)]^2, na.rm = TRUE)
}

.rdm_decorrelation_vectors_to_matrices <- function(V, n, labels, names) {
  out <- vector("list", ncol(V))
  for (j in seq_len(ncol(V))) {
    M <- matrix(0, n, n, dimnames = list(labels, labels))
    M[lower.tri(M)] <- V[, j]
    M <- M + t(M)
    out[[j]] <- M
  }
  names(out) <- names
  out
}

.rdm_decorrelation_matrix_list_to_vectors <- function(mats) {
  V <- vapply(mats, function(M) M[lower.tri(M)], numeric(nrow(mats[[1]]) * (nrow(mats[[1]]) - 1L) / 2L))
  colnames(V) <- names(mats)
  V
}

.rdm_decorrelation_project_correlation <- function(D) {
  K <- 1 - D
  K <- (K + t(K)) / 2
  diag(K) <- 1
  projected <- as.matrix(Matrix::nearPD(K, corr = TRUE, keepDiag = TRUE)$mat)
  Dp <- 1 - projected
  diag(Dp) <- 0
  dimnames(Dp) <- dimnames(D)
  Dp
}
