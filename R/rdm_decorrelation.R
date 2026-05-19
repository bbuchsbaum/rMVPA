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
#' @param output Character string. `"dist"` returns adjusted RDMs as `dist`
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
                            output = c("dist", "matrix", "vector"),
                            project_correlation = FALSE,
                            zca_ridge = 1e-8,
                            max_combinations = 2e5,
                            tol = 1e-10) {
  method <- match.arg(method)
  similarity <- match.arg(similarity)
  objective <- match.arg(objective)
  target <- match.arg(target)
  output <- match.arg(output)

  parsed <- .rdm_decorrelation_parse_rdms(rdms)
  Y <- .rdm_standardize_columns(parsed$vectors, similarity = similarity,
                                tol = tol, label = "RDM vectors")

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

  out_rdms <- switch(output,
    vector = X,
    matrix = mats,
    dist = lapply(mats, stats::as.dist)
  )
  if (output != "vector") {
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

#' Parse a list of RDMs into a lower-triangular vector matrix
#'
#' Shared by [rdm_decorrelate()] and the model-space connectivity helpers.
#'
#' @keywords internal
#' @noRd
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

#' Rank-optionally standardize the columns of an RDM-vector matrix
#'
#' Shared standardization used by [rdm_decorrelate()],
#' [rdm_model_space_connectivity()], and the RSA model-space basis. When
#' `similarity == "spearman"` each column is rank-transformed (average ties)
#' before centering and unit scaling. Errors if any column is constant after
#' the transform.
#'
#' @keywords internal
#' @noRd
.rdm_standardize_columns <- function(mat, similarity, tol, label = "RDM vectors") {
  mat <- as.matrix(mat)
  dn <- dimnames(mat)
  vals <- if (identical(similarity, "spearman")) {
    apply(mat, 2L, rank, ties.method = "average")
  } else {
    mat
  }
  out <- apply(vals, 2L, function(x) {
    x <- as.numeric(x)
    x <- x - mean(x)
    s <- stats::sd(x)
    if (!is.finite(s) || s <= tol) {
      stop(label, " must be non-constant after the requested transformation.", call. = FALSE)
    }
    x / s
  })
  out <- as.matrix(out)
  dimnames(out) <- dn
  out
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
    # qr.fitted/qr.resid honor the QR-detected rank, so rank-deficient
    # innovations (duplicate/collinear RDMs) are handled correctly without
    # manually slicing qr.Q() columns under base R's column pivoting.
    qr_prev <- qr(prev, tol = tol)
    P[, j] <- qr.fitted(qr_prev, Y[, j])
    E[, j] <- qr.resid(qr_prev, Y[, j])
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
  # Materialize the gamma grid as a numeric matrix; row extraction from a
  # data.frame (expand.grid's default) is slow in the per-combination loop.
  combos <- as.matrix(expand.grid(grids, KEEP.OUT.ATTRS = FALSE))
  best <- NULL
  best_score <- Inf
  best_feasible <- NULL
  best_feasible_pres <- -Inf
  best_feasible_dec <- Inf

  for (i in seq_len(nrow(combos))) {
    gamma <- combos[i, ]
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
  # Columns are always either ~unit-sd standardized vectors or exactly zero
  # (degenerate columns are zeroed upstream), so stats::cor() yields NA only
  # for the degenerate columns; treat those as zero correlation.
  C <- suppressWarnings(stats::cor(X))
  C[!is.finite(C)] <- 0
  diag(C) <- 1
  dimnames(C) <- dimnames(out)
  C
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
