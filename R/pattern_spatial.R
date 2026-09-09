# Feature-aligned spatial graphs
#
# A spatial_graph is the anatomical adjacency structure used by spatially
# regularized estimators. Its single invariant is that graph vertex j refers
# to column j of the feature matrix produced by get_feature_matrix() for the
# dataset it was built from. Everything else (Laplacian, degree, edge lists)
# is derived from the adjacency matrix.

#' @keywords internal
#' @noRd
new_spatial_graph <- function(A, feature_ids, domain_type, geometry_id,
                              basis = NULL, weighted = FALSE,
                              parent_index = NULL) {
  if (!inherits(A, "sparseMatrix")) {
    A <- Matrix::Matrix(as.matrix(A), sparse = TRUE)
  }
  # Normalize to a numeric general column-sparse matrix (pattern matrices such
  # as the edgeless ngCMatrix returned for isolated voxels have no @x slot).
  A <- methods::as(methods::as(methods::as(A, "generalMatrix"), "CsparseMatrix"), "dMatrix")
  if (nrow(A) != ncol(A)) {
    stop("spatial_graph: adjacency matrix must be square.", call. = FALSE)
  }
  n <- nrow(A)
  feature_ids <- as.integer(feature_ids)
  if (length(feature_ids) != n) {
    stop(sprintf("spatial_graph: %d feature ids supplied for a %d-vertex graph.",
                 length(feature_ids), n), call. = FALSE)
  }
  if (!is.null(basis) && length(basis) != n) {
    stop("spatial_graph: 'basis' must have one entry per vertex.", call. = FALSE)
  }

  # Validate raw weights before any symmetrization can mask them.
  A <- Matrix::drop0(A)
  if (length(A@x) && (anyNA(A@x) || any(!is.finite(A@x)) || any(A@x < 0))) {
    stop("spatial_graph: edge weights must be finite and non-negative.", call. = FALSE)
  }
  # Symmetrize by taking the larger of the two directions. base::pmax has no
  # sparse method and would coerce to a dense p x p matrix, which is fatal at
  # whole-brain feature counts, so use the sparse-preserving identity
  # max(a, b) = (a + b + |a - b|) / 2.
  At <- Matrix::t(A)
  A <- Matrix::drop0((A + At + abs(A - At)) / 2)
  Matrix::diag(A) <- 0
  A <- Matrix::drop0(A)
  if (!isTRUE(weighted) && length(A@x)) A@x[] <- 1
  A <- methods::as(methods::as(A, "generalMatrix"), "CsparseMatrix")

  degree <- as.numeric(Matrix::rowSums(A))
  L <- Matrix::Diagonal(n = n, x = degree) - A

  structure(
    list(
      A = A,
      degree = degree,
      L = L,
      L_norm = if (n > 0L) .pattern_L_norm(L) else 0,
      feature_ids = feature_ids,
      n_features = n,
      domain_type = domain_type,
      geometry_id = geometry_id,
      basis = basis,
      weighted = isTRUE(weighted),
      n_edges = as.integer(Matrix::nnzero(A) / 2),
      parent_index = parent_index
    ),
    class = c("spatial_graph", "list")
  )
}

#' @export
print.spatial_graph <- function(x, ...) {
  cat(sprintf("spatial_graph [%s]: %d features, %d edges%s\n",
              x$domain_type, x$n_features, x$n_edges,
              if (x$weighted) " (weighted)" else ""))
  if (!is.null(x$basis)) {
    cat(sprintf("  basis channels: %d\n", length(unique(x$basis))))
  }
  if (!is.null(x$parent_index)) {
    cat("  restricted from a parent graph\n")
  }
  cat("  geometry:", x$geometry_id, "\n")
  invisible(x)
}

#' @keywords internal
#' @noRd
.mask_active_vec <- function(mask) {
  vals <- if (inherits(mask, c("NeuroVol", "NeuroVec"))) as.numeric(neuroim2::values(mask)) else as.numeric(mask)
  vals[is.na(vals)] <- 0
  vals > 0
}

#' @keywords internal
#' @noRd
.volume_geometry_id <- function(mask, dims, neighbors) {
  spacing <- tryCatch(
    paste(signif(neuroim2::spacing(neuroim2::space(mask)), 6), collapse = "x"),
    error = function(e) "NA"
  )
  sprintf("volume:%s:spacing=%s:nbr=%d", paste(dims, collapse = "x"), spacing, as.integer(neighbors))
}

#' @keywords internal
#' @noRd
.volume_adjacency_from_mask <- function(mask, neighbors) {
  dims <- .infer_spatial_dims(mask, NULL)
  active <- .mask_active_vec(mask)
  if (length(active) != prod(dims)) {
    stop("spatial_graph: mask length does not match its spatial dimensions.", call. = FALSE)
  }
  A <- build_voxel_adjacency(active, dims = dims, neighbors = neighbors)
  list(A = A, feature_ids = which(active), dims = dims)
}

#' @rdname spatial_graph
#' @param neighbors Voxel neighbourhood for volumetric grids: 6 (faces),
#'   18 (faces and edges), or 26 (faces, edges, and corners).
#' @export
spatial_graph.mvpa_image_dataset <- function(x, neighbors = 6, ...) {
  vol <- .volume_adjacency_from_mask(x$mask, neighbors)
  new_spatial_graph(
    vol$A,
    feature_ids = vol$feature_ids,
    domain_type = "volume",
    geometry_id = .volume_geometry_id(x$mask, vol$dims, neighbors)
  )
}

#' @rdname spatial_graph
#' @param connect_basis Logical; for multibasis datasets, also connect each
#'   voxel to the same voxel in every other basis channel (default
#'   \code{FALSE}: channels form disconnected components, so smoothing never
#'   crosses channels).
#' @export
spatial_graph.mvpa_multibasis_image_dataset <- function(x, neighbors = 6,
                                                        connect_basis = FALSE, ...) {
  vol <- .volume_adjacency_from_mask(x$mask, neighbors)
  k <- as.integer(x$basis_count)
  n_vox <- length(vol$feature_ids)
  A <- Matrix::bdiag(replicate(k, vol$A, simplify = FALSE))
  if (isTRUE(connect_basis) && k > 1L) {
    pairs <- utils::combn(k, 2)
    ii <- unlist(lapply(seq_len(ncol(pairs)), function(p) (pairs[1, p] - 1L) * n_vox + seq_len(n_vox)))
    jj <- unlist(lapply(seq_len(ncol(pairs)), function(p) (pairs[2, p] - 1L) * n_vox + seq_len(n_vox)))
    cross <- Matrix::sparseMatrix(i = ii, j = jj, x = 1, dims = c(k * n_vox, k * n_vox))
    A <- A + cross + Matrix::t(cross)
  }
  new_spatial_graph(
    A,
    feature_ids = rep(vol$feature_ids, times = k),
    domain_type = "multibasis_volume",
    geometry_id = paste0(.volume_geometry_id(x$mask, vol$dims, neighbors), ":k=", k),
    basis = rep(seq_len(k), each = n_vox)
  )
}

#' @rdname spatial_graph
#' @export
spatial_graph.mvpa_surface_dataset <- function(x, ...) {
  require_package("neurosurf", "for surface adjacency graphs")
  geom <- x$train_data@geometry
  A_full <- neurosurf::adjacency(geom)
  idx <- which(as.numeric(x$mask) > 0)
  if (max(idx) > nrow(A_full)) {
    stop("spatial_graph: surface mask indexes nodes beyond the mesh.", call. = FALSE)
  }
  new_spatial_graph(
    A_full[idx, idx, drop = FALSE],
    feature_ids = idx,
    domain_type = "surface",
    geometry_id = sprintf("surface:nodes=%d:edges=%d", nrow(A_full),
                          as.integer(Matrix::nnzero(A_full) / 2))
  )
}

#' @rdname spatial_graph
#' @details
#' For clustered datasets the graph vertices are the parcels, in the column
#' order of the cluster time series, and \code{feature_ids} holds the actual
#' cluster labels (which need not be \code{1..K}). Column positions, not
#' labels, index the feature matrix.
#' @export
spatial_graph.mvpa_clustered_dataset <- function(x, neighbors = 6, ...) {
  cvol <- x$train_data@cvol
  vol <- .volume_adjacency_from_mask(cvol@mask, neighbors)
  clusters <- as.integer(cvol@clusters)
  if (length(clusters) != length(vol$feature_ids)) {
    stop("spatial_graph: cluster labels do not align with the clustered mask.", call. = FALSE)
  }
  cl_ids <- sort(unique(clusters))
  K <- ncol(x$train_data@ts)
  if (length(cl_ids) != K) {
    stop(sprintf("spatial_graph: %d cluster labels but %d cluster time series.",
                 length(cl_ids), K), call. = FALSE)
  }
  M <- Matrix::sparseMatrix(
    i = seq_along(clusters),
    j = match(clusters, cl_ids),
    x = 1,
    dims = c(length(clusters), K)
  )
  A_c <- Matrix::crossprod(M, vol$A %*% M)
  new_spatial_graph(
    A_c,
    feature_ids = cl_ids,
    domain_type = "cluster",
    geometry_id = paste0(.volume_geometry_id(cvol@mask, vol$dims, neighbors), ":clusters=", K)
  )
}

#' @rdname spatial_graph
#' @param feature_ids Integer identifiers mapping graph vertices to dataset
#'   locations when a raw adjacency is supplied. Defaults to \code{seq_len(n)}.
#' @param weighted Logical; keep edge weights of a raw adjacency instead of
#'   binarizing them.
#' @param domain_type Label for a raw adjacency's domain (default
#'   \code{"custom"}).
#' @export
spatial_graph.default <- function(x, feature_ids = NULL, weighted = FALSE,
                                  domain_type = "custom", ...) {
  A <- if (is.list(x) && !is.data.frame(x) && !is.null(x$A)) {
    # Fields carried by the list are defaults; explicit arguments win.
    if (is.null(feature_ids) && !is.null(x$feature_ids)) feature_ids <- x$feature_ids
    if (missing(weighted) && !is.null(x$weighted)) weighted <- isTRUE(x$weighted)
    x$A
  } else if (is.matrix(x) || inherits(x, "Matrix")) {
    x
  } else {
    stop(sprintf("spatial_graph: no method for class '%s'.", paste(class(x), collapse = "/")),
         call. = FALSE)
  }
  n <- nrow(A)
  if (is.null(feature_ids)) feature_ids <- seq_len(n)
  new_spatial_graph(
    A,
    feature_ids = feature_ids,
    domain_type = domain_type,
    geometry_id = sprintf("%s:n=%d", domain_type, n),
    weighted = weighted
  )
}

#' Restrict a Spatial Graph to a Subset of Features
#'
#' Returns the induced subgraph on the selected feature columns, keeping the
#' feature-to-column alignment. Use it after dropping invalid or constant
#' columns from the feature matrix so the graph continues to describe exactly
#' the columns the estimator sees.
#'
#' @param graph A \code{spatial_graph}.
#' @param keep Either a logical vector with one entry per feature column or an
#'   integer vector of column positions to keep (in the order they should
#'   appear in the restricted matrix).
#' @return A \code{spatial_graph} over the kept columns. Its
#'   \code{parent_index} field records the position of each kept column in
#'   the parent graph.
#' @examples
#' ds <- gen_sample_dataset(c(4, 4, 4), 8)
#' g <- spatial_graph(ds$dataset)
#' g2 <- restrict_graph(g, seq_len(10))
#' g2$n_features
#' @seealso \code{\link{spatial_graph}}
#' @export
restrict_graph <- function(graph, keep) {
  stopifnot(inherits(graph, "spatial_graph"))
  if (is.logical(keep)) {
    if (length(keep) != graph$n_features) {
      stop("restrict_graph: logical 'keep' must have one entry per feature.", call. = FALSE)
    }
    keep <- which(keep)
  }
  keep <- as.integer(keep)
  if (any(keep < 1L | keep > graph$n_features) || anyDuplicated(keep)) {
    stop("restrict_graph: 'keep' must be unique column positions within the graph.", call. = FALSE)
  }
  parent <- if (is.null(graph$parent_index)) keep else graph$parent_index[keep]
  new_spatial_graph(
    graph$A[keep, keep, drop = FALSE],
    feature_ids = graph$feature_ids[keep],
    domain_type = graph$domain_type,
    geometry_id = graph$geometry_id,
    basis = if (!is.null(graph$basis)) graph$basis[keep] else NULL,
    weighted = graph$weighted,
    parent_index = parent
  )
}

#' Edge List of a Spatial Graph
#'
#' @param graph A \code{spatial_graph}.
#' @return A data frame with integer columns \code{from} and \code{to}
#'   (column positions, \code{from < to}) and a numeric \code{weight}.
#' @examples
#' ds <- gen_sample_dataset(c(4, 4, 4), 8)
#' head(graph_edges(spatial_graph(ds$dataset)))
#' @seealso \code{\link{spatial_graph}}
#' @export
graph_edges <- function(graph) {
  stopifnot(inherits(graph, "spatial_graph"))
  s <- Matrix::summary(graph$A)
  s <- s[s$i < s$j, , drop = FALSE]
  data.frame(from = as.integer(s$i), to = as.integer(s$j), weight = as.numeric(s$x))
}

#' @keywords internal
#' @noRd
.assert_graph_aligned <- function(graph, n_columns, what = "feature matrix") {
  if (!inherits(graph, "spatial_graph")) {
    stop("A 'spatial_graph' is required.", call. = FALSE)
  }
  if (!identical(as.integer(graph$n_features), as.integer(n_columns))) {
    stop(sprintf("spatial_graph has %d vertices but the %s has %d columns.",
                 graph$n_features, what, as.integer(n_columns)), call. = FALSE)
  }
  invisible(TRUE)
}


# ---------------------------------------------------------------------------
# Spatial penalties for the pattern model
#
# The A-step minimises, for fixed target directions C and residual covariance
# Psi,
#
#   g(A) + h(A),
#   g(A) = 1/(2n) tr[(X - T A') Psi^{-1} (X - T A')'] + lambda_l/2 tr(A' L A)
#   h(A) = lambda_s sum_v ||A_v||_2
#
# with T = Y_w C. g is smooth and convex, h is separable across features, so
# the step is a proximal gradient problem. Expanding g removes every n x p
# operation from the inner loop: with XtT = X'T and the constant tr(X Psi^{-1} X')
# precomputed once per outer iteration, each FISTA iteration costs
# O(p r^2 + nnz(L) r).
# ---------------------------------------------------------------------------

#' Row-wise group soft-threshold (the prox of the group lasso).
#' @keywords internal
#' @noRd
.prox_group_lasso <- function(B, thr) {
  if (thr <= 0) return(B)
  nrm <- sqrt(rowSums(B^2))
  scale <- pmax(0, 1 - thr / pmax(nrm, .Machine$double.eps))
  B * scale
}

#' A deterministic starting vector for power iteration.
#'
#' Fits must be reproducible, so the power iterations that set the step size
#' cannot draw from the global RNG stream: two identical calls would otherwise
#' return slightly different step sizes and therefore slightly different
#' iterates. This is a fixed low-discrepancy pattern, unlikely to be orthogonal
#' to a leading eigenvector.
#' @keywords internal
#' @noRd
.power_start <- function(m) {
  i <- seq_len(m)
  v <- sin(i * 0.7391) + cos(i * 1.2153)
  v / sqrt(sum(v^2))
}

#' Spectral norm of the A-step's linear operator, by power iteration.
#'
#' The gradient of the smooth part is linear in A:
#'   A -> (1/n) Psi^{-1} A (T'T) + lambda_l L A
#' and its spectral norm is the Lipschitz constant that sets the step size. The
#' 1/n is part of the operator: dropping it would make the step roughly n times
#' too short and inflate any penalty calibrated against this norm by the same
#' factor. Power iteration converges from below, so `safety` inflates the
#' result when it is used as a Lipschitz bound; leave it at 1 when the value is
#' wanted as a scale reference instead.
#' @keywords internal
#' @noRd
.pattern_grad_norm <- function(noise, TtT, L, lambda_l, n, p, r, iter = 60L, tol = 1e-8,
                               safety = 1.02) {
  V <- matrix(.power_start(p * r), p, r)
  lam <- 0
  for (i in seq_len(iter)) {
    W <- .noise_apply_precision(noise, V %*% TtT) / n
    if (lambda_l > 0 && !is.null(L)) W <- W + lambda_l * as.matrix(L %*% V)
    nw <- sqrt(sum(W^2))
    if (!is.finite(nw) || nw <= 0) return(1)
    V <- W / nw
    if (abs(nw - lam) <= tol * nw) { lam <- nw; break }
    lam <- nw
  }
  lam * safety
}

#' Smallest group-lasso penalty that zeroes every feature at A = 0.
#'
#' At A = 0 the smooth gradient is -(1/n) Psi^{-1} X'T, so feature v survives
#' the prox only when its gradient row norm exceeds lambda_s.
#' @keywords internal
#' @noRd
.pattern_lambda_max <- function(PXtT, n) {
  max(sqrt(rowSums((PXtT / n)^2)))
}

#' Objective of the A-step, without any n x p operation.
#' @keywords internal
#' @noRd
.pattern_astep_objective <- function(A, prob) {
  PA <- .noise_apply_precision(prob$noise, A)
  quad <- sum((A %*% prob$TtT) * PA)
  cross <- sum(prob$XtT * PA)
  g <- (prob$const - 2 * cross + quad) / (2 * prob$n)
  if (prob$lambda_l > 0) g <- g + 0.5 * prob$lambda_l * sum(A * as.matrix(prob$L %*% A))
  g + prob$lambda_s * sum(sqrt(rowSums(A^2)))
}

#' Gradient of the smooth part of the A-step.
#' @keywords internal
#' @noRd
.pattern_astep_gradient <- function(A, prob) {
  G <- .noise_apply_precision(prob$noise, A %*% prob$TtT - prob$XtT) / prob$n
  if (prob$lambda_l > 0) G <- G + prob$lambda_l * as.matrix(prob$L %*% A)
  G
}

#' Monotone FISTA for the penalized A-step.
#'
#' Uses the monotone variant (Beck & Teboulle 2009), which accepts a candidate
#' only when it does not increase the objective, so the returned trace is
#' non-increasing by construction.
#' @keywords internal
#' @noRd
.pattern_astep_fista <- function(A_init, prob, max_iter = 300L, tol = 1e-7,
                                 tol_iterate = 1e-5) {
  A <- A_init
  step <- 1 / prob$lipschitz
  y <- A
  tk <- 1
  restarted <- FALSE
  obj <- .pattern_astep_objective(A, prob)
  trace <- obj
  converged <- FALSE
  iters <- 0L
  for (k in seq_len(max_iter)) {
    iters <- k
    z <- .prox_group_lasso(y - step * .pattern_astep_gradient(y, prob), step * prob$lambda_s)
    obj_z <- .pattern_astep_objective(z, prob)

    if (obj_z > obj) {
      # No progress. Never report this as convergence: a rejected candidate
      # leaves the objective unchanged, which would otherwise look like a zero
      # relative change. Momentum overshoot is the common cause, so drop the
      # momentum first; if a plain proximal-gradient step from the current
      # iterate also fails, the step length itself is too long.
      if (restarted) step <- step / 2 else restarted <- TRUE
      tk <- 1
      y <- A
      next
    }

    A_prev <- A
    A <- z
    restarted <- FALSE
    tk_new <- (1 + sqrt(1 + 4 * tk^2)) / 2
    y <- A + ((tk - 1) / tk_new) * (A - A_prev)
    tk <- tk_new
    rel <- abs(obj - obj_z) / max(abs(obj), .Machine$double.eps)
    # The objective is flat near the optimum, so a small change in it does not
    # imply a settled iterate. Require both: otherwise the patterns, which are
    # the scientific output, can still be a percent away when the loop stops.
    rel_A <- sqrt(sum((A - A_prev)^2)) / max(sqrt(sum(A^2)), .Machine$double.eps)
    trace <- c(trace, obj_z)
    obj <- obj_z
    if (rel < tol && rel_A < tol_iterate) { converged <- TRUE; break }
  }
  list(A = A, objective = obj, trace = trace, iterations = iters, converged = converged)
}

#' Assemble the fixed quantities of one A-step subproblem.
#'
#' `XtT` and the constant term are the only n x p work; everything the solver
#' then does is O(p r^2 + nnz(L) r) per iteration.
#' @keywords internal
#' @noRd
.pattern_astep_problem <- function(X, Tm, noise, L, lambda_s, lambda_l, const, XtT = NULL) {
  n <- nrow(X)
  prob <- list(
    n = n, TtT = crossprod(Tm),
    XtT = if (is.null(XtT)) crossprod(X, Tm) else XtT,
    noise = noise,
    L = L, lambda_s = lambda_s, lambda_l = lambda_l, const = const
  )
  prob$lipschitz <- .pattern_grad_norm(noise, prob$TtT, L, lambda_l, n, ncol(X), ncol(Tm))
  prob
}

#' Resolve a penalty specification into concrete lambda values.
#'
#' \code{sparse} is a fraction of the smallest penalty that zeroes every
#' feature, so it means the same thing across folds and feature domains.
#' \code{signed_smooth} is the weight of the graph-Laplacian term relative to
#' the curvature of the data-fit term, again dimensionless: a value of 1 makes
#' smoothing as influential as the data fit. Both are resolved on training
#' rows only.
#' @keywords internal
#' @noRd
.pattern_resolve_penalty <- function(penalty, X, Tm, noise, graph, XtT = NULL) {
  out <- list(lambda_s = 0, lambda_l = 0, alpha = 0, rho = 0,
              lambda_max = NA_real_, L = NULL, active = FALSE)
  if (is.null(penalty)) return(out)
  alpha <- penalty$sparse %||% 0
  rho <- penalty$signed_smooth %||% 0
  if (!is.numeric(alpha) || length(alpha) != 1L || alpha < 0 || alpha >= 1) {
    stop("pattern_model: penalty$sparse must be a single number in [0, 1).", call. = FALSE)
  }
  if (!is.numeric(rho) || length(rho) != 1L || rho < 0) {
    stop("pattern_model: penalty$signed_smooth must be a single number >= 0.", call. = FALSE)
  }
  if (alpha == 0 && rho == 0) return(out)

  n <- nrow(X)
  PXtT <- .noise_apply_precision(noise, if (is.null(XtT)) crossprod(X, Tm) else XtT)
  out$lambda_max <- .pattern_lambda_max(PXtT, n)
  out$alpha <- alpha
  out$lambda_s <- alpha * out$lambda_max

  if (rho > 0) {
    if (is.null(graph)) {
      stop("pattern_model: penalty$signed_smooth requires a spatial graph; ",
           "pass one via spatial_graph(dataset).", call. = FALSE)
    }
    .assert_graph_aligned(graph, ncol(X), "feature matrix")
    # Normalise the Laplacian to unit spectral norm so rho is comparable across
    # domains, then scale it to the data-fit curvature.
    L <- graph$L
    Lnorm <- graph$L_norm %||% .pattern_L_norm(L)
    if (!is.finite(Lnorm) || Lnorm <= 0) {
      out$L <- NULL
    } else {
      out$L <- L / Lnorm
      # scale reference, not a step-size bound, so no safety inflation
      dat <- .pattern_grad_norm(noise, crossprod(Tm), NULL, 0, n, ncol(X), ncol(Tm),
                                safety = 1)
      out$lambda_l <- rho * dat
    }
  }
  out$rho <- rho
  out$active <- out$lambda_s > 0 || out$lambda_l > 0
  out
}

#' Spectral norm of a graph Laplacian (power iteration; bounded by 2 max degree).
#' @keywords internal
#' @noRd
.pattern_L_norm <- function(L, iter = 100L, tol = 1e-9) {
  p <- nrow(L)
  if (p == 0L) return(0)
  v <- .power_start(p)
  lam <- 0
  for (i in seq_len(iter)) {
    w <- as.numeric(L %*% v)
    nw <- sqrt(sum(w^2))
    if (!is.finite(nw) || nw <= 0) return(0)
    v <- w / nw
    if (abs(nw - lam) <= tol * nw) { lam <- nw; break }
    lam <- nw
  }
  lam
}
