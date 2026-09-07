# Numerical core of pattern_model
#
# Forward model:  X = T A' + E,  T = Y_w C,  C' C = I,  E ~ (0, Psi)
#   X    n x p  brain measurements (column-centred on training rows)
#   Y_w  n x q  whitened target coding (see .pattern_encode_targets)
#   C    q x r  target directions (orthonormal columns)
#   A    p x r  spatial patterns
#
# This file knows nothing about rMVPA datasets, designs, folds, or futures.
# .pattern_fit() takes matrices in and returns a serializable `pattern_fit`.

# ---------------------------------------------------------------------------
# Feature transform (columnwise only, so regional restriction stays local)
# ---------------------------------------------------------------------------

.pattern_x_transform_fit <- function(X, scale = c("none", "sd")) {
  scale <- match.arg(scale)
  mu <- colMeans(X)
  sd <- if (scale == "sd") {
    s <- sqrt(colSums(sweep(X, 2L, mu, "-")^2) / max(nrow(X) - 1, 1))
    s[!is.finite(s) | s <= 0] <- 1
    s
  } else {
    NULL
  }
  list(mu = mu, sd = sd, scale = scale)
}

.pattern_x_transform_apply <- function(xt, X, cols = NULL) {
  X <- as.matrix(X)
  mu <- xt$mu; sd <- xt$sd
  if (!is.null(cols)) { mu <- mu[cols]; if (!is.null(sd)) sd <- sd[cols] }
  X <- sweep(X, 2L, mu, "-")
  if (!is.null(sd)) X <- sweep(X, 2L, sd, "/")
  X
}

# Columns that are finite and non-constant on the training rows. Positions are
# kept (feature_index) so dropped columns stay distinguishable in maps.
.pattern_screen_columns <- function(X) {
  finite <- colSums(!is.finite(X)) == 0
  v <- colSums(sweep(X, 2L, colMeans(X), "-")^2)
  keep <- which(finite & v > 0)
  keep
}

# ---------------------------------------------------------------------------
# Target coding
# ---------------------------------------------------------------------------

#' Encode targets for the pattern model.
#'
#' Categorical targets become centred one-hot codes; continuous targets are
#' centred (and optionally scaled). Both are then whitened on the training
#' rows: Y_w = Y_c S^{-1/2} with S = Y_c' Y_c / n restricted to its
#' numerically non-null eigenspace. Whitening (i) removes the rank deficiency
#' of centred one-hot codes (K classes -> K - 1 columns), (ii) makes the
#' C-step an exact orthogonal Procrustes problem, and (iii) gives the working
#' prior y_w ~ N(0, I) used for decoding.
#' @keywords internal
#' @noRd
.pattern_encode_targets <- function(values, scale = c("none", "sd"), tol = 1e-8) {
  scale <- match.arg(scale)
  if (is.character(values) || is.logical(values)) values <- factor(values)
  if (anyNA(values)) {
    stop("pattern_model: targets contain missing values; drop or impute those ",
         "observations before fitting.", call. = FALSE)
  }

  if (is.factor(values)) {
    values <- droplevels(values)
    lev <- levels(values)
    K <- length(lev)
    if (K < 2L) stop("pattern_model: categorical targets need at least two classes.", call. = FALSE)
    Y <- matrix(0, length(values), K, dimnames = list(NULL, lev))
    Y[cbind(seq_along(values), as.integer(values))] <- 1
    priors <- colMeans(Y)
    mu <- priors
    sdv <- NULL
    type <- "categorical"
  } else {
    Y <- as.matrix(values)
    if (!is.numeric(Y)) stop("pattern_model: continuous targets must be numeric.", call. = FALSE)
    if (is.null(colnames(Y))) colnames(Y) <- paste0("y", seq_len(ncol(Y)))
    lev <- NULL; priors <- NULL; K <- NA_integer_
    mu <- colMeans(Y)
    sdv <- if (scale == "sd") {
      s <- sqrt(colSums(sweep(Y, 2L, mu, "-")^2) / max(nrow(Y) - 1, 1))
      s[!is.finite(s) | s <= 0] <- 1
      s
    } else {
      NULL
    }
    type <- "continuous"
  }

  Yc <- sweep(Y, 2L, mu, "-")
  if (!is.null(sdv)) Yc <- sweep(Yc, 2L, sdv, "/")
  n <- nrow(Yc)

  # Whitening on the non-null eigenspace of the target covariance.
  S <- crossprod(Yc) / n
  eg <- eigen(S, symmetric = TRUE)
  keep <- eg$values > tol * max(eg$values[1], .Machine$double.eps)
  if (!any(keep)) stop("pattern_model: targets are constant on the training rows.", call. = FALSE)
  V <- eg$vectors[, keep, drop = FALSE]
  lam <- eg$values[keep]
  Wy <- V %*% diag(1 / sqrt(lam), nrow = length(lam))     # q x q_eff
  Wy_inv <- diag(sqrt(lam), nrow = length(lam)) %*% t(V)  # q_eff x q  (S^{1/2} restricted)

  yt <- list(
    type = type, levels = lev, priors = priors, K = K,
    mu = mu, sd = sdv, scale = scale,
    Wy = Wy, Wy_inv = Wy_inv, q = ncol(Y), q_eff = ncol(Wy),
    response_ids = colnames(Y)
  )
  # class representatives in whitened coordinates (K x q_eff)
  if (type == "categorical") {
    yt$class_codes <- .pattern_targets_apply(yt, factor(lev, levels = lev))
  }
  list(Yw = Yc %*% Wy, transform = yt)
}

# Raw targets -> whitened coordinates using training statistics.
.pattern_targets_apply <- function(yt, values) {
  if (yt$type == "categorical") {
    if (!is.factor(values)) values <- factor(values, levels = yt$levels)
    values <- factor(as.character(values), levels = yt$levels)
    if (anyNA(values)) stop("pattern_model: test labels contain classes absent from training.", call. = FALSE)
    Y <- matrix(0, length(values), yt$K)
    Y[cbind(seq_along(values), as.integer(values))] <- 1
  } else {
    Y <- as.matrix(values)
    if (ncol(Y) != yt$q) stop("pattern_model: test targets have the wrong number of responses.", call. = FALSE)
  }
  Yc <- sweep(Y, 2L, yt$mu, "-")
  if (!is.null(yt$sd)) Yc <- sweep(Yc, 2L, yt$sd, "/")
  Yc %*% yt$Wy
}

# Whitened coordinates -> original target scale (continuous targets).
.pattern_targets_invert <- function(yt, Yw) {
  Yc <- as.matrix(Yw) %*% yt$Wy_inv
  if (!is.null(yt$sd)) Yc <- sweep(Yc, 2L, yt$sd, "*")
  out <- sweep(Yc, 2L, yt$mu, "+")
  colnames(out) <- yt$response_ids
  out
}

# ---------------------------------------------------------------------------
# Objective and updates (Phase 2: no spatial penalty)
# ---------------------------------------------------------------------------

# f(A, C) = 1/(2n) tr[(X - Yw C A') Psi^{-1} (X - Yw C A')']
.pattern_objective <- function(X, Yw, A, C, noise, penalty = NULL) {
  R <- X - (Yw %*% C) %*% t(A)
  val <- 0.5 * .noise_quadform(noise, R) / nrow(X)
  if (!is.null(penalty) && isTRUE(penalty$active)) {
    if (penalty$lambda_s > 0) val <- val + penalty$lambda_s * sum(sqrt(rowSums(A^2)))
    if (penalty$lambda_l > 0) val <- val + 0.5 * penalty$lambda_l * sum(A * as.matrix(penalty$L %*% A))
  }
  val
}

# A-step. With no penalty the stationarity condition is
# Psi^{-1}(A T'T - X'T)/n = 0, so Psi cancels and A = X'T (T'T)^{-1} is the
# exact conditional minimiser for any Psi.
#
# A ridge term would make the condition Psi^{-1}(A T'T - X'T)/n + lambda A = 0,
# a Sylvester equation whose solution is NOT X'T(T'T + n lambda I)^{-1} unless
# Psi is a multiple of the identity. Rather than converge to a point that does
# not minimise its own objective, pattern_control() rejects lambda_2 != 0; the
# penalized solver arrives with the spatial penalties.
.pattern_step_A <- function(X, Tm) {
  t(solve(crossprod(Tm), crossprod(Tm, X)))
}

# C-step: with Yw' Yw = n I the objective in C is const - 2 tr(C' M),
# M = Yw' X Psi^{-1} A, so the constrained minimiser is the polar factor of M.
.pattern_step_C <- function(X, Yw, A, noise) {
  M <- crossprod(Yw, X %*% .noise_apply_precision(noise, A))
  .polar_factor(M)
}

.polar_factor <- function(M) {
  sv <- svd(M)
  sv$u %*% t(sv$v)
}

# Supervised spectral initialization = exact reduced-rank regression of the
# whitened data X W' on Yw. Returns A, C for every rank up to r_max (nested).
.pattern_init_rrr <- function(X, Yw, noise, r_max) {
  n <- nrow(X)
  Xw <- .noise_whiten_rows(noise, X)                  # X W'      (n x p)
  B <- crossprod(Yw, Xw) / n                          # (Yw'Yw)^{-1} Yw' Xw = Yw'Xw / n  (q_eff x p)
  k <- min(r_max, ncol(Yw), ncol(X), n - 1L)
  # The fitted values are Yw %*% B, and Yw / sqrt(n) has orthonormal columns, so
  # (Yw B)'(Yw B) = n B'B: the right singular vectors of the n x p fitted matrix
  # are exactly those of the q_eff x p matrix B, with singular values scaled by
  # sqrt(n). Decomposing B instead avoids an n x p factorization (and its
  # min(n, p) x p workspace) at whole-brain feature counts.
  sv <- svd(B, nu = 0, nv = k)
  V <- sv$v[, seq_len(k), drop = FALSE]                # p x k, whitened-space directions
  A_all <- .noise_unwhiten(noise, V)                   # W^{-1} V
  C_raw <- B %*% V                                     # q_eff x k
  list(A_all = A_all, C_raw = C_raw,
       singular_values = sqrt(n) * sv$d[seq_len(k)], k = k)
}

# Truncate a nested initialization to rank r with C orthonormal and all scale in A.
.pattern_init_rank <- function(init, r) {
  A <- init$A_all[, seq_len(r), drop = FALSE]
  C <- init$C_raw[, seq_len(r), drop = FALSE]
  sv <- svd(C)
  C_orth <- sv$u %*% t(sv$v)
  # C = U S V'  =>  C A' = (U V') (V S V' A')  =>  A <- A V S V'
  A <- A %*% (sv$v %*% (sv$d * t(sv$v)))
  list(A = A, C = C_orth)
}

# ---------------------------------------------------------------------------
# The fit
# ---------------------------------------------------------------------------

#' Fit the pattern model on numeric matrices (single training set).
#'
#' @param X n x p numeric matrix (raw; centred internally).
#' @param targets factor, numeric vector, or n x q numeric matrix.
#' @param rank integer rank (<= eligible rank) or "path" to return nested
#'   unpenalized solutions for every rank up to `control$max_rank`.
#' @param control list from `pattern_control()`.
#' @param graph optional spatial_graph aligned to the columns of X (unused in
#'   Phase 2; stored for Phase 3).
#' @return A `pattern_fit`, or for rank = "path" a list of `pattern_fit`s.
#' @keywords internal
#' @noRd
.pattern_fit <- function(X, targets, rank = 1L, control = pattern_control(), graph = NULL,
                         cap_rank = FALSE, penalty = NULL, start = NULL) {
  X <- as.matrix(X)
  n <- nrow(X); p_input <- ncol(X)
  if (n < 3L) stop("pattern_model: at least three training observations are required.", call. = FALSE)

  # --- feature screening and transform (training rows only) ---
  keep <- .pattern_screen_columns(X)
  if (length(keep) == 0L) stop("pattern_model: no usable (finite, non-constant) features.", call. = FALSE)
  # avoid a full copy when screening drops nothing
  Xk <- if (length(keep) == ncol(X)) X else X[, keep, drop = FALSE]
  xt <- .pattern_x_transform_fit(Xk, scale = control$x_scale)
  Xc <- .pattern_x_transform_apply(xt, Xk)
  p <- ncol(Xc)

  # The graph must describe exactly the columns the estimator sees, so restrict
  # it whenever screening dropped features.
  if (!is.null(graph)) {
    if (identical(as.integer(graph$n_features), as.integer(ncol(X)))) {
      if (length(keep) != ncol(X)) graph <- restrict_graph(graph, keep)
    } else {
      .assert_graph_aligned(graph, p, "screened feature matrix")
    }
  }

  # --- target coding ---
  enc <- .pattern_encode_targets(targets, scale = control$y_scale)
  Yw <- enc$Yw; yt <- enc$transform
  r_elig <- min(yt$q_eff, p, n - 1L)
  r_max <- min(control$max_rank, r_elig)
  if (r_max < 1L) stop("pattern_model: eligible rank is zero.", call. = FALSE)

  # --- pilot fit (identity noise) -> residual covariance, held fixed ---
  pilot <- .pattern_init_rrr(Xc, Yw, .estimate_pattern_noise(Xc, type = "identity"), r_max)
  pil <- .pattern_init_rank(pilot, pilot$k)
  E <- Xc - (Yw %*% pil$C) %*% t(pil$A)
  noise <- .estimate_pattern_noise(
    E, type = control$noise$type, rank = control$noise$rank,
    max_rank = control$noise$max_rank, shrink = control$noise$shrink,
    df = max(n - pilot$k - 1L, 1L)
  )

  # --- supervised spectral initialization under Psi ---
  init <- .pattern_init_rrr(Xc, Yw, noise, r_max)

  # One O(n p q) pass supplies every quantity the alternation needs: X'Yw gives
  # both X'T (= X'Yw C) and the C-step's cross-product, and the constant term
  # completes the expanded objective. After this the whole alternation is
  # O(p q r + p r^2 + nnz(L) r) per iteration, with no n x p work at all.
  XtYw <- crossprod(Xc, Yw)
  quad_const <- .noise_quadform(noise, Xc)

  make_fit <- function(r, refine = TRUE, A_start = NULL) {
    st <- .pattern_init_rank(init, r)
    A <- if (!is.null(A_start) && identical(dim(A_start), c(p, r))) A_start else st$A
    C <- st$C
    # The penalty scale is resolved at this rank's initial C, so `sparse` means
    # the same fraction of "the penalty that empties the model" at every rank.
    pen <- .pattern_resolve_penalty(penalty, Xc, Yw %*% C, noise, graph,
                                    XtT = XtYw %*% C)
    build_prob <- function(Cmat) {
      .pattern_astep_problem(Xc, Yw %*% Cmat, noise, pen$L, pen$lambda_s, pen$lambda_l,
                             quad_const, XtT = XtYw %*% Cmat)
    }
    inner <- list()
    prob <- build_prob(C)
    obj <- .pattern_astep_objective(A, prob)
    trace <- obj
    converged <- FALSE
    iter <- 0L
    if (refine) {
      for (it in seq_len(control$max_outer)) {
        iter <- it
        C <- .polar_factor(crossprod(XtYw, .noise_apply_precision(noise, A)))
        prob <- build_prob(C)
        if (pen$active) {
          sol <- .pattern_astep_fista(A, prob, max_iter = control$max_inner,
                                      tol = control$tol_inner,
                                      tol_iterate = control$tol_iterate)
          A <- sol$A
          inner[[length(inner) + 1L]] <- list(iterations = sol$iterations,
                                              converged = sol$converged)
        } else {
          A <- t(solve(prob$TtT, t(prob$XtT)))
        }
        obj_new <- .pattern_astep_objective(A, prob)
        trace <- c(trace, obj_new)
        rel <- abs(obj - obj_new) / max(abs(obj), .Machine$double.eps)
        obj <- obj_new
        if (rel < control$tol) { converged <- TRUE; break }
      }
    }
    .new_pattern_fit(A, C, noise, xt, yt, Yw, keep, p_input, r, control,
                     diagnostics = list(objective = trace, outer_iterations = iter,
                                        converged = converged,
                                        inner = inner,
                                        init_singular_values = init$singular_values,
                                        noise_spectrum = noise$meta$spectrum,
                                        noise_rank = noise$h,
                                        rank_eligible = r_elig,
                                        n_nonzero = sum(rowSums(A^2) > 0)),
                     graph = graph, penalty = pen)
  }

  if (identical(rank, "path")) {
    # Each rank starts from its own spectral initialization. Warm starting rank
    # r from rank r-1 is not the free lunch it looks like: the objective is
    # invariant under (A, C) -> (A Q, C Q) for orthogonal Q, so a padded
    # lower-rank solution sits in an arbitrary rotation of the new coordinates
    # and was observed to leave the alternation short of convergence.
    refine <- control$refine_path || !is.null(penalty)
    fits <- lapply(seq_len(r_max), function(r) make_fit(r, refine = refine))
    names(fits) <- paste0("rank", seq_len(r_max))
    return(fits)
  }
  rank <- as.integer(rank)
  if (isTRUE(cap_rank)) rank <- min(rank, r_max)
  if (rank < 1L || rank > r_max) {
    stop(sprintf("pattern_model: requested rank %d exceeds the eligible rank %d.", rank, r_max), call. = FALSE)
  }
  make_fit(rank, refine = TRUE, A_start = start)
}

.new_pattern_fit <- function(A, C, noise, xt, yt, Yw, keep, p_input, rank, control,
                             diagnostics = list(), graph = NULL, penalty = NULL) {
  PA <- .noise_apply_precision(noise, A)          # Psi^{-1} A  (p x r)
  G <- crossprod(A, PA)                          # A' Psi^{-1} A (r x r)
  Tm <- Yw %*% C
  Phi <- crossprod(Tm) / nrow(Tm)                # Cov(T); identity up to numerical error
  structure(
    list(
      A = A, C = C, noise = noise,
      precision_A = PA, G = G,
      target_cov = Phi,
      x_transform = xt, y_transform = yt,
      feature_index = keep, p_input = p_input,
      rank = as.integer(rank),
      control = control,
      penalty = penalty,
      diagnostics = diagnostics,
      graph = graph,
      n_train = nrow(Tm)
    ),
    class = c("pattern_fit", "list")
  )
}

#' Control parameters for the pattern model estimator
#'
#' @param max_rank Maximum rank considered (capped at the eligible rank:
#'   number of classes minus one for categorical targets, the effective
#'   number of target dimensions otherwise, and never above the number of
#'   features or observations).
#' @param x_scale Feature scaling on training rows: \code{"none"} (centre
#'   only) or \code{"sd"}.
#' @param y_scale Continuous-target scaling before whitening: \code{"none"}
#'   or \code{"sd"}.
#' @param noise Residual covariance specification: a list with \code{type}
#'   (\code{"diag_lowrank"}, \code{"diag"}, or \code{"identity"}),
#'   \code{rank} (\code{"auto"} selects components above the
#'   Marchenko-Pastur edge, or an integer), \code{max_rank}, and
#'   \code{shrink} (shrinkage of the diagonal toward its median).
#' @param lambda_2 Reserved for the penalized solver; must be 0. A ridge on
#'   the patterns turns the A-step into a Sylvester equation under a general
#'   residual covariance, so it is rejected rather than solved approximately.
#' @param max_outer Maximum number of alternating (C-step, A-step) updates.
#' @param tol Relative objective change that declares convergence.
#' @param refine_path Logical; when fitting a rank path, run the alternating
#'   refinement for every rank (default \code{FALSE}: unpenalized path
#'   solutions are the exact reduced-rank optima, so refinement changes
#'   nothing). Refinement is always used when a spatial penalty is active.
#' @param max_inner Maximum proximal-gradient iterations per penalized A-step.
#' @param tol_inner Relative objective change that stops the A-step solver.
#' @param tol_iterate Relative change in the patterns required alongside
#'   \code{tol_inner}. The objective is flat near the optimum, so it can settle
#'   while the patterns are still moving; both must be small.
#' @return A list of class \code{pattern_control}.
#' @examples
#' pattern_control(max_rank = 3)
#' @export
pattern_control <- function(max_rank = 8L, x_scale = c("none", "sd"), y_scale = c("none", "sd"),
                            noise = list(type = "diag_lowrank", rank = "auto", max_rank = 10L, shrink = 0.1),
                            lambda_2 = 0, max_outer = 50L, tol = 1e-8, refine_path = FALSE,
                            max_inner = 500L, tol_inner = 1e-9, tol_iterate = 1e-6) {
  x_scale <- match.arg(x_scale)
  y_scale <- match.arg(y_scale)
  noise_default <- list(type = "diag_lowrank", rank = "auto", max_rank = 10L, shrink = 0.1)
  noise <- utils::modifyList(noise_default, as.list(noise))
  noise$type <- match.arg(noise$type, c("diag_lowrank", "diag", "identity"))
  if (!is.numeric(max_rank) || length(max_rank) != 1L || max_rank < 1) {
    stop("pattern_control: max_rank must be a positive integer.", call. = FALSE)
  }
  if (!is.numeric(lambda_2) || length(lambda_2) != 1L || !isTRUE(lambda_2 == 0)) {
    stop("pattern_control: lambda_2 must be 0; the penalized A-step is not ",
         "implemented in this version.", call. = FALSE)
  }
  structure(
    list(max_rank = as.integer(max_rank), x_scale = x_scale, y_scale = y_scale, noise = noise,
         lambda_2 = lambda_2, max_outer = as.integer(max_outer), tol = tol,
         refine_path = isTRUE(refine_path),
         max_inner = as.integer(max_inner), tol_inner = tol_inner,
         tol_iterate = tol_iterate),
    class = c("pattern_control", "list")
  )
}

#' @export
print.pattern_fit <- function(x, ...) {
  cat(sprintf("pattern_fit: rank %d, %d features (of %d input), %s targets (%d coded dims)\n",
              x$rank, nrow(x$A), x$p_input, x$y_transform$type, x$y_transform$q_eff))
  if (!is.null(x$penalty) && isTRUE(x$penalty$active)) {
    cat(sprintf("  penalty: sparse alpha = %.3g (lambda_s = %.4g), signed_smooth rho = %.3g; %d of %d features non-zero\n",
                x$penalty$alpha, x$penalty$lambda_s, x$penalty$rho,
                x$diagnostics$n_nonzero %||% NA_integer_, nrow(x$A)))
  }
  cat(sprintf("  noise: %s (h = %d), n_train = %d, converged = %s after %d outer iterations\n",
              x$noise$type, x$noise$h, x$n_train,
              if (isTRUE(x$diagnostics$converged)) "TRUE" else "FALSE",
              x$diagnostics$outer_iterations %||% 0L))
  invisible(x)
}
