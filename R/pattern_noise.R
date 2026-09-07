# Structured residual covariance for pattern_model
#
# Psi = D + U U'  with D positive diagonal (p) and U p x h (h small).
# Nothing here ever forms a dense p x p matrix. All operations are expressed
# through D, U, and a few h x h matrices so they restrict cheaply to any subset
# of features (Psi_RR = D_R + U_R U_R').
#
# The object is a plain list (no closures) so it serializes without dragging
# training data along.

#' @keywords internal
#' @noRd
new_pattern_noise <- function(D, U = NULL, type = "diag_lowrank", meta = list()) {
  D <- as.numeric(D)
  if (any(!is.finite(D)) || any(D <= 0)) {
    stop("pattern noise: diagonal variances must be finite and positive.", call. = FALSE)
  }
  p <- length(D)
  if (!is.null(U)) {
    U <- as.matrix(U)
    if (nrow(U) != p) stop("pattern noise: U must have one row per feature.", call. = FALSE)
    if (ncol(U) == 0L) U <- NULL
  }
  obj <- list(type = type, D = D, U = U, h = if (is.null(U)) 0L else ncol(U), p = p, meta = meta)
  .noise_precompute(obj)
}

# Cache the small matrices used by Woodbury and the symmetric square roots.
.noise_precompute <- function(obj) {
  if (is.null(obj$U)) {
    obj$Kinv <- NULL; obj$P <- NULL; obj$s <- NULL
    return(obj)
  }
  Dinv_U <- obj$U / obj$D                       # D^{-1} U   (p x h)
  K <- diag(obj$h) + crossprod(obj$U, Dinv_U)   # I + U' D^{-1} U
  obj$Kinv <- solve(K)
  # V = D^{-1/2} U = P diag(s) Q'  ->  (I + V V')^{+-1/2} = I + P (f(s) - 1) P'
  V <- obj$U / sqrt(obj$D)
  sv <- svd(V, nu = obj$h, nv = 0)
  obj$P <- sv$u
  obj$s <- sv$d
  obj
}

#' Apply Psi^{-1} to a matrix with p rows (Woodbury).
#' @keywords internal
#' @noRd
.noise_apply_precision <- function(noise, M) {
  M <- as.matrix(M)
  DinvM <- M / noise$D
  if (is.null(noise$U)) return(DinvM)
  Dinv_U <- noise$U / noise$D
  DinvM - Dinv_U %*% (noise$Kinv %*% crossprod(Dinv_U, M))
}

#' Apply Psi to a matrix with p rows.
#' @keywords internal
#' @noRd
.noise_apply <- function(noise, M) {
  M <- as.matrix(M)
  out <- M * noise$D
  if (!is.null(noise$U)) out <- out + noise$U %*% crossprod(noise$U, M)
  out
}

#' diag(Psi^{-1}) without forming Psi^{-1}.
#' @keywords internal
#' @noRd
.noise_precision_diag <- function(noise) {
  d <- 1 / noise$D
  if (is.null(noise$U)) return(d)
  Dinv_U <- noise$U / noise$D
  d - rowSums((Dinv_U %*% noise$Kinv) * Dinv_U)
}

#' Whitening operator W with W' W = Psi^{-1}:  W M = (I + V V')^{-1/2} D^{-1/2} M.
#' Any such W makes ||R W'||_F^2 = tr(R Psi^{-1} R'), which is all the
#' estimator needs; W need not be the symmetric square root.
#' @keywords internal
#' @noRd
.noise_whiten <- function(noise, M) {
  M <- as.matrix(M) / sqrt(noise$D)
  if (is.null(noise$U)) return(M)
  f <- 1 / sqrt(1 + noise$s^2) - 1
  M + noise$P %*% (f * crossprod(noise$P, M))
}

#' Inverse of .noise_whiten:  W^{-1} M = D^{1/2} (I + V V')^{1/2} M.
#' @keywords internal
#' @noRd
.noise_unwhiten <- function(noise, M) {
  M <- as.matrix(M)
  if (!is.null(noise$U)) {
    f <- sqrt(1 + noise$s^2) - 1
    M <- M + noise$P %*% (f * crossprod(noise$P, M))
  }
  M * sqrt(noise$D)
}

#' Restrict the noise model to a feature subset (exact marginal covariance).
#' @keywords internal
#' @noRd
.noise_restrict <- function(noise, keep) {
  if (is.logical(keep)) keep <- which(keep)
  keep <- as.integer(keep)
  U <- if (is.null(noise$U)) NULL else noise$U[keep, , drop = FALSE]
  new_pattern_noise(noise$D[keep], U, type = noise$type, meta = noise$meta)
}

#' Quadratic form tr(R Psi^{-1} R') for R with p columns (rows = observations).
#' @keywords internal
#' @noRd
.noise_quadform <- function(noise, R) {
  Rt <- t(as.matrix(R))
  sum(Rt * .noise_apply_precision(noise, Rt))
}

#' Estimate Psi = D + U U' from residuals (training rows only).
#'
#' D is the shrunken residual variance; U comes from the leading principal
#' components of the standardized residuals whose eigenvalues exceed the
#' Marchenko-Pastur noise edge, capped at max_rank. The variance carried by U
#' is removed from D (with a floor) so Psi does not double count it.
#' @keywords internal
#' @noRd
.estimate_pattern_noise <- function(E, type = c("diag_lowrank", "diag", "identity"),
                                    rank = "auto", max_rank = 10L, shrink = 0.1,
                                    df = NULL) {
  type <- match.arg(type)
  E <- as.matrix(E)
  n <- nrow(E); p <- ncol(E)
  if (is.null(df)) df <- max(n - 1, 1)

  if (type == "identity") {
    return(new_pattern_noise(rep(1, p), NULL, type = "identity",
                             meta = list(df = df, spectrum = numeric(0))))
  }

  v <- colSums(E^2) / df
  v[!is.finite(v)] <- 0
  med <- stats::median(v[v > 0])
  if (!is.finite(med) || med <= 0) med <- 1
  shrink <- min(max(shrink, 0), 1)
  D <- (1 - shrink) * v + shrink * med
  D <- pmax(D, 1e-8 * med)

  if (type == "diag" || max_rank < 1L || n < 3L) {
    return(new_pattern_noise(D, NULL, type = "diag",
                             meta = list(df = df, spectrum = numeric(0), h_selected = 0L)))
  }

  # Standardized residuals: their covariance is I + (signal beyond D).
  Es <- sweep(E, 2L, sqrt(D), "/")
  k <- min(max_rank, n - 1L, p)
  # With p >> n, decompose the n x n Gram matrix instead of the n x p residual
  # matrix: Es Es' = U D^2 U', and the right singular vectors follow as
  # V = Es' U / d. Only the leading components (well separated from the noise
  # bulk) are used, so squaring the condition number is not a concern here.
  sv <- if (p > 2L * n) {
    eg <- eigen(tcrossprod(Es), symmetric = TRUE)
    dvals <- sqrt(pmax(eg$values[seq_len(k)], 0))
    pos <- dvals > max(dvals[1], .Machine$double.eps) * 1e-8
    V <- matrix(0, p, k)
    if (any(pos)) {
      V[, pos] <- crossprod(Es, eg$vectors[, seq_len(k), drop = FALSE][, pos, drop = FALSE]) %*%
        diag(1 / dvals[pos], nrow = sum(pos))
    }
    list(d = dvals, v = V)
  } else {
    svd(Es, nu = 0, nv = k)
  }
  ev <- sv$d[seq_len(k)]^2 / df                # eigenvalues of standardized residual cov
  mp_edge <- (1 + sqrt(p / df))^2              # Marchenko-Pastur upper edge for unit variance
  h <- if (identical(rank, "auto")) {
    sum(ev > mp_edge)
  } else {
    min(as.integer(rank), k)
  }
  h <- max(min(h, k), 0L)

  U <- NULL
  if (h > 0L) {
    excess <- pmax(ev[seq_len(h)] - 1, 0)
    U <- sweep(sv$v[, seq_len(h), drop = FALSE], 2L, sqrt(excess), "*")
    U <- U * sqrt(D)                           # back to the original scale
    D <- pmax(D - rowSums(U^2), 0.1 * D)       # avoid double counting; keep a floor
    if (all(excess == 0)) U <- NULL
  }
  new_pattern_noise(D, U, type = "diag_lowrank",
                    meta = list(df = df, spectrum = ev, mp_edge = mp_edge,
                                h_selected = if (is.null(U)) 0L else ncol(U)))
}
