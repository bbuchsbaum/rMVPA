# Native shrinkage discriminant analysis: the estimator of
# sda::sda(Xtrain, L, diagonal = FALSE) with all shrinkage intensities
# estimated (Ahdesmaki & Strimmer 2010; Schafer & Strimmer 2005), without the
# sda/corpcor/entropy dependency chain.
#
# sda applies the shrunk inverse correlation matrix through two SVDs of the
# centred data plus an eigendecomposition. Both uses reduce to the n x n Gram
# matrix G = Xs Xs' of the standardised, centred training data:
#   * the correlation shrinkage intensity needs sum(d^4) over the singular
#     values, which is ||G||_F^2;
#   * (lambda I + (1 - lambda) R)^-1 with R = Xs'Xs / (n - 1) of rank <= n is,
#     by the Woodbury identity,
#     (1/lambda) [I - Xs' (lambda (n - 1) / (1 - lambda) I + G)^-1 Xs],
#     an n x n symmetric positive definite solve (one Cholesky). When there
#     are fewer columns than rows, the p x p primal system is smaller and is
#     solved directly instead.
# Moments follow corpcor::wt.moments() exactly (one-pass variance, values
# below .Machine$double.eps set to zero), so degenerate columns are treated as
# sda treats them. The only case not covered, an estimated lambda of exactly
# zero (sda then uses a pseudoinverse), is delegated to sda::sda().

#' corpcor::wt.moments() with uniform weights
#' @keywords internal
#' @noRd
.sda_moments <- function(x) {
  n <- nrow(x)
  w <- 1 / n
  m <- colSums(w * x)
  v <- (1 / (1 - n * w * w)) * (colSums(w * x^2) - m^2)
  v[v < .Machine$double.eps] <- 0
  list(mean = m, var = v)
}

#' Per-column statistics of the sda fit
#'
#' Everything sda computes column by column (class means, pooled mean,
#' within-class centring, moments, standardisation, the lambda.var
#' ingredients). Because each is per column, the statistics of any column
#' subset are the subset of these, bit for bit. A searchlight engine can
#' therefore compute them once for the whole brain per fold.
#' @keywords internal
#' @noRd
.sda_column_stats <- function(X, L) {
  X <- as.matrix(X)
  y <- factor(L)
  lev <- levels(y)
  K <- length(lev)
  n <- nrow(X)
  yi <- as.integer(y)
  nk <- tabulate(yi, K)

  # Class frequencies shrunk toward uniform (entropy::freqs.shrink).
  u <- nk / n
  tf <- 1 / K
  msp <- sum((u - tf)^2)
  lf <- if (n <= 1 || msp == 0) 1 else min(1, max(0, sum(u * (1 - u) / (n - 1)) / msp))
  freqs <- lf * tf + (1 - lf) * u

  # Class means, pooled mean, within-class centred data (sda::centroids).
  mu <- matrix(0, K, ncol(X))
  for (k in seq_len(K)) mu[k, ] <- colMeans(X[yi == k, , drop = FALSE])
  mup <- colSums(mu * freqs)
  xc <- X - mu[yi, , drop = FALSE]

  # corpcor moments and the lambda.var ingredients.
  mom <- .sda_moments(xc)
  xc0 <- sweep(xc, 2, mom$mean)
  zz <- xc0^2
  h1 <- n / (n - 1)
  q1 <- colSums(zz / n)
  sdv <- sqrt(mom$var)
  zero <- sdv == 0
  xs <- sweep(xc0, 2, ifelse(zero, 1, sdv), "/")
  xs[, zero] <- 0
  list(
    lev = lev, K = K, n = n, h1 = h1, freqs = freqs, lf = lf,
    mu = mu, mup = mup, v = mom$var, v_lv = h1 * colSums(zz / n),
    q1 = q1, q2 = colSums(zz^2 / n) - q1^2,
    zero = zero, xs = xs, cs2 = colSums(xs^2 / n), a = xs^2 / sqrt(n),
    colnames = colnames(X)
  )
}

#' Finish the sda fit on a subset of columns
#'
#' The cross-column steps: variance and correlation shrinkage intensities and
#' the Woodbury solve. Returns NULL when the fit must be delegated to
#' sda::sda() (estimated correlation shrinkage of exactly zero).
#' @keywords internal
#' @noRd
.sda_fit_columns <- function(st, cols = seq_along(st$v)) {
  n <- st$n
  K <- st$K
  h1 <- st$h1
  p <- length(cols)
  v <- st$v[cols]

  # Pooled variances shrunk toward their median (corpcor::var.shrink).
  target <- stats::median(v)
  # estimate.lambda.var() uses its own two-pass variance for its target.
  target_lv <- stats::median(st$v_lv[cols])
  den_v <- sum((st$q1[cols] - target_lv / h1)^2)
  lv <- if (den_v == 0) 1 else max(0, min(1, sum(st$q2[cols]) / den_v / (n - 1)))
  sc <- sqrt((lv * target + (1 - lv) * v) * (n - 1) / (n - K))

  mu <- st$mu[, cols, drop = FALSE]
  mup <- st$mup[cols]
  pw <- t((mu - rep(mup, each = K)) / rep(sc, each = K))  # p x K

  # Correlation shrinkage intensity (corpcor estimate.lambda) from the Gram.
  zero <- st$zero[cols]
  xs <- st$xs[, cols, drop = FALSE]
  # Work in whichever Gram is smaller: ||Xs Xs'||_F = ||Xs' Xs||_F, and the
  # shrunk solve is primal (p x p) for p < n, Woodbury (n x n) otherwise.
  primal <- p < n
  lam <- 1
  G <- NULL
  if (p > 1L) {
    G <- if (primal) crossprod(xs) else tcrossprod(xs)
    sE2R <- sum(G^2) / n^2 - sum(st$cs2[cols]^2)
    a <- st$a[, cols, drop = FALSE]
    sER2 <- sum(rowSums(a)^2 - rowSums(a^2))
    lam <- if (sE2R == 0) 1 else max(0, min(1, (sER2 - sE2R) / sE2R / (n - 1)))
  }
  if (lam == 0) return(NULL)

  chol_solve <- function(A, B) {
    R <- chol(A)
    backsolve(R, backsolve(R, B, transpose = TRUE))
  }
  cp <- if (lam == 1) {
    pw
  } else if (primal) {
    # (lambda I + (1 - lambda) Xs'Xs / (n - 1)) cp = pw
    A <- G * ((1 - lam) / (n - 1))
    diag(A) <- diag(A) + lam
    chol_solve(A, pw)
  } else {
    # Woodbury on the n x n Gram.
    A <- G
    diag(A) <- diag(A) + lam * (n - 1) / (1 - lam)
    (pw - crossprod(xs, chol_solve(A, xs %*% pw))) / lam
  }
  cp[zero, ] <- pw[zero, ]
  cp <- cp / sc

  alpha <- log(st$freqs) - colSums(cp * (t(mu) + mup)) / 2
  names(alpha) <- st$lev
  beta <- t(cp)
  rownames(beta) <- st$lev
  colnames(beta) <- st$colnames[cols]
  structure(
    list(regularization = c(lambda = lam, lambda.var = lv, lambda.freqs = st$lf),
         freqs = stats::setNames(st$freqs, st$lev), alpha = alpha, beta = beta),
    class = "sda_native"
  )
}

#' Fit the sda estimator natively
#'
#' Returns NULL when the fit must be delegated to sda::sda().
#' @keywords internal
#' @noRd
.sda_native_fit <- function(X, L) {
  X <- as.matrix(X)
  if (nrow(X) < 3L || nlevels(factor(L)) < 2L) return(NULL)
  .sda_fit_columns(.sda_column_stats(X, L))
}

#' Posterior probabilities, computed as sda's predict.sda() does
#' @keywords internal
#' @noRd
.sda_native_posterior <- function(fit, Xtest) {
  Xtest <- as.matrix(Xtest)
  probs <- t(tcrossprod(fit$beta, Xtest) + fit$alpha)
  probs <- exp(probs - probs[cbind(seq_len(nrow(probs)), max.col(probs, ties.method = "first"))])
  probs <- zapsmall(probs / rowSums(probs))
  colnames(probs) <- names(fit$alpha)
  rownames(probs) <- rownames(Xtest)
  probs
}
