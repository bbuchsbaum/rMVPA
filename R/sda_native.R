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
#     an n x n symmetric positive definite solve (one Cholesky).
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

#' Fit the sda estimator natively
#'
#' Returns NULL when the fit must be delegated to sda::sda().
#' @keywords internal
#' @noRd
.sda_native_fit <- function(X, L) {
  X <- as.matrix(X)
  y <- factor(L)
  lev <- levels(y)
  K <- length(lev)
  n <- nrow(X)
  p <- ncol(X)
  if (n < 3L || K < 2L) return(NULL)
  yi <- as.integer(y)
  nk <- tabulate(yi, K)

  # Class frequencies shrunk toward uniform (entropy::freqs.shrink).
  u <- nk / n
  tf <- 1 / K
  msp <- sum((u - tf)^2)
  lf <- if (n <= 1 || msp == 0) 1 else min(1, max(0, sum(u * (1 - u) / (n - 1)) / msp))
  freqs <- lf * tf + (1 - lf) * u

  # Class means, pooled mean, within-class centred data (sda::centroids).
  mu <- matrix(0, K, p)
  for (k in seq_len(K)) mu[k, ] <- colMeans(X[yi == k, , drop = FALSE])
  mup <- colSums(mu * freqs)
  xc <- X - mu[yi, , drop = FALSE]

  # Pooled variances shrunk toward their median (corpcor::var.shrink).
  mom <- .sda_moments(xc)
  v <- mom$var
  target <- stats::median(v)
  xc0 <- sweep(xc, 2, mom$mean)
  zz <- xc0^2
  h1 <- n / (n - 1)
  # estimate.lambda.var() uses its own two-pass variance for its target.
  target_lv <- stats::median(h1 * colSums(zz / n))
  q1 <- colSums(zz / n)
  q2 <- colSums(zz^2 / n) - q1^2
  den_v <- sum((q1 - target_lv / h1)^2)
  lv <- if (den_v == 0) 1 else max(0, min(1, sum(q2) / den_v / (n - 1)))
  sc <- sqrt((lv * target + (1 - lv) * v) * (n - 1) / (n - K))

  pw <- t((mu - rep(mup, each = K)) / rep(sc, each = K))  # p x K

  # Correlation shrinkage intensity (corpcor estimate.lambda) from the Gram.
  sdv <- sqrt(v)
  zero <- sdv == 0
  xs <- sweep(xc0, 2, ifelse(zero, 1, sdv), "/")
  xs[, zero] <- 0
  lam <- 1
  G <- NULL
  if (p > 1L) {
    G <- tcrossprod(xs)
    sE2R <- sum(G^2) / n^2 - sum(colSums(xs^2 / n)^2)
    a <- xs^2 / sqrt(n)
    sER2 <- sum(rowSums(a)^2 - rowSums(a^2))
    lam <- if (sE2R == 0) 1 else max(0, min(1, (sER2 - sE2R) / sE2R / (n - 1)))
  }
  if (lam == 0) return(NULL)

  # Shrunk inverse correlation times pw, by Woodbury on the n x n Gram.
  cp <- if (lam == 1) {
    pw
  } else {
    A <- G
    diag(A) <- diag(A) + lam * (n - 1) / (1 - lam)
    R <- chol(A)
    (pw - crossprod(xs, backsolve(R, forwardsolve(t(R), xs %*% pw)))) / lam
  }
  cp[zero, ] <- pw[zero, ]
  cp <- cp / sc

  alpha <- log(freqs) - colSums(cp * (t(mu) + mup)) / 2
  names(alpha) <- lev
  beta <- t(cp)
  rownames(beta) <- lev
  colnames(beta) <- colnames(X)
  structure(
    list(regularization = c(lambda = lam, lambda.var = lv, lambda.freqs = lf),
         freqs = stats::setNames(freqs, lev), alpha = alpha, beta = beta),
    class = "sda_native"
  )
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
