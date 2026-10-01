# The vectorised Gaussian naive Bayes must agree with the per-feature
# implementation it replaced.
ref_nb_fit <- function(x, y) {
  classes <- levels(y)
  mus <- vars <- matrix(NA, length(classes), ncol(x), dimnames = list(classes, NULL))
  for (k in seq_along(classes)) {
    s <- x[y == classes[k], , drop = FALSE]
    nk <- nrow(s)
    mus[k, ] <- colMeans(s)
    v <- apply(s, 2, var) * (nk - 1) / nk
    z <- v <= .Machine$double.eps
    if (any(z)) {
      nz <- v[v > .Machine$double.eps]
      v[z] <- max(if (length(nz)) min(nz) * 1e-6 else 1e-10, 1e-10)
    }
    vars[k, ] <- v
  }
  list(mus = mus, vars = vars, log_class_priors = as.numeric(log(table(y) / length(y))),
       classes = classes)
}
ref_log_post <- function(fit, newdata) {
  out <- matrix(NA, nrow(newdata), length(fit$classes))
  for (k in seq_along(fit$classes)) {
    ll <- sapply(seq_len(ncol(newdata)), function(j)
      dnorm(newdata[, j], fit$mus[k, j], sqrt(fit$vars[k, j]), log = TRUE))
    ll <- matrix(ll, nrow(newdata))
    ll[!is.finite(ll)] <- -1e100
    out[, k] <- rowSums(ll) + fit$log_class_priors[k]
  }
  out
}

test_that("vectorised naive Bayes matches the per-feature reference", {
  nb <- rMVPA:::MVPAModels$naive_bayes
  set.seed(5501)
  for (rep in 1:25) {
    n <- 30; p <- sample(2:40, 1)
    y <- factor(sample(letters[1:3], n, replace = TRUE), levels = letters[1:3])
    if (any(table(y) < 2)) next
    x <- matrix(rnorm(n * p, sd = runif(1, 0.1, 50)), n, p)
    if (p > 3) x[y == "a", 2] <- 7  # zero variance within one class
    newdata <- matrix(rnorm(10 * p), 10, p)
    newdata[sample(length(newdata), 3)] <- NA

    fit <- suppressWarnings(nb$fit(x, y, NULL, NULL, levels(y), NULL, NULL, TRUE))
    ref <- ref_nb_fit(x, y)
    # Bit-identical, not merely close: downstream rank metrics (AUC) can
    # reorder near-tied tiny probabilities on last-bit differences.
    expect_identical(fit$mus, ref$mus)
    expect_identical(fit$vars, ref$vars)

    lp <- suppressWarnings(rMVPA:::calculate_log_posteriors(fit, newdata))
    expect_identical(unname(lp), ref_log_post(ref, newdata))

    probs <- suppressWarnings(nb$prob(fit, newdata))
    ref_probs <- t(apply(lp, 1, function(r) exp(r - max(r)) / sum(exp(r - max(r)))))
    expect_equal(unname(probs), unname(ref_probs), tolerance = 1e-14)
    expect_identical(colnames(probs), levels(y))
  }
})
