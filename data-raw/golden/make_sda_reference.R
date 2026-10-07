# Frozen sda::sda() outputs for tests/testthat/test_sda_native.R, so parity
# with the reference implementation is checked even where sda is absent.
# Run from the package root: Rscript data-raw/golden/make_sda_reference.R
suppressMessages(library(sda))
make_case <- function(seed, n, p, K, prob = NULL, scale = 1, offset = 0, zero_col = FALSE) {
  set.seed(seed)
  y <- factor(sample(letters[seq_len(K)], n, replace = TRUE, prob = prob))
  X <- offset + scale * matrix(rnorm(n * p), n) + outer(as.integer(y), rnorm(p, sd = 0.5))
  if (zero_col) X[, 3] <- as.integer(y)  # varies across classes only
  Xt <- offset + scale * matrix(rnorm(12 * p), 12)
  fit <- sda(X, y, verbose = FALSE)
  list(X = X, y = y, Xt = Xt, regularization = fit$regularization,
       alpha = fit$alpha, beta = fit$beta,
       posterior = predict(fit, Xt, verbose = FALSE)$posterior)
}
cases <- list(
  balanced = make_case(1, 40, 30, 4),
  unbalanced = make_case(2, 45, 60, 3, prob = c(.6, .3, .1)),
  n_gt_p = make_case(3, 80, 12, 2),
  raw_scale = make_case(4, 40, 30, 3, scale = 20, offset = 1500),
  within_class_constant = make_case(5, 40, 25, 4, zero_col = TRUE)
)
attr(cases, "sda_version") <- as.character(packageVersion("sda"))
saveRDS(cases, file.path("tests", "testthat", "fixtures", "sda_reference.rds"), version = 2)
