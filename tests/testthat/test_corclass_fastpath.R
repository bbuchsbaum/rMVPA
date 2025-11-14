context("corclass fast path (pearson)")

test_that("corclass pearson matches base cor() for predictions", {
  set.seed(123)
  n <- 36; p <- 50; k <- 3
  y <- factor(rep(letters[1:k], each = n/k))
  x <- matrix(rnorm(n * p), n, p)
  newx <- matrix(rnorm(10 * p), 10, p)

  mod <- load_model("corclass")
  param <- data.frame(method = "pearson", robust = FALSE)
  fit <- mod$fit(x, y, wts = NULL, param = param, lev = levels(y), last = FALSE, weights = NULL, classProbs = TRUE)

  # Model predictions
  pred <- mod$predict(fit, newx)
  probs <- mod$prob(fit, newx)

  # Manual baseline via cor()
  cond_means <- group_means(x, 1, y)
  scores <- cor(t(newx), t(cond_means), method = "pearson")
  pred_ref <- factor(colnames(scores)[max.col(scores)], levels = levels(y))

  expect_equal(as.character(pred), as.character(pred_ref))
  expect_equal(nrow(probs), nrow(newx))
  expect_equal(colnames(probs), levels(y))
  rs <- rowSums(probs)
  expect_true(all(abs(rs - 1) < 1e-6))
})

test_that("corclass spearman path executes and returns sane shapes", {
  set.seed(456)
  n <- 24; p <- 30; k <- 3
  y <- factor(rep(letters[1:k], length.out = n))
  x <- matrix(rnorm(n * p), n, p)
  newx <- matrix(rnorm(7 * p), 7, p)
  mod <- load_model("corclass")
  param <- data.frame(method = "spearman", robust = FALSE)
  fit <- mod$fit(x, y, wts = NULL, param = param, lev = levels(y), last = FALSE, weights = NULL, classProbs = TRUE)
  pred <- mod$predict(fit, newx)
  probs <- mod$prob(fit, newx)
  expect_equal(nrow(probs), nrow(newx))
  expect_equal(colnames(probs), levels(y))
  expect_equal(length(pred), nrow(newx))
})
