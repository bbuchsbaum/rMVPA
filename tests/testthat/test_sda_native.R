# The native sda fit (R/sda_native.R) must reproduce sda::sda()'s estimator:
# the same shrinkage intensities, discriminant coefficients and posteriors.

ref_cases <- readRDS(test_path("fixtures", "sda_reference.rds"))

for (nm in names(ref_cases)) {
  local({
    case <- ref_cases[[nm]]
    test_that(sprintf("native sda matches frozen sda::sda output (%s)", nm), {
      fit <- rMVPA:::.sda_native_fit(case$X, case$y)
      expect_s3_class(fit, "sda_native")
      expect_equal(as.numeric(fit$regularization), as.numeric(case$regularization), tolerance = 1e-10)
      expect_equal(as.numeric(fit$alpha), as.numeric(case$alpha), tolerance = 1e-9)
      expect_equal(as.numeric(fit$beta), as.numeric(case$beta), tolerance = 1e-9)
      expect_equal(unname(rMVPA:::.sda_native_posterior(fit, case$Xt)),
                   unname(case$posterior), tolerance = 1e-9)
    })
  })
}

test_that("native sda matches live sda::sda on random problems", {
  skip_if_not_installed("sda")
  set.seed(77)
  for (rep in 1:10) {
    n <- sample(20:90, 1); p <- sample(5:150, 1); K <- sample(2:5, 1)
    y <- factor(sample(seq_len(K), n, replace = TRUE))
    if (any(table(y) < 2)) next
    X <- matrix(rnorm(n * p), n) * rep(runif(p, 0.2, 5), each = n) + outer(as.integer(y), rnorm(p))
    ref <- suppressWarnings(sda::sda(X, y, verbose = FALSE))
    fit <- suppressWarnings(rMVPA:::.sda_native_fit(X, y))
    if (is.null(fit)) next  # delegated case
    Xt <- matrix(rnorm(8 * p), 8)
    expect_equal(unname(fit$regularization), unname(ref$regularization), tolerance = 1e-10)
    expect_equal(unname(rMVPA:::.sda_native_posterior(fit, Xt)),
                 unname(predict(ref, Xt, verbose = FALSE)$posterior), tolerance = 1e-9)
  }
})

test_that("sda_notune uses the native fit and needs no sda package for it", {
  set.seed(78)
  y <- factor(rep(c("a", "b", "c"), 15))
  X <- matrix(rnorm(45 * 20), 45) + outer(as.integer(y), rnorm(20))
  model <- load_model("sda_notune")
  fit <- model$fit(X, y, NULL, NULL, levels(y), NULL, NULL, TRUE)
  expect_s3_class(fit, "sda_native")
  probs <- model$prob(fit, X[1:5, ])
  expect_equal(dim(probs), c(5L, 3L))
  expect_equal(unname(rowSums(probs)), rep(1, 5), tolerance = 1e-6)
  expect_identical(levels(model$predict(fit, X[1:5, ])), levels(y))
})


test_that("native SDA retains weight extraction and training-data importance", {
  skip_if_not_installed("sda")
  set.seed(79)
  y <- factor(rep(c("a", "b"), 30))
  X <- matrix(rnorm(60 * 12), 60, 12,
              dimnames = list(NULL, paste0("v", 1:12)))
  fit <- rMVPA:::.sda_native_fit(X, y)
  ref <- sda::sda(X, y, verbose = FALSE)
  # The backend adds a shrinkage class to its coefficient matrix.
  expect_equal(as.numeric(extract_weights(fit)), as.numeric(extract_weights(ref)),
               tolerance = 1e-9)
  expect_equal(dim(extract_weights(fit)), dim(extract_weights(ref)))
  expect_equal(model_importance(fit, X), model_importance(ref, X), tolerance = 1e-9)
  expect_equal(rownames(extract_weights(fit)), colnames(X))
  expect_equal(colnames(extract_weights(fit)), levels(y))
})
