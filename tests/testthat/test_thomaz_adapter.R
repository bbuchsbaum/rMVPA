test_that("Thomaz adapter preserves the backend prediction and probability contracts", {
  skip_if_not_installed("sparsediscrim")
  set.seed(20261006)
  x <- matrix(rnorm(48 * 10), 48, 10,
              dimnames = list(NULL, paste0("v", 1:10)))
  y <- factor(rep(c("a", "b", "c"), 16))
  model <- load_model("lda_thomaz")
  fit <- model$fit(x, y, NULL, data.frame(parameter = "none"), levels(y), TRUE, NULL, TRUE)
  newdata <- x[1:7, , drop = FALSE]
  expect_identical(model$predict(fit, newdata),
                   factor(predict(fit, newdata, type = "class"), levels = levels(y)))
  expect_equal(model$prob(fit, newdata),
               as.matrix(predict(fit, newdata, type = "prob")))
  expect_identical(colnames(model$prob(fit, newdata)), levels(y))

  boot <- load_model("lda_thomaz_boot")
  bfit <- boot$fit(x, y, NULL, data.frame(reps = 2, frac = 0.8), levels(y), TRUE, NULL, TRUE)
  probabilities <- boot$prob(bfit, newdata)
  expect_equal(dim(probabilities), c(7L, 3L))
  expect_true(all(is.finite(probabilities)))
  # The adapter rounds each bootstrap probability with zapsmall (7 digits).
  expect_equal(unname(rowSums(probabilities)), rep(1, 7), tolerance = 1e-6)
  expect_equal(as.character(boot$predict(bfit, newdata)),
               levels(y)[max.col(probabilities, ties.method = "first")])
})


test_that("Thomaz adapters accept unnamed ROI matrices without losing supplied names", {
  skip_if_not_installed("sparsediscrim")
  set.seed(20261007)
  x <- matrix(rnorm(48 * 10), 48, 10)
  y <- factor(rep(c("a", "b", "c"), 16))
  for (name in c("lda_thomaz", "lda_thomaz_boot")) {
    model <- load_model(name)
    fit <- model$fit(x, y, NULL, data.frame(reps = 2, frac = 0.8),
                     levels(y), TRUE, NULL, TRUE)
    expect_equal(dim(model$prob(fit, x[1:7, , drop = FALSE])), c(7L, 3L))
    expect_length(model$predict(fit, x[1:7, , drop = FALSE]), 7L)
  }
  colnames(x) <- paste0("voxel_", 11:20)
  expect_identical(rMVPA:::.thomaz_matrix(x), x)
})
