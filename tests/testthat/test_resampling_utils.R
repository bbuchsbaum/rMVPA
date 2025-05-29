context("resampling_utils")

test_that("create_mvpa_folds returns correct number of folds", {
  y <- rep(1:3, each = 10)
  folds <- rMVPA:::create_mvpa_folds(y, k = 5, list = TRUE)
  
  expect_equal(length(folds), 5)
  expect_equal(sum(lengths(folds)), length(y))
  
  # Check all indices are included
  all_indices <- sort(unlist(folds))
  expect_equal(all_indices, seq_len(length(y)))
})

test_that("create_mvpa_folds with list=FALSE returns fold assignments", {
  y <- rep(1:2, each = 10)
  fold_vec <- rMVPA:::create_mvpa_folds(y, k = 5, list = FALSE, seed = 123)
  
  expect_type(fold_vec, "integer")
  expect_equal(length(fold_vec), length(y))
  expect_equal(sort(unique(fold_vec)), 1:5)
  
  # Check each fold has observations
  for (i in 1:5) {
    expect_true(sum(fold_vec == i) > 0)
  }
})

test_that("create_mvpa_folds with stratification works for factors", {
  y <- factor(rep(c("A", "B", "C"), each = 10))
  folds <- rMVPA:::create_mvpa_folds(y, k = 5, list = TRUE, seed = 456)
  
  # Check stratification - each fold should have roughly equal class distribution
  for (fold in folds) {
    fold_classes <- table(y[fold])
    # Each class should be represented (approximately)
    expect_true(all(fold_classes > 0))
  }
})

test_that("create_mvpa_folds seed produces reproducible results", {
  y <- 1:20
  
  folds1 <- rMVPA:::create_mvpa_folds(y, k = 4, seed = 789)
  folds2 <- rMVPA:::create_mvpa_folds(y, k = 4, seed = 789)
  folds3 <- rMVPA:::create_mvpa_folds(y, k = 4, seed = 999)
  
  expect_identical(folds1, folds2)
  expect_false(identical(folds1, folds3))
})

test_that("create_mvpa_folds handles edge cases", {
  # k equals n
  y <- 1:5
  folds <- rMVPA:::create_mvpa_folds(y, k = 5, list = TRUE)
  expect_equal(length(folds), 5)
  expect_true(all(lengths(folds) == 1))
  
  # Very small n
  y <- 1:3
  folds <- rMVPA:::create_mvpa_folds(y, k = 3, list = TRUE)
  expect_equal(length(folds), 3)
})