library(testthat)
library(rMVPA)
library(modelr)

context("cross-validation helper functions")

# Helper data
set.seed(1)
df <- data.frame(x = rnorm(20), y = rnorm(20))
y <- df$y

# ---- crossv_k ----

test_that("crossv_k creates k distinct folds", {
  res <- crossv_k(df, y, k = 5)
  expect_s3_class(res, "tbl_df")
  expect_equal(nrow(res), 5)
  expect_equal(res$.id, sprintf("%02d", 1:5))

  for (i in seq_len(nrow(res))) {
    train_idx <- res$train[[i]]$idx
    test_idx  <- res$test[[i]]$idx

    expect_length(intersect(train_idx, test_idx), 0)
    expect_equal(length(train_idx) + length(test_idx), nrow(df))
    expect_equal(length(res$ytrain[[i]]) + length(res$ytest[[i]]), length(y))
  }
})


test_that("crossv_k validates k argument", {
  expect_error(crossv_k(df, y, k = 1), "at least 2")
  expect_error(crossv_k(df, y, k = c(2,3)), "single integer")
})

# ---- crossv_twofold ----

test_that("crossv_twofold splits blocks by half", {
  block_var <- rep(1:4, each = 5)
  res <- crossv_twofold(df, y, block_var, nreps = 3)
  expect_s3_class(res, "tbl_df")
  expect_equal(nrow(res), 3)
  expect_equal(res$.id, sprintf("%02d", 1:3))

  half_blocks <- length(unique(block_var)) / 2
  for (i in seq_len(nrow(res))) {
    test_idx  <- res$test[[i]]$idx
    train_idx <- res$train[[i]]$idx
    expect_length(intersect(train_idx, test_idx), 0)
    expect_equal(length(unique(block_var[test_idx])), half_blocks)
    expect_equal(length(train_idx) + length(test_idx), nrow(df))
  }
})


test_that("crossv_twofold validates inputs", {
  block_var <- rep(1, 10)
  expect_error(crossv_twofold(df, y, block_var, nreps = 2), "at least two unique blocks")
  expect_error(crossv_twofold(df, y, rep(1:2, each=5), nreps = 1), "at least 2")
})

