test_that("feature_sets supports block, by-set, and list specifications", {
  X <- matrix(seq_len(24), nrow = 4)

  by_block <- feature_sets(X, blocks(low = 2, high = 4),
                           row_weights = c(1, 0.5, 1, 0.25))
  expect_s3_class(by_block, "feature_sets")
  expect_equal(by_block$dims, c(low = 2L, high = 4L))
  expect_equal(by_block$indices$low, 1:2)
  expect_equal(by_block$indices$high, 3:6)
  expect_equal(by_block$row_weights, c(1, 0.5, 1, 0.25))
  expect_output(print(by_block), "feature_sets")

  interleaved <- feature_sets(
    X,
    by_set(c("a", "b", "a", "b", "b", "a"), order = c("b", "a"))
  )
  expect_equal(interleaved$set_order, c("b", "a"))
  expect_equal(interleaved$indices$b, c(2L, 4L, 5L))
  expect_equal(interleaved$indices$a, c(1L, 3L, 6L))

  from_list <- feature_sets(
    list(sem = X[, 1:2], low = X[, 3:6]),
    set_order = c("low", "sem")
  )
  expect_equal(from_list$set_order, c("low", "sem"))
  expect_equal(dim(from_list$X), c(4L, 6L))
  expect_equal(from_list$dims, c(low = 4L, sem = 2L))
})

test_that("feature_sets validates specifications and weights", {
  X <- matrix(1, nrow = 3, ncol = 4)

  expect_error(blocks(), "provide one or more named")
  expect_error(blocks(a = 0), "positive")
  expect_error(by_set(1:3), "character/factor")
  expect_error(by_set(c("a", "b"), order = ""), "order")

  expect_error(feature_sets(X), "spec")
  expect_error(feature_sets(X, blocks(a = 2)), "sizes sum")
  expect_error(feature_sets(X, by_set(c("a", "b"))), "provided 2 set labels")
  expect_error(feature_sets(X, by_set(c("a", "b", "c", "d"), order = c("a", "b"))),
               "not present in 'order'")
  expect_error(feature_sets(X, "bad"), "unknown 'spec'")
  expect_error(feature_sets(X, blocks(a = 4), row_weights = c(1, 2)), "row_weights")
  expect_error(feature_sets(X, blocks(a = 4), row_weights = c(1, 1, -1)), "finite and non-negative")

  expect_error(feature_sets(list()), "at least one")
  expect_error(feature_sets(list(matrix(1, 2, 2))), "named")
  expect_error(feature_sets(list(a = matrix(1, 2, 2)), set_order = "b"), "not present")
  expect_error(feature_sets(list(a = data.frame(x = 1:2))), "numeric matrix")
  expect_error(feature_sets(list(a = matrix(1, 2, 2), b = matrix(1, 3, 2))), "same number of rows")
  expect_error(feature_sets(1:3), "matrix or a named list")
})

test_that("expected_features handles null mass and validation", {
  X <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 3, byrow = TRUE)
  fs <- feature_sets(X, blocks(a = 1, b = 1))

  gamma <- matrix(c(
    0.2, 0.3, 0.5, 0.0,
    0.4, 0.1, 0.1, 0.4
  ), nrow = 2, byrow = TRUE)

  weighted <- expected_features(fs, gamma, drop_null = TRUE, renormalize = FALSE)
  expect_equal(weighted$row_weights, c(0.8, 0.6))
  expect_equal(weighted$X, gamma[, -1, drop = FALSE] %*% fs$X)
  expect_equal(weighted$indices, fs$indices)

  renorm <- expected_features(fs, gamma, drop_null = TRUE, renormalize = TRUE)
  G <- gamma[, -1, drop = FALSE]
  expect_equal(renorm$row_weights, c(1, 1))
  expect_equal(renorm$X, sweep(G, 1, rowSums(G), "/") %*% fs$X)

  expect_error(expected_features(X, gamma), "feature_sets object")
  expect_error(expected_features(fs, matrix(character(), nrow = 1)), "numeric matrix")
  expect_error(expected_features(fs, matrix(1, nrow = 2, ncol = 5), drop_null = FALSE),
               "gamma has 5 columns")
})

test_that("feature_sets_design validates target layouts and builders", {
  X <- matrix(seq_len(24), nrow = 4)
  train <- feature_sets(X, blocks(a = 2, b = 4))
  target <- feature_sets(matrix(seq_len(18), nrow = 3), blocks(a = 2, b = 4),
                         row_weights = c(1, 0.5, 1))

  des <- feature_sets_design(train, target, block_var_test = c(1, 1, 2))
  expect_s3_class(des, "feature_sets_design")
  expect_true(has_test_set(des))
  expect_equal(des$n_test, 3L)
  expect_equal(des$block_var_test, c(1, 1, 2))
  expect_output(print(des), "Recall blocks")

  builder <- function(X_train, n_test, builder_data) {
    list(
      X = matrix(builder_data$value, nrow = n_test, ncol = ncol(X_train$X)),
      row_weights = seq_len(n_test)
    )
  }
  builder_des <- feature_sets_design(
    train,
    n_test = 3,
    block_var_test = c(1, 2, 2),
    target_builder = builder,
    target_builder_data = list(value = 2)
  )
  expect_output(print(builder_des), "built per fold")
  built <- rMVPA:::.feature_sets_build_target_fold(builder_des, train_idx = 1:2, test_idx = 3, fold_id = 1)
  expect_s3_class(built, "feature_sets")
  expect_equal(dim(built$X), c(3L, 6L))
  expect_equal(built$row_weights, 1:3)

  folds <- list(list(train = 1:2, test = 3), list(train = 2:3, test = 1))
  precomputed <- rMVPA:::.feature_sets_precompute_fold_targets(builder_des, folds)
  expect_equal(length(precomputed), 2)
  expect_null(rMVPA:::.feature_sets_precompute_fold_targets(des, folds))

  expect_error(feature_sets_design(train, target_builder = "bad", n_test = 3), "must be a function")
  expect_error(feature_sets_design(train, target_builder = function(...) train), "supply X_test, n_test")
  expect_error(feature_sets_design(train, n_test = 0), "finite positive scalar")
  expect_error(feature_sets_design(train, block_var_test = 1:2, n_test = 3), "block_var_test")
})

test_that("feature_sets_design materializes and rejects target-builder outputs", {
  train <- feature_sets(matrix(1:24, nrow = 4), blocks(a = 2, b = 4))

  matrix_target <- rMVPA:::.feature_sets_materialize_target(
    train,
    matrix(1, nrow = 3, ncol = 6),
    n_test = 3
  )
  expect_s3_class(matrix_target, "feature_sets")
  expect_equal(matrix_target$row_weights, rep(1, 3))

  list_target <- rMVPA:::.feature_sets_materialize_target(
    train,
    list(X_test = matrix(2, nrow = 3, ncol = 6), row_weights = c(1, 0, 1)),
    n_test = 3
  )
  expect_equal(list_target$row_weights, c(1, 0, 1))

  gamma_target <- rMVPA:::.feature_sets_materialize_target(
    train,
    list(gamma = diag(4)[1:3, , drop = FALSE], renormalize = TRUE),
    n_test = 3
  )
  expect_equal(dim(gamma_target$X), c(3L, 6L))

  expect_error(rMVPA:::.feature_sets_validate_target_fs(train, matrix(1, 3, 6)), "feature_sets object")
  expect_error(
    rMVPA:::.feature_sets_validate_target_fs(train, feature_sets(matrix(1, 3, 5), blocks(a = 2, b = 3))),
    "same number of columns"
  )
  expect_error(
    rMVPA:::.feature_sets_validate_target_fs(train, feature_sets(matrix(1, 3, 6), blocks(b = 4, a = 2))),
    "identical set names/order"
  )
  bad_weights <- feature_sets(matrix(1, 3, 6), blocks(a = 2, b = 4))
  bad_weights$row_weights <- c(1, -1, 1)
  expect_error(rMVPA:::.feature_sets_validate_target_fs(train, bad_weights), "finite and non-negative")
  expect_error(rMVPA:::.feature_sets_materialize_target(train, matrix(1, 2, 5), n_test = 2), "same number of columns")
  expect_error(rMVPA:::.feature_sets_materialize_target(train, list(foo = 1), n_test = 2), "must return")
  expect_error(rMVPA:::.feature_sets_materialize_target(train, 1:3, n_test = 2), "compatible list")
  expect_error(rMVPA:::.feature_sets_build_target_fold(feature_sets_design(train, n_test = 2), 1, 2),
               "does not have a target_builder")
})
