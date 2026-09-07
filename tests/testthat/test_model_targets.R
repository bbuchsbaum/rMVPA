library(testthat)

test_that("model_targets on a classification design mirrors cv_labels and leaves y_train alone", {
  df <- data.frame(cond = rep(c("a", "b"), 10), block = rep(1:4, each = 5))
  des <- mvpa_design(df, y_train = ~ cond, block_var = ~ block)

  mt <- model_targets(des)
  expect_s3_class(mt, "model_targets")
  expect_equal(mt$type, "categorical")
  expect_equal(mt$partition, "train")
  expect_equal(mt$values, des$cv_labels)
  expect_equal(mt$observation_ids, 1:20)
  expect_equal(mt$response_ids, c("a", "b"))
  expect_null(mt$response_groups)
  expect_null(mt$row_weights)

  # backward compatibility: y_train still returns cv_labels
  expect_equal(y_train(des), des$cv_labels)
  expect_null(model_targets(des, "test"))
})

test_that("model_targets returns matrix-valued training targets row-aligned to the design", {
  set.seed(1)
  n <- 24
  feats <- matrix(rnorm(n * 3), n, 3, dimnames = list(NULL, c("f1", "f2", "f3")))
  des <- mvpa_design(data.frame(id = seq_len(n)), cv_labels = seq_len(n), targets = feats,
                     block_var = rep(1:4, each = 6))

  mt <- model_targets(des)
  expect_equal(mt$type, "matrix")
  expect_equal(dim(mt$values), c(n, 3))
  expect_equal(mt$values, feats)
  expect_equal(mt$response_ids, c("f1", "f2", "f3"))
  expect_equal(mt$observation_ids, seq_len(n))

  # cv_labels untouched by matrix targets
  expect_equal(y_train(des), seq_len(n))
  expect_equal(nobs(des), n)
})

test_that("mvpa_design rejects misaligned targets and targets_test", {
  df <- data.frame(id = 1:10)
  expect_error(mvpa_design(df, cv_labels = 1:10, targets = matrix(0, 9, 2)), "9 rows")
  expect_error(mvpa_design(df, cv_labels = 1:10, targets = rnorm(11)), "11 rows")
  expect_error(mvpa_design(df, cv_labels = 1:10, targets_test = matrix(0, 5, 2)),
               "requires a 'test_design'")
  expect_error(
    mvpa_design(df, test_design = data.frame(id = 1:5), cv_labels = 1:10,
                y_test = 1:5, targets_test = matrix(0, 4, 2)),
    "4 rows"
  )
})

test_that("model_targets test partition uses targets_test and falls back to y_test", {
  train_df <- data.frame(cond = rep(c("a", "b"), 6))
  test_df <- data.frame(cond = rep(c("a", "b"), 3))
  tt <- matrix(seq_len(12), 6, 2, dimnames = list(NULL, c("u", "v")))

  des <- mvpa_design(train_df, test_df, y_train = ~ cond, y_test = ~ cond, targets_test = tt)
  mt_test <- model_targets(des, "test")
  expect_equal(mt_test$partition, "test")
  expect_equal(mt_test$values, tt)
  expect_equal(mt_test$observation_ids, 1:6)

  # y_test semantics preserved
  expect_equal(y_test(des), rep(c("a", "b"), 3))

  des2 <- mvpa_design(train_df, test_df, y_train = ~ cond, y_test = ~ cond)
  mt2 <- model_targets(des2, "test")
  expect_equal(mt2$type, "categorical")
  expect_equal(mt2$values, factor(rep(c("a", "b"), 3)))
})

test_that("model_targets on a feature_sets_design exposes sets and row weights", {
  set.seed(2)
  X <- matrix(rnorm(20 * 8), 20, 8)
  w <- runif(20)
  fs <- feature_sets(X, blocks(low = 3, sem = 5), row_weights = w)
  fs_test <- feature_sets(matrix(rnorm(10 * 8), 10, 8), blocks(low = 3, sem = 5))
  des <- feature_sets_design(fs, X_test = fs_test, block_var_train = rep(1:4, each = 5))

  mt <- model_targets(des)
  expect_equal(mt$type, "matrix")
  expect_equal(dim(mt$values), c(20, 8))
  expect_equal(unname(mt$values), unname(X))
  expect_equal(as.character(mt$response_groups), rep(c("low", "sem"), c(3, 5)))
  expect_equal(mt$row_weights, w)
  expect_equal(length(mt$response_ids), 8)

  mt_test <- model_targets(des, "test")
  expect_equal(dim(mt_test$values), c(10, 8))
  expect_equal(mt_test$partition, "test")

  # the CV side still sees integer bookkeeping labels
  expect_equal(y_train(des), 1:20)
})

test_that("model_targets on a feature_rsa_design returns the feature matrix", {
  set.seed(3)
  F <- matrix(rnorm(15 * 4), 15, 4)
  des <- feature_rsa_design(F = F, labels = paste0("s", 1:15))
  mt <- model_targets(des)
  expect_equal(mt$type, "matrix")
  expect_equal(unname(mt$values), F)
  expect_equal(mt$response_ids, paste0("F", 1:4))
  expect_null(model_targets(des, "test"))
})

test_that("model_targets rejects non-numeric matrix targets and types logical targets", {
  df <- data.frame(id = 1:10)
  bad <- data.frame(a = rnorm(10), b = rep(c("x", "y"), 5))
  des <- mvpa_design(df, cv_labels = 1:10, targets = bad)
  expect_error(model_targets(des), "numeric columns only")

  chr_mat <- matrix(letters[1:20], 10, 2)
  des2 <- mvpa_design(df, cv_labels = 1:10, targets = chr_mat)
  expect_error(model_targets(des2), "must be numeric")

  num_df <- data.frame(a = rnorm(10), b = rnorm(10))
  des3 <- mvpa_design(df, cv_labels = 1:10, targets = num_df)
  mt3 <- model_targets(des3)
  expect_equal(mt3$type, "matrix")
  expect_equal(dim(mt3$values), c(10, 2))
  expect_equal(mt3$response_ids, c("a", "b"))

  des4 <- mvpa_design(df, cv_labels = 1:10, targets = rep(c(TRUE, FALSE), 5))
  mt4 <- model_targets(des4)
  expect_equal(mt4$type, "categorical")
  expect_true(is.factor(mt4$values))
  expect_equal(mt4$response_ids, c("FALSE", "TRUE"))
})

test_that("mvpa_design names the offending argument and rejects formula targets", {
  df <- data.frame(cond = rep(c("a", "b"), 5))
  expect_error(mvpa_design(df, y_train = rep(c("a", "b"), 6)), "'y_train' has 12")
  expect_error(mvpa_design(df, cv_labels = 1:10, targets = ~ cond), "do not accept formulas")
})

test_that("model_targets prints", {
  des <- mvpa_design(data.frame(cond = rep(c("a", "b"), 5)), y_train = ~ cond)
  expect_output(print(model_targets(des)), "categorical targets, 10 observations x 2 responses")
})
