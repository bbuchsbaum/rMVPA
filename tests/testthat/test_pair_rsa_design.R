test_that("pair_rsa_design within-mode reproduces rsa_design pair vectors", {
  set.seed(20260425)
  n <- 8L
  R1 <- as.matrix(stats::dist(matrix(rnorm(n * 5), n, 5)))
  R2 <- as.matrix(stats::dist(matrix(rnorm(n * 5), n, 5)))
  rownames(R1) <- colnames(R1) <- paste0("it", seq_len(n))
  rownames(R2) <- colnames(R2) <- paste0("it", seq_len(n))

  classic <- rsa_design(~ R1 + R2, list(R1 = R1, R2 = R2))
  paired  <- pair_rsa_design(items_a = paste0("it", seq_len(n)),
                             model = list(R1 = R1, R2 = R2))

  expect_equal(paired$model_mat$R1, classic$model_mat$R1, tolerance = 1e-12)
  expect_equal(paired$model_mat$R2, classic$model_mat$R2, tolerance = 1e-12)
  expect_equal(nrow(paired$pair_index), n * (n - 1L) / 2L)
  expect_s3_class(paired, "pair_rsa_design")
  expect_s3_class(paired, "rsa_design")
})

test_that("pair_rsa_design within-mode block masking matches rsa_design include", {
  set.seed(2026)
  n <- 10L
  R1 <- as.matrix(stats::dist(matrix(rnorm(n * 4), n, 4)))
  blocks <- rep(1:2, length.out = n)

  classic <- rsa_design(~ R1, list(R1 = R1, blocks = blocks),
                        block_var = ~ blocks)
  paired  <- pair_rsa_design(items_a = paste0("it", seq_len(n)),
                             model = list(R1 = R1),
                             block_var_a = blocks)

  expect_equal(sum(paired$include), sum(classic$include))
  expect_equal(length(paired$model_mat$R1), sum(classic$include))
  expect_equal(paired$model_mat$R1, classic$model_mat$R1, tolerance = 1e-12)
})

test_that("pair_rsa_design between-mode produces rectangular pair vectors", {
  set.seed(99)
  n_a <- 5L; n_b <- 7L
  M <- matrix(rnorm(n_a * n_b), n_a, n_b)
  des <- pair_rsa_design(
    items_a   = paste0("a", seq_len(n_a)),
    items_b   = paste0("b", seq_len(n_b)),
    model     = list(M = M),
    pairs     = "between",
    row_idx_a = seq_len(n_a),
    row_idx_b = n_a + seq_len(n_b)
  )

  expect_equal(des$pair_kind, "between")
  expect_equal(nrow(des$pair_index), n_a * n_b)
  expect_equal(des$model_mat$M, as.vector(M), tolerance = 1e-12)
  expect_identical(des$row_idx_a, seq_len(n_a))
  expect_identical(des$row_idx_b, n_a + seq_len(n_b))
})

test_that("pair_rsa_design accepts function-valued model entries", {
  items <- c("apple", "pear", "plum", "fig")
  semantic <- function(a, b) as.numeric(a == b)
  des <- pair_rsa_design(items_a = items, model = list(same = semantic))
  expect_equal(des$model_mat$same, rep(0, length(des$model_mat$same)))
  expect_equal(length(des$model_mat$same), length(items) * (length(items) - 1L) / 2L)
})

test_that("pair_rsa_design passes feature rows to function-valued entries", {
  items <- c("apple", "pear", "plum", "fig")
  features <- data.frame(length = c(5, 4, 4, 3), vowel = c(2, 2, 1, 1))
  fdiff <- function(a, b, fa, fb) abs(fa$length - fb$length) + abs(fa$vowel - fb$vowel)

  des <- pair_rsa_design(items_a = items,
                         features_a = features,
                         model = list(feature_distance = fdiff))

  idx <- des$pair_index
  expected <- abs(features$length[idx$i] - features$length[idx$j]) +
              abs(features$vowel[idx$i] - features$vowel[idx$j])
  expect_equal(des$model_mat$feature_distance, expected)
})

test_that("pair_rsa_design feature functions work for rectangular between-domain pairs", {
  items_a <- paste0("movie1_", 1:3)
  items_b <- paste0("movie2_", 1:2)
  fa <- data.frame(semantic = c(0.1, 0.4, 0.9), visual = c(2, 1, 3))
  fb <- data.frame(semantic = c(0.2, 0.8), visual = c(1, 4))
  pair_fun <- function(a, b, xa, xb) {
    abs(xa$semantic - xb$semantic) + 0.25 * abs(xa$visual - xb$visual)
  }

  des <- pair_rsa_design(
    items_a = items_a,
    items_b = items_b,
    features_a = fa,
    features_b = fb,
    model = list(feature_pair = pair_fun),
    pairs = "between",
    row_idx_a = 1:3,
    row_idx_b = 4:5
  )

  expected <- as.vector(outer(seq_len(nrow(fa)), seq_len(nrow(fb)),
                              Vectorize(function(i, j) {
                                abs(fa$semantic[i] - fb$semantic[j]) +
                                  0.25 * abs(fa$visual[i] - fb$visual[j])
                              })))
  expect_equal(des$model_mat$feature_pair, expected)
  expect_equal(nrow(des$pair_index), length(items_a) * length(items_b))
})

test_that("rsa_model end-to-end reproduces classical RSA when fed pair_rsa_design", {
  skip_if_not_installed("neuroim2")
  set.seed(7)
  n <- 12L
  v <- 30L
  X <- matrix(rnorm(n * v), n, v)
  R <- as.matrix(stats::dist(matrix(rnorm(n * 4), n, 4)))

  classic_des <- rsa_design(~ R, list(R = R))
  paired_des  <- pair_rsa_design(items_a = paste0("it", seq_len(n)),
                                 model = list(R = R))

  arr  <- array(rnorm(prod(c(2, 2, 2, n))), c(2, 2, 2, n))
  sp   <- neuroim2::NeuroSpace(c(2, 2, 2, n))
  vec  <- neuroim2::NeuroVec(arr, sp)
  mask <- neuroim2::LogicalNeuroVol(array(1, c(2, 2, 2)),
                                    neuroim2::NeuroSpace(c(2, 2, 2)))
  ds   <- mvpa_dataset(train_data = vec, mask = mask)

  m_classic <- rsa_model(ds, classic_des, distmethod = "pearson", regtype = "pearson")
  m_paired  <- rsa_model(ds, paired_des,  distmethod = "pearson", regtype = "pearson")

  out_classic <- train_model(m_classic, X, y = NULL, indices = NULL)
  out_paired  <- train_model(m_paired,  X, y = NULL, indices = NULL)

  expect_equal(unname(out_classic), unname(out_paired), tolerance = 1e-10)
})

test_that("rsa_model with within pair_rsa_design can target row-indexed subsets", {
  skip_if_not_installed("neuroim2")
  set.seed(77)
  n <- 6L
  row_idx <- 2:5
  v <- 18L
  X <- matrix(rnorm(n * v), n, v)
  R <- as.matrix(stats::dist(matrix(rnorm(length(row_idx) * 3), length(row_idx), 3)))

  arr  <- array(rnorm(prod(c(2, 2, 2, n))), c(2, 2, 2, n))
  sp   <- neuroim2::NeuroSpace(c(2, 2, 2, n))
  vec  <- neuroim2::NeuroVec(arr, sp)
  mask <- neuroim2::LogicalNeuroVol(array(1, c(2, 2, 2)),
                                    neuroim2::NeuroSpace(c(2, 2, 2)))
  ds   <- mvpa_dataset(train_data = vec, mask = mask)

  des <- pair_rsa_design(items_a = paste0("it", seq_along(row_idx)),
                         model = list(R = R),
                         row_idx_a = row_idx)
  m <- rsa_model(ds, des, distmethod = "pearson", regtype = "pearson")
  out <- train_model(m, X, y = NULL, indices = NULL)

  dsub <- 1 - stats::cor(t(X[row_idx, , drop = FALSE]), method = "pearson")
  expect_equal(out[["R"]], stats::cor(dsub[lower.tri(dsub)], R[lower.tri(R)]),
               tolerance = 1e-10)
})

test_that("rsa_model with pair_rsa_design between-mode runs and respects rectangular geometry", {
  skip_if_not_installed("neuroim2")
  set.seed(101)
  n_a <- 4L; n_b <- 6L
  n   <- n_a + n_b
  v   <- 20L
  X <- matrix(rnorm(n * v), n, v)

  M <- matrix(rnorm(n_a * n_b), n_a, n_b)
  des <- pair_rsa_design(
    items_a   = paste0("a", seq_len(n_a)),
    items_b   = paste0("b", seq_len(n_b)),
    model     = list(M = M),
    pairs     = "between",
    row_idx_a = seq_len(n_a),
    row_idx_b = n_a + seq_len(n_b)
  )

  arr  <- array(rnorm(prod(c(2, 2, 2, n))), c(2, 2, 2, n))
  sp   <- neuroim2::NeuroSpace(c(2, 2, 2, n))
  vec  <- neuroim2::NeuroVec(arr, sp)
  mask <- neuroim2::LogicalNeuroVol(array(1, c(2, 2, 2)),
                                    neuroim2::NeuroSpace(c(2, 2, 2)))
  ds   <- mvpa_dataset(train_data = vec, mask = mask)

  m <- rsa_model(ds, des, distmethod = "pearson", regtype = "pearson")
  out <- train_model(m, X, y = NULL, indices = NULL)
  expect_named(out, "M")
  expect_true(is.finite(out[["M"]]))

  # Verify the engine actually computed a rectangular block: cor of A rows vs B rows
  M_neural <- 1 - stats::cor(t(X[seq_len(n_a), ]), t(X[n_a + seq_len(n_b), ]),
                              method = "pearson")
  expect_equal(out[["M"]],
               stats::cor(as.vector(M_neural), as.vector(M)),
               tolerance = 1e-10)
})

test_that("pair_rsa_design rejects malformed inputs", {
  expect_error(pair_rsa_design(items_a = "x", model = list()),
               "at least two items")
  expect_error(pair_rsa_design(items_a = c("a", "b", "c"),
                               pairs = "between",
                               model = list(M = matrix(1, 3, 3))),
               "items_b")
  expect_error(pair_rsa_design(items_a = c("a", "b", "c"),
                               items_b = c("x", "y"),
                               pairs = "between",
                               model = list(M = matrix(1, 3, 3))),
               "row_idx_a")
  expect_error(pair_rsa_design(items_a = c("a", "b"),
                               items_b = c("x", "y"),
                               pairs = "between",
                               row_idx_a = 1:2,
                               row_idx_b = 3:4,
                               model = list(D = stats::dist(matrix(1:4, 2, 2)))),
               "cannot represent rectangular")
  expect_error(pair_rsa_design(items_a = c("a", "b"),
                               features_a = data.frame(x = 1:3),
                               model = list(x = function(a, b, xa, xb) xa$x)),
               "one row per")
})

test_that("rsa_model rejects within pair designs that do not match dataset rows", {
  skip_if_not_installed("neuroim2")
  set.seed(19)
  n <- 6L
  R <- as.matrix(stats::dist(matrix(rnorm(4 * 2), 4, 2)))
  arr  <- array(rnorm(prod(c(2, 2, 2, n))), c(2, 2, 2, n))
  sp   <- neuroim2::NeuroSpace(c(2, 2, 2, n))
  vec  <- neuroim2::NeuroVec(arr, sp)
  mask <- neuroim2::LogicalNeuroVol(array(1, c(2, 2, 2)),
                                    neuroim2::NeuroSpace(c(2, 2, 2)))
  ds   <- mvpa_dataset(train_data = vec, mask = mask)
  des <- pair_rsa_design(items_a = paste0("it", 1:4), model = list(R = R))

  expect_error(rsa_model(ds, des), "Use `row_idx_a`")
})
