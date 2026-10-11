library(testthat)

# Regression tests for two CV correctness bugs:
#  1. crossv_seq_block wrote per-block fold labels back in concatenated
#     (block-by-block) order, so fold j tested the wrong rows whenever blocks
#     were not contiguous in row order.
#  2. kfold_cross_validation stored a fixed fold assignment in block_var, but
#     crossval_samples() redrew random folds on every call.

# ---- helpers -----------------------------------------------------------------

# Canonical string per set, so that a collection of index sets can be compared
# independently of list order.
canon_sets <- function(sets) {
  unname(sort(vapply(sets, function(s) paste(sort(as.integer(s)), collapse = ","), character(1))))
}

# Oracle: the contiguous segments of a block's own row positions, as produced by
# cut() over the block's positions (the intended semantics of crossv_seq_block).
expected_segments <- function(pos, nfolds) {
  segs <- unname(split(pos, cut(pos, nfolds)))
  segs[lengths(segs) > 0]
}

# Checks a sequential-blocked CV result against the oracle, rep by rep.
# Folds are stored rep-major: rep r, fold j is at index (r - 1) * nfolds + j.
check_seq_block_partition <- function(res, block_var, nfolds, nreps) {
  n <- length(block_var)
  tests <- lapply(res$test, function(t) t$idx)
  expect_equal(length(tests), nfolds * nreps)

  for (r in seq_len(nreps)) {
    rep_tests <- tests[((r - 1) * nfolds + 1):(r * nfolds)]

    # every row is tested exactly once per repetition
    expect_equal(sort(as.integer(unlist(rep_tests))), seq_len(n), info = paste("rep", r))

    for (b in sort(unique(block_var))) {
      pos <- which(block_var == b)
      expected <- expected_segments(pos, nfolds)

      # the rows of block b tested in each fold, one entry per fold
      got <- lapply(rep_tests, function(tt) intersect(tt, pos))
      got <- got[lengths(got) > 0]

      # each fold holds exactly one contiguous segment of the block's own rows,
      # and the segments form the oracle partition (labels may be permuted)
      expect_equal(canon_sets(got), canon_sets(expected),
                   info = paste("rep", r, "block", b))
    }

    # when blocks are at least as large as nfolds, every fold has rows from every block
    if (all(table(block_var) >= nfolds)) {
      for (tt in rep_tests) {
        for (b in sort(unique(block_var))) {
          expect_gt(length(intersect(tt, which(block_var == b))), 0)
        }
      }
    }
  }
  invisible(TRUE)
}

# ---- sequential blocked CV: non-contiguous blocks ----------------------------

test_that("seq-block CV with non-contiguous blocks tests each block's own contiguous segments", {
  n <- 12
  X <- data.frame(x = seq_len(n))
  y <- rep(letters[1:3], length.out = n)
  block_var <- rep(1:3, times = 4)   # blocks interleaved in row order

  set.seed(11)
  res <- crossv_seq_block(X, y, nfolds = 2, block_var = block_var, nreps = 3)
  check_seq_block_partition(res, block_var, nfolds = 2, nreps = 3)

  # the ytest column must follow the test rows
  for (i in seq_len(nrow(res))) {
    tt <- res$test[[i]]$idx
    expect_equal(res$ytest[[i]], y[tt])
  }
})

test_that("seq-block CV with non-contiguous blocks, three folds", {
  n <- 18
  X <- data.frame(x = seq_len(n))
  y <- rep(letters[1:3], length.out = n)
  block_var <- rep(1:3, times = 6)

  set.seed(12)
  res <- crossv_seq_block(X, y, nfolds = 3, block_var = block_var, nreps = 2)
  check_seq_block_partition(res, block_var, nfolds = 3, nreps = 2)
})

test_that("seq-block CV with non-contiguous blocks gives train = complement of test", {
  n <- 12
  X <- data.frame(x = seq_len(n))
  y <- rep(letters[1:3], length.out = n)
  block_var <- rep(1:3, times = 4)

  set.seed(13)
  res <- crossv_seq_block(X, y, nfolds = 2, block_var = block_var, nreps = 2)
  for (i in seq_len(nrow(res))) {
    tt <- res$test[[i]]$idx
    tr <- res$train[[i]]$idx
    expect_length(intersect(tt, tr), 0)
    expect_equal(sort(as.integer(c(tt, tr))), seq_len(n))
  }
})

test_that("seq-block CV with contiguous blocks keeps the same partition structure", {
  n <- 12
  X <- data.frame(x = seq_len(n))
  y <- rep(letters[1:3], length.out = n)
  block_var <- rep(1:3, each = 4)   # contiguous blocks

  set.seed(14)
  res <- crossv_seq_block(X, y, nfolds = 2, block_var = block_var, nreps = 3)
  check_seq_block_partition(res, block_var, nfolds = 2, nreps = 3)

  # for contiguous blocks, block 1 (rows 1:4) is split into the consecutive runs 1:2 and 3:4
  tests <- lapply(res$test, function(t) t$idx)
  expect_equal(canon_sets(lapply(tests[1:2], function(tt) intersect(tt, 1:4))),
               canon_sets(list(1:2, 3:4)))
})

# ---- k-fold CV: folds fixed at construction ----------------------------------

kfold_fixture <- function() {
  n <- 20
  set.seed(21)
  list(
    cv = kfold_cross_validation(len = n, nfolds = 4),
    X = data.frame(x1 = rnorm(n), x2 = rnorm(n)),
    y = rep(letters[1:4], length.out = n),
    n = n
  )
}

test_that("kfold crossval_samples returns identical folds across calls with different RNG state", {
  fx <- kfold_fixture()

  set.seed(1)
  a <- crossval_samples(fx$cv, fx$X, fx$y)
  set.seed(2)
  b <- crossval_samples(fx$cv, fx$X, fx$y)

  expect_equal(nrow(a), fx$cv$nfolds)
  expect_equal(nrow(b), fx$cv$nfolds)
  for (j in seq_len(nrow(a))) {
    expect_identical(a$test[[j]]$idx, b$test[[j]]$idx)
    expect_identical(a$train[[j]]$idx, b$train[[j]]$idx)
    expect_identical(a$ytest[[j]], b$ytest[[j]])
    expect_identical(a$ytrain[[j]], b$ytrain[[j]])
  }
  expect_identical(a$.id, b$.id)
})

test_that("kfold test rows of fold j equal which(block_var == j)", {
  fx <- kfold_fixture()
  set.seed(3)
  samp <- crossval_samples(fx$cv, fx$X, fx$y)

  for (j in seq_len(fx$cv$nfolds)) {
    expect_equal(as.integer(samp$test[[j]]$idx), which(fx$cv$block_var == j))
    expect_equal(samp$ytest[[j]], fx$y[which(fx$cv$block_var == j)])
  }
})

test_that("train_indices is consistent with crossval_samples for kfold", {
  fx <- kfold_fixture()
  set.seed(4)
  samp <- crossval_samples(fx$cv, fx$X, fx$y)

  for (j in seq_len(fx$cv$nfolds)) {
    expect_equal(as.integer(train_indices(fx$cv, j)), as.integer(samp$train[[j]]$idx))
    expect_equal(as.integer(partition_indices(fx$cv, j)), as.integer(samp$test[[j]]$idx))
    expect_length(intersect(train_indices(fx$cv, j), partition_indices(fx$cv, j)), 0)
  }
})

test_that("kfold folds partition 1..n", {
  fx <- kfold_fixture()
  set.seed(5)
  samp <- crossval_samples(fx$cv, fx$X, fx$y)

  all_test <- unlist(lapply(samp$test, function(t) t$idx))
  expect_equal(sort(as.integer(all_test)), seq_len(fx$n))
  expect_equal(sort(unlist(lapply(seq_len(fx$cv$nfolds), partition_indices, obj = fx$cv))),
               seq_len(fx$n))
})

test_that("kfold crossval_samples folds are stable across repeated calls on one object", {
  fx <- kfold_fixture()
  folds <- lapply(1:5, function(i) {
    set.seed(100 + i)
    crossval_samples(fx$cv, fx$X, fx$y)$test
  })
  ref <- lapply(folds[[1]], function(t) t$idx)
  for (f in folds[-1]) {
    expect_equal(lapply(f, function(t) t$idx), ref)
  }
})

# ---- end to end: run_regional with kfold CV ----------------------------------

test_that("run_regional with kfold CV is reproducible across differing seeds", {
  dset <- gen_sample_dataset(c(5, 5, 5), nobs = 48, nlevels = 2, data_mode = "image",
                             response_type = "categorical")
  cval <- kfold_cross_validation(len = 48, nfolds = 4)
  region_mask <- NeuroVol(sample(1:4, size = length(dset$dataset$mask), replace = TRUE),
                          space(dset$dataset$mask))
  model <- load_model("corclass")
  mspec <- mvpa_model(model, dset$dataset, dset$design, model_type = "classification",
                      crossval = cval)

  set.seed(1)
  res1 <- run_regional(mspec, region_mask)
  set.seed(2)
  res2 <- run_regional(mspec, region_mask)

  expect_equal(as.data.frame(res1$performance_table),
               as.data.frame(res2$performance_table))
})

test_that("kfold crossval_samples rejects data whose row count differs from the fold assignment", {
  set.seed(7)
  cv <- kfold_cross_validation(len = 20, nfolds = 4)
  y <- factor(rep(c("a", "b"), 10))
  expect_error(crossval_samples(cv, data.frame(x = rnorm(24)), y), "built for 20 observations")
  expect_error(crossval_samples(cv, data.frame(x = rnorm(16)), y), "built for 20 observations")
  expect_silent(crossval_samples(cv, data.frame(x = rnorm(20)), y))
})
