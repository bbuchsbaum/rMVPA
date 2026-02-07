context("Temporal confound RDMs")

test_that("temporal_rdm returns dist and matrix forms correctly", {
  on <- 1:10
  run <- rep(1:2, each=5)
  d1 <- temporal_rdm(on, block=run, kernel="exp", units="TR", TR=1, metric="distance")
  expect_true(inherits(d1, "dist"))
  expect_equal(length(d1), 45)

  m1 <- temporal_rdm(on, block=run, kernel="exp", units="TR", TR=1, metric="distance", as_dist=FALSE)
  expect_true(is.matrix(m1))
  expect_equal(dim(m1), c(10,10))
  expect_true(all(diag(m1) == 0))
  expect_true(isTRUE(all.equal(m1, t(m1))))
})

test_that("within_blocks_only masks cross-run pairs to zero nuisance (distance)", {
  on <- 1:10
  run <- rep(1:2, each=5)
  m <- temporal_rdm(on, block=run, kernel="exp", units="TR", TR=1, metric="distance", as_dist=FALSE)
  # cross-run block should be zeros due to masking
  expect_true(all(m[1:5, 6:10] == 0))
  expect_true(all(m[6:10, 1:5] == 0))
})

test_that("temporal_hrf_overlap structure and masking", {
  set.seed(1)
  on <- c(0, 2, 4, 6,  20, 22, 24, 26) # two runs separated in time
  run <- rep(1:2, each=4)
  d <- temporal_hrf_overlap(on, run=run, TR=2, similarity="overlap", metric="distance")
  expect_true(inherits(d, "dist"))
  m <- as.matrix(d)
  diag(m) <- 0
  # within-run distances should be smaller on average than far pairs
  expect_true(mean(m[1:4,1:4][lower.tri(m[1:4,1:4])]) < mean(m))
  # cross-run masked to zero nuisance
  expect_true(all(m[1:4,5:8] == 0))
})

test_that("temporal_confounds creates a named list of dists", {
  on <- 1:12
  run <- rep(1:3, each=4)
  spec <- list(
    adj = list(kernel="adjacent", width=1),
    exp3 = list(kernel="exp", lambda=3),
    hrf  = list(kind="hrf", TR=1.0)
  )
  lst <- temporal_confounds(spec, on, run=run, units="TR", TR=1)
  expect_true(is.list(lst))
  expect_true(all(c("adj","exp3","hrf") %in% names(lst)))
  expect_true(all(vapply(lst, function(x) inherits(x, "dist"), logical(1))))
})

test_that("temporal_nuisance_for_msreve returns KxK matrix", {
  # Build a tiny mvpa_design with 3 conditions across runs
  train_df <- data.frame(cond = factor(rep(letters[1:3], each=4)), run = rep(1:2, 6))
  mvdes <- mvpa_design(train_df, y_train = ~ cond, block_var = ~ run)
  time_idx <- seq_len(nrow(train_df))
  Kmat <- temporal_nuisance_for_msreve(mvpa_design = mvdes, time_idx = time_idx, kernel = "exp", units = "index", metric = "distance")
  expect_true(is.matrix(Kmat))
  expect_equal(dim(Kmat), c(3,3))
  expect_true(all(diag(Kmat) == 0))
  expect_true(isTRUE(all.equal(Kmat, t(Kmat))))
})

test_that("msreve_temporal_confounds builds multiple KxK nuisances", {
  train_df <- data.frame(cond = factor(rep(letters[1:3], each=4)), run = rep(1:2, 6))
  mvdes <- mvpa_design(train_df, y_train = ~ cond, block_var = ~ run)
  time_idx <- seq_len(nrow(train_df))
  spec <- list(
    expk = list(kernel="exp", lambda=2),
    hrf  = list(kind="hrf", TR=1.0)
  )
  lst <- msreve_temporal_confounds(mvpa_design = mvdes, time_idx = time_idx, spec = spec)
  expect_true(is.list(lst))
  expect_true(all(names(spec) %in% names(lst)))
  expect_true(all(vapply(lst, function(M) is.matrix(M) && all(dim(M) == c(3,3)), logical(1))))
})

