test_that("center_patterns subtracts stimulus mean across rows", {
  set.seed(1)
  X <- matrix(rnorm(20), nrow = 5, ncol = 4)
  Xc <- rMVPA:::center_patterns(X, method = "stimulus_mean")
  expect_true(all(abs(colMeans(Xc)) < 1e-12))
})

test_that("center_patterns_train_test uses train mean for test centering", {
  Xtr <- matrix(1, nrow = 3, ncol = 4)
  Xte <- matrix(3, nrow = 2, ncol = 4)
  res <- rMVPA:::center_patterns_train_test(Xtr, Xte, method = "stimulus_mean")
  expect_true(all(abs(colMeans(res$train)) < 1e-12))
  expect_true(all(abs(colMeans(res$test) - 2) < 1e-12))
})

test_that("pairwise_dist honors stimulus-mean centering", {
  set.seed(2)
  X <- matrix(rnorm(30), nrow = 6, ncol = 5)
  Xc <- rMVPA:::center_patterns(X, method = "stimulus_mean")
  dist_obj <- eucdist(center = "stimulus_mean")
  D1 <- pairwise_dist(dist_obj, X)
  D2 <- as.matrix(dist(Xc))
  expect_true(isTRUE(all.equal(D1, D2, tolerance = 1e-12)))
})

test_that("rsa_model pattern_center removes shared stimulus-invariant pattern", {
  set.seed(3)
  n <- 6
  p <- 5
  base <- matrix(rnorm(n * p), nrow = n, ncol = p)
  shared <- matrix(rep(seq_len(p), each = n), nrow = n, ncol = p)
  X1 <- base
  X2 <- base + shared

  d1 <- as.dist(1 - cor(t(X1), method = "pearson"))
  rdes <- rsa_design(~ d1, list(d1 = d1))

  dims <- c(2, 2, 1)
  mask <- neuroim2::NeuroVol(array(1, dim = dims), neuroim2::NeuroSpace(dims))
  arr <- array(rnorm(prod(dims) * n), dim = c(dims, n))
  vec <- neuroim2::NeuroVec(arr, neuroim2::NeuroSpace(c(dims, n)))
  dset <- mvpa_dataset(vec, mask = mask)

  mod_center <- rsa_model(dset, rdes, regtype = "pearson", distmethod = "pearson", pattern_center = "stimulus_mean")
  res_center <- train_model(mod_center, X2, y = NULL, indices = NULL)
  res_base_center <- train_model(mod_center, X1, y = NULL, indices = NULL)

  expect_true(isTRUE(all.equal(res_center, res_base_center, tolerance = 1e-8)))
})
