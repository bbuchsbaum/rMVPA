context("vector RSA permutations respect blocks")

test_that("block permutations never move a trial out of its block", {
  set.seed(1)
  block <- sample(rep(c("r1", "r2", "r3", "r4"), times = c(3, 5, 4, 6)))
  n <- length(block)
  for (i in seq_len(500)) {
    perm <- rMVPA:::.vector_rsa_block_perm(block, n)
    expect_identical(sort(perm), seq_len(n))
    expect_identical(block[perm], block)
  }
})

test_that("block permutations reach every within-block arrangement", {
  set.seed(2)
  block <- rep(1:2, each = 3)
  seen <- unique(t(replicate(400, rMVPA:::.vector_rsa_block_perm(block, 6L))))
  # 3! arrangements in each of two blocks
  expect_equal(nrow(seen), 36L)
})

test_that("without block structure trials are permuted freely", {
  set.seed(3)
  free <- t(replicate(200, rMVPA:::.vector_rsa_block_perm(NULL, 6L)))
  expect_true(any(free[, 1] > 3))
  singletons <- t(replicate(200, rMVPA:::.vector_rsa_block_perm(1:6, 6L)))
  expect_true(any(singletons[, 1] != 1L))
})

.vrsa_block_fixture <- function(n_blocks = 4L, per_block = 3L, seed = 9L) {
  set.seed(seed)
  n <- n_blocks * per_block
  block <- rep(seq_len(n_blocks), each = per_block)
  # Every trial in a block has the same pattern, so any within-block
  # permutation leaves the ROI distance matrix unchanged.
  base <- matrix(stats::rnorm(n_blocks * 27L), n_blocks)
  X <- base[block, , drop = FALSE]
  labels <- paste0("L", seq_len(n))
  D <- as.matrix(stats::dist(matrix(stats::rnorm(n * 3L), n)))
  rownames(D) <- colnames(D) <- labels
  vec <- neuroim2::NeuroVec(
    array(as.numeric(t(X)), dim = c(3L, 3L, 3L, n)),
    neuroim2::NeuroSpace(c(3L, 3L, 3L, n))
  )
  mask <- neuroim2::NeuroVol(array(1, dim = c(3L, 3L, 3L)),
                             neuroim2::NeuroSpace(c(3L, 3L, 3L)))
  list(X = X, block = block,
       dataset = mvpa_dataset(vec, mask = mask),
       design = vector_rsa_design(D, factor(labels), block))
}

test_that("the permutation null keeps each trial in its block", {
  fx <- .vrsa_block_fixture()
  vm <- vector_rsa_model(fx$dataset, fx$design, distfun = euclidean(),
                         rsa_simfun = "pearson", nperm = 30L,
                         save_distributions = TRUE)
  scores <- rMVPA:::compute_trial_scores(vm, fx$X)
  expect_true(is.finite(mean(scores)))
  set.seed(4)
  perm <- evaluate_model.vector_rsa_model(
    vm, predicted = NULL, observed = scores, roi_data_for_perm = fx$X,
    nperm = 30L, save_distributions = TRUE
  )$permutation_results
  # Within-block shuffles cannot change the data, so every null draw equals
  # the observed score. A cross-block shuffle would not.
  expect_equal(unname(perm$permutation_distributions$rsa_score),
               rep(mean(scores), 30L), tolerance = 1e-12)
  expect_equal(unname(perm$p_values[["rsa_score"]]), 1)
})

test_that("regional vector RSA reports the within-block p-value", {
  fx <- .vrsa_block_fixture(seed = 10L)
  vm <- vector_rsa_model(fx$dataset, fx$design, distfun = euclidean(),
                         rsa_simfun = "pearson", nperm = 20L)
  region <- neuroim2::NeuroVol(array(1L, dim = c(3L, 3L, 3L)),
                               neuroim2::NeuroSpace(c(3L, 3L, 3L)))
  set.seed(5)
  res <- run_regional(vm, region)
  perf <- res$performance_table
  expect_equal(perf$p_rsa_score, 1)
})
