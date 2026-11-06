# Tests for R/permutation_regional.R: run_permutation_regional() and its scoring.

regional_perm_fixture <- function(signal = 1.5, seed = 11, model = "corclass") {
  set.seed(seed)
  dset <- gen_sample_dataset(c(6, 6, 6), nobs = 48, blocks = 4, nlevels = 2)
  y <- dset$design$y_train
  arr <- neuroim2::as.array(dset$dataset$train_data)
  # Region 1 (x = 1:2) carries the label as a spatial pattern (correlation
  # classifiers ignore a uniform shift); regions 2 and 3 are noise.
  sgn <- ifelse(y == levels(y)[1], signal, -signal)
  pat <- array(sample(c(-1, 1), 2 * 36, replace = TRUE), c(2, 6, 6))
  for (i in seq_along(y)) arr[1:2, , , i] <- arr[1:2, , , i] + sgn[i] * pat
  dset$dataset <- mvpa_dataset(neuroim2::NeuroVec(arr, neuroim2::space(dset$dataset$train_data)), mask = dset$dataset$mask)
  lab <- array(0L, c(6, 6, 6)); lab[1:2, , ] <- 1L; lab[3:4, , ] <- 2L; lab[5:6, , ] <- 3L
  regions <- neuroim2::NeuroVol(lab, neuroim2::space(dset$dataset$mask))
  cv <- blocked_cross_validation(dset$design$block_var)
  ms <- mvpa_model(load_model(model), dset$dataset, dset$design, "classification", crossval = cv)
  list(ms = ms, regions = regions)
}

test_that(".score_regional_null matches a brute-force count", {
  set.seed(1)
  B <- 50; R <- 4
  null <- matrix(stats::rnorm(B * R), B, R)
  null[3, 2] <- NA
  obs <- c(3, 0, -1, NA)
  sc <- rMVPA:::.score_regional_null(obs, null)
  n <- colSums(!is.na(null))
  for (r in 1:3) expect_equal(sc$p[r], (1 + sum(null[, r] >= obs[r], na.rm = TRUE)) / (1 + n[r]))
  expect_true(is.na(sc$p[4]))
  expect_equal(sc$n, as.integer(n))
  # minP by brute force
  null_p <- sapply(1:R, function(r) sapply(1:B, function(b) {
    v <- null[b, r]; if (is.na(v)) return(NA)
    pool <- c(null[, r], obs[r]); (sum(pool >= v, na.rm = TRUE)) / (n[r] + 1)
  }))
  mins <- apply(null_p, 1, min, na.rm = TRUE)
  for (r in 1:3) expect_equal(sc$p_fwer[r], (1 + sum(mins <= sc$p[r])) / (1 + B))
  expect_true(all(sc$p_fwer[1:3] >= sc$p[1:3]))
})

test_that("run_permutation_regional finds the signal region and not the noise", {
  skip_on_cran()
  fx <- regional_perm_fixture()
  pc <- permutation_control(n_perm = 99, seed = 7)
  res <- suppressMessages(run_permutation_regional(fx$ms, fx$regions, perm_ctrl = pc, metric = "Accuracy"))
  expect_s3_class(res, "regional_permutation_result")
  t <- res$table
  expect_equal(t$roinum, 1:3)
  expect_equal(res$p_floor, 1 / 100)
  expect_true(all(t$p >= 1 / 100 & t$p <= 1))
  expect_equal(t$p[1], 1 / 100)                     # nothing in the null beats the real signal
  expect_true(all(t$p[2:3] > 0.05))
  expect_true(all(t$p_adj >= t$p) && all(t$p_fwer >= t$p))
  expect_lt(t$p_fwer[1], 0.05)
  expect_equal(t$n_perm, rep(99L, 3))
  expect_s4_class(res$vol_results$Accuracy$p, "NeuroVol")
  expect_output(print(res), "Regional permutation result")
})

test_that("permutations are reproducible and independent of batching", {
  skip_on_cran()
  fx <- regional_perm_fixture(signal = 0.4)
  pc <- permutation_control(n_perm = 12, seed = 3)
  obs <- suppressMessages(run_regional(fx$ms, fx$regions))
  a <- suppressMessages(run_permutation_regional(fx$ms, fx$regions, observed = obs, perm_ctrl = pc, keep_null = TRUE))
  b <- suppressMessages(run_permutation_regional(fx$ms, fx$regions, observed = obs, perm_ctrl = pc, keep_null = TRUE,
                                                 batch_size = 1))
  expect_identical(a$table, b$table)
  expect_identical(a$null, b$null)
  expect_setequal(a$metrics, c("Accuracy", "AUC"))
  # The null actually varies: relabelling reached the fitted models.
  expect_gt(stats::sd(a$null[, "2", "Accuracy"]), 0)
})

test_that("run_permutation_regional checks its inputs", {
  fx <- regional_perm_fixture()
  expect_error(run_permutation_regional(fx$ms, fx$regions, perm_ctrl = list()), "permutation_control")
  obs <- suppressMessages(run_regional(fx$ms, fx$regions))
  expect_error(run_permutation_regional(fx$ms, fx$regions, observed = obs, metric = "nope",
                                        perm_ctrl = permutation_control(n_perm = 2)), "not found")
})

test_that("both engines give the same null as refitting with run_regional", {
  skip_on_cran()
  skip_if_not_installed("sda")
  for (model in c("corclass", "sda_notune")) {
  fx <- regional_perm_fixture(signal = 0.4, model = model)
  pc <- permutation_control(n_perm = 4, seed = 21)
  res <- suppressMessages(run_permutation_regional(fx$ms, fx$regions, perm_ctrl = pc, keep_null = TRUE))
  for (i in 1:4) {
    spec <- fx$ms   # corclass mvpa_model: the permuted spec is the spec with the permuted design
    spec$design <- permute_labels(fx$ms$design, method = pc$shuffle, seed = pc$seed + i)
    naive <- as.data.frame(suppressMessages(run_regional(spec, fx$regions))$performance_table)
    expect_equal(unname(res$null[i, , "Accuracy"]), naive$Accuracy)
    expect_equal(unname(res$null[i, , "AUC"]), naive$AUC)
  }
  expect_equal(res$engine, if (model == "corclass") "aggregate_fast" else "per_region")
  }
})
