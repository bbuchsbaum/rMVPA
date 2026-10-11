context("filter_roi test-set missing-value policy and centre preservation")

# Build a list(train_roi, test_roi) ROI from explicit matrices on a 1-D grid.
.na_policy_roi <- function(train, test = NULL) {
  p <- ncol(train)
  sp <- neuroim2::NeuroSpace(c(p, 1L, 1L))
  coords <- cbind(seq_len(p), 1L, 1L)
  roi <- list(
    train_roi = neuroim2::ROIVec(sp, coords = coords, data = train),
    test_roi = if (is.null(test)) NULL else
      neuroim2::ROIVec(sp, coords = coords, data = test)
  )
  roi
}

.kept_ids <- function(res) as.integer(neuroim2::indices(res$train_roi))

.na_policy_mats <- function(n = 8L, p = 5L, seed = 11L) {
  set.seed(seed)
  list(
    train = matrix(rnorm(n * p), n, p),
    test = matrix(rnorm(n * p), n, p)
  )
}

test_that("a feature with NA in the test data is dropped by default", {
  m <- .na_policy_mats()
  m$test[3L, 2L] <- NA_real_
  res <- rMVPA:::filter_roi(.na_policy_roi(m$train, m$test))
  expect_identical(.kept_ids(res), c(1L, 3L, 4L, 5L))
  expect_false(anyNA(neuroim2::values(res$test_roi)))
  # train and test stay column-aligned
  expect_equal(neuroim2::values(res$test_roi), m$test[, -2L],
               ignore_attr = TRUE)
  expect_equal(neuroim2::values(res$train_roi), m$train[, -2L],
               ignore_attr = TRUE)
})

test_that("test_na = 'keep' retains test NA but still drops infinite test values", {
  m <- .na_policy_mats()
  m$test[3L, 2L] <- NA_real_
  m$test[5L, 4L] <- Inf
  res <- rMVPA:::filter_roi(.na_policy_roi(m$train, m$test), test_na = "keep")
  expect_identical(.kept_ids(res), c(1L, 2L, 3L, 5L))
  expect_true(is.na(neuroim2::values(res$test_roi)[3L, 2L]))

  res_drop <- rMVPA:::filter_roi(.na_policy_roi(m$train, m$test))
  expect_identical(.kept_ids(res_drop), c(1L, 3L, 5L))
})

test_that("an invalid test_na policy is rejected", {
  m <- .na_policy_mats()
  expect_error(
    rMVPA:::filter_roi(.na_policy_roi(m$train, m$test), test_na = "bogus"),
    "should be one of"
  )
})

test_that("a centre with any missing training value is never preserved", {
  m <- .na_policy_mats()
  m$train[2L, 3L] <- NA_real_   # partly-NA centre (was preserved before)
  res <- rMVPA:::filter_roi(.na_policy_roi(m$train, m$test), preserve = 3L)
  expect_false(3L %in% .kept_ids(res))
  expect_false(anyNA(neuroim2::values(res$train_roi)))

  m$train[, 3L] <- NA_real_     # all-NA centre
  res <- rMVPA:::filter_roi(.na_policy_roi(m$train, m$test), preserve = 3L)
  expect_false(3L %in% .kept_ids(res))
})

test_that("a zero-variance centre is still preserved", {
  m <- .na_policy_mats()
  m$train[, 3L] <- 7
  res <- rMVPA:::filter_roi(.na_policy_roi(m$train, m$test), preserve = 3L)
  expect_true(3L %in% .kept_ids(res))
  # and the same column without preserve is removed
  res_np <- rMVPA:::filter_roi(.na_policy_roi(m$train, m$test))
  expect_false(3L %in% .kept_ids(res_np))
  # train-only ROIs follow the same rule
  res_tr <- rMVPA:::filter_roi(.na_policy_roi(m$train), preserve = 3L)
  expect_true(3L %in% .kept_ids(res_tr))
})

test_that("a zero-variance centre with test NA follows the test policy", {
  m <- .na_policy_mats()
  m$train[, 3L] <- 7
  m$test[4L, 3L] <- NA_real_
  roi <- .na_policy_roi(m$train, m$test)
  expect_false(3L %in% .kept_ids(rMVPA:::filter_roi(roi, preserve = 3L)))
  expect_true(3L %in% .kept_ids(
    rMVPA:::filter_roi(roi, preserve = 3L, test_na = "keep")
  ))
  m$test[4L, 3L] <- -Inf
  roi <- .na_policy_roi(m$train, m$test)
  expect_false(3L %in% .kept_ids(
    rMVPA:::filter_roi(roi, preserve = 3L, test_na = "keep")
  ))
})

test_that("multibasis groups apply the same test-NA and centre rules", {
  set.seed(5)
  n <- 8L; V <- 4L; k <- 2L
  train <- matrix(rnorm(n * V * k), n, V * k)
  test <- matrix(rnorm(n * V * k), n, V * k)
  test[2L, V + 2L] <- NA_real_         # basis 2 of voxel 2
  sp <- neuroim2::NeuroSpace(c(V, 1L, 1L))
  coords <- cbind(rep(seq_len(V), k), 1L, 1L)
  roi <- list(
    train_roi = neuroim2::ROIVec(sp, coords = coords, data = train),
    test_roi = neuroim2::ROIVec(sp, coords = coords, data = test),
    basis_count = k
  )
  res <- rMVPA:::filter_roi(roi)
  expect_identical(res$column_positions, c(1L, 3L, 4L, 5L, 7L, 8L))
  res_keep <- rMVPA:::filter_roi(roi, test_na = "keep")
  expect_identical(res_keep, roi)

  train[3L, 2L] <- NA_real_            # basis 1 of voxel 2 in train
  test[2L, V + 2L] <- 0
  roi$train_roi <- neuroim2::ROIVec(sp, coords = coords, data = train)
  roi$test_roi <- neuroim2::ROIVec(sp, coords = coords, data = test)
  res <- rMVPA:::filter_roi(roi, preserve = 2L)
  expect_identical(res$column_positions, c(1L, 3L, 4L, 5L, 7L, 8L))
})

test_that("era_rsa_model declares that it handles missing test values", {
  expect_identical(rMVPA:::.model_test_na_policy(list()), "drop")
  expect_identical(
    rMVPA:::.model_test_na_policy(list(handles_missing_test = TRUE)),
    "keep"
  )
  expect_true(any(grepl("handles_missing_test = TRUE",
                        deparse(body(rMVPA::era_rsa_model)), fixed = TRUE)))
})

test_that("extract_roi applies the requested test policy", {
  ds <- gen_sample_dataset(c(4, 4, 4), 20, blocks = 4, nlevels = 2,
                           external_test = TRUE)
  test_arr <- neuroim2::as.array(ds$dataset$test_data)
  test_arr[1, 1, 1, 3] <- NA_real_
  ds$dataset$test_data <- neuroim2::NeuroVec(test_arr,
                                             neuroim2::space(ds$dataset$test_data))
  vox <- which(as.logical(ds$dataset$mask))[1:10]
  sample <- rMVPA:::data_sample(ds$dataset, vox)
  drop_roi <- rMVPA:::extract_roi(sample, ds$dataset, min_voxels = 2)
  keep_roi <- rMVPA:::extract_roi(sample, ds$dataset, min_voxels = 2,
                                  test_na = "keep")
  expect_false(1L %in% neuroim2::indices(drop_roi$test_roi))
  expect_true(1L %in% neuroim2::indices(keep_roi$test_roi))
  expect_false(anyNA(neuroim2::values(drop_roi$test_roi)))
})

test_that("a regional external-test analysis survives a test-set NA voxel", {
  set.seed(42)
  ds <- gen_sample_dataset(c(5, 5, 5), 40, blocks = 4, nlevels = 2,
                           external_test = TRUE)
  test_arr <- neuroim2::as.array(ds$dataset$test_data)
  test_arr[2, 2, 2, 5] <- NA_real_
  ds$dataset$test_data <- neuroim2::NeuroVec(test_arr,
                                             neuroim2::space(ds$dataset$test_data))
  region_mask <- neuroim2::NeuroVol(
    array(rep(1:2, length.out = 125), c(5, 5, 5)),
    neuroim2::space(ds$dataset$mask)
  )
  cv <- blocked_cross_validation(ds$design$block_var)
  mspec <- mvpa_model(load_model("sda_notune"), ds$dataset, ds$design,
                      "classification", crossval = cv)
  res <- run_regional(mspec, region_mask)
  perf <- res$performance_table
  expect_equal(nrow(perf), 2L)
  expect_true(all(is.finite(perf$Accuracy)))
  pred <- res$prediction_table
  expect_false(anyNA(pred$predicted))
})

test_that("shard extraction applies the same test policy as extract_roi", {
  skip_if_not_installed("shard")
  ds <- gen_sample_dataset(c(4, 4, 4), 20, blocks = 4, nlevels = 2,
                           external_test = TRUE)
  test_arr <- neuroim2::as.array(ds$dataset$test_data)
  test_arr[1, 1, 1, 3] <- NA_real_
  ds$dataset$test_data <- neuroim2::NeuroVec(test_arr,
                                             neuroim2::space(ds$dataset$test_data))
  mspec <- mvpa_model(load_model("sda_notune"), ds$dataset, ds$design,
                      "classification",
                      crossval = blocked_cross_validation(ds$design$block_var))
  mspec_shard <- use_shard(mspec)
  on.exit(rMVPA:::shard_cleanup(mspec_shard$shard_data), add = TRUE)
  vox <- which(as.logical(ds$dataset$mask))[1:10]
  sample <- rMVPA:::data_sample(ds$dataset, vox)
  for (policy in c("drop", "keep")) {
    roi_default <- rMVPA:::extract_roi(sample, ds$dataset, min_voxels = 2,
                                       test_na = policy)
    roi_shard <- rMVPA:::shard_extract_roi(vox, mspec_shard$shard_data,
                                           min_voxels = 2, test_na = policy)
    expect_equal(neuroim2::indices(roi_shard$test_roi),
                 neuroim2::indices(roi_default$test_roi), info = policy)
  }
})
