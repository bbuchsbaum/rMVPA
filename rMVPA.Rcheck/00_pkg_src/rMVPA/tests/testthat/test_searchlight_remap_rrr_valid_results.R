test_that("remap_rrr_model searchlight yields non-empty metrics (standard & randomized)", {
  skip_if_not_installed("rrpack")

  set.seed(2025)
  D <- c(5, 5, 5)
  nobs <- 18

  # Volumetric train/test with full mask to avoid tiny ROIs
  train <- neuroim2::NeuroVec(array(rnorm(prod(D) * nobs), c(D, nobs)),
                              neuroim2::NeuroSpace(c(D, nobs)))
  test  <- neuroim2::NeuroVec(array(rnorm(prod(D) * nobs), c(D, nobs)),
                              neuroim2::NeuroSpace(c(D, nobs)))
  mask  <- neuroim2::NeuroVol(array(1, D), neuroim2::NeuroSpace(D))

  train_des <- data.frame(Y = factor(rep(letters[1:3], length.out = nobs)),
                          block_var = rep(1, nobs))
  test_des  <- data.frame(Ytest = factor(rep(letters[1:3], length.out = nobs)))
  design <- mvpa_design(train_des, test_design = test_des,
                        y_train = ~Y, y_test = ~Ytest, block_var = "block_var")

  mspec <- remap_rrr_model(
    mvpa_dataset(train, test, mask),
    design,
    rank = 2,
    leave_one_key_out = FALSE,
    return_adapter = FALSE
  )

  old_opts <- options(futile.logger.threshold = "ERROR")
  on.exit(options(old_opts), add = TRUE)

  # Standard searchlight
  res_std <- suppressWarnings(run_searchlight(mspec, radius = 1, method = "standard"))
  expect_s3_class(res_std, "searchlight_result")
  expect_equal(res_std$active_voxels, sum(neuroim2::values(mask) != 0))
  expect_gt(sum(!is.na(neuroim2::values(res_std$results[[1]]))), 0)

  # Randomized searchlight (2 iterations)
  res_rand <- suppressWarnings(run_searchlight(mspec, radius = 1, method = "randomized", niter = 2))
  expect_s3_class(res_rand, "searchlight_result")
  expect_equal(res_rand$active_voxels, sum(neuroim2::values(mask) != 0))
  expect_gt(sum(!is.na(neuroim2::values(res_rand$results[[1]]))), 0)
})
