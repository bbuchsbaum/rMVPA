test_that("randomized searchlight respects the provided mask", {
  set.seed(42)

  # Small synthetic 3D dataset with an external test set
  D <- c(5, 5, 5)
  nobs <- 6
  train <- neuroim2::NeuroVec(array(rnorm(prod(D) * nobs), c(D, nobs)),
                              neuroim2::NeuroSpace(c(D, nobs)))
  test  <- neuroim2::NeuroVec(array(rnorm(prod(D) * nobs), c(D, nobs)),
                              neuroim2::NeuroSpace(c(D, nobs)))

  # Two masks: full and a smaller version with some active voxels zeroed out
  full_mask_vals <- sample(c(0, 1), prod(D), replace = TRUE, prob = c(0.2, 0.8))
  full_mask  <- neuroim2::NeuroVol(array(full_mask_vals, D), neuroim2::NeuroSpace(D))

  small_mask_vals <- full_mask_vals
  nz <- which(small_mask_vals != 0)
  to_zero <- head(nz, min(15, length(nz)))  # deterministically drop some active voxels
  small_mask_vals[to_zero] <- 0
  small_mask <- neuroim2::NeuroVol(array(small_mask_vals, D), neuroim2::NeuroSpace(D))

  # Simple design with external test labels
  train_des <- data.frame(
    Y = factor(rep(letters[1:3], length.out = nobs)),
    block_var = rep(1, nobs)
  )
  test_des <- data.frame(
    Ytest = factor(rep(letters[1:3], length.out = nobs))
  )
  design <- mvpa_design(train_des, test_design = test_des,
                        y_train = ~Y, y_test = ~Ytest, block_var = "block_var")

  make_model <- function(mask) {
    remap_rrr_model(
      mvpa_dataset(train, test, mask),
      design,
      rank = 0,                    # fast path; avoids rrpack dependency
      leave_one_key_out = FALSE    # keep runtime low
    )
  }

  ms_full  <- make_model(full_mask)
  ms_small <- make_model(small_mask)

  withr::local_options(list(futile.logger.threshold = "ERROR"))

  res_full  <- suppressWarnings(run_searchlight(ms_full,  radius = 1, method = "randomized", niter = 1))
  res_small <- suppressWarnings(run_searchlight(ms_small, radius = 1, method = "randomized", niter = 1))

  full_active  <- sum(neuroim2::values(full_mask) != 0)
  small_active <- sum(neuroim2::values(small_mask) != 0)

  # Reported active voxels should match mask size and shrink with the smaller mask
  expect_equal(res_full$active_voxels, full_active)
  expect_equal(res_small$active_voxels, small_active)
  expect_lt(res_small$active_voxels, res_full$active_voxels)

  # Output maps should be sized to the corresponding mask space
  expect_equal(length(neuroim2::values(res_full$results[[1]])),
               length(neuroim2::values(full_mask)))
  expect_equal(length(neuroim2::values(res_small$results[[1]])),
               length(neuroim2::values(small_mask)))
})
