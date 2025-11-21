test_that("resampled searchlight respects mask size and uses total samples", {
  set.seed(123)

  D <- c(4, 4, 4)
  nobs <- 18  # Increased from 5 to ensure sufficient paired items (3 levels * 6 reps)
  space_tr <- neuroim2::NeuroSpace(c(D, nobs))
  space_te <- neuroim2::NeuroSpace(c(D, nobs))
  train <- neuroim2::NeuroVec(array(rnorm(prod(D) * nobs), c(D, nobs)), space_tr)
  test  <- neuroim2::NeuroVec(array(rnorm(prod(D) * nobs), c(D, nobs)), space_te)

  mask_vals <- sample(c(0, 1), prod(D), replace = TRUE, prob = c(0.3, 0.7))
  mask <- neuroim2::NeuroVol(array(mask_vals, D), neuroim2::NeuroSpace(D))

  train_des <- data.frame(Y = factor(rep(letters[1:3], length.out = nobs)), block_var = rep(1, nobs))
  test_des  <- data.frame(Ytest = factor(rep(letters[1:3], length.out = nobs)))
  design <- mvpa_design(train_des, test_design = test_des, y_train = ~Y, y_test = ~Ytest, block_var = "block_var")

  mspec <- remap_rrr_model(mvpa_dataset(train, test, mask), design, rank = 0, leave_one_key_out = FALSE)

  niter <- 15
  res <- suppressWarnings(run_searchlight(mspec, radius = 1, method = "resampled", niter = niter))

  # Active voxels should equal non-zero mask voxels
  expect_equal(res$active_voxels, sum(neuroim2::values(mask) != 0))

  # Output map length should match full mask size
  expect_equal(length(neuroim2::values(res$results[[1]])), length(neuroim2::values(mask)))
})

test_that("resampled searchlight supports vector radii", {
  set.seed(321)

  D <- c(3, 3, 3)
  nobs <- 12  # Increased from 4 to ensure sufficient paired items (2 levels * 6 reps)
  train <- neuroim2::NeuroVec(array(rnorm(prod(D) * nobs), c(D, nobs)), neuroim2::NeuroSpace(c(D, nobs)))
  test  <- neuroim2::NeuroVec(array(rnorm(prod(D) * nobs), c(D, nobs)), neuroim2::NeuroSpace(c(D, nobs)))

  mask <- neuroim2::NeuroVol(array(1, D), neuroim2::NeuroSpace(D))

  train_des <- data.frame(Y = factor(rep(letters[1:2], length.out = nobs)), block_var = rep(1, nobs))
  test_des  <- data.frame(Ytest = factor(rep(letters[1:2], length.out = nobs)))
  design <- mvpa_design(train_des, test_design = test_des, y_train = ~Y, y_test = ~Ytest, block_var = "block_var")

  mspec <- remap_rrr_model(mvpa_dataset(train, test, mask), design, rank = 0, leave_one_key_out = FALSE)

  radii <- c(1, 2)
  res <- suppressWarnings(run_searchlight(mspec, radius = radii, method = "resampled", niter = 10))

  expect_s3_class(res, "searchlight_result")
  expect_equal(res$active_voxels, sum(neuroim2::values(mask) != 0))
})
