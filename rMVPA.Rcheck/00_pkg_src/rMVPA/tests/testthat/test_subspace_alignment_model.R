context("subspace_alignment_model")

test_that("subspace_alignment_model runs regional pipeline", {
  skip_on_cran()
  set.seed(123)

  # Small synthetic dataset with external test set
  toy <- gen_sample_dataset(D = c(4,4,4), nobs = 24, nlevels = 3, blocks = 3, external_test = TRUE)

  # Single-region mask so run_regional is lightweight
  regionMask <- neuroim2::NeuroVol(array(1, c(4,4,4)), neuroim2::space(toy$dataset$mask))

  ms <- subspace_alignment_model(toy$dataset, toy$design, d = 10)
  res <- run_regional(ms, regionMask)

  expect_s3_class(res, "regional_mvpa_result")
  expect_true(is.data.frame(res$performance_table))
  expect_true(all(c("Accuracy", "AUC", "d_used", "alignment_frob") %in% names(res$performance_table)))
})


test_that("subspace_alignment caps subspace dimension appropriately", {
  skip_on_cran()
  set.seed(321)

  toy <- gen_sample_dataset(D = c(3,3,3), nobs = 12, nlevels = 2, blocks = 2, external_test = TRUE)
  regionMask <- neuroim2::NeuroVol(array(1, c(3,3,3)), neuroim2::space(toy$dataset$mask))

  # Request a very large d; should cap at min(p, n_train-1, n_test-1)
  ms <- subspace_alignment_model(toy$dataset, toy$design, d = 50)
  res <- run_regional(ms, regionMask)

  p <- prod(dim(toy$dataset$mask))  # number of voxels/features in ROI
  n_train <- nobs(toy$dataset)
  n_test <- dim(toy$dataset$test_data)[length(dim(toy$dataset$test_data))]
  cap_expected <- min(p, n_train - 1L, n_test - 1L)

  expect_equal(unique(res$performance_table$d_used), cap_expected)
})
