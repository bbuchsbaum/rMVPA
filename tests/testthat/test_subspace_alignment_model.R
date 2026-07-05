context("subspace_alignment_model")

test_that("subspace_alignment_model fits a direct ROI and reports diagnostics", {
  set.seed(123)

  toy <- gen_sample_dataset(D = c(3, 3, 3), nobs = 12, nlevels = 2,
                            blocks = 2, external_test = TRUE)
  model <- subspace_alignment_model(toy$dataset, toy$design, d = 50,
                                    return_predictions = TRUE)
  roi_data <- list(
    train_data = matrix(rnorm(12 * 5), nrow = 12),
    test_data = matrix(rnorm(12 * 5), nrow = 12),
    indices = seq_len(5)
  )

  out <- suppressWarnings(fit_roi(model, roi_data, context = list(id = 7L)))

  expect_false(out$error, info = out$error_message)
  expect_true(all(c("Accuracy", "AUC", "d_used", "alignment_frob") %in% names(out$metrics)))
  expect_equal(unname(out$metrics[["d_used"]]), 5)
  expect_true(is.finite(out$metrics[["alignment_frob"]]))
  expect_s3_class(out$result, "classification_result")
  expect_null(output_schema(model))
  expect_named(compute_performance(model, out$result), c("Accuracy", "AUC"))
})

test_that("subspace_alignment_model validates external test and ROI dimensions", {
  toy <- gen_sample_dataset(D = c(3, 3, 3), nobs = 8, nlevels = 2, blocks = 2)
  expect_error(
    subspace_alignment_model(toy$dataset, toy$design),
    "requires an external test set"
  )

  toy_ext <- gen_sample_dataset(D = c(3, 3, 3), nobs = 8, nlevels = 2,
                                blocks = 2, external_test = TRUE)
  model <- subspace_alignment_model(toy_ext$dataset, toy_ext$design, d = 2)
  bad <- fit_roi(
    model,
    roi_data = list(
      train_data = matrix(1, nrow = 1, ncol = 2),
      test_data = matrix(1, nrow = 1, ncol = 2),
      indices = 1:2
    ),
    context = list(id = 1L)
  )

  expect_true(bad$error)
  expect_match(bad$error_message, "insufficient samples or features")
})

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
