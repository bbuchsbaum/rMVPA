testthat::skip_if_not_installed("neuroim2")
library(neuroim2)

test_that("mvpa_multibasis_dataset accepts basis-wise NeuroVec list", {
  fx <- make_multibasis_fixture()
  dset <- mvpa_multibasis_dataset(train_data = fx$basis_vecs, mask = fx$mask)

  expect_true(inherits(dset, "mvpa_multibasis_image_dataset"))
  expect_equal(dset$basis_count, fx$k)
  expect_equal(nobs(dset), fx$n_events)

  X <- get_feature_matrix(dset)
  expect_true(is.matrix(X))
  expect_equal(nrow(X), fx$n_events)
  expect_equal(ncol(X), sum(fx$mask > 0) * fx$k)
})

test_that("mvpa_multibasis_dataset splits event-major concatenated series", {
  fx <- make_multibasis_fixture(n_events = 5, k = 3)
  concat_vec <- make_concat_vec(fx, ordering = "event_major")
  dset <- mvpa_multibasis_dataset(
    train_data = concat_vec,
    mask = fx$mask,
    basis_count = fx$k,
    ordering = "event_major"
  )

  vox <- which(fx$mask > 0)
  for (b in seq_len(fx$k)) {
    got <- neuroim2::series(dset$train_data[[b]], vox)
    expect <- neuroim2::series(fx$basis_vecs[[b]], vox)
    expect_equal(got, expect)
  }
})

test_that("mvpa_multibasis_dataset splits basis-major concatenated series", {
  fx <- make_multibasis_fixture(n_events = 5, k = 3)
  concat_vec <- make_concat_vec(fx, ordering = "basis_major")
  dset <- mvpa_multibasis_dataset(
    train_data = concat_vec,
    mask = fx$mask,
    basis_count = fx$k,
    ordering = "basis_major"
  )

  vox <- which(fx$mask > 0)
  for (b in seq_len(fx$k)) {
    got <- neuroim2::series(dset$train_data[[b]], vox)
    expect <- neuroim2::series(fx$basis_vecs[[b]], vox)
    expect_equal(got, expect)
  }
})

test_that("mvpa_multibasis_dataset reads one file per basis", {
  fx <- make_multibasis_fixture(n_events = 4, k = 2)
  files <- write_basis_files(fx$basis_vecs)
  on.exit(unlink(files), add = TRUE)

  dset <- mvpa_multibasis_dataset(train_data = files, mask = fx$mask)
  expect_equal(dset$basis_count, fx$k)
  expect_equal(nobs(dset), fx$n_events)

  vox <- which(fx$mask > 0)
  for (b in seq_len(fx$k)) {
    got <- neuroim2::series(dset$train_data[[b]], vox)
    expect <- neuroim2::series(fx$basis_vecs[[b]], vox)
    expect_equal(got, expect)
  }
})

test_that("mvpa_multibasis_dataset reads concatenated file with event-major ordering", {
  fx <- make_multibasis_fixture(n_events = 5, k = 2)
  concat_vec <- make_concat_vec(fx, ordering = "event_major")
  f <- tempfile(fileext = ".nii.gz")
  neuroim2::write_vec(concat_vec, f)
  on.exit(unlink(f), add = TRUE)

  dset <- mvpa_multibasis_dataset(
    train_data = f,
    mask = fx$mask,
    basis_count = fx$k,
    ordering = "event_major"
  )
  expect_equal(dset$basis_count, fx$k)
  expect_equal(nobs(dset), fx$n_events)
})

test_that("mvpa_multibasis_dataset reads concatenated file with basis-major ordering", {
  fx <- make_multibasis_fixture(n_events = 5, k = 2)
  concat_vec <- make_concat_vec(fx, ordering = "basis_major")
  f <- tempfile(fileext = ".nii.gz")
  neuroim2::write_vec(concat_vec, f)
  on.exit(unlink(f), add = TRUE)

  dset <- mvpa_multibasis_dataset(
    train_data = f,
    mask = fx$mask,
    basis_count = fx$k,
    ordering = "basis_major"
  )
  expect_equal(dset$basis_count, fx$k)
  expect_equal(nobs(dset), fx$n_events)
})

test_that("mvpa_multibasis_dataset handles test_data and returns test ROI", {
  fx_train <- make_multibasis_fixture(n_events = 6, k = 2)
  fx_test <- make_multibasis_fixture(n_events = 3, k = 2)
  dset <- mvpa_multibasis_dataset(
    train_data = fx_train$basis_vecs,
    test_data = fx_test$basis_vecs,
    mask = fx_train$mask
  )

  expect_true(has_test_set(dset))
  expect_equal(nobs(dset), fx_train$n_events)

  vox <- which(dset$mask > 0)[1:5]
  roi <- as_roi(data_sample(dset, vox), dset)
  expect_true(inherits(roi$train_roi, "ROIVec"))
  expect_true(inherits(roi$test_roi, "ROIVec"))
  expect_equal(nrow(neuroim2::values(roi$train_roi)), fx_train$n_events)
  expect_equal(nrow(neuroim2::values(roi$test_roi)), fx_test$n_events)
  expect_equal(ncol(neuroim2::values(roi$train_roi)), length(vox) * fx_train$k)
})

test_that("as_roi on multibasis samples concatenates basis features", {
  fx <- make_multibasis_fixture(n_events = 6, k = 2)
  dset <- mvpa_multibasis_dataset(train_data = fx$basis_vecs, mask = fx$mask)

  vox <- which(fx$mask > 0)[1:4]
  sam <- data_sample(dset, vox)
  roi <- as_roi(sam, dset)

  expect_true(inherits(roi$train_roi, "ROIVec"))
  vals <- neuroim2::values(roi$train_roi)
  expect_equal(nrow(vals), fx$n_events)
  expect_equal(ncol(vals), length(vox) * fx$k)

  expected <- do.call(cbind, lapply(fx$basis_vecs, function(v) neuroim2::series(v, vox)))
  expect_equal(vals, expected)
})

test_that("build_output_map aggregates duplicated voxel ids for multibasis features", {
  fx <- make_multibasis_fixture(n_events = 5, k = 3)
  dset <- mvpa_multibasis_dataset(train_data = fx$basis_vecs, mask = fx$mask)

  ids <- rep(which(dset$mask > 0), times = dset$basis_count)
  V <- sum(dset$mask > 0)
  metric_vec <- c(rep(1, V), rep(2, V), rep(4, V))

  out_map <- rMVPA:::build_output_map(dset, metric_vec, ids)
  mask_ids <- which(dset$mask > 0)
  expect_true(all(abs(out_map[mask_ids] - (7 / 3)) < 1e-10))
})

test_that("mvpa_multibasis_dataset validates ambiguous and inconsistent inputs", {
  fx <- make_multibasis_fixture(n_events = 4, k = 3)

  expect_error(
    mvpa_multibasis_dataset(train_data = make_concat_vec(fx, "event_major"), mask = fx$mask),
    "basis_count.*must be supplied"
  )

  expect_error(
    mvpa_multibasis_dataset(
      train_data = make_concat_vec(fx, "event_major"),
      mask = fx$mask,
      basis_count = 5,
      ordering = "event_major"
    ),
    "not divisible by basis_count"
  )

  expect_error(
    mvpa_multibasis_dataset(train_data = fx$basis_vecs, mask = fx$mask, basis_count = 2),
    "resolved 3 basis series but basis_count=2"
  )

  bad <- fx$basis_vecs
  bad[[2]] <- neuroim2::sub_vector(bad[[2]], 1:3)
  expect_error(
    mvpa_multibasis_dataset(train_data = bad, mask = fx$mask),
    "identical dimensions"
  )

  expect_error(
    mvpa_multibasis_dataset(train_data = c("/tmp/definitely_missing_1.nii.gz",
                                           "/tmp/definitely_missing_2.nii.gz"),
                            mask = fx$mask),
    "file not found"
  )
})

test_that("run_global works with multibasis dataset", {
  skip_if_not_installed("sda")
  fx <- make_multibasis_fixture(D = c(4, 4, 4), n_events = 24, k = 2)

  y <- factor(rep(c("a", "b"), length.out = fx$n_events))
  blocks <- rep(1:4, length.out = fx$n_events)
  design <- mvpa_design(data.frame(y = y, block = blocks), y_train = ~y, block_var = ~block)

  dset <- mvpa_multibasis_dataset(train_data = fx$basis_vecs, mask = fx$mask)
  cv <- blocked_cross_validation(design$block_var)
  model <- load_model("sda_notune")
  mspec <- mvpa_model(model, dset, design, model_type = "classification", crossval = cv)

  res <- run_global(mspec)
  expect_s3_class(res, "global_mvpa_result")
  expect_true(!is.null(res$importance_vector))
  expect_equal(length(res$importance_vector), sum(fx$mask > 0) * fx$k)
  expect_true(inherits(res$importance_map, "NeuroVol"))
})

test_that("run_searchlight works with multibasis dataset", {
  skip_if_not_installed("sda")
  fx <- make_multibasis_fixture(D = c(4, 4, 4), n_events = 20, k = 2)

  y <- factor(rep(c("a", "b"), length.out = fx$n_events))
  blocks <- rep(1:4, length.out = fx$n_events)
  design <- mvpa_design(data.frame(y = y, block = blocks), y_train = ~y, block_var = ~block)

  dset <- mvpa_multibasis_dataset(train_data = fx$basis_vecs, mask = fx$mask)
  cv <- blocked_cross_validation(design$block_var)
  model <- load_model("sda_notune")
  mspec <- mvpa_model(model, dset, design, model_type = "classification", crossval = cv)

  res <- suppressWarnings(run_searchlight(mspec, radius = 2, method = "standard"))
  expect_s3_class(res, "searchlight_result")
  expect_true(length(res$results) > 0)
})

# ===========================================================================
# REQ-1 Tests: filter_roi grouped filtering for multibasis datasets
# ===========================================================================

test_that("filter_roi removes entire voxel group when one basis has zero variance", {
  fx <- make_multibasis_fixture(D = c(3, 3, 3), n_events = 6, k = 3)
  dset <- mvpa_multibasis_dataset(train_data = fx$basis_vecs, mask = fx$mask)
  vox <- which(fx$mask > 0)
  sam <- data_sample(dset, vox)
  roi <- as_roi(sam, dset)

  # Inject zero-variance into one basis for one voxel
  vals <- neuroim2::values(roi$train_roi)
  V_phys <- length(vox)
  # Column layout is basis-major: [b1_v1..b1_vV, b2_v1..b2_vV, b3_v1..b3_vV]
  # Voxel 3, basis 2 -> column index V_phys + 3
  vals[, V_phys + 3] <- 999  # constant value = zero variance
  # Rebuild ROI with modified data
  roi$train_roi <- neuroim2::ROIVec(neuroim2::space(roi$train_roi),
                                     neuroim2::coords(roi$train_roi),
                                     data = vals)

  filtered <- rMVPA:::filter_roi(roi)
  fvals <- neuroim2::values(filtered$train_roi)
  # Should have (V_phys - 1) * 3 columns
  expect_equal(ncol(fvals), (V_phys - 1) * 3)
  expect_equal(filtered$basis_count, 3L)
})

test_that("filter_roi all-pass multibasis ROI is unchanged", {
  fx <- make_multibasis_fixture(D = c(3, 3, 3), n_events = 6, k = 2)
  dset <- mvpa_multibasis_dataset(train_data = fx$basis_vecs, mask = fx$mask)
  vox <- which(fx$mask > 0)
  sam <- data_sample(dset, vox)
  roi <- as_roi(sam, dset)

  filtered <- rMVPA:::filter_roi(roi)
  expect_equal(ncol(neuroim2::values(filtered$train_roi)), length(vox) * 2)
  expect_equal(filtered$basis_count, 2L)
})

test_that("filter_roi rejects multibasis ROI when too few physical voxels survive", {
  fx <- make_multibasis_fixture(D = c(2, 2, 2), n_events = 6, k = 2)
  dset <- mvpa_multibasis_dataset(train_data = fx$basis_vecs, mask = fx$mask)
  vox <- which(fx$mask > 0)[1:3]  # only 3 voxels
  sam <- data_sample(dset, vox)
  roi <- as_roi(sam, dset)

  # Make 2 voxels have zero variance -> only 1 survives
  vals <- neuroim2::values(roi$train_roi)
  vals[, 1] <- 0  # voxel 1 basis 1
  vals[, 2] <- 0  # voxel 2 basis 1
  roi$train_roi <- neuroim2::ROIVec(neuroim2::space(roi$train_roi),
                                     neuroim2::coords(roi$train_roi),
                                     data = vals)

  expect_error(rMVPA:::filter_roi(roi, min_voxels = 2), "physical voxels")
})

test_that("filter_roi standard (non-multibasis) path is unchanged", {
  ds <- gen_sample_dataset(c(3, 3, 3), 10)
  vox <- which(ds$dataset$mask > 0)[1:10]
  sam <- data_sample(ds$dataset, vox)
  roi <- as_roi(sam, ds$dataset)

  # Standard ROI has no basis_count field
  expect_null(roi$basis_count)
  filtered <- rMVPA:::filter_roi(roi)
  expect_null(filtered$basis_count)
})

test_that("filter_roi preserves center voxel group in multibasis", {
  fx <- make_multibasis_fixture(D = c(3, 3, 3), n_events = 6, k = 2)
  dset <- mvpa_multibasis_dataset(train_data = fx$basis_vecs, mask = fx$mask)
  vox <- which(fx$mask > 0)
  sam <- data_sample(dset, vox)
  roi <- as_roi(sam, dset)

  # Make voxel 1 zero-variance in basis 1
  vals <- neuroim2::values(roi$train_roi)
  vals[, 1] <- 42  # constant
  roi$train_roi <- neuroim2::ROIVec(neuroim2::space(roi$train_roi),
                                     neuroim2::coords(roi$train_roi),
                                     data = vals)

  V_phys <- length(vox)

  # Without preserve: voxel 1 is removed
  filtered_no_preserve <- rMVPA:::filter_roi(roi)
  expect_equal(ncol(neuroim2::values(filtered_no_preserve$train_roi)), (V_phys - 1) * 2)

  # With preserve: voxel 1 is kept
  center_id <- neuroim2::indices(roi$train_roi)[1]
  filtered_preserve <- rMVPA:::filter_roi(roi, preserve = center_id)
  expect_equal(ncol(neuroim2::values(filtered_preserve$train_roi)), V_phys * 2)
})

test_that("filter_roi k=1 multibasis behaves like standard", {
  fx <- make_multibasis_fixture(D = c(3, 3, 3), n_events = 6, k = 1)
  dset <- mvpa_multibasis_dataset(train_data = fx$basis_vecs, mask = fx$mask)
  vox <- which(fx$mask > 0)
  sam <- data_sample(dset, vox)
  roi <- as_roi(sam, dset)

  filtered <- rMVPA:::filter_roi(roi)
  expect_equal(ncol(neuroim2::values(filtered$train_roi)), length(vox))
})

# ===========================================================================
# REQ-2 Tests: contrast_rsa center voxel aggregation index math
# ===========================================================================

test_that("multibasis center voxel index math is correct", {
  # Simulate the computation from train_model.contrast_rsa_model
  basis_count <- 3L
  V_phys <- 10L
  center_idx <- 5L  # 5th physical voxel

  center_col_indices <- center_idx + (seq_len(basis_count) - 1L) * V_phys
  expect_equal(center_col_indices, c(5, 15, 25))

  # Simulate Delta_sl matrix (V_total x Q)
  V_total <- V_phys * basis_count
  Q <- 4  # contrasts
  Delta_sl <- matrix(seq_len(V_total * Q), nrow = V_total, ncol = Q)

  # Aggregation should be colMeans of 3 rows
  aggregated <- colMeans(Delta_sl[center_col_indices, , drop = FALSE])
  expected <- colMeans(Delta_sl[c(5, 15, 25), , drop = FALSE])
  expect_equal(aggregated, expected)
})

# ===========================================================================
# REQ-5: Regional analysis with multibasis dataset
# ===========================================================================

test_that("run_regional works with multibasis dataset", {
  skip_if_not_installed("sda")
  mb <- gen_multibasis_sample_dataset(D = c(5, 5, 5), n_events = 24, k = 2,
                                       n_classes = 2, n_blocks = 4)

  # Create region mask with 2 regions from the active voxels
  mask_arr <- array(0L, dim(mb$dataset$mask))
  mask_vals <- which(mb$dataset$mask > 0)
  half <- length(mask_vals) %/% 2
  mask_arr[mask_vals[1:half]] <- 1L
  mask_arr[mask_vals[(half + 1):length(mask_vals)]] <- 2L
  region_mask <- neuroim2::NeuroVol(mask_arr,
                                     neuroim2::NeuroSpace(dim(mb$dataset$mask)))

  model <- load_model("sda_notune")
  mspec <- mvpa_model(model, mb$dataset, mb$design,
                      model_type = "classification", crossval = mb$crossval)

  res <- run_regional(mspec, region_mask)
  expect_s3_class(res, "regional_mvpa_result")
  expect_true(nrow(res$performance_table) >= 2)
})

# ===========================================================================
# REQ-6: Vector RSA model with multibasis dataset
# ===========================================================================

test_that("run_regional with vector_rsa works on multibasis dataset", {
  mb <- gen_multibasis_sample_dataset(D = c(4, 4, 4), n_events = 12, k = 2,
                                       n_classes = 3, n_blocks = 3)

  # Create a reference dissimilarity matrix (n_unique_classes x n_unique_classes)
  # with rownames matching the class labels used in the design
  n_classes <- 3
  Dref <- as.matrix(dist(matrix(rnorm(n_classes * 3), n_classes, 3)))
  class_labels <- levels(mb$design$y_train)
  rownames(Dref) <- colnames(Dref) <- class_labels

  # vector_rsa_design needs labels (per observation) and block_var
  labels <- mb$design$y_train
  block_var <- mb$design$block_var

  vdes <- vector_rsa_design(Dref, labels, block_var)
  mspec <- vector_rsa_model(dataset = mb$dataset, design = vdes)

  # Create region mask with 2 regions
  mask_arr <- array(0L, dim(mb$dataset$mask))
  mask_vals <- which(mb$dataset$mask > 0)
  half <- length(mask_vals) %/% 2
  mask_arr[mask_vals[1:half]] <- 1L
  mask_arr[mask_vals[(half + 1):length(mask_vals)]] <- 2L
  region_mask <- neuroim2::NeuroVol(mask_arr,
                                     neuroim2::NeuroSpace(dim(mb$dataset$mask)))

  res <- suppressWarnings(run_regional(mspec, region_mask))
  expect_s3_class(res, "regional_mvpa_result")
  expect_true(nrow(res$performance_table) >= 2)
})

# ===========================================================================
# REQ-7: Naive cross-decoding with multibasis dataset
# ===========================================================================

test_that("run_searchlight with naive_xdec works on multibasis dataset", {
  # naive_xdec needs train and test data from different domains
  fx_train <- make_multibasis_fixture(D = c(4, 4, 4), n_events = 20, k = 2)
  fx_test  <- make_multibasis_fixture(D = c(4, 4, 4), n_events = 20, k = 2)

  dset <- mvpa_multibasis_dataset(
    train_data = fx_train$basis_vecs,
    test_data  = fx_test$basis_vecs,
    mask       = fx_train$mask
  )

  y_train <- factor(rep(c("a", "b"), length.out = 20))
  y_test  <- factor(rep(c("a", "b"), length.out = 20))
  blocks  <- rep(1:4, length.out = 20)

  design <- mvpa_design(
    train_design = data.frame(y = y_train, block = blocks),
    test_design  = data.frame(y = y_test),
    y_train = ~ y, y_test = ~ y, block_var = ~ block
  )

  mspec <- naive_xdec_model(dataset = dset, design = design)

  res <- suppressWarnings(run_searchlight(mspec, radius = 2, method = "standard"))
  expect_s3_class(res, "searchlight_result")
  expect_true(length(res$results) > 0)
})

# ===========================================================================
# REQ-8: MANOVA model with multibasis dataset
# ===========================================================================

test_that("run_searchlight with manova_model works on multibasis dataset", {
  skip_if_not_installed("ffmanova")
  mb <- gen_multibasis_sample_dataset(D = c(4, 4, 4), n_events = 20, k = 2,
                                       n_classes = 2, n_blocks = 4)

  # manova_design needs a formula and a named data list
  # The response is "response" (assigned internally by train_model.manova_model)
  # so the formula LHS is always "y" and we supply a predictor variable
  group <- factor(rep(c("a", "b"), length.out = 20))
  mdes <- manova_design(~ group, data = list(group = group))
  mspec <- manova_model(dataset = mb$dataset, design = mdes)

  res <- suppressWarnings(run_searchlight(mspec, radius = 2, method = "standard"))
  expect_s3_class(res, "searchlight_result")
  expect_true(length(res$results) > 0)
})
