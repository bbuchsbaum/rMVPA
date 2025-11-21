library(testthat)
library(rMVPA)

make_small_dataset <- function(D = c(3, 3, 3), nobs = 6) {
  train <- neuroim2::NeuroVec(array(rnorm(prod(D) * nobs), c(D, nobs)), neuroim2::NeuroSpace(c(D, nobs)))
  mask  <- neuroim2::NeuroVol(array(1, D), neuroim2::NeuroSpace(D))
  block_var <- rep(seq_len(3), length.out = nobs)
  design <- mvpa_design(
    data.frame(Y = factor(rep(letters[1:2], length.out = nobs)), block = block_var),
    y_train = ~Y,
    block_var = "block"
  )
  list(train = train, mask = mask, design = design)
}

test_that("get_searchlight + get_samples yields multi-voxel ROI", {
  ds <- make_small_dataset()
  dset <- mvpa_dataset(ds$train, mask = ds$mask)
  sl <- get_searchlight(dset, type = "randomized", radius = 1)
  expect_gt(length(sl), 0)

  sf <- rMVPA:::get_samples(dset, sl)
  vox_first <- sf$sample[[1]]$vox
  expect_gt(length(vox_first), 1)
})

test_that("mvpa_iterate keeps ROIs with >=2 voxels", {
  ds <- make_small_dataset()
  dset <- mvpa_dataset(ds$train, mask = ds$mask)
  base_mod <- load_model("corclass")
  model <- mvpa_model(base_mod, dset, ds$design)
  sl <- get_searchlight(dset, type = "randomized", radius = 1)
  ids <- vapply(sl, function(x) x@parent_index, integer(1))
  res <- suppressWarnings(
    mvpa_iterate(model, sl, ids = ids, analysis_type = "searchlight", drop_probs = TRUE)
  )
  expect_true(any(!res$error))  # at least one ROI succeeds
})
