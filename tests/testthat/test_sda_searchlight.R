# The sda_notune searchlight engine computes per-voxel statistics once per
# fold and runs only the cross-column steps per sphere, in the general path's
# voxel order, so its maps must equal the general path's.
testthat::skip_if_not_installed("neuroim2")

sda_sl_spec <- function(seed, nlevels = 3, nobs = 48, class_metrics = FALSE,
                        scale = 1, offset = 0) {
  set.seed(seed)
  ds <- gen_sample_dataset(D = c(5, 5, 5), nobs = nobs, nlevels = nlevels, blocks = 4)
  if (scale != 1 || offset != 0) {
    arr <- as.array(ds$dataset$train_data) * scale + offset
    ds$dataset$train_data <- neuroim2::NeuroVec(arr, neuroim2::space(ds$dataset$train_data))
  }
  mvpa_model(load_model("sda_notune"), ds$dataset, ds$design, "classification",
             crossval = blocked_cross_validation(ds$design$block_var),
             class_metrics = class_metrics)
}

expect_sda_engine_parity <- function(ms, radius = 2) {
  fast <- run_searchlight(ms, radius = radius, engine = "sda_fast", backend = "default")
  ref <- run_searchlight(ms, radius = radius, engine = "legacy", backend = "default")
  expect_identical(attr(fast, "searchlight_engine"), "sda_fast")
  expect_identical(sort(names(fast$results)), sort(names(ref$results)))
  for (m in names(ref$results)) {
    a <- as.numeric(neuroim2::values(fast$results[[m]]))
    b <- as.numeric(neuroim2::values(ref$results[[m]]))
    expect_identical(is.na(a), is.na(b), info = m)
    expect_equal(a, b, tolerance = 1e-12, info = m)
  }
}

test_that("sda engine matches the general path (multiclass, binary, per-class AUC)", {
  expect_sda_engine_parity(sda_sl_spec(1201))
  expect_sda_engine_parity(sda_sl_spec(1202, nlevels = 2, nobs = 40))
  expect_sda_engine_parity(sda_sl_spec(1203, class_metrics = TRUE))
})

test_that("sda engine matches the general path on raw-BOLD-scale data", {
  expect_sda_engine_parity(sda_sl_spec(1204, scale = 20, offset = 1500))
})

test_that("auto selects the sda engine and falls back on missing values", {
  ms <- sda_sl_spec(1205)
  expect_identical(rMVPA:::.resolve_searchlight_engine(ms, "standard", "auto"), "sda_fast")
  arr <- as.array(ms$dataset$train_data)
  arr[2, 2, 2, 5] <- NA
  ms$dataset$train_data <- neuroim2::NeuroVec(arr, neuroim2::space(ms$dataset$train_data))
  res <- expect_no_warning(run_searchlight(ms, radius = 2, backend = "default"))
  expect_identical(attr(res, "searchlight_engine"), "legacy")
})

test_that("a subset fit from shared column statistics is bit-identical to a direct fit", {
  set.seed(1206)
  y <- factor(rep(1:4, 20))
  X <- matrix(rnorm(80 * 60), 80) + outer(as.integer(y), rnorm(60))
  st <- rMVPA:::.sda_column_stats(X, y)
  for (cols in list(c(5, 2, 40, 33, 17, 9), 1:60, sample(60, 45))) {
    expect_identical(unclass(rMVPA:::.sda_fit_columns(st, cols)),
                     unclass(rMVPA:::.sda_native_fit(X[, cols], y)))
  }
})
