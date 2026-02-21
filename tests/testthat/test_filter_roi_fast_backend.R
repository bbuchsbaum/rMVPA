testthat::skip_if_not_installed("neuroim2")

run_filter_backend <- function(roi, fast_enabled, preserve = NULL, min_voxels = 2) {
  old <- options(
    rMVPA.searchlight_mode = if (isTRUE(fast_enabled)) "fast" else "legacy",
    rMVPA.fast_filter_roi = NULL
  )
  on.exit(options(old), add = TRUE)
  try(rMVPA:::filter_roi(roi, preserve = preserve, min_voxels = min_voxels), silent = TRUE)
}

expect_filtered_equivalent <- function(reference, candidate) {
  expect_identical(inherits(reference, "try-error"), inherits(candidate, "try-error"))
  if (inherits(reference, "try-error")) {
    expect_match(as.character(reference), "filter_roi")
    expect_match(as.character(candidate), "filter_roi")
    return(invisible(TRUE))
  }

  expect_equal(neuroim2::indices(reference$train_roi), neuroim2::indices(candidate$train_roi))
  expect_equal(
    neuroim2::values(reference$train_roi),
    neuroim2::values(candidate$train_roi),
    tolerance = 1e-12
  )
  expect_identical(reference$basis_count, candidate$basis_count)

  if (is.null(reference$test_roi)) {
    expect_null(candidate$test_roi)
  } else {
    expect_equal(neuroim2::indices(reference$test_roi), neuroim2::indices(candidate$test_roi))
    expect_equal(
      neuroim2::values(reference$test_roi),
      neuroim2::values(candidate$test_roi),
      tolerance = 1e-12
    )
  }

  invisible(TRUE)
}

test_that("fast filter backend matches slow backend on standard ROIVec", {
  ds <- gen_sample_dataset(D = c(4, 4, 4), nobs = 18, nlevels = 2, blocks = 3)
  vox <- which(ds$dataset$mask > 0)[1:24]
  roi <- as_roi(data_sample(ds$dataset, vox), ds$dataset)

  vals <- neuroim2::values(roi$train_roi)
  vals[, 2] <- NA_real_
  vals[, 5] <- 7
  vals[nrow(vals), 6] <- vals[nrow(vals), 6] + 1e-12
  roi$train_roi <- neuroim2::ROIVec(
    neuroim2::space(roi$train_roi),
    neuroim2::coords(roi$train_roi),
    data = vals
  )

  ref <- run_filter_backend(roi, fast_enabled = FALSE, preserve = neuroim2::indices(roi$train_roi)[1])
  cand <- run_filter_backend(roi, fast_enabled = TRUE, preserve = neuroim2::indices(roi$train_roi)[1])
  expect_filtered_equivalent(ref, cand)
})

test_that("fast filter backend matches slow backend for min_voxels error behavior", {
  ds <- gen_sample_dataset(D = c(3, 3, 3), nobs = 12, nlevels = 2, blocks = 3)
  vox <- which(ds$dataset$mask > 0)[1:6]
  roi <- as_roi(data_sample(ds$dataset, vox), ds$dataset)

  vals <- neuroim2::values(roi$train_roi)
  vals[, 1:5] <- 0
  roi$train_roi <- neuroim2::ROIVec(
    neuroim2::space(roi$train_roi),
    neuroim2::coords(roi$train_roi),
    data = vals
  )

  ref <- run_filter_backend(roi, fast_enabled = FALSE, min_voxels = 2)
  cand <- run_filter_backend(roi, fast_enabled = TRUE, min_voxels = 2)
  expect_filtered_equivalent(ref, cand)
})

test_that("fast filter backend matches slow backend for multibasis grouped semantics", {
  fx <- make_multibasis_fixture(D = c(3, 3, 3), n_events = 8, k = 3)
  dset <- mvpa_multibasis_dataset(train_data = fx$basis_vecs, mask = fx$mask)
  vox <- which(fx$mask > 0)[1:12]
  roi <- as_roi(data_sample(dset, vox), dset)

  vals <- neuroim2::values(roi$train_roi)
  V_phys <- length(vox)
  vals[, V_phys + 2] <- 3.14
  vals[, (2 * V_phys) + 4] <- NA_real_
  roi$train_roi <- neuroim2::ROIVec(
    neuroim2::space(roi$train_roi),
    neuroim2::coords(roi$train_roi),
    data = vals
  )

  center_id <- neuroim2::indices(roi$train_roi)[1]
  ref <- run_filter_backend(roi, fast_enabled = FALSE, preserve = center_id)
  cand <- run_filter_backend(roi, fast_enabled = TRUE, preserve = center_id)
  expect_filtered_equivalent(ref, cand)
})

build_fast_filter_mvpa_spec <- function() {
  ds <- gen_sample_dataset(D = c(4, 4, 4), nobs = 24, nlevels = 2, blocks = 4)
  cv <- blocked_cross_validation(ds$design$block_var)
  mvpa_model(
    model = load_model("sda_notune"),
    dataset = ds$dataset,
    design = ds$design,
    model_type = "classification",
    crossval = cv
  )
}

build_fast_filter_rsa_spec <- function() {
  ds <- gen_sample_dataset(D = c(4, 4, 4), nobs = 24, nlevels = 3, blocks = 4)
  dmat <- dist(matrix(rnorm(24 * 10), nrow = 24, ncol = 10))
  rdes <- rsa_design(~ dmat, list(dmat = dmat), block_var = ds$design$block_var)
  rsa_model(ds$dataset, rdes, regtype = "pearson")
}

build_fast_filter_naive_spec <- function() {
  ds <- gen_sample_dataset(
    D = c(4, 4, 4),
    nobs = 24,
    nlevels = 3,
    blocks = 4,
    external_test = TRUE,
    ntest_obs = 24
  )
  naive_xdec_model(ds$dataset, ds$design)
}

run_searchlight_fast_filter <- function(mspec, radius, fast_enabled) {
  old <- options(
    rMVPA.searchlight_mode = if (isTRUE(fast_enabled)) "fast" else "legacy",
    rMVPA.fast_filter_roi = NULL
  )
  on.exit(options(old), add = TRUE)
  run_searchlight(mspec, radius = radius, method = "standard", backend = "default")
}

test_that("fast filter backend preserves searchlight outputs across model families", {
  specs <- list(
    mvpa = build_fast_filter_mvpa_spec(),
    rsa = build_fast_filter_rsa_spec(),
    naive = build_fast_filter_naive_spec()
  )

  for (nm in names(specs)) {
    set.seed(4300)
    base <- run_searchlight_fast_filter(specs[[nm]], radius = 2, fast_enabled = FALSE)
    set.seed(4300)
    fast <- run_searchlight_fast_filter(specs[[nm]], radius = 2, fast_enabled = TRUE)

    expect_searchlight_parity(
      reference = base,
      candidate = fast,
      atol = 1e-10,
      rtol = 1e-8
    )
  }
})

test_that("fast filter backend is not slower on medium ROI matrices (guardrail)", {
  skip_on_cran()
  skip_if_not_perf_tests()

  ds <- gen_sample_dataset(D = c(8, 8, 8), nobs = 80, nlevels = 2, blocks = 5)
  vox <- which(ds$dataset$mask > 0)[1:400]
  roi <- as_roi(data_sample(ds$dataset, vox), ds$dataset)

  vals <- neuroim2::values(roi$train_roi)
  vals[, 5] <- 1
  vals[, 9] <- NA_real_
  roi$train_roi <- neuroim2::ROIVec(
    neuroim2::space(roi$train_roi),
    neuroim2::coords(roi$train_roi),
    data = vals
  )

  nrep <- 3L
  slow_t <- numeric(nrep)
  fast_t <- numeric(nrep)

  for (i in seq_len(nrep)) {
    old_slow <- options(
      rMVPA.searchlight_mode = "legacy",
      rMVPA.fast_filter_roi = NULL
    )
    on.exit(options(old_slow), add = TRUE)
    slow_t[i] <- as.numeric(system.time(
      rMVPA:::filter_roi(roi, preserve = neuroim2::indices(roi$train_roi)[1], min_voxels = 2)
    )["elapsed"])

    old_fast <- options(
      rMVPA.searchlight_mode = "fast",
      rMVPA.fast_filter_roi = NULL
    )
    on.exit(options(old_fast), add = TRUE)
    fast_t[i] <- as.numeric(system.time(
      rMVPA:::filter_roi(roi, preserve = neuroim2::indices(roi$train_roi)[1], min_voxels = 2)
    )["elapsed"])
  }

  med_slow <- stats::median(slow_t)
  med_fast <- stats::median(fast_t)
  message(sprintf(
    "fast_filter guardrail: median fast=%.4fs slow=%.4fs speedup=%.2fx",
    med_fast, med_slow, med_slow / med_fast
  ))
  expect_lte(med_fast, med_slow * 1.20)
})
