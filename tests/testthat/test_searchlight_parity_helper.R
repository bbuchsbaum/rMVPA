test_that("searchlight_allclose handles absolute tolerance around zero", {
  ref <- c(0, 1e-10, -2e-10, NA_real_, Inf, -Inf)
  cand <- c(5e-10, -2e-10, 4e-10, NA_real_, Inf, -Inf)

  cmp <- searchlight_allclose(ref, cand, atol = 1e-9, rtol = 0)
  expect_true(cmp$ok)
  expect_equal(cmp$n_bad_finite, 0L)
  expect_equal(length(cmp$non_comparable_indices), 0L)
})

test_that("searchlight_allclose handles relative tolerance at scale", {
  ref <- c(1000, -2500, 40)
  cand <- c(1001, -2501, 40.02)

  cmp <- searchlight_allclose(ref, cand, atol = 1e-8, rtol = 1e-3)
  expect_true(cmp$ok)
  expect_gt(cmp$finite_overlap, 0L)
})

test_that("expect_searchlight_parity compares fake metric maps and catches mismatches", {
  ref <- list(
    metrics = c("Accuracy"),
    results = list(Accuracy = c(0.25, NA_real_, 0.75, Inf)),
    active_voxels = 4L
  )
  cand_ok <- list(
    metrics = c("Accuracy"),
    results = list(Accuracy = c(0.2500001, NA_real_, 0.7500002, Inf)),
    active_voxels = 4L
  )
  cand_bad <- list(
    metrics = c("Accuracy"),
    results = list(Accuracy = c(0.35, NA_real_, 0.75, Inf)),
    active_voxels = 4L
  )

  expect_no_error(
    expect_searchlight_parity(ref, cand_ok, atol = 1e-7, rtol = 1e-5)
  )
  expect_error(
    expect_searchlight_parity(ref, cand_bad, atol = 1e-7, rtol = 1e-5),
    "parity failed"
  )
})
