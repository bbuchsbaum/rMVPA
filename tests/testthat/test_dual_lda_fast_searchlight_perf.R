test_that("dual_lda incremental searchlight is not slower than full recompute (guardrail)", {
  skip_on_cran()
  skip_if_not_perf_tests()

  set.seed(5501)
  built <- build_dual_lda_mspec(
    D = c(8, 8, 8),
    nobs = 80,
    nlevels = 3,
    blocks = 5,
    gamma = 1e-2
  )

  radius <- 3

  # Warm-up to reduce one-time load effects.
  invisible(run_searchlight(built$mspec, radius = radius, method = "standard", incremental = TRUE))
  invisible(run_searchlight(built$mspec, radius = radius, method = "standard", incremental = FALSE))

  nrep <- 3L
  inc_t <- numeric(nrep)
  full_t <- numeric(nrep)

  for (i in seq_len(nrep)) {
    inc_t[i] <- as.numeric(system.time(
      run_searchlight(built$mspec, radius = radius, method = "standard", incremental = TRUE)
    )["elapsed"])
    full_t[i] <- as.numeric(system.time(
      run_searchlight(built$mspec, radius = radius, method = "standard", incremental = FALSE)
    )["elapsed"])
  }

  med_inc <- stats::median(inc_t)
  med_full <- stats::median(full_t)
  speedup <- med_full / med_inc

  message(sprintf(
    "dual_lda perf guardrail: median incremental=%.3fs, full=%.3fs, speedup=%.2fx",
    med_inc, med_full, speedup
  ))

  # Soft performance guardrail for heterogeneous CI hardware.
  # Allows small run-to-run noise while catching clear regressions.
  expect_lte(med_inc, med_full * 1.10)
})

test_that("dual_lda incremental searchlight scales on larger nightly workload", {
  skip_on_cran()
  skip_if_not_nightly_perf_tests()

  set.seed(6601)
  built <- build_dual_lda_mspec(
    D = c(10, 10, 10),
    nobs = 120,
    nlevels = 3,
    blocks = 6,
    gamma = 1e-2
  )

  radius <- 4

  # Warm-up to reduce one-time load effects.
  invisible(run_searchlight(built$mspec, radius = radius, method = "standard", incremental = TRUE))
  invisible(run_searchlight(built$mspec, radius = radius, method = "standard", incremental = FALSE))

  nrep <- 2L
  inc_t <- numeric(nrep)
  full_t <- numeric(nrep)

  for (i in seq_len(nrep)) {
    inc_t[i] <- as.numeric(system.time(
      run_searchlight(built$mspec, radius = radius, method = "standard", incremental = TRUE)
    )["elapsed"])
    full_t[i] <- as.numeric(system.time(
      run_searchlight(built$mspec, radius = radius, method = "standard", incremental = FALSE)
    )["elapsed"])
  }

  med_inc <- stats::median(inc_t)
  med_full <- stats::median(full_t)
  speedup <- med_full / med_inc

  message(sprintf(
    "dual_lda nightly perf: median incremental=%.3fs, full=%.3fs, speedup=%.2fx",
    med_inc, med_full, speedup
  ))

  # Nightly runs on larger geometry where noise is still possible.
  # Allow moderate variance while detecting substantial regressions.
  expect_lte(med_inc, med_full * 1.20)
})
