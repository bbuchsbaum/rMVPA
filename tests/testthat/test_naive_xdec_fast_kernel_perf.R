testthat::skip_if_not_installed("neuroim2")

naive_fast_perf_env_int <- function(name, default) {
  raw <- Sys.getenv(name, "")
  if (!nzchar(raw)) {
    return(as.integer(default))
  }
  val <- suppressWarnings(as.integer(raw))
  if (!is.finite(val) || is.na(val) || val < 1L) {
    as.integer(default)
  } else {
    as.integer(val)
  }
}

naive_fast_perf_env_numeric <- function(name, default) {
  raw <- Sys.getenv(name, "")
  if (!nzchar(raw)) {
    return(as.numeric(default))
  }
  val <- suppressWarnings(as.numeric(raw))
  if (!is.finite(val) || is.na(val) || val < 0) {
    as.numeric(default)
  } else {
    as.numeric(val)
  }
}

naive_fast_benchmark_pair <- function(run_baseline, run_candidate, nrep = 3L) {
  invisible(run_baseline())
  invisible(run_candidate())

  baseline <- numeric(nrep)
  candidate <- numeric(nrep)

  for (i in seq_len(nrep)) {
    order <- if ((i %% 2L) == 1L) c("baseline", "candidate") else c("candidate", "baseline")
    for (path in order) {
      set.seed(9500 + i)
      gc(verbose = FALSE)
      elapsed <- as.numeric(system.time({
        if (identical(path, "baseline")) {
          run_baseline()
        } else {
          run_candidate()
        }
      })["elapsed"])
      if (identical(path, "baseline")) {
        baseline[i] <- elapsed
      } else {
        candidate[i] <- elapsed
      }
    }
  }

  med_base <- stats::median(baseline)
  med_cand <- stats::median(candidate)
  list(
    baseline = baseline,
    candidate = candidate,
    median_baseline = med_base,
    median_candidate = med_cand,
    speedup = med_base / med_cand,
    slowdown = (med_cand / med_base) - 1
  )
}

build_naive_fast_perf_model <- function(fast_enabled, D = c(8, 8, 8), nobs = 54, blocks = 6, radius = 2) {
  old_opt <- options(
    rMVPA.searchlight_mode = if (isTRUE(fast_enabled)) "fast" else "legacy",
    rMVPA.naive_xdec_fast_kernel = NULL,
    rMVPA.warn_legacy_options = FALSE
  )
  on.exit(options(old_opt), add = TRUE)

  ds <- gen_sample_dataset(
    D = D,
    nobs = nobs,
    nlevels = 3,
    blocks = blocks,
    external_test = TRUE,
    ntest_obs = nobs
  )
  mspec <- naive_xdec_model(ds$dataset, ds$design, return_predictions = FALSE)
  list(mspec = mspec, radius = radius)
}

run_naive_fast_searchlight <- function(fast_enabled, D = c(8, 8, 8), nobs = 54, blocks = 6, radius = 2) {
  built <- build_naive_fast_perf_model(
    fast_enabled = fast_enabled,
    D = D,
    nobs = nobs,
    blocks = blocks,
    radius = radius
  )
  invisible(run_searchlight(built$mspec, radius = built$radius, method = "standard"))
}

test_that("naive_xdec fast kernel does not regress searchlight runtime", {
  skip_on_cran()
  skip_if_not_perf_tests()

  nrep <- naive_fast_perf_env_int("RMVPA_PERF_REP", 2L)
  max_slowdown <- naive_fast_perf_env_numeric("RMVPA_NAIVE_FAST_MAX_SLOWDOWN", 0.20)

  perf <- naive_fast_benchmark_pair(
    run_baseline = function() run_naive_fast_searchlight(fast_enabled = FALSE),
    run_candidate = function() run_naive_fast_searchlight(fast_enabled = TRUE),
    nrep = nrep
  )

  message(sprintf(
    "naive_xdec fast kernel guardrail: median fast=%.3fs baseline=%.3fs speedup=%.2fx slowdown=%.1f%%",
    perf$median_candidate,
    perf$median_baseline,
    perf$speedup,
    100 * perf$slowdown
  ))

  expect_lte(perf$median_candidate, perf$median_baseline * (1 + max_slowdown))
})
