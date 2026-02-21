testthat::skip_if_not_installed("neuroim2")

clustered_perf_env_int <- function(name, default) {
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

clustered_perf_env_numeric <- function(name, default) {
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

clustered_benchmark_pair <- function(run_baseline, run_candidate, nrep = 3L) {
  invisible(run_baseline())
  invisible(run_candidate())

  baseline <- numeric(nrep)
  candidate <- numeric(nrep)

  for (i in seq_len(nrep)) {
    order <- if ((i %% 2L) == 1L) c("baseline", "candidate") else c("candidate", "baseline")
    for (path in order) {
      set.seed(9700 + i)
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

build_clustered_perf_dataset <- function(D = c(35, 35, 35), nobs = 12, K = 300, blocks = 4) {
  ds <- gen_clustered_sample_dataset(D = D, nobs = nobs, K = K, nlevels = 2, blocks = blocks)
  ds$dataset
}

run_clustered_neighbor_build <- function(dset, radius, k, fast_enabled) {
  old_opt <- options(
    rMVPA.searchlight_mode = if (isTRUE(fast_enabled)) "fast" else "legacy",
    rMVPA.clustered_nn_fastpath = NULL,
    rMVPA.warn_legacy_options = FALSE
  )
  on.exit(options(old_opt), add = TRUE)
  invisible(get_searchlight(dset, type = "standard", radius = radius, k = k))
}

test_that("clustered fastpath does not regress neighbor construction at larger K", {
  skip_on_cran()
  skip_if_not_perf_tests()

  set.seed(7601)
  dset <- build_clustered_perf_dataset()
  nrep <- clustered_perf_env_int("RMVPA_PERF_REP", 2L)
  max_slowdown <- clustered_perf_env_numeric("RMVPA_CLUSTERED_FASTPATH_MAX_SLOWDOWN", 0.20)

  perf <- clustered_benchmark_pair(
    run_baseline = function() run_clustered_neighbor_build(dset, radius = NULL, k = 6, fast_enabled = FALSE),
    run_candidate = function() run_clustered_neighbor_build(dset, radius = NULL, k = 6, fast_enabled = TRUE),
    nrep = nrep
  )

  message(sprintf(
    "clustered fastpath guardrail: median fast=%.3fs dense=%.3fs speedup=%.2fx slowdown=%.1f%%",
    perf$median_candidate,
    perf$median_baseline,
    perf$speedup,
    100 * perf$slowdown
  ))

  expect_lte(perf$median_candidate, perf$median_baseline * (1 + max_slowdown))
})

test_that("clustered fastpath scales on nightly larger-k workload", {
  skip_on_cran()
  skip_if_not_nightly_perf_tests()

  set.seed(7602)
  dset <- build_clustered_perf_dataset(D = c(45, 45, 45), nobs = 10, K = 450, blocks = 5)
  nrep <- clustered_perf_env_int("RMVPA_NIGHTLY_PERF_REP", 2L)
  max_slowdown <- clustered_perf_env_numeric("RMVPA_NIGHTLY_CLUSTERED_FASTPATH_MAX_SLOWDOWN", 0.15)

  perf <- clustered_benchmark_pair(
    run_baseline = function() run_clustered_neighbor_build(dset, radius = NULL, k = 6, fast_enabled = FALSE),
    run_candidate = function() run_clustered_neighbor_build(dset, radius = NULL, k = 6, fast_enabled = TRUE),
    nrep = nrep
  )

  message(sprintf(
    "clustered fastpath nightly: median fast=%.3fs dense=%.3fs speedup=%.2fx slowdown=%.1f%%",
    perf$median_candidate,
    perf$median_baseline,
    perf$speedup,
    100 * perf$slowdown
  ))

  expect_lte(perf$median_candidate, perf$median_baseline * (1 + max_slowdown))
})
