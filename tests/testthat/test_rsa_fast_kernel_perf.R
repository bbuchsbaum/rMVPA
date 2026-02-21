testthat::skip_if_not_installed("neuroim2")

rsa_fast_perf_env_int <- function(name, default) {
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

rsa_fast_perf_env_numeric <- function(name, default) {
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

rsa_fast_benchmark_pair <- function(run_baseline, run_candidate, nrep = 3L) {
  invisible(run_baseline())
  invisible(run_candidate())

  baseline <- numeric(nrep)
  candidate <- numeric(nrep)

  for (i in seq_len(nrep)) {
    order <- if ((i %% 2L) == 1L) c("baseline", "candidate") else c("candidate", "baseline")
    for (path in order) {
      set.seed(9100 + i)
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

build_rsa_fast_perf_model <- function(fast_enabled, D = c(10, 10, 10), nobs = 36, radius = 2) {
  old_opt <- options(
    rMVPA.searchlight_mode = if (isTRUE(fast_enabled)) "fast" else "legacy",
    rMVPA.rsa_fast_kernel = NULL,
    rMVPA.warn_legacy_options = FALSE
  )
  on.exit(options(old_opt), add = TRUE)

  ds <- gen_sample_dataset(D = D, nobs = nobs, blocks = 3)
  D1 <- dist(matrix(rnorm(nobs * nobs), nobs, nobs))
  D2 <- dist(matrix(rnorm(nobs * nobs), nobs, nobs))
  rdes <- rsa_design(~ D1 + D2, list(D1 = D1, D2 = D2, block = ds$design$block_var), block_var = "block")
  mspec <- rsa_model(ds$dataset, rdes, regtype = "lm", distmethod = "pearson", check_collinearity = FALSE)

  list(mspec = mspec, radius = radius)
}

run_rsa_fast_searchlight <- function(fast_enabled, D = c(10, 10, 10), nobs = 36, radius = 2) {
  built <- build_rsa_fast_perf_model(
    fast_enabled = fast_enabled,
    D = D,
    nobs = nobs,
    radius = radius
  )
  invisible(run_searchlight(built$mspec, radius = built$radius, method = "standard"))
}

test_that("rsa fast kernel does not regress searchlight runtime", {
  skip_on_cran()
  skip_if_not_perf_tests()

  nrep <- rsa_fast_perf_env_int("RMVPA_PERF_REP", 2L)
  max_slowdown <- rsa_fast_perf_env_numeric("RMVPA_RSA_FAST_MAX_SLOWDOWN", 0.20)

  perf <- rsa_fast_benchmark_pair(
    run_baseline = function() run_rsa_fast_searchlight(fast_enabled = FALSE),
    run_candidate = function() run_rsa_fast_searchlight(fast_enabled = TRUE),
    nrep = nrep
  )

  message(sprintf(
    "rsa fast kernel guardrail: median fast=%.3fs baseline=%.3fs speedup=%.2fx slowdown=%.1f%%",
    perf$median_candidate,
    perf$median_baseline,
    perf$speedup,
    100 * perf$slowdown
  ))

  expect_lte(perf$median_candidate, perf$median_baseline * (1 + max_slowdown))
})
