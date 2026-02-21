testthat::skip_if_not_installed("neuroim2")

perf_env_int <- function(name, default) {
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

perf_env_numeric <- function(name, default) {
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

benchmark_pair <- function(run_baseline, run_candidate, nrep = 3L) {
  invisible(run_baseline())
  invisible(run_candidate())

  baseline <- numeric(nrep)
  candidate <- numeric(nrep)

  for (i in seq_len(nrep)) {
    order <- if ((i %% 2L) == 1L) {
      c("baseline", "candidate")
    } else {
      c("candidate", "baseline")
    }
    for (path in order) {
      set.seed(8100 + i)
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

build_framework_perf_mvpa_spec <- function(D = c(6, 6, 6), nobs = 54, nlevels = 3, blocks = 6) {
  ds <- gen_sample_dataset(D = D, nobs = nobs, nlevels = nlevels, blocks = blocks)
  cv <- blocked_cross_validation(ds$design$block_var)
  mvpa_model(
    model = load_model("sda_notune"),
    dataset = ds$dataset,
    design = ds$design,
    model_type = "classification",
    crossval = cv
  )
}

build_framework_perf_rsa_spec <- function(D = c(6, 6, 6), nobs = 42, blocks = 6) {
  ds <- gen_sample_dataset(D = D, nobs = nobs, nlevels = 3, blocks = blocks)
  dmat <- dist(matrix(rnorm(nobs * 16), nrow = nobs, ncol = 16))
  rdes <- rsa_design(~ dmat, list(dmat = dmat), block_var = ds$design$block_var)
  rsa_model(ds$dataset, rdes, regtype = "pearson")
}

build_framework_perf_naive_xdec_spec <- function(D = c(6, 6, 6), nobs = 48, blocks = 6) {
  ds <- gen_sample_dataset(
    D = D,
    nobs = nobs,
    nlevels = 3,
    blocks = blocks,
    external_test = TRUE,
    ntest_obs = nobs
  )
  naive_xdec_model(ds$dataset, ds$design, return_predictions = FALSE)
}

build_framework_perf_geometry_dataset <- function(D = c(80, 80, 80), nobs = 2, blocks = 2) {
  ds <- gen_sample_dataset(D = D, nobs = nobs, nlevels = 2, blocks = blocks)
  ds$dataset
}

resolve_perf_mode <- function(mode) {
  match.arg(mode, c("legacy", "fast"))
}

perf_mode_options <- function(mode) {
  mode <- resolve_perf_mode(mode)
  list(
    rMVPA.searchlight_mode = mode,
    rMVPA.searchlight_profile = NULL,
    rMVPA.warn_legacy_options = FALSE
  )
}

run_with_perf_mode <- function(mspec, radius, mode = c("legacy", "fast")) {
  old_opt <- options(perf_mode_options(mode))
  on.exit(options(old_opt), add = TRUE)
  if (requireNamespace("future", quietly = TRUE)) {
    old_plan <- future::plan()
    on.exit(future::plan(old_plan), add = TRUE)
    future::plan(future::sequential)
  }
  run_searchlight(mspec, radius = radius, method = "standard", backend = "default")
}

run_with_fold_cache <- function(mspec, radius, fold_cache_enabled) {
  mode <- if (isTRUE(fold_cache_enabled)) "fast" else "legacy"
  run_with_perf_mode(mspec, radius = radius, mode = mode)
}

run_with_geometry_cache <- function(dset, radius, cache_enabled) {
  old_opt <- options(
    rMVPA.searchlight_mode = "fast",
    rMVPA.searchlight_profile = NULL,
    rMVPA.searchlight_geometry_cache = isTRUE(cache_enabled),
    rMVPA.searchlight_geometry_cache_max_entries = 8L,
    rMVPA.warn_legacy_options = FALSE
  )
  on.exit(options(old_opt), add = TRUE)
  rMVPA:::.searchlight_geometry_cache_clear()

  invisible(get_searchlight(dset, type = "standard", radius = radius))
  invisible(get_searchlight(dset, type = "standard", radius = radius))
}

test_that("mvpa_model fold cache path does not regress runtime (guardrail)", {
  skip_on_cran()
  skip_if_not_perf_tests()

  set.seed(7101)
  mspec <- build_framework_perf_mvpa_spec()
  radius <- 2
  nrep <- perf_env_int("RMVPA_PERF_REP", 3L)
  max_slowdown <- perf_env_numeric("RMVPA_FRAMEWORK_PERF_MAX_SLOWDOWN", 0.20)

  perf <- benchmark_pair(
    run_baseline = function() run_with_fold_cache(mspec, radius = radius, fold_cache_enabled = FALSE),
    run_candidate = function() run_with_fold_cache(mspec, radius = radius, fold_cache_enabled = TRUE),
    nrep = nrep
  )

  message(sprintf(
    "mvpa_model perf guardrail: median cached=%.3fs, uncached=%.3fs, speedup=%.2fx, slowdown=%.1f%%",
    perf$median_candidate,
    perf$median_baseline,
    perf$speedup,
    100 * perf$slowdown
  ))

  expect_lte(perf$median_candidate, perf$median_baseline * (1 + max_slowdown))
})

test_that("rsa_model path does not regress when fold cache option is enabled (guardrail)", {
  skip_on_cran()
  skip_if_not_perf_tests()

  set.seed(7102)
  mspec <- build_framework_perf_rsa_spec()
  radius <- 2
  nrep <- perf_env_int("RMVPA_PERF_REP", 3L)
  max_slowdown <- perf_env_numeric("RMVPA_FRAMEWORK_RSA_PERF_MAX_SLOWDOWN", 0.25)

  perf <- benchmark_pair(
    run_baseline = function() run_with_fold_cache(mspec, radius = radius, fold_cache_enabled = FALSE),
    run_candidate = function() run_with_fold_cache(mspec, radius = radius, fold_cache_enabled = TRUE),
    nrep = nrep
  )

  message(sprintf(
    "rsa_model perf guardrail: median cache_on=%.3fs, cache_off=%.3fs, speedup=%.2fx, slowdown=%.1f%%",
    perf$median_candidate,
    perf$median_baseline,
    perf$speedup,
    100 * perf$slowdown
  ))

  expect_lte(perf$median_candidate, perf$median_baseline * (1 + max_slowdown))
})

test_that("geometry cache improves repeated searchlight setup path (guardrail)", {
  skip_on_cran()
  skip_if_not_perf_tests()

  set.seed(7150)
  dset <- build_framework_perf_geometry_dataset()
  radius <- 8
  nrep <- perf_env_int("RMVPA_PERF_REP", 3L)
  max_slowdown <- perf_env_numeric("RMVPA_GEOMETRY_CACHE_MAX_SLOWDOWN", 0.10)

  perf <- benchmark_pair(
    run_baseline = function() run_with_geometry_cache(dset, radius = radius, cache_enabled = FALSE),
    run_candidate = function() run_with_geometry_cache(dset, radius = radius, cache_enabled = TRUE),
    nrep = nrep
  )

  message(sprintf(
    "geometry cache guardrail: median cache_on=%.3fs cache_off=%.3fs speedup=%.2fx slowdown=%.1f%%",
    perf$median_candidate,
    perf$median_baseline,
    perf$speedup,
    100 * perf$slowdown
  ))

  expect_lte(perf$median_candidate, perf$median_baseline * (1 + max_slowdown))
})

test_that("mvpa_model fast mode does not regress vs legacy baseline (guardrail)", {
  skip_on_cran()
  skip_if_not_perf_tests()

  set.seed(7301)
  mspec <- build_framework_perf_mvpa_spec()
  radius <- 2
  nrep <- perf_env_int("RMVPA_PERF_REP", 3L)
  max_slowdown <- perf_env_numeric("RMVPA_PHASE2_MVPA_MAX_SLOWDOWN", 0.20)

  perf <- benchmark_pair(
    run_baseline = function() run_with_perf_mode(mspec, radius = radius, mode = "legacy"),
    run_candidate = function() run_with_perf_mode(mspec, radius = radius, mode = "fast"),
    nrep = nrep
  )

  message(sprintf(
    "mvpa_model fast-mode guardrail: median fast=%.3fs legacy=%.3fs speedup=%.2fx slowdown=%.1f%%",
    perf$median_candidate,
    perf$median_baseline,
    perf$speedup,
    100 * perf$slowdown
  ))

  expect_lte(perf$median_candidate, perf$median_baseline * (1 + max_slowdown))
})

test_that("rsa_model fast mode does not regress vs legacy baseline (guardrail)", {
  skip_on_cran()
  skip_if_not_perf_tests()

  set.seed(7302)
  mspec <- build_framework_perf_rsa_spec()
  radius <- 2
  nrep <- perf_env_int("RMVPA_PERF_REP", 3L)
  max_slowdown <- perf_env_numeric("RMVPA_PHASE2_RSA_MAX_SLOWDOWN", 0.25)

  perf <- benchmark_pair(
    run_baseline = function() run_with_perf_mode(mspec, radius = radius, mode = "legacy"),
    run_candidate = function() run_with_perf_mode(mspec, radius = radius, mode = "fast"),
    nrep = nrep
  )

  message(sprintf(
    "rsa_model fast-mode guardrail: median fast=%.3fs legacy=%.3fs speedup=%.2fx slowdown=%.1f%%",
    perf$median_candidate,
    perf$median_baseline,
    perf$speedup,
    100 * perf$slowdown
  ))

  expect_lte(perf$median_candidate, perf$median_baseline * (1 + max_slowdown))
})

test_that("naive_xdec fast mode does not regress vs legacy baseline (guardrail)", {
  skip_on_cran()
  skip_if_not_perf_tests()

  set.seed(7303)
  mspec <- build_framework_perf_naive_xdec_spec()
  radius <- 2
  nrep <- perf_env_int("RMVPA_PERF_REP", 2L)
  max_slowdown <- perf_env_numeric("RMVPA_PHASE2_NAIVE_XDEC_MAX_SLOWDOWN", 0.30)

  perf <- benchmark_pair(
    run_baseline = function() run_with_perf_mode(mspec, radius = radius, mode = "legacy"),
    run_candidate = function() run_with_perf_mode(mspec, radius = radius, mode = "fast"),
    nrep = nrep
  )

  message(sprintf(
    "naive_xdec fast-mode guardrail: median fast=%.3fs legacy=%.3fs speedup=%.2fx slowdown=%.1f%%",
    perf$median_candidate,
    perf$median_baseline,
    perf$speedup,
    100 * perf$slowdown
  ))

  expect_lte(perf$median_candidate, perf$median_baseline * (1 + max_slowdown))
})

test_that("mvpa_model fold cache path scales on nightly workload", {
  skip_on_cran()
  skip_if_not_nightly_perf_tests()

  set.seed(7201)
  mspec <- build_framework_perf_mvpa_spec(D = c(8, 8, 8), nobs = 72, blocks = 6)
  radius <- 3
  nrep <- perf_env_int("RMVPA_NIGHTLY_PERF_REP", 2L)
  max_slowdown <- perf_env_numeric("RMVPA_NIGHTLY_FRAMEWORK_PERF_MAX_SLOWDOWN", 0.20)

  perf <- benchmark_pair(
    run_baseline = function() run_with_fold_cache(mspec, radius = radius, fold_cache_enabled = FALSE),
    run_candidate = function() run_with_fold_cache(mspec, radius = radius, fold_cache_enabled = TRUE),
    nrep = nrep
  )

  message(sprintf(
    "mvpa_model nightly perf: median cached=%.3fs, uncached=%.3fs, speedup=%.2fx, slowdown=%.1f%%",
    perf$median_candidate,
    perf$median_baseline,
    perf$speedup,
    100 * perf$slowdown
  ))

  expect_lte(perf$median_candidate, perf$median_baseline * (1 + max_slowdown))
})

test_that("mvpa_model fast mode stays within nightly slowdown guardrail", {
  skip_on_cran()
  skip_if_not_nightly_perf_tests()

  set.seed(7401)
  mspec <- build_framework_perf_mvpa_spec(D = c(8, 8, 8), nobs = 72, blocks = 6)
  radius <- 3
  nrep <- perf_env_int("RMVPA_NIGHTLY_PERF_REP", 2L)
  max_slowdown <- perf_env_numeric("RMVPA_NIGHTLY_PHASE2_MVPA_MAX_SLOWDOWN", 0.20)

  perf <- benchmark_pair(
    run_baseline = function() run_with_perf_mode(mspec, radius = radius, mode = "legacy"),
    run_candidate = function() run_with_perf_mode(mspec, radius = radius, mode = "fast"),
    nrep = nrep
  )

  message(sprintf(
    "mvpa_model fast-mode nightly perf: median fast=%.3fs legacy=%.3fs speedup=%.2fx slowdown=%.1f%%",
    perf$median_candidate,
    perf$median_baseline,
    perf$speedup,
    100 * perf$slowdown
  ))

  expect_lte(perf$median_candidate, perf$median_baseline * (1 + max_slowdown))
})

test_that("rsa_model fast mode stays within nightly slowdown guardrail", {
  skip_on_cran()
  skip_if_not_nightly_perf_tests()

  set.seed(7402)
  mspec <- build_framework_perf_rsa_spec(D = c(8, 8, 8), nobs = 54, blocks = 6)
  radius <- 3
  nrep <- perf_env_int("RMVPA_NIGHTLY_PERF_REP", 2L)
  max_slowdown <- perf_env_numeric("RMVPA_NIGHTLY_PHASE2_RSA_MAX_SLOWDOWN", 0.25)

  perf <- benchmark_pair(
    run_baseline = function() run_with_perf_mode(mspec, radius = radius, mode = "legacy"),
    run_candidate = function() run_with_perf_mode(mspec, radius = radius, mode = "fast"),
    nrep = nrep
  )

  message(sprintf(
    "rsa_model fast-mode nightly perf: median fast=%.3fs legacy=%.3fs speedup=%.2fx slowdown=%.1f%%",
    perf$median_candidate,
    perf$median_baseline,
    perf$speedup,
    100 * perf$slowdown
  ))

  expect_lte(perf$median_candidate, perf$median_baseline * (1 + max_slowdown))
})

test_that("rsa_model path stays within nightly slowdown guardrail", {
  skip_on_cran()
  skip_if_not_nightly_perf_tests()

  set.seed(7202)
  mspec <- build_framework_perf_rsa_spec(D = c(8, 8, 8), nobs = 54, blocks = 6)
  radius <- 3
  nrep <- perf_env_int("RMVPA_NIGHTLY_PERF_REP", 2L)
  max_slowdown <- perf_env_numeric("RMVPA_NIGHTLY_FRAMEWORK_RSA_PERF_MAX_SLOWDOWN", 0.25)

  perf <- benchmark_pair(
    run_baseline = function() run_with_fold_cache(mspec, radius = radius, fold_cache_enabled = FALSE),
    run_candidate = function() run_with_fold_cache(mspec, radius = radius, fold_cache_enabled = TRUE),
    nrep = nrep
  )

  message(sprintf(
    "rsa_model nightly perf: median cache_on=%.3fs, cache_off=%.3fs, speedup=%.2fx, slowdown=%.1f%%",
    perf$median_candidate,
    perf$median_baseline,
    perf$speedup,
    100 * perf$slowdown
  ))

  expect_lte(perf$median_candidate, perf$median_baseline * (1 + max_slowdown))
})
