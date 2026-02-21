testthat::skip_if_not_installed("neuroim2")

swift_perf_env_int <- function(name, default) {
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

swift_perf_env_numeric <- function(name, default) {
  raw <- Sys.getenv(name, "")
  if (!nzchar(raw)) {
    return(as.numeric(default))
  }
  val <- suppressWarnings(as.numeric(raw))
  if (!is.finite(val) || is.na(val) || val <= 0) {
    as.numeric(default)
  } else {
    as.numeric(val)
  }
}

swift_benchmark_pair <- function(run_baseline, run_candidate, nrep = 2L, warmup = TRUE) {
  if (isTRUE(warmup)) {
    invisible(run_baseline())
    invisible(run_candidate())
  }

  baseline <- numeric(nrep)
  candidate <- numeric(nrep)

  for (i in seq_len(nrep)) {
    order <- if ((i %% 2L) == 1L) {
      c("baseline", "candidate")
    } else {
      c("candidate", "baseline")
    }

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
    speedup = med_base / med_cand
  )
}

build_swift_perf_mspec <- function(D = c(8, 8, 8), nobs = 90, blocks = 6) {
  ds <- gen_sample_dataset(D = D, nobs = nobs, nlevels = 3, blocks = blocks)

  classes <- letters[1:3]
  stopifnot(nobs %% blocks == 0L)
  per_block <- nobs %/% blocks
  stopifnot(per_block %% length(classes) == 0L)

  y_block <- rep(classes, each = per_block %/% length(classes))
  y <- factor(rep(y_block, times = blocks), levels = classes)
  block <- rep(seq_len(blocks), each = per_block)

  des <- mvpa_design(
    data.frame(y = y, block = block),
    y_train = ~ y,
    block_var = ~ block
  )

  cval <- blocked_cross_validation(des$block_var)
  mvpa_model(
    model = load_model("sda_notune"),
    dataset = ds$dataset,
    design = des,
    model_type = "classification",
    crossval = cval
  )
}

run_swift_engine <- function(mspec, radius) {
  old_opt <- options(
    rMVPA.searchlight_mode = "fast",
    rMVPA.searchlight_backend_default = NULL,
    rMVPA.swift_searchlight = TRUE
  )
  on.exit(options(old_opt), add = TRUE)

  if (requireNamespace("future", quietly = TRUE)) {
    old_plan <- future::plan()
    on.exit(future::plan(old_plan), add = TRUE)
    future::plan(future::sequential)
  }

  run_searchlight(
    mspec,
    radius = radius,
    method = "standard",
    backend = "default",
    engine = "swift"
  )
}

run_legacy_engine <- function(mspec, radius) {
  old_opt <- options(
    rMVPA.searchlight_mode = "fast",
    rMVPA.searchlight_backend_default = NULL,
    rMVPA.swift_searchlight = TRUE
  )
  on.exit(options(old_opt), add = TRUE)

  if (requireNamespace("future", quietly = TRUE)) {
    old_plan <- future::plan()
    on.exit(future::plan(old_plan), add = TRUE)
    future::plan(future::sequential)
  }

  run_searchlight(
    mspec,
    radius = radius,
    method = "standard",
    backend = "default",
    engine = "legacy"
  )
}

test_that("swift searchlight engine is faster than legacy iterator (guardrail)", {
  skip_on_cran()
  skip_if_not_perf_tests()

  set.seed(9201)
  mspec <- build_swift_perf_mspec(D = c(5, 5, 5), nobs = 36, blocks = 6)
  radius <- 2
  nrep <- swift_perf_env_int("RMVPA_PERF_REP", 2L)
  min_speedup <- swift_perf_env_numeric("RMVPA_SWIFT_MIN_SPEEDUP", 1.00)

  perf <- swift_benchmark_pair(
    run_baseline = function() run_legacy_engine(mspec, radius = radius),
    run_candidate = function() run_swift_engine(mspec, radius = radius),
    nrep = nrep,
    warmup = nrep > 1L
  )

  message(sprintf(
    "swift perf guardrail: median swift=%.3fs, legacy=%.3fs, speedup=%.2fx",
    perf$median_candidate,
    perf$median_baseline,
    perf$speedup
  ))

  expect_gte(perf$speedup, min_speedup)
})

test_that("swift searchlight engine remains faster on nightly workload", {
  skip_on_cran()
  skip_if_not_nightly_perf_tests()

  set.seed(9202)
  mspec <- build_swift_perf_mspec(D = c(8, 8, 8), nobs = 90, blocks = 6)
  radius <- 3
  nrep <- swift_perf_env_int("RMVPA_NIGHTLY_PERF_REP", 2L)
  min_speedup <- swift_perf_env_numeric("RMVPA_NIGHTLY_SWIFT_MIN_SPEEDUP", 1.05)

  perf <- swift_benchmark_pair(
    run_baseline = function() run_legacy_engine(mspec, radius = radius),
    run_candidate = function() run_swift_engine(mspec, radius = radius),
    nrep = nrep,
    warmup = TRUE
  )

  message(sprintf(
    "swift nightly perf: median swift=%.3fs, legacy=%.3fs, speedup=%.2fx",
    perf$median_candidate,
    perf$median_baseline,
    perf$speedup
  ))

  expect_gte(perf$speedup, min_speedup)
})
