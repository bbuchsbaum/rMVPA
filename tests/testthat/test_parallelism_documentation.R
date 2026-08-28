test_that("parallelism vignette separates the three runtime layers", {
  vignette_path <- testthat::test_path("..", "..", "vignettes", "Parallelism.Rmd")
  skip_if_not(file.exists(vignette_path),
              "source vignette is not installed in the package tarball")
  text <- paste(readLines(vignette_path, warn = FALSE), collapse = "\n")

  expect_match(text, "VignetteIndexEntry\\{Choose a Safe Parallel Runtime for rMVPA\\}")
  expect_match(text, "Process topology", fixed = TRUE)
  expect_match(text, "rMVPA data-transport")
  expect_match(text, "OpenMP and BLAS")
  expect_match(text, "does **not** export the complete rMVPA dataset",
               fixed = TRUE)
  expect_match(text, "Copy-on-write is a memory\\s+mechanism")
  expect_match(text, "max_abs_result_error = 0", fixed = TRUE)
  expect_match(text, "123-fold reduction", fixed = TRUE)
  expect_false(grepl("rMVPA:::", text, fixed = TRUE))
})

test_that("parallel runtime receipt proves parity and smaller shard frames", {
  receipt_path <- testthat::test_path(
    "..", "..", "inst", "extdata", "parallel_runtime_benchmark_results.csv"
  )
  if (!file.exists(receipt_path)) {
    receipt_path <- system.file(
      "extdata", "parallel_runtime_benchmark_results.csv", package = "rMVPA"
    )
  }
  expect_true(file.exists(receipt_path))
  receipt <- utils::read.csv(receipt_path, stringsAsFactors = FALSE)

  required <- c(
    "future_backend", "data_backend", "repetitions", "requested_workers",
    "configured_workers", "omp_threads", "blas_threads", "batch_size",
    "n_batches", "analysis_median_seconds", "analysis_min_seconds",
    "analysis_max_seconds", "frame_bytes_sum_median",
    "frame_bytes_max_median", "max_abs_result_error", "fixture_data_bytes",
    "package_version", "load_mode", "future_version", "furrr_version",
    "shard_version", "source_fingerprint",
    "default_to_shard_frame_sum_ratio", "rss_interpretation"
  )
  expect_true(all(required %in% names(receipt)))
  expect_equal(nrow(receipt), 8L)
  expect_setequal(
    receipt$future_backend,
    rep(c("sequential", "multisession", "multicore", "mirai_multisession"),
        each = 2L)
  )
  expect_setequal(receipt$data_backend, rep(c("default", "shard"), 4L))
  expect_true(all(receipt$repetitions == 3L))
  expect_true(all(receipt$n_batches == 5L))
  expect_true(all(receipt$omp_threads == 1L))
  expect_true(all(receipt$blas_threads == 1L))
  expect_true(all(receipt$load_mode == "installed"))
  expect_lte(max(receipt$max_abs_result_error), 1e-12)
  expect_true(all(is.finite(receipt$analysis_median_seconds)))
  expect_true(all(receipt$analysis_min_seconds <= receipt$analysis_median_seconds))
  expect_true(all(receipt$analysis_median_seconds <= receipt$analysis_max_seconds))

  defaults <- receipt[receipt$data_backend == "default", ]
  shards <- receipt[receipt$data_backend == "shard", ]
  defaults <- defaults[order(defaults$future_backend), ]
  shards <- shards[order(shards$future_backend), ]
  expect_identical(defaults$future_backend, shards$future_backend)
  expect_true(all(shards$frame_bytes_sum_median < defaults$frame_bytes_sum_median))
  expect_true(all(shards$frame_bytes_max_median < defaults$frame_bytes_max_median))
  expect_true(all(shards$default_to_shard_frame_sum_ratio > 100))
  expect_true(all(grepl("not published", receipt$rss_interpretation,
                        fixed = TRUE)))
})

test_that("benchmark receipt fingerprints its executable evidence", {
  script_path <- testthat::test_path(
    "..", "..", "inst", "benchmarks", "parallel_runtime_benchmark.R"
  )
  if (!file.exists(script_path)) {
    script_path <- system.file(
      "benchmarks", "parallel_runtime_benchmark.R", package = "rMVPA"
    )
  }
  expect_true(file.exists(script_path))
  env <- new.env(parent = globalenv())
  sys.source(script_path, envir = env)

  expect_true(is.function(env$parallel_benchmark_validate))
  expect_true(is.function(env$parallel_benchmark_summarise))
  expect_true(is.function(env$parallel_benchmark_run))

  receipt_path <- testthat::test_path(
    "..", "..", "inst", "extdata", "parallel_runtime_benchmark_results.csv"
  )
  if (!file.exists(receipt_path)) {
    receipt_path <- system.file(
      "extdata", "parallel_runtime_benchmark_results.csv", package = "rMVPA"
    )
  }
  receipt <- utils::read.csv(receipt_path, stringsAsFactors = FALSE)
  expected <- unique(receipt$source_fingerprint)
  expect_length(expected, 1L)
  expect_true(!is.na(expected) && nzchar(expected))

  source_root <- testthat::test_path("..", "..")
  if (file.exists(file.path(source_root, "scripts",
                            "sweep_parallel_runtime_grid.R"))) {
    expect_identical(
      env$parallel_benchmark_source_fingerprint(source_root),
      expected
    )
  }
})

test_that("parallelism vignette is published in the site navigation", {
  pkgdown_path <- testthat::test_path("..", "..", "_pkgdown.yml")
  skip_if_not(file.exists(pkgdown_path),
              "source pkgdown configuration is excluded from the tarball")
  pkgdown <- paste(
    readLines(pkgdown_path, warn = FALSE),
    collapse = "\n"
  )
  expect_match(pkgdown, "title: Performance and Deployment", fixed = TRUE)
  expect_match(pkgdown, "- Parallelism", fixed = TRUE)
})
