testthat::skip_if_not_installed("neuroim2")

.era_perf_env_int <- function(name, default) {
  value <- suppressWarnings(as.integer(Sys.getenv(name, "")))
  if (length(value) != 1L || is.na(value) || value < 1L) {
    as.integer(default)
  } else {
    value
  }
}

.era_reference_specificity <- function(E, R) {
  S <- rMVPA:::.era_cross_similarity(
    E, R, method = "pearson", min_voxels = 2L
  )
  diag(S) - vapply(
    seq_len(nrow(S)),
    function(j) mean(S[-j, j]),
    numeric(1L)
  )
}

.era_benchmark <- function(fun, nrep, inner) {
  invisible(fun())
  elapsed <- numeric(nrep)
  for (i in seq_len(nrep)) {
    gc(verbose = FALSE)
    elapsed[[i]] <- unname(system.time({
      for (j in seq_len(inner)) invisible(fun())
    })[["elapsed"]])
  }
  elapsed
}

test_that("era item-specificity fast kernel improves representative runtime", {
  skip_on_cran()
  skip_if_not_perf_tests()

  K <- .era_perf_env_int("RMVPA_ERA_PERF_ITEMS", 240L)
  P <- .era_perf_env_int("RMVPA_ERA_PERF_VOXELS", 80L)
  nrep <- .era_perf_env_int("RMVPA_PERF_REP", 3L)
  inner <- .era_perf_env_int("RMVPA_ERA_PERF_INNER", 5L)

  set.seed(20260812)
  E <- matrix(stats::rnorm(K * P), nrow = K)
  R <- 0.25 * E + matrix(stats::rnorm(K * P), nrow = K)

  baseline_fun <- function() .era_reference_specificity(E, R)
  candidate_fun <- function() {
    rMVPA:::.era_item_similarity_scores(
      E,
      R,
      method = "pearson",
      min_voxels = 2L,
      need_matrix = FALSE
    )$specificity
  }

  expect_equal(candidate_fun(), baseline_fun(), tolerance = 1e-12)

  # Alternate measurement order to reduce warm-cache/order bias.
  baseline <- candidate <- numeric(nrep)
  for (i in seq_len(nrep)) {
    if (i %% 2L == 1L) {
      baseline[[i]] <- .era_benchmark(baseline_fun, 1L, inner)
      candidate[[i]] <- .era_benchmark(candidate_fun, 1L, inner)
    } else {
      candidate[[i]] <- .era_benchmark(candidate_fun, 1L, inner)
      baseline[[i]] <- .era_benchmark(baseline_fun, 1L, inner)
    }
  }

  median_baseline <- stats::median(baseline)
  median_candidate <- stats::median(candidate)
  speedup <- median_baseline / median_candidate
  message(sprintf(
    paste0(
      "ERA item kernel guardrail: K=%d P=%d inner=%d ",
      "fast=%.3fs reference=%.3fs speedup=%.2fx"
    ),
    K, P, inner, median_candidate, median_baseline, speedup
  ))

  expect_lte(median_candidate, median_baseline * 1.20)
})
