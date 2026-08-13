testthat::skip_if_not_installed("neuroim2")

.era_engine_perf_int <- function(name, default) {
  value <- suppressWarnings(as.integer(Sys.getenv(name, "")))
  if (length(value) != 1L || is.na(value) || value < 1L) {
    as.integer(default)
  } else {
    value
  }
}

.era_engine_perf_fixture <- function(K = 40L, dims = c(5L, 5L, 5L)) {
  set.seed(20260812)
  keys <- sprintf("item_%03d", seq_len(K))
  p <- prod(dims)
  E <- matrix(stats::rnorm(K * p), nrow = K)
  R <- 0.3 * E + matrix(stats::rnorm(K * p), nrow = K)
  mask <- neuroim2::NeuroVol(array(1, dims), neuroim2::NeuroSpace(dims))
  enc <- neuroim2::NeuroVec(
    array(as.numeric(t(E)), c(dims, K)),
    neuroim2::NeuroSpace(c(dims, K))
  )
  ret <- neuroim2::NeuroVec(
    array(as.numeric(t(R)), c(dims, K)),
    neuroim2::NeuroSpace(c(dims, K))
  )
  tab <- data.frame(item = factor(keys, levels = keys))
  dataset <- mvpa_dataset(enc, test_data = ret, mask = mask)
  design <- mvpa_design(
    tab,
    test_design = tab,
    y_train = ~ item,
    y_test = ~ item
  )
  suppressWarnings(era_rsa_model(
    dataset,
    design,
    key_var = ~ item,
    pairing = "one_to_one",
    era_components = "item"
  ))
}

.era_engine_elapsed <- function(fun, nrep) {
  elapsed <- numeric(nrep)
  for (i in seq_len(nrep)) {
    gc(verbose = FALSE)
    elapsed[[i]] <- unname(system.time(invisible(fun()))[["elapsed"]])
  }
  elapsed
}

test_that("ERA-RSA direct searchlight engine improves representative runtime", {
  skip_on_cran()
  skip_if_not_perf_tests()

  K <- .era_engine_perf_int("RMVPA_ERA_ENGINE_PERF_ITEMS", 40L)
  side <- .era_engine_perf_int("RMVPA_ERA_ENGINE_PERF_SIDE", 5L)
  nrep <- .era_engine_perf_int("RMVPA_PERF_REP", 3L)
  model <- .era_engine_perf_fixture(K = K, dims = rep(side, 3L))

  legacy_fun <- function() run_searchlight(
    model,
    radius = 2,
    method = "standard",
    engine = "legacy",
    backend = "default",
    preflight = "off"
  )
  fast_fun <- function() run_searchlight(
    model,
    radius = 2,
    method = "standard",
    engine = "era_rsa_fast",
    backend = "default",
    preflight = "off"
  )

  legacy_reference <- legacy_fun()
  fast_reference <- fast_fun()
  expect_identical(fast_reference$metrics, legacy_reference$metrics)
  for (nm in legacy_reference$metrics) {
    expect_equal(
      as.numeric(neuroim2::values(fast_reference$results[[nm]])),
      as.numeric(neuroim2::values(legacy_reference$results[[nm]])),
      tolerance = 1e-12,
      info = nm
    )
  }

  legacy <- fast <- numeric(nrep)
  for (i in seq_len(nrep)) {
    if (i %% 2L == 1L) {
      legacy[[i]] <- .era_engine_elapsed(legacy_fun, 1L)
      fast[[i]] <- .era_engine_elapsed(fast_fun, 1L)
    } else {
      fast[[i]] <- .era_engine_elapsed(fast_fun, 1L)
      legacy[[i]] <- .era_engine_elapsed(legacy_fun, 1L)
    }
  }

  median_legacy <- stats::median(legacy)
  median_fast <- stats::median(fast)
  speedup <- median_legacy / median_fast
  message(sprintf(
    paste0(
      "ERA-RSA engine guardrail: K=%d side=%d fast=%.3fs ",
      "legacy=%.3fs speedup=%.2fx"
    ),
    K, side, median_fast, median_legacy, speedup
  ))

  # This is a broad regression bound rather than a brittle benchmark target.
  expect_lte(median_fast, median_legacy * 1.10)
})
