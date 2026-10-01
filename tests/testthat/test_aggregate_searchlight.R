# The exact sphere-aggregation engine must reproduce the general per-sphere
# path for corclass, and must decline anything outside its proven regime.
testthat::skip_if_not_installed("neuroim2")

agg_maps <- function(res) lapply(res$results, function(m) as.numeric(neuroim2::values(m)))

agg_spec <- function(seed, D = c(5, 5, 5), nobs = 48, nlevels = 3, blocks = 4,
                     class_metrics = FALSE, scale = 1, offset = 0, model = "corclass") {
  set.seed(seed)
  ds <- gen_sample_dataset(D = D, nobs = nobs, nlevels = nlevels, blocks = blocks)
  if (scale != 1 || offset != 0) {
    arr <- as.array(ds$dataset$train_data) * scale + offset
    ds$dataset$train_data <- neuroim2::NeuroVec(arr, neuroim2::space(ds$dataset$train_data))
  }
  mvpa_model(load_model(model), ds$dataset, ds$design, "classification",
             crossval = blocked_cross_validation(ds$design$block_var),
             class_metrics = class_metrics)
}

expect_engine_parity <- function(mspec, radius = 2, tol = 1e-10) {
  fast <- run_searchlight(mspec, radius = radius, engine = "aggregate_fast", backend = "default")
  ref <- run_searchlight(mspec, radius = radius, engine = "legacy", backend = "default")
  expect_identical(attr(fast, "searchlight_engine"), "aggregate_fast")
  expect_identical(sort(names(fast$results)), sort(names(ref$results)))
  a <- agg_maps(fast)[names(ref$results)]
  b <- agg_maps(ref)
  for (m in names(b)) {
    expect_identical(is.na(a[[m]]), is.na(b[[m]]), info = m)
    expect_equal(a[[m]], b[[m]], tolerance = tol, info = m)
  }
  invisible(fast)
}

test_that("aggregation engine matches the general path (multiclass)", {
  expect_engine_parity(agg_spec(901))
})

test_that("aggregation engine matches the general path (binary)", {
  expect_engine_parity(agg_spec(902, nlevels = 2, nobs = 40))
})

test_that("aggregation engine matches the general path with per-class AUC maps", {
  expect_engine_parity(agg_spec(903, class_metrics = TRUE))
})

test_that("aggregation engine matches on raw-BOLD-scale data (large mean, small spread)", {
  # Large offsets stress the aggregated centred sums; flagged centres are
  # recomputed exactly.
  expect_engine_parity(agg_spec(904, scale = 5, offset = 1500))
})

test_that("auto selects the aggregation engine only for eligible corclass specs", {
  ms <- agg_spec(905)
  expect_identical(rMVPA:::.resolve_searchlight_engine(ms, "standard", "auto"), "aggregate_fast")
  expect_identical(rMVPA:::.resolve_searchlight_engine(ms, "randomized", "auto"), "legacy")

  robust <- ms
  robust$tune_grid <- data.frame(method = "pearson", robust = TRUE)
  expect_false(rMVPA:::.is_aggregate_fast_path(robust, "standard"))

  spearman <- ms
  spearman$tune_grid <- data.frame(method = "spearman", robust = FALSE)
  expect_false(rMVPA:::.is_aggregate_fast_path(spearman, "standard"))
})

test_that("data outside the regime fall back quietly under auto and error when requested", {
  ms <- agg_spec(906)
  arr <- as.array(ms$dataset$train_data)
  arr[1, 1, 1, 3] <- NA
  ms$dataset$train_data <- neuroim2::NeuroVec(arr, neuroim2::space(ms$dataset$train_data))

  res <- expect_no_warning(run_searchlight(ms, radius = 2, backend = "default"))
  expect_identical(attr(res, "searchlight_engine"), "legacy")
  expect_error(run_searchlight(ms, radius = 2, engine = "aggregate_fast", backend = "default"),
               "missing or non-finite")
})

test_that("zapsmall emulation matches base::zapsmall on this R", {
  rule <- rMVPA:::.aggregate_zapsmall_rule()
  skip_if(isFALSE(rule), "no vectorised rule for this R version; per-block fallback is used")
  set.seed(907)
  for (scale in c(1e-4, 0.01, 0.12, 0.5, 0.99, 7)) {
    m <- matrix(stats::runif(40) * scale, 8)
    expect_identical(round(m, digits = rule(max(abs(m)), getOption("digits"))), base::zapsmall(m))
  }
})

test_that("aggregation engine matches the general path for naive_bayes", {
  expect_engine_parity(agg_spec(911, model = "naive_bayes"))
  expect_engine_parity(agg_spec(912, model = "naive_bayes", nlevels = 2, nobs = 40))
  expect_engine_parity(agg_spec(913, model = "naive_bayes", class_metrics = TRUE))
  expect_engine_parity(agg_spec(914, model = "naive_bayes", scale = 5, offset = 1500))
})

test_that("naive_bayes engine repairs near-tied centres exactly", {
  # Strong signal saturates probabilities, so one-vs-rest scores tie at the
  # last bit; flagged centres are recomputed in the general path's voxel order.
  set.seed(915)
  ds <- gen_sample_dataset(D = c(5, 5, 5), nobs = 48, nlevels = 3, blocks = 4)
  arr <- as.array(ds$dataset$train_data)
  y <- as.integer(ds$design$y_train)
  for (i in seq_along(y)) arr[, , , i] <- arr[, , , i] + 3 * y[i]
  ds$dataset$train_data <- neuroim2::NeuroVec(arr, neuroim2::space(ds$dataset$train_data))
  ms <- mvpa_model(load_model("naive_bayes"), ds$dataset, ds$design, "classification",
                   crossval = blocked_cross_validation(ds$design$block_var))
  fast <- expect_engine_parity(ms)
  expect_true(attr(fast, "aggregate_recomputed_centres") >= 0L)
})

test_that("naive_bayes engine declines zero within-class variance", {
  ms <- agg_spec(916, model = "naive_bayes")
  arr <- as.array(ms$dataset$train_data)
  y <- ms$design$y_train
  arr[2, 2, 2, y == levels(y)[1]] <- 4.2  # constant within one class, varies overall
  ms$dataset$train_data <- neuroim2::NeuroVec(arr, neuroim2::space(ms$dataset$train_data))
  expect_error(run_searchlight(ms, radius = 2, engine = "aggregate_fast", backend = "default"),
               "zero within-class variance")
  res <- suppressWarnings(run_searchlight(ms, radius = 2, backend = "default"))
  expect_identical(attr(res, "searchlight_engine"), "legacy")
})
