testthat::skip_if_not_installed("neuroim2")

test_that("fold cache is enabled by default in fast mode", {
  old_opt <- options(
    rMVPA.searchlight_mode = "fast",
    rMVPA.fold_cache_enabled = NULL
  )
  on.exit(options(old_opt), add = TRUE)

  expect_true(rMVPA:::.fold_cache_enabled())
})

build_fold_cache_mspec <- function(crossval = NULL) {
  ds <- gen_sample_dataset(D = c(4, 4, 4), nobs = 24, nlevels = 2, blocks = 4)
  cv <- if (is.null(crossval)) blocked_cross_validation(ds$design$block_var) else crossval
  mspec <- mvpa_model(
    model = load_model("sda_notune"),
    dataset = ds$dataset,
    design = ds$design,
    model_type = "classification",
    crossval = cv
  )
  list(mspec = mspec, dataset = ds)
}

test_that("fold cache builder reproduces blocked CV train/test indices", {
  built <- build_fold_cache_mspec()
  mspec <- built$mspec

  old_opt <- options(
    rMVPA.searchlight_mode = "fast",
    rMVPA.fold_cache_enabled = NULL
  )
  on.exit(options(old_opt), add = TRUE)

  cache <- rMVPA:::.build_cv_fold_cache(mspec)
  expect_type(cache, "list")
  expect_true(length(cache$train_idx) > 0)
  expect_equal(length(cache$train_idx), length(cache$test_idx))

  y <- y_train(mspec)
  baseline <- crossval_samples(mspec$crossval, tibble::tibble(.row = seq_along(y)), y)
  expect_equal(length(cache$train_idx), nrow(baseline))

  for (i in seq_len(nrow(baseline))) {
    expect_equal(as.integer(cache$train_idx[[i]]), as.integer(baseline$train[[i]]$idx))
    expect_equal(as.integer(cache$test_idx[[i]]), as.integer(baseline$test[[i]]$idx))
    expect_equal(cache$ytrain[[i]], baseline$ytrain[[i]])
    expect_equal(cache$ytest[[i]], baseline$ytest[[i]])
  }
})

test_that("fold cache supports custom_cross_validation sample sets", {
  built <- build_fold_cache_mspec()
  n <- length(y_train(built$mspec))
  sample_set <- list(
    list(train = c(1:12, 19:24), test = 13:18),
    list(train = c(7:24), test = 1:6)
  )
  cv <- custom_cross_validation(sample_set)
  built_custom <- build_fold_cache_mspec(crossval = cv)

  old_opt <- options(
    rMVPA.searchlight_mode = "fast",
    rMVPA.fold_cache_enabled = NULL
  )
  on.exit(options(old_opt), add = TRUE)

  cache <- rMVPA:::.build_cv_fold_cache(built_custom$mspec)
  expect_equal(cache$n_obs, n)
  expect_equal(length(cache$train_idx), length(sample_set))
  expect_equal(cache$train_idx, lapply(sample_set, function(x) as.integer(x$train)))
  expect_equal(cache$test_idx, lapply(sample_set, function(x) as.integer(x$test)))
})

test_that("fold cache is disabled for nondeterministic CV specifications", {
  built <- build_fold_cache_mspec(
    crossval = twofold_blocked_cross_validation(block_var = rep(1:4, each = 6), nreps = 3)
  )

  old_opt <- options(
    rMVPA.searchlight_mode = "fast",
    rMVPA.fold_cache_enabled = NULL
  )
  on.exit(options(old_opt), add = TRUE)

  cache <- rMVPA:::.build_cv_fold_cache(built$mspec)
  expect_null(cache)
})

test_that("fold cache resolve invalidates when labels or dimensions change", {
  built <- build_fold_cache_mspec()
  mspec <- built$mspec

  old_opt <- options(
    rMVPA.searchlight_mode = "fast",
    rMVPA.fold_cache_enabled = NULL
  )
  on.exit(options(old_opt), add = TRUE)

  cache <- rMVPA:::.build_cv_fold_cache(mspec)
  mspec_cached <- mspec
  mspec_cached$.cv_fold_cache <- cache

  n_obs <- length(y_train(mspec_cached))
  expect_type(rMVPA:::.resolve_cv_fold_cache(mspec_cached, n_obs = n_obs), "list")
  expect_null(rMVPA:::.resolve_cv_fold_cache(mspec_cached, n_obs = n_obs - 1L))

  mspec_changed <- mspec_cached
  mspec_changed$design$cv_labels <- rev(mspec_changed$design$cv_labels)
  expect_null(rMVPA:::.resolve_cv_fold_cache(mspec_changed, n_obs = n_obs))
})

test_that("internal_crossval skips crossval_samples when valid fold cache is attached", {
  built <- build_fold_cache_mspec()
  mspec <- built$mspec
  ds <- built$dataset

  vox <- sample(which(ds$dataset$mask > 0), 24)
  roi <- as_roi(data_sample(ds$dataset, vox), ds$dataset)

  old_opt <- options(
    rMVPA.searchlight_mode = "fast",
    rMVPA.fold_cache_enabled = NULL
  )
  on.exit(options(old_opt), add = TRUE)

  original_crossval_samples <- get("crossval_samples", envir = asNamespace("rMVPA"))
  call_count <- 0L
  counting_crossval_samples <- function(...) {
    call_count <<- call_count + 1L
    original_crossval_samples(...)
  }

  mspec_uncached <- mspec
  res_uncached <- testthat::with_mocked_bindings(
    rMVPA:::internal_crossval(mspec_uncached, roi, id = 1),
    crossval_samples = counting_crossval_samples,
    .package = "rMVPA"
  )
  expect_gt(call_count, 0L)

  call_count <- 0L
  mspec_cached <- mspec
  mspec_cached$.cv_fold_cache <- rMVPA:::.build_cv_fold_cache(mspec_cached)
  res_cached <- testthat::with_mocked_bindings(
    rMVPA:::internal_crossval(mspec_cached, roi, id = 1),
    crossval_samples = counting_crossval_samples,
    .package = "rMVPA"
  )
  expect_equal(call_count, 0L)

  expect_equal(res_uncached$error, res_cached$error)
  expect_equal(res_uncached$performance[[1]], res_cached$performance[[1]], tolerance = 1e-10)
})

test_that("searchlight fold cache preserves output maps when enabled", {
  built <- build_fold_cache_mspec()
  mspec <- built$mspec

  set.seed(2201)
  old_off <- options(
    rMVPA.searchlight_mode = "legacy",
    rMVPA.fold_cache_enabled = NULL
  )
  on.exit(options(old_off), add = TRUE)
  res_base <- run_searchlight(mspec, radius = 2, method = "standard", backend = "default")

  set.seed(2201)
  old_on <- options(
    rMVPA.searchlight_mode = "fast",
    rMVPA.fold_cache_enabled = NULL
  )
  on.exit(options(old_on), add = TRUE)
  res_cache <- run_searchlight(mspec, radius = 2, method = "standard", backend = "default")

  expect_searchlight_parity(
    reference = res_base,
    candidate = res_cache,
    atol = 1e-10,
    rtol = 1e-8
  )
})
