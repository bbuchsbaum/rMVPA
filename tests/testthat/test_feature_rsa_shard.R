library(testthat)
library(rMVPA)

# Skip entire file if shard is not installed
skip_if_not_installed("shard")

# Muffle environment-specific package version warnings emitted on some workers
muffle_worker_version_warnings <- function(expr) {
  withCallingHandlers(
    expr,
    warning = function(w) {
      if (grepl("was built under R version", conditionMessage(w), fixed = TRUE)) {
        invokeRestart("muffleWarning")
      }
    }
  )
}

# ---- Regional: shard backend via run_regional(..., backend = "shard") ----

test_that("feature_rsa_model regional with shard backend runs without error (PLS)", {
  skip_on_cran()
  set.seed(123)

  dset <- gen_sample_dataset(c(5, 5, 5), 80, blocks = 3)
  Fmat <- matrix(rnorm(80 * 12), 80, 12)
  fdes <- feature_rsa_design(F = Fmat, labels = paste0("obs", 1:80), max_comps = 5)
  mspec <- feature_rsa_model(dset$dataset, fdes, method = "pls",
                              crossval = blocked_cross_validation(dset$design$block_var))

  region_mask <- neuroim2::NeuroVol(
    sample(1:3, length(dset$dataset$mask), replace = TRUE),
    neuroim2::space(dset$dataset$mask)
  )

  res <- muffle_worker_version_warnings(
    run_regional(mspec, region_mask, backend = "shard")
  )

  expect_s3_class(res, "regional_mvpa_result")
  expect_true(is.data.frame(res$performance_table))
  expect_true(nrow(res$performance_table) > 0)
  expect_true("pattern_correlation" %in% colnames(res$performance_table))
  expect_true("rdm_correlation" %in% colnames(res$performance_table))
})

test_that("feature_rsa_model regional with shard backend runs without error (PCA)", {
  skip_on_cran()
  set.seed(42)

  dset <- gen_sample_dataset(c(5, 5, 5), 60, blocks = 3)
  Fmat <- matrix(rnorm(60 * 10), 60, 10)
  fdes <- feature_rsa_design(F = Fmat, labels = paste0("obs", 1:60), max_comps = 5)
  mspec <- feature_rsa_model(dset$dataset, fdes, method = "pca",
                              crossval = blocked_cross_validation(dset$design$block_var))

  region_mask <- neuroim2::NeuroVol(
    sample(1:3, length(dset$dataset$mask), replace = TRUE),
    neuroim2::space(dset$dataset$mask)
  )

  res <- muffle_worker_version_warnings(
    run_regional(mspec, region_mask, backend = "shard")
  )

  expect_s3_class(res, "regional_mvpa_result")
  expect_true(is.data.frame(res$performance_table))
  expect_true("pattern_correlation" %in% colnames(res$performance_table))
  expect_true("pattern_discrimination" %in% colnames(res$performance_table))
  expect_true("voxel_correlation" %in% colnames(res$performance_table))
  expect_true("mse" %in% colnames(res$performance_table))
  expect_true("r_squared" %in% colnames(res$performance_table))
})

# ---- Regional: shard vs default parity ----

test_that("feature_rsa_model regional shard parity with default backend", {
  skip_on_cran()
  set.seed(99)

  dset <- gen_sample_dataset(c(5, 5, 5), 60, blocks = 3)
  Fmat <- matrix(rnorm(60 * 8), 60, 8)
  fdes <- feature_rsa_design(F = Fmat, labels = paste0("obs", 1:60), max_comps = 5)
  mspec <- feature_rsa_model(dset$dataset, fdes, method = "pca",
                              crossval = blocked_cross_validation(dset$design$block_var))

  region_mask <- neuroim2::NeuroVol(
    sample(1:2, length(dset$dataset$mask), replace = TRUE),
    neuroim2::space(dset$dataset$mask)
  )

  res_default <- muffle_worker_version_warnings(
    run_regional(mspec, region_mask, backend = "default")
  )

  res_shard <- muffle_worker_version_warnings(
    run_regional(mspec, region_mask, backend = "shard")
  )

  expect_s3_class(res_default, "regional_mvpa_result")
  expect_s3_class(res_shard, "regional_mvpa_result")

  # Same number of ROIs processed
  expect_equal(nrow(res_default$performance_table), nrow(res_shard$performance_table))

  # Same metric columns
  expect_equal(
    sort(colnames(res_default$performance_table)),
    sort(colnames(res_shard$performance_table))
  )
})

# ---- Regional: shard with S-based design ----

test_that("feature_rsa_model with S-based design works with shard backend", {
  skip_on_cran()
  set.seed(77)

  dset <- gen_sample_dataset(c(5, 5, 5), 60, blocks = 3)
  obs_features <- matrix(rnorm(60 * 10), 60, 10)
  S <- tcrossprod(base::scale(obs_features))
  fdes <- feature_rsa_design(S = S, labels = paste0("obs", 1:60), k = 8)
  mspec <- feature_rsa_model(dset$dataset, fdes, method = "pls",
                              crossval = blocked_cross_validation(dset$design$block_var))

  region_mask <- neuroim2::NeuroVol(
    sample(1:3, length(dset$dataset$mask), replace = TRUE),
    neuroim2::space(dset$dataset$mask)
  )

  res <- muffle_worker_version_warnings(
    run_regional(mspec, region_mask, backend = "shard")
  )

  expect_s3_class(res, "regional_mvpa_result")
  expect_true(is.data.frame(res$performance_table))
  expect_true("rdm_correlation" %in% colnames(res$performance_table))
})

# ---- Regional: shard with permutation testing ----

test_that("feature_rsa_model with permutation testing works under shard backend", {
  skip_on_cran()
  set.seed(88)

  dset <- gen_sample_dataset(c(3, 3, 3), 40, blocks = 2)
  Fmat <- matrix(rnorm(40 * 6), 40, 6)
  fdes <- feature_rsa_design(F = Fmat, labels = paste0("obs", 1:40))
  mspec <- feature_rsa_model(
    dset$dataset, fdes, method = "pca",
    crossval = blocked_cross_validation(dset$design$block_var),
    nperm = 5,
    permute_by = "observations"
  )

  region_mask <- neuroim2::NeuroVol(
    rep(1L, length(dset$dataset$mask)),
    neuroim2::space(dset$dataset$mask)
  )

  res <- muffle_worker_version_warnings(
    run_regional(mspec, region_mask, backend = "shard")
  )

  expect_s3_class(res, "regional_mvpa_result")
  perf_cols <- colnames(res$performance_table)
  expect_true(any(grepl("^p_", perf_cols)))
  expect_true(any(grepl("^z_", perf_cols)))
  expect_true("p_pattern_correlation" %in% perf_cols)
  expect_true("z_pattern_correlation" %in% perf_cols)
})

# ---- Regional: shard with ncomp_selection variants ----

test_that("feature_rsa_model ncomp_selection='loo' works with shard backend", {
  skip_on_cran()
  set.seed(42)

  dset <- gen_sample_dataset(c(4, 4, 4), 60, blocks = 3)
  Fmat <- matrix(rnorm(60 * 8), 60, 8)
  fdes <- feature_rsa_design(F = Fmat, labels = paste0("o", 1:60), max_comps = 8)
  mspec <- feature_rsa_model(dset$dataset, fdes, method = "pls",
                              ncomp_selection = "loo",
                              crossval = blocked_cross_validation(dset$design$block_var))

  region_mask <- neuroim2::NeuroVol(
    sample(1:2, length(dset$dataset$mask), replace = TRUE),
    neuroim2::space(dset$dataset$mask)
  )

  res <- muffle_worker_version_warnings(
    run_regional(mspec, region_mask, backend = "shard")
  )

  expect_s3_class(res, "regional_mvpa_result")
  expect_true("ncomp" %in% colnames(res$performance_table))
  ncomps <- res$performance_table$ncomp
  expect_true(all(ncomps >= 1 & ncomps <= 8))
})

# ---- Shard with use_shard() + mvpa_iterate (low-level) ----

test_that("feature_rsa_model shard via use_shard + mvpa_iterate produces valid results", {
  skip_on_cran()
  set.seed(42)

  dset <- gen_sample_dataset(c(5, 5, 5), 60, blocks = 3)
  Fmat <- matrix(rnorm(60 * 10), 60, 10)
  fdes <- feature_rsa_design(F = Fmat, labels = paste0("obs", 1:60), max_comps = 5)
  mspec <- feature_rsa_model(dset$dataset, fdes, method = "pca",
                              crossval = blocked_cross_validation(dset$design$block_var))

  # Tag with shard
  mspec_shard <- use_shard(mspec)
  expect_true(inherits(mspec_shard, "shard_model_spec"))
  expect_false(is.null(mspec_shard$shard_data))
  expect_equal(mspec_shard$shard_data$roi_type, "volumetric")

  # Get a few searchlight neighborhoods
 sl <- get_searchlight(dset$dataset, radius = 3)
  vox_iter <- lapply(sl, function(x) x)
  vox_iter <- vox_iter[1:min(5, length(vox_iter))]

  res <- muffle_worker_version_warnings(
    mvpa_iterate(mspec_shard, vox_iter, ids = seq_along(vox_iter))
  )

  expect_equal(nrow(res), length(vox_iter))
  # feature_rsa may legitimately error on tiny ROIs; check structure
  expect_true("id" %in% colnames(res))
  expect_true("error" %in% colnames(res))
  expect_true("performance" %in% colnames(res))
})

# ---- Shard parity at mvpa_iterate level ----

test_that("feature_rsa_model shard vs default parity at mvpa_iterate level", {
  skip_on_cran()
  set.seed(42)

  dset <- gen_sample_dataset(c(5, 5, 5), 60, blocks = 3)
  Fmat <- matrix(rnorm(60 * 10), 60, 10)
  fdes <- feature_rsa_design(F = Fmat, labels = paste0("obs", 1:60), max_comps = 5)
  mspec <- feature_rsa_model(dset$dataset, fdes, method = "pca",
                              crossval = blocked_cross_validation(dset$design$block_var))

  sl <- get_searchlight(dset$dataset, radius = 3)
  vox_iter <- lapply(sl, function(x) x)
  vox_iter <- vox_iter[1:min(5, length(vox_iter))]

  # Default path
  set.seed(42)
  res_default <- muffle_worker_version_warnings(
    mvpa_iterate(mspec, vox_iter, ids = seq_along(vox_iter))
  )

  # Shard path
  mspec_shard <- use_shard(mspec)
  set.seed(42)
  res_shard <- muffle_worker_version_warnings(
    mvpa_iterate(mspec_shard, vox_iter, ids = seq_along(vox_iter))
  )

  expect_equal(nrow(res_default), nrow(res_shard))
  expect_equal(sort(res_default$id), sort(res_shard$id))
  expect_equal(res_default$error, res_shard$error)

  # Performance parity for non-error rows
  if (!all(res_default$error)) {
    ok <- !res_default$error
    perf_default <- lapply(res_default$performance[ok], function(p) {
      if (is.null(p)) NULL else unlist(p)
    })
    perf_shard <- lapply(res_shard$performance[ok], function(p) {
      if (is.null(p)) NULL else unlist(p)
    })
    for (i in seq_along(perf_default)) {
      if (!is.null(perf_default[[i]]) && !is.null(perf_shard[[i]])) {
        expect_equal(perf_default[[i]], perf_shard[[i]], tolerance = 1e-8)
      }
    }
  }
})

# ---- Shard with multisession workers ----

test_that("feature_rsa_model shard works across multisession workers", {
  skip_on_cran()
  skip_if_not_installed("future")
  set.seed(42)

  dset <- gen_sample_dataset(c(5, 5, 5), 60, blocks = 3)
  Fmat <- matrix(rnorm(60 * 8), 60, 8)
  fdes <- feature_rsa_design(F = Fmat, labels = paste0("obs", 1:60), max_comps = 5)
  mspec <- feature_rsa_model(dset$dataset, fdes, method = "pca",
                              crossval = blocked_cross_validation(dset$design$block_var))

  region_mask <- neuroim2::NeuroVol(
    sample(1:2, length(dset$dataset$mask), replace = TRUE),
    neuroim2::space(dset$dataset$mask)
  )

  old_plan <- future::plan(future::multisession, workers = 2)
  on.exit(future::plan(old_plan), add = TRUE)

  res <- muffle_worker_version_warnings(
    run_regional(mspec, region_mask, backend = "shard")
  )

  expect_s3_class(res, "regional_mvpa_result")
  expect_true(is.data.frame(res$performance_table))
  expect_true(nrow(res$performance_table) > 0)
})

# ---- backend = "auto" works for feature_rsa_model ----

test_that("feature_rsa_model with backend='auto' works", {
  skip_on_cran()
  set.seed(42)

  dset <- gen_sample_dataset(c(5, 5, 5), 60, blocks = 3)
  Fmat <- matrix(rnorm(60 * 8), 60, 8)
  fdes <- feature_rsa_design(F = Fmat, labels = paste0("obs", 1:60), max_comps = 5)
  mspec <- feature_rsa_model(dset$dataset, fdes, method = "pca",
                              crossval = blocked_cross_validation(dset$design$block_var))

  region_mask <- neuroim2::NeuroVol(
    sample(1:2, length(dset$dataset$mask), replace = TRUE),
    neuroim2::space(dset$dataset$mask)
  )

  res <- muffle_worker_version_warnings(
    run_regional(mspec, region_mask, backend = "auto")
  )

  expect_s3_class(res, "regional_mvpa_result")
  expect_true(is.data.frame(res$performance_table))
})

# ---- Shard with multicore workers (typical HPC scenario) ----

test_that("feature_rsa_model shard works across multicore workers", {
  skip_on_cran()
  skip_if_not_installed("future")
  skip_if_not(isTRUE(future::supportsMulticore()),
              "future::multicore not supported in this environment")
  set.seed(42)

  dset <- gen_sample_dataset(c(5, 5, 5), 60, blocks = 3)
  Fmat <- matrix(rnorm(60 * 8), 60, 8)
  fdes <- feature_rsa_design(F = Fmat, labels = paste0("obs", 1:60), max_comps = 5)
  mspec <- feature_rsa_model(dset$dataset, fdes, method = "pca",
                              crossval = blocked_cross_validation(dset$design$block_var))

  region_mask <- neuroim2::NeuroVol(
    sample(1:3, length(dset$dataset$mask), replace = TRUE),
    neuroim2::space(dset$dataset$mask)
  )

  old_plan <- future::plan(future::multicore, workers = 2)
  on.exit(future::plan(old_plan), add = TRUE)

  res <- muffle_worker_version_warnings(
    run_regional(mspec, region_mask, backend = "shard")
  )

  expect_s3_class(res, "regional_mvpa_result")
  expect_true(is.data.frame(res$performance_table))
  expect_true(nrow(res$performance_table) > 0)
  expect_true("pattern_correlation" %in% colnames(res$performance_table))
  expect_true("rdm_correlation" %in% colnames(res$performance_table))
})

test_that("feature_rsa_model shard multicore parity with sequential", {
  skip_on_cran()
  skip_if_not_installed("future")
  skip_if_not(isTRUE(future::supportsMulticore()),
              "future::multicore not supported in this environment")
  set.seed(42)

  dset <- gen_sample_dataset(c(5, 5, 5), 60, blocks = 3)
  Fmat <- matrix(rnorm(60 * 8), 60, 8)
  fdes <- feature_rsa_design(F = Fmat, labels = paste0("obs", 1:60), max_comps = 5)
  mspec <- feature_rsa_model(dset$dataset, fdes, method = "pls",
                              crossval = blocked_cross_validation(dset$design$block_var))

  sl <- get_searchlight(dset$dataset, radius = 3)
  vox_iter <- lapply(sl, function(x) x)
  vox_iter <- vox_iter[1:min(5, length(vox_iter))]

  # Sequential baseline with shard
  old_plan <- future::plan(future::sequential)
  on.exit(future::plan(old_plan), add = TRUE)

  mspec_shard <- use_shard(mspec)
  set.seed(42)
  res_seq <- muffle_worker_version_warnings(
    mvpa_iterate(mspec_shard, vox_iter, ids = seq_along(vox_iter))
  )

  # Multicore with shard
  future::plan(future::multicore, workers = 2)
  mspec_shard2 <- use_shard(mspec)
  set.seed(42)
  res_mc <- muffle_worker_version_warnings(
    mvpa_iterate(mspec_shard2, vox_iter, ids = seq_along(vox_iter))
  )

  expect_equal(nrow(res_seq), nrow(res_mc))
  expect_equal(sort(res_seq$id), sort(res_mc$id))
  expect_equal(res_seq$error, res_mc$error)

  # Performance parity for non-error rows
  if (!all(res_seq$error)) {
    ok <- !res_seq$error
    perf_seq <- lapply(res_seq$performance[ok], function(p) {
      if (is.null(p)) NULL else unlist(p)
    })
    perf_mc <- lapply(res_mc$performance[ok], function(p) {
      if (is.null(p)) NULL else unlist(p)
    })
    for (i in seq_along(perf_seq)) {
      if (!is.null(perf_seq[[i]]) && !is.null(perf_mc[[i]])) {
        expect_equal(perf_seq[[i]], perf_mc[[i]], tolerance = 1e-8)
      }
    }
  }
})
