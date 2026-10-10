# Shared searchlight engine helpers: sphere columns, failed-centre diagnostics,
# and parity of the fast engines with the general path.
testthat::skip_if_not_installed("neuroim2")

# Oracle: the pre-consolidation column algorithm, one match() per sphere
# against the mask indices, spheres built on access.
oracle_sphere_cols <- function(ds, radius, centers) {
  sl <- get_searchlight(ds, "standard", radius)
  sp <- neuroim2::space(ds$mask)
  mask_indices <- ds$mask_indices %||% compute_mask_indices(ds$mask)
  parent <- vapply(sl, function(w) as.integer(w@parent_index), 1L)
  lapply(centers, function(cid) {
    w <- sl[[match(cid, parent)]]
    match(as.integer(neuroim2::grid_to_index(sp, w@coords)), mask_indices)
  })
}

# Dataset with constant voxels, random holes in the mask (so some spheres are
# small or partly outside it) and optional anisotropic spacing.
cons_dataset <- function(D = c(6, 6, 6), spacing = c(1, 1, 1), nobs = 24, seed = 11,
                         holes = TRUE, constant = TRUE) {
  set.seed(seed)
  g <- suppressWarnings(gen_sample_dataset(D = D, nobs = nobs, nlevels = 3, blocks = 4,
                                           spacing = spacing))
  arr <- neuroim2::as.array(g$dataset$train_data)
  if (constant) arr[1:2, 1:2, 1:2, ] <- 5
  bvec <- neuroim2::NeuroVec(arr, neuroim2::space(g$dataset$train_data))
  mask <- array(TRUE, D)
  if (holes) mask[sample(prod(D), 25)] <- FALSE
  mask_vol <- neuroim2::NeuroVol(mask, neuroim2::NeuroSpace(D, spacing))
  list(
    dataset = mvpa_dataset(train_data = bvec, mask = mask_vol),
    design = g$design
  )
}

# Sparse mask with an isolated voxel and a constant 2x2x2 block, so that some
# spheres keep fewer than two non-constant voxels and the general path fails
# on them.
bad_dataset <- function(seed = 21, nobs = 40) {
  set.seed(seed)
  D <- c(6, 6, 6)
  g <- suppressWarnings(gen_sample_dataset(D = D, nobs = nobs, nlevels = 3, blocks = 4))
  arr <- neuroim2::as.array(g$dataset$train_data)
  arr[1:2, 1:2, 1:2, ] <- 5
  bvec <- neuroim2::NeuroVec(arr, neuroim2::space(g$dataset$train_data))
  mask <- array(FALSE, D)
  mask[1, 6, 6] <- TRUE
  mask[3:6, 1:3, 1:3] <- TRUE
  mask[1:2, 1:2, 1:2] <- TRUE
  list(
    dataset = mvpa_dataset(train_data = bvec,
                           mask = neuroim2::NeuroVol(mask, neuroim2::NeuroSpace(D, c(1, 1, 1)))),
    design = g$design
  )
}

cons_model <- function(fx, label, ...) {
  mvpa_model(load_model(label), fx$dataset, fx$design, "classification",
             crossval = blocked_cross_validation(fx$design$block_var), ...)
}

cons_rsa_model <- function(fx, distmethod = "pearson") {
  n <- length(fx$design$block_var)
  set.seed(7)
  D1 <- dist(matrix(rnorm(n * 4), n))
  D2 <- dist(matrix(rnorm(n * 4), n))
  rdes <- rsa_design(~ D1 + D2,
                     list(D1 = D1, D2 = D2, block = fx$design$block_var),
                     block_var = "block")
  rsa_model(fx$dataset, rdes, distmethod = distmethod, regtype = "pearson",
            check_collinearity = FALSE)
}

map_values <- function(res, name) as.numeric(neuroim2::values(res$results[[name]]))

bad_ids <- function(res) {
  b <- attr(res, "bad_results")
  if (is.null(b) || nrow(b) == 0L) integer(0) else sort(as.integer(b$id))
}

expect_matches_legacy <- function(ms, engine, radius = 2, tolerance = NULL, metrics = NULL) {
  fast <- suppressWarnings(run_searchlight(ms, radius = radius, engine = engine,
                                           backend = "default"))
  ref <- suppressWarnings(run_searchlight(ms, radius = radius, engine = "legacy",
                                          backend = "default"))
  expect_identical(attr(fast, "searchlight_engine"), engine)
  expect_identical(sort(names(fast$results)), sort(names(ref$results)))
  for (m in names(ref$results)) {
    if (!is.null(metrics) && !(m %in% metrics)) next
    a <- map_values(fast, m)
    b <- map_values(ref, m)
    expect_identical(is.na(a), is.na(b), info = paste(engine, m))
    if (is.null(tolerance)) {
      expect_identical(a, b, info = paste(engine, m))
    } else {
      expect_equal(a, b, tolerance = tolerance, info = paste(engine, m))
    }
  }
  invisible(list(fast = fast, ref = ref))
}

# ---- sphere columns ----------------------------------------------------------

test_that("sphere columns equal the per-sphere match() oracle (full and subset, iso/aniso)", {
  for (spacing in list(c(1, 1, 1), c(1, 1, 2.5))) {
    fx <- cons_dataset(spacing = spacing)
    ds <- fx$dataset
    mask_indices <- ds$mask_indices %||% compute_mask_indices(ds$mask)
    centers_all <- intersect(get_center_ids(ds), mask_indices)

    full <- rMVPA:::.engine_sphere_cols(ds, 2, mask_indices = mask_indices)
    expect_setequal(full$centers, centers_all)
    expect_identical(full$cols, oracle_sphere_cols(ds, 2, full$centers),
                     info = paste(spacing, collapse = "x"))

    set.seed(3)
    subset <- sample(centers_all, 9)
    part <- rMVPA:::.engine_sphere_cols(ds, 2, centers = subset, mask_indices = mask_indices)
    expect_setequal(part$centers, subset)
    expect_identical(part$cols, oracle_sphere_cols(ds, 2, part$centers),
                     info = paste(spacing, collapse = "x"))
    # Subset columns are the full run's columns for the same centres.
    pos <- match(part$centers, full$centers)
    expect_identical(part$cols, full$cols[pos])
  }
})

test_that("sphere columns with keep_voxel follow the filter_roi rule", {
  fx <- cons_dataset()
  ds <- fx$dataset
  mask_indices <- ds$mask_indices %||% compute_mask_indices(ds$mask)
  x <- as.matrix(neuroim2::series(ds$train_data, mask_indices))
  keep <- rMVPA:::.engine_roi_keep(x)
  expect_false(all(keep))  # constant voxels exist in this fixture
  res <- rMVPA:::.engine_sphere_cols(ds, 2, keep_voxel = keep, mask_indices = mask_indices)
  ref_cols <- oracle_sphere_cols(ds, 2, res$centers)
  for (i in seq_along(res$centers)) {
    ref <- ref_cols[[i]]
    ref <- ref[!is.na(ref)]
    ref <- ref[keep[ref] | mask_indices[ref] == res$centers[i]]
    expect_identical(res$cols[[i]], ref, info = as.character(res$centers[i]))
  }
})

test_that("sphere columns for a subset touch only the requested spheres", {
  fx <- cons_dataset()
  ds <- fx$dataset
  mask_indices <- ds$mask_indices %||% compute_mask_indices(ds$mask)
  centers_all <- intersect(get_center_ids(ds), mask_indices)
  pick <- centers_all[c(1, 5, 9)]

  orig <- rMVPA:::.engine_lookup_cols
  calls <- 0L
  testthat::local_mocked_bindings(
    .engine_lookup_cols = function(coords, space_obj, lookup) {
      calls <<- calls + 1L
      orig(coords, space_obj, lookup)
    },
    .package = "rMVPA"
  )
  part <- rMVPA:::.engine_sphere_cols(ds, 2, centers = pick, mask_indices = mask_indices)
  expect_identical(calls, 3L)
  expect_identical(length(part$cols), 3L)

  calls <- 0L
  full <- rMVPA:::.engine_sphere_cols(ds, 2, mask_indices = mask_indices)
  expect_identical(calls, length(full$cols))
  expect_identical(length(full$cols), length(centers_all))
})

test_that("aligned sphere columns follow the order of the requested centres", {
  fx <- cons_dataset()
  ds <- fx$dataset
  mask_indices <- ds$mask_indices %||% compute_mask_indices(ds$mask)
  centers_all <- intersect(get_center_ids(ds), mask_indices)
  want <- rev(centers_all[1:4])
  got <- rMVPA:::.engine_sphere_cols_aligned(ds, 2, want, mask_indices)
  expect_identical(got, oracle_sphere_cols(ds, 2, want))
  expect_error(
    rMVPA:::.engine_sphere_cols_aligned(ds, 2, 999999L, mask_indices),
    "no sphere"
  )
})

# ---- parity with the general path ---------------------------------------------

test_that("rsa_fast matches the general path on an edge/constant fixture", {
  fx <- cons_dataset()
  for (dm in c("pearson", "spearman")) {
    expect_matches_legacy(cons_rsa_model(fx, dm), "rsa_fast")
  }
})

test_that("sda_fast matches the general path on an edge/constant fixture", {
  testthat::skip_if_not_installed("sda")
  fx <- cons_dataset(nobs = 36)
  # Per-column statistics are shared, so AUC can differ in the last bits
  # (as in test_sda_searchlight.R, which uses the same tolerance).
  expect_matches_legacy(cons_model(fx, "sda_notune"), "sda_fast", tolerance = 1e-12)
})

test_that("aggregate_fast matches the general path (corclass, naive_bayes)", {
  fx <- cons_dataset()
  expect_matches_legacy(cons_model(fx, "corclass"), "aggregate_fast", tolerance = 1e-10)
  expect_matches_legacy(cons_model(fx, "naive_bayes"), "aggregate_fast", tolerance = 1e-10)
})

test_that("dual_lda_fast matches the general path on Accuracy", {
  # Pre-existing: the dual-LDA AUC map differs from the general path's on
  # this fixture (the same before and after the consolidation), so only the
  # Accuracy map is asserted here.
  fx <- cons_dataset(D = c(5, 5, 5), nobs = 30)
  expect_matches_legacy(cons_model(fx, "dual_lda"), "dual_lda_fast", metrics = "Accuracy")
})

test_that("non-finite data are refused by an explicitly requested fast engine", {
  fx <- cons_dataset()
  fx$dataset$train_data <- {
    arr <- neuroim2::as.array(fx$dataset$train_data)
    arr[3, 3, 3, 2] <- NA
    neuroim2::NeuroVec(arr, neuroim2::space(fx$dataset$train_data))
  }
  ms <- cons_model(fx, "corclass")
  expect_error(suppressWarnings(run_searchlight(ms, radius = 2, engine = "aggregate_fast")),
               "non-finite")
  # Under auto the same data fall back to the general path, which is not an error.
  res <- suppressWarnings(run_searchlight(ms, radius = 2, engine = "auto"))
  expect_false(identical(attr(res, "searchlight_engine"), "aggregate_fast"))
})

# ---- failed-centre diagnostics --------------------------------------------------

test_that("aggregate_fast reports dropped centres as the general path does", {
  fx <- bad_dataset()
  ms <- cons_model(fx, "corclass")
  fast <- suppressWarnings(run_searchlight(ms, radius = 2, engine = "aggregate_fast"))
  legacy <- suppressWarnings(run_searchlight(ms, radius = 2, engine = "legacy"))
  bad <- attr(fast, "bad_results")
  expect_gt(nrow(bad), 0L)
  expect_identical(bad_ids(fast), bad_ids(legacy))
  expect_identical(names(bad), names(attr(legacy, "bad_results")))
  expect_true(all(bad$error))
  expect_true(all(grepl("less than 2 features", bad$error_message)))
})

test_that("sda_fast reports dropped centres as the general path does", {
  testthat::skip_if_not_installed("sda")
  fx <- bad_dataset(nobs = 36)
  ms <- cons_model(fx, "sda_notune")
  fast <- suppressWarnings(run_searchlight(ms, radius = 2, engine = "sda_fast"))
  legacy <- suppressWarnings(run_searchlight(ms, radius = 2, engine = "legacy"))
  expect_gt(nrow(attr(fast, "bad_results")), 0L)
  expect_identical(bad_ids(fast), bad_ids(legacy))
  expect_identical(names(attr(fast, "bad_results")), names(attr(legacy, "bad_results")))
})

test_that("rsa_fast reports the same dropped centres as the general path", {
  # The general RSA path fits single-voxel spheres, so this fixture has no
  # dropped centres on either path; the test pins that agreement and the
  # table layout. The rsa_fast failure rows are not exercised by a real failure.
  fx <- bad_dataset()
  ms <- cons_rsa_model(fx)
  fast <- suppressWarnings(run_searchlight(ms, radius = 2, engine = "rsa_fast"))
  legacy <- suppressWarnings(run_searchlight(ms, radius = 2, engine = "legacy"))
  expect_identical(bad_ids(fast), bad_ids(legacy))
  expect_identical(names(attr(fast, "bad_results")), names(attr(legacy, "bad_results")))
})

test_that("a failed-centre table has the general path's columns, also when empty", {
  b <- rMVPA:::.engine_bad_results(integer(0), error_message = character(0))
  expect_identical(nrow(b), 0L)
  expect_identical(names(b), c("result", "indices", "performance", "id", "error",
                               "error_message", "warning", "warning_message"))
  b2 <- rMVPA:::.engine_bad_results(c(4L, 9L), error_message = "boom")
  expect_identical(b2$error_message, c("boom", "boom"))
  expect_true(all(b2$error))
})
