testthat::skip_if_not_installed("neuroim2")

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

iter_core_dataset <- function(dims = c(5, 5, 5), nobs = 8, test_nobs = NULL,
                              seed = 1, train_array = NULL, test_array = NULL) {
  set.seed(seed)
  n <- prod(dims)
  tr_arr <- if (is.null(train_array)) array(rnorm(n * nobs), c(dims, nobs)) else train_array
  tr <- neuroim2::NeuroVec(tr_arr, neuroim2::NeuroSpace(dim(tr_arr)))
  te <- NULL
  if (!is.null(test_nobs) || !is.null(test_array)) {
    te_arr <- if (is.null(test_array)) array(rnorm(n * test_nobs), c(dims, test_nobs)) else test_array
    te <- neuroim2::NeuroVec(te_arr, neuroim2::NeuroSpace(dim(te_arr)))
  }
  mask <- neuroim2::NeuroVol(array(1, dims), neuroim2::NeuroSpace(dims))
  mvpa_dataset(tr, test_data = te, mask = mask)
}

iter_core_spec <- function() {
  ds <- gen_sample_dataset(D = c(4, 4, 4), nobs = 18, nlevels = 2, blocks = 3)
  cv <- blocked_cross_validation(ds$design$block_var)
  mspec <- mvpa_model(
    model = load_model("sda_notune"),
    dataset = ds$dataset,
    design = ds$design,
    model_type = "classification",
    crossval = cv
  )
  list(mspec = mspec, dataset = ds)
}

iter_core_passthrough <- function() {
  function(obj, roi, rnum, center_global_id = NA) {
    tibble::tibble(
      result = list(NULL),
      indices = list(NULL),
      performance = list(NULL),
      id = rnum,
      error = FALSE,
      error_message = "~",
      warning = FALSE,
      warning_message = "~"
    )
  }
}

# ---------------------------------------------------------------------------
# Oracles: copies of the code as it stood before this change
# ---------------------------------------------------------------------------

# Old as_roi.data_sample voxel restriction (setdiff / intersect per ROI).
oracle_restrict_vox <- function(vox, mask_indices) {
  invalid_vox <- setdiff(vox, mask_indices)
  if (length(invalid_vox) > 0) {
    vox <- intersect(vox, mask_indices)
  }
  vox
}

# Old .filter_roi_stats_fast (row loop of pmin/pmax).
oracle_filter_roi_stats <- function(trdat) {
  p <- ncol(trdat)
  if (is.null(p) || p == 0L) {
    return(list(nas = logical(0), sdnonzero = logical(0)))
  }
  nas <- colSums(is.na(trdat)) > 0L
  sdnonzero <- rep(FALSE, p)
  valid_cols <- which(!nas)
  if (length(valid_cols) > 0L && nrow(trdat) >= 2L) {
    x <- trdat[, valid_cols, drop = FALSE]
    mins <- x[1, , drop = TRUE]
    maxs <- mins
    if (nrow(x) > 1L) {
      for (i in 2:nrow(x)) {
        row_i <- x[i, , drop = TRUE]
        mins <- pmin(mins, row_i)
        maxs <- pmax(maxs, row_i)
      }
    }
    span <- maxs - mins
    sdnonzero[valid_cols] <- is.finite(span) & (span > 0)
  }
  list(nas = nas, sdnonzero = sdnonzero)
}

# Old prep_regional per-region scans.
oracle_prep_regional <- function(region_vec, mask_vec) {
  region_set <- sort(as.integer(unique(region_vec[region_vec > 0])))
  vox_iter <- lapply(region_set, function(rnum) which(region_vec == rnum & mask_vec > 0))
  lens <- sapply(vox_iter, length)
  keep <- lens > 1
  if (any(lens < 2)) {
    vox_iter <- vox_iter[keep]
    region_set <- region_set[keep]
  }
  list(region_set = region_set, vox_iter = vox_iter, lens = lens, keep = keep)
}

# ---------------------------------------------------------------------------
# as_roi: voxel restriction
# ---------------------------------------------------------------------------

test_that("as_roi restriction matches the setdiff/intersect oracle exactly", {
  dims <- c(5, 5, 5)
  N <- prod(dims)
  dset <- iter_core_dataset(dims = dims, nobs = 8, seed = 11)
  set.seed(7)
  pm <- numeric(N)
  pm[sample(N, 70)] <- 1
  dset$mask <- neuroim2::NeuroVol(array(pm, dims), neuroim2::NeuroSpace(dims))
  mask_idx <- which(pm != 0)
  out_idx <- which(pm == 0)

  # Lookup variant mirrors what mvpa_iterate stores on its local dataset copy.
  dset_lk <- dset
  dset_lk$mask_indices <- mask_idx
  lk <- logical(max(mask_idx))
  lk[mask_idx] <- TRUE
  dset_lk$mask_lookup <- lk

  for (trial in 1:60) {
    n_in <- if (trial %% 9 == 0) 0L else sample(0:20, 1)
    n_out <- if (trial %% 5 == 0) 0L else sample(0:6, 1)
    if (n_in + n_out == 0L) n_in <- 1L
    vox <- c(sample(mask_idx, n_in, replace = TRUE),
             if (n_out > 0L) sample(out_idx, n_out, replace = TRUE))
    if (trial %% 7 == 0) vox <- c(vox, N + 5L)      # beyond the grid
    if (trial %% 11 == 0) vox <- c(vox, NA_integer_)
    vox <- sample(vox)                              # unsorted, duplicated
    if (trial %% 2 == 0) vox <- as.double(vox)      # double-typed indices

    ov <- oracle_restrict_vox(vox, mask_idx)
    # Old as_roi returns a try-error (no series) when restriction empties the ROI.
    emptied <- length(setdiff(vox, mask_idx)) > 0 && length(ov) == 0
    exp_train <- if (emptied) {
      structure(list(message = "empty"), class = "try-error")
    } else {
      try(neuroim2::series_roi(dset$train_data, ov), silent = TRUE)
    }

    for (d in list(dset, dset_lk)) {
      got <- as_roi(data_sample(d, vox), d)
      if (inherits(exp_train, "try-error")) {
        expect_true(inherits(got$train_roi, "try-error"))
      } else {
        expect_identical(neuroim2::indices(got$train_roi), neuroim2::indices(exp_train))
        expect_identical(neuroim2::values(got$train_roi), neuroim2::values(exp_train))
      }
    }
  }
})

test_that("as_roi keeps an all-in-mask ROI unchanged, duplicates and order included", {
  dset <- iter_core_dataset(nobs = 6, seed = 5)
  vox <- c(9L, 3L, 9L, 20L, 3L)
  got <- as_roi(data_sample(dset, vox), dset)
  expected <- neuroim2::series_roi(dset$train_data, vox)
  expect_identical(neuroim2::indices(got$train_roi), neuroim2::indices(expected))
  expect_identical(neuroim2::values(got$train_roi), neuroim2::values(expected))
})

# ---------------------------------------------------------------------------
# filter_roi: keep decisions, unchanged return, test-set and centre rules
# ---------------------------------------------------------------------------

test_that("filter_roi stats match the old row-loop decisions (NA, constant, single-row, Inf)", {
  set.seed(99)
  for (trial in 1:150) {
    n <- sample(1:6, 1)
    p <- sample(0:9, 1)
    if (p == 0L) {
      m <- matrix(numeric(0), n, 0)
    } else if (trial %% 3 == 0) {
      # integer-typed trial: constants, NA and a few distinct values
      m <- matrix(sample(c(-2L, 0L, 1L, 2L, NA_integer_), n * p, replace = TRUE), n, p)
      m[, sample(p, 1)] <- 4L
    } else {
      m <- matrix(rnorm(n * p), n, p)
      for (j in seq_len(p)) {
        r <- runif(1)
        if (r < 0.2) {
          m[, j] <- 2.5                                         # constant
        } else if (r < 0.35) {
          m[sample(n, 1), j] <- NA                              # NA
        } else if (r < 0.45) {
          m[sample(n, 1), j] <- Inf                             # infinite
        } else if (r < 0.5) {
          m[, j] <- c(1e308, -1e308, rep(0, max(0L, n - 2L)))[seq_len(n)] # overflowing span
        } else if (r < 0.55) {
          m[, j] <- NaN                                         # NaN counts as NA
        }
      }
    }
    new <- .filter_roi_stats_fast(m)
    old <- oracle_filter_roi_stats(m)
    expect_identical(new$nas, old$nas)
    expect_identical(new$sdnonzero, old$sdnonzero)
  }
})

test_that("filter_roi returns the input unchanged when nothing is filtered", {
  dset <- iter_core_dataset(nobs = 10, seed = 21)
  roi <- as_roi(data_sample(dset, c(1L, 2L, 7L, 30L)), dset)
  expect_identical(filter_roi(roi, min_voxels = 2), roi)
})

test_that("filter_roi drops a voxel with an infinite test-set value", {
  dims <- c(4, 4, 4)
  set.seed(31)
  tr <- array(rnorm(prod(dims) * 10), c(dims, 10))
  te <- array(rnorm(prod(dims) * 6), c(dims, 6))
  te[2, 1, 1, 3] <- Inf               # voxel 2: one infinite test value
  dset <- iter_core_dataset(dims = dims, train_array = tr, test_array = te)
  roi <- as_roi(data_sample(dset, c(1L, 2L, 5L, 9L, 13L)), dset)
  filtered <- filter_roi(roi, min_voxels = 2)
  expect_setequal(neuroim2::indices(filtered$train_roi), c(1L, 5L, 9L, 13L))
  expect_setequal(neuroim2::indices(filtered$test_roi), c(1L, 5L, 9L, 13L))
  expect_true(all(!is.infinite(neuroim2::values(filtered$test_roi))))
})

test_that("filter_roi drops test-set NA by default and keeps it under test_na = 'keep'", {
  # Only models that handle missing test values (era_rsa_model's
  # pairwise-complete fallback) opt in to keeping such columns.
  dims <- c(4, 4, 4)
  set.seed(32)
  tr <- array(rnorm(prod(dims) * 10), c(dims, 10))
  te <- array(rnorm(prod(dims) * 6), c(dims, 6))
  te[1, 1, 1, ] <- NA
  dset <- iter_core_dataset(dims = dims, train_array = tr, test_array = te)
  roi <- as_roi(data_sample(dset, c(1L, 2L, 5L, 9L)), dset)
  expect_identical(filter_roi(roi, min_voxels = 2, test_na = "keep"), roi)
  filtered <- filter_roi(roi, min_voxels = 2)
  expect_setequal(neuroim2::indices(filtered$test_roi), c(2L, 5L, 9L))
})

test_that("filter_roi does not preserve an all-NA centre voxel", {
  dims <- c(4, 4, 4)
  set.seed(41)
  tr <- array(rnorm(prod(dims) * 10), c(dims, 10))
  tr[1, 1, 1, ] <- NA                 # centre voxel 1 is all-NA in training data
  dset <- iter_core_dataset(dims = dims, train_array = tr)
  roi <- as_roi(data_sample(dset, c(1L, 2L, 5L, 9L)), dset)
  filtered <- filter_roi(roi, preserve = 1L, min_voxels = 2)
  expect_false(1L %in% neuroim2::indices(filtered$train_roi))
  expect_setequal(neuroim2::indices(filtered$train_roi), c(2L, 5L, 9L))
})

# ---------------------------------------------------------------------------
# prep_regional: one-pass grouping
# ---------------------------------------------------------------------------

test_that("prep_regional matches the per-region scan on a mask with gaps and zero labels", {
  dims <- c(6, 7, 5)
  N <- prod(dims)
  set.seed(51)
  labels <- sample(c(0, 0, 0, 2, 5, 9, 11, 14), N, replace = TRUE)
  labels[sample(N, 3)] <- 0          # extra background
  region <- neuroim2::NeuroVol(array(labels, dims), neuroim2::NeuroSpace(dims))
  dmask_vals <- as.numeric(runif(N) < 0.7)
  dmask <- neuroim2::NeuroVol(array(dmask_vals, dims), neuroim2::NeuroSpace(dims))
  model_spec <- list(dataset = list(mask = dmask))

  res <- suppressWarnings(prep_regional(model_spec, region))
  oracle <- oracle_prep_regional(as.vector(region), as.vector(dmask))

  expect_identical(res$region_set, oracle$region_set)
  expect_identical(res$vox_iter, oracle$vox_iter)
  expect_identical(res$lens, oracle$lens)
  expect_identical(res$keep, oracle$keep)
})

# ---------------------------------------------------------------------------
# mvpa_iterate: dropped ROIs, argument errors
# ---------------------------------------------------------------------------

test_that("mvpa_iterate emits an error row for every size-dropped ROI and counts them", {
  built <- iter_core_spec()
  vox_list <- list(1:10, 20L, 30:40, integer(0))
  old_opt <- options(rMVPA.profile_searchlight = TRUE)
  on.exit(options(old_opt), add = TRUE)

  res <- mvpa_iterate(
    mod_spec = built$mspec,
    vox_list = vox_list,
    ids = 11:14,
    batch_size = 2,
    verbose = FALSE,
    analysis_type = "regional",
    processor = iter_core_passthrough()
  )

  expect_identical(sort(as.integer(res$id)), 11:14)
  expect_identical(sort(res$id[res$error %in% TRUE]), c(12L, 14L))
  dropped <- res[res$id %in% c(12L, 14L), , drop = FALSE]
  expect_true(all(dropped$error_message == "ROI filtered out (size < 2 voxels)"))
  expect_true(all(dropped$warning))
  timing <- attr(res, "timing")
  expect_identical(as.integer(timing$skipped_rois), 2L)
})

test_that("mvpa_iterate searchlight threshold in dropped-ROI messages is 1 voxel", {
  built <- iter_core_spec()
  vox_list <- list(1:10, 20L, 30:40, integer(0))
  res <- mvpa_iterate(
    mod_spec = built$mspec,
    vox_list = vox_list,
    ids = 11:14,
    batch_size = 2,
    verbose = FALSE,
    analysis_type = "searchlight",
    processor = iter_core_passthrough()
  )
  expect_identical(res$id[res$error %in% TRUE], 14L)
  expect_identical(res$error_message[res$id == 14L], "ROI filtered out (size < 1 voxels)")
})

test_that("mvpa_iterate lets argument errors propagate instead of returning a row", {
  built <- iter_core_spec()
  expect_error(
    mvpa_iterate(
      mod_spec = built$mspec,
      vox_list = list(1:10, 30:40),
      ids = 1:3,
      verbose = FALSE,
      analysis_type = "regional",
      processor = iter_core_passthrough()
    ),
    "length\\(ids\\)"
  )
})

# ---------------------------------------------------------------------------
# summarize_errors: grouping and NA-robustness
# ---------------------------------------------------------------------------

test_that("summarize_errors groups identical failures across ROIs into one cause", {
  built <- iter_core_spec()
  N <- length(built$dataset$dataset$mask)
  # Voxels outside the grid: as_roi leaves nothing, so each ROI fails validation.
  vox_list <- list(N + 10L:11L, N + 20L:21L, N + 30L:31L)
  res <- mvpa_iterate(
    mod_spec = built$mspec,
    vox_list = vox_list,
    ids = 1:3,
    batch_size = 3,
    verbose = FALSE,
    analysis_type = "searchlight",
    processor = iter_core_passthrough()
  )
  expect_identical(sum(res$error %in% TRUE), 3L)

  lines <- character(0)
  old_app <- futile.logger::flog.appender()
  old_thr <- futile.logger::flog.threshold()
  futile.logger::flog.appender(function(line) lines <<- c(lines, line))
  futile.logger::flog.threshold("INFO")
  on.exit({
    futile.logger::flog.appender(old_app)
    futile.logger::flog.threshold(old_thr)
  }, add = TRUE)

  n <- summarize_errors(res, "core")
  expect_identical(as.integer(n), 3L)
  expect_true(any(grepl("\\[3 ROIs\\]", lines)))
  expect_false(any(grepl("\\[1 ROIs\\]", lines)))
})

test_that("summarize_errors tolerates NA in the error column", {
  res <- tibble::tibble(
    id = 1:3,
    error = c(TRUE, NA, FALSE),
    error_message = c("boom", "boom", "~")
  )
  n <- NULL
  expect_no_error(n <- summarize_errors(res, "na-check"))
  expect_identical(as.integer(n), 1L)
})

# ---------------------------------------------------------------------------
# cv_run_global: folds from indices, results identical to the old code
# ---------------------------------------------------------------------------

test_that("cv_run_global folds are identical whether built from the matrix or row indices", {
  set.seed(3)
  X <- matrix(rnorm(24 * 5), 24, 5)
  y <- factor(rep(c("a", "b"), length.out = 24))
  cv <- blocked_cross_validation(rep(1:4, each = 6))
  old_folds <- crossval_samples(cv, as.data.frame(X), y)
  new_folds <- crossval_samples(cv, tibble::tibble(.row = seq_len(nrow(X))), y)
  expect_identical(lapply(old_folds$train, function(r) as.integer(r$idx)),
                   lapply(new_folds$train, function(r) as.integer(r$idx)))
  expect_identical(lapply(old_folds$test, function(r) as.integer(r$idx)),
                   lapply(new_folds$test, function(r) as.integer(r$idx)))
})

test_that("cv_run_global gives the same results as the pre-change code", {
  # Reference values recorded from the code before this change, with the same
  # seed, data and model as below.
  set.seed(42)
  T_obs <- 24L
  P <- 7L
  X <- matrix(rnorm(T_obs * P), T_obs, P)
  y <- factor(rep(c("a", "b"), length.out = T_obs))
  X[y == "a", 1:2] <- X[y == "a", 1:2] + 1.5
  cv <- blocked_cross_validation(rep(1:4, each = 6))
  ds <- gen_sample_dataset(D = c(4, 4, 4), nobs = T_obs, nlevels = 2, blocks = 4)
  mspec <- mvpa_model(model = load_model("sda_notune"), dataset = ds$dataset,
                      design = ds$design, model_type = "classification", crossval = cv)
  res <- cv_run_global(mspec, X, y, cv, feature_ids = seq_len(P), extract_weights_fn = TRUE)
  tab <- res$result_table

  expect_identical(nrow(tab), 4L)
  expect_identical(tab$error, c(FALSE, FALSE, FALSE, FALSE))
  expect_identical(
    vapply(tab$test_ind, paste, character(1), collapse = " "),
    c("1 2 3 4 5 6", "7 8 9 10 11 12", "13 14 15 16 17 18", "19 20 21 22 23 24")
  )
  expect_identical(
    vapply(tab$y_true, function(v) paste(as.character(v), collapse = ""), character(1)),
    c("ababab", "ababab", "ababab", "ababab")
  )
  expect_identical(
    vapply(tab$class, function(v) paste(as.character(unlist(v)), collapse = ""), character(1)),
    c("ababab", "ababab", "bbbbab", "abbbaa")
  )
  expect_equal(
    sum(vapply(res$fold_weights, function(w) if (is.null(w)) 0 else sum(abs(w)), numeric(1))),
    22.5248501851, tolerance = 1e-8
  )
  expect_identical(vapply(res$fold_weights, is.null, logical(1)), rep(FALSE, 4))
})

# ---------------------------------------------------------------------------
# external_crossval: standard error-row schema and the ncol < 2 guard
# ---------------------------------------------------------------------------

test_that("external_crossval failure rows use the standard CV error schema", {
  built <- iter_core_spec()
  std_names <- c("result", "indices", "performance", "id", "error",
                 "error_message", "warning", "warning_message")

  # One feature: the ncol < 2 guard internal_crossval also has.
  row1 <- external_crossval(
    built$mspec, roi = list(), id = 7L,
    xtrain_mat = matrix(rnorm(18), 18, 1),
    xtest_mat = matrix(rnorm(6), 6, 1),
    ind = 1L
  )
  expect_named(row1, std_names)
  expect_true(row1$error[[1]])
  expect_identical(row1$id[[1]], 7L)
  expect_identical(row1$error_message[[1]], "error: less than 2 features")

  # Centre voxel removed by filtering: same schema.
  row2 <- external_crossval(
    built$mspec, roi = list(), id = 8L, center_global_id = 99L,
    xtrain_mat = matrix(rnorm(36), 18, 2),
    xtest_mat = matrix(rnorm(12), 6, 2),
    ind = c(1L, 2L)
  )
  expect_named(row2, std_names)
  expect_true(row2$error[[1]])
  expect_identical(row2$id[[1]], 8L)
})
