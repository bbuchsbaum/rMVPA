context("one error-row schema for every failed ROI")

.err_cols <- c("result", "indices", "performance", "id", "error",
               "error_message", "error_trace", "warning", "warning_message")

.expect_error_schema <- function(rows, info = NULL) {
  expect_identical(names(rows), .err_cols, info = info)
  expect_true(all(rows$error), info = info)
  expect_true(all(rows$warning), info = info)
  expect_type(rows$error_message, "character")
  expect_type(rows$error_trace, "character")
  expect_type(rows$warning_message, "character")
  expect_true(is.list(rows$result) && is.list(rows$indices) &&
                is.list(rows$performance), info = info)
  # Messages carry no ROI id or traceback, so identical failures group.
  expect_false(any(grepl("trace:", rows$error_message, fixed = TRUE)), info = info)
}

.err_spec <- function() {
  ds <- gen_sample_dataset(D = c(4, 4, 4), nobs = 18, nlevels = 2, blocks = 3)
  mspec <- mvpa_model(load_model("sda_notune"), ds$dataset, ds$design,
                      "classification",
                      crossval = blocked_cross_validation(ds$design$block_var))
  list(mspec = mspec, ds = ds)
}

.err_failing_processor <- function() {
  function(obj, roi, rnum, center_global_id = NA) {
    if (rnum %% 2L == 0L) stop("deliberate failure")
    tibble::tibble(
      result = list(NULL), indices = list(NULL), performance = list(NULL),
      id = rnum, error = FALSE, error_message = "~",
      warning = FALSE, warning_message = "~"
    )
  }
}

test_that("the constructor gives the schema, recycles, and handles zero rows", {
  r <- rMVPA:::.roi_error_rows(c(3L, 5L), "boom", indices = list(1:2))
  .expect_error_schema(r)
  expect_identical(r$indices, list(1:2, 1:2))
  expect_identical(r$error_trace, c(NA_character_, NA_character_))
  expect_identical(nrow(rMVPA:::.roi_error_rows(integer(0), character(0))), 0L)
  expect_named(rMVPA:::.roi_error_rows(integer(0), character(0)), .err_cols)
})

test_that("model-reported failures (roi_result) use the same schema", {
  rr <- rMVPA:::roi_result(metrics = NULL, indices = 1:3, id = 4L,
                           error = TRUE, error_message = "too few")
  row <- rMVPA:::roi_result_to_tibble(rr)
  .expect_error_schema(row)
  expect_identical(row$indices, list(1:3))
  expect_identical(row$warning_message, "too few")
})

test_that("processor errors, invalid ROIs and size drops share one schema (default backend)", {
  fx <- .err_spec()
  N <- length(fx$ds$dataset$mask)
  vox_list <- list(1:10, 11:20, N + 1:2, 21L)
  res <- mvpa_iterate(fx$mspec, vox_list, ids = 1:4, batch_size = 4,
                      verbose = FALSE, analysis_type = "regional",
                      processor = .err_failing_processor())
  bad <- res[res$error %in% TRUE, , drop = FALSE]
  expect_identical(sort(as.integer(bad$id)), 2:4)
  .expect_error_schema(bad[, .err_cols])
  thrown <- bad[bad$id == 2L, ]
  expect_identical(thrown$error_message, "Error processing ROI: deliberate failure")
  expect_false(is.na(thrown$error_trace))
  expect_identical(
    bad$error_message[bad$id == 3L],
    "ROI failed validation (e.g., <2 voxels after filtering or other extract_roi issue)"
  )
  expect_identical(bad$error_message[bad$id == 4L], "ROI filtered out (size < 2 voxels)")
})

test_that("the shard backend writes the same error rows as the default backend", {
  skip_if_not_installed("shard")
  fx <- .err_spec()
  N <- length(fx$ds$dataset$mask)
  vox_list <- list(1:10, 11:20, N + 1:2, 21L)
  run <- function(spec) {
    res <- mvpa_iterate(spec, vox_list, ids = 1:4, batch_size = 4,
                        verbose = FALSE, analysis_type = "regional",
                        processor = .err_failing_processor())
    res <- res[res$error %in% TRUE, , drop = FALSE]
    res[order(res$id), , drop = FALSE]
  }
  default <- run(fx$mspec)
  mspec_shard <- use_shard(fx$mspec)
  on.exit(rMVPA:::shard_cleanup(mspec_shard$shard_data), add = TRUE)
  shard <- run(mspec_shard)
  .expect_error_schema(shard[, .err_cols])
  expect_identical(shard$id, default$id)
  expect_identical(shard$error_message, default$error_message)
  expect_identical(shard$warning_message, default$warning_message)
  expect_identical(is.na(shard$error_trace), is.na(default$error_trace))
})

test_that("internal and external cross-validation failure rows use the schema", {
  fx <- .err_spec()
  ext <- rMVPA:::external_crossval(
    fx$mspec, roi = list(), id = 7L,
    xtrain_mat = matrix(rnorm(18), 18, 1), xtest_mat = matrix(rnorm(6), 6, 1),
    ind = 1L
  )
  .expect_error_schema(ext, "external <2 features")

  vox <- which(as.logical(fx$ds$dataset$mask))[1:6]
  roi <- rMVPA:::extract_roi(rMVPA:::data_sample(fx$ds$dataset, vox),
                             fx$ds$dataset)
  int <- rMVPA:::internal_crossval(fx$mspec, roi, id = 9L,
                                   center_global_id = 10^6)
  .expect_error_schema(int, "internal centre missing")
  expect_identical(int$warning_message, "Center voxel removed by filter_roi")
})

test_that("fast-engine failed centres use the schema", {
  b <- rMVPA:::.engine_bad_results(c(4L, 9L), error_message = "boom")
  .expect_error_schema(b)
  expect_identical(b$indices, list(integer(0), integer(0)))
})

# ---------------------------------------------------------------------------
# Row alignment: an invalid ROI must not shift or recycle other ROIs' data
# ---------------------------------------------------------------------------

.err_index_processor <- function() {
  function(obj, roi, rnum, center_global_id = NA) {
    tibble::tibble(
      result = list(NULL), indices = list(NULL),
      performance = list(c(first = min(neuroim2::indices(roi$train_roi)))),
      id = rnum, error = FALSE, error_message = "~",
      warning = FALSE, warning_message = "~"
    )
  }
}

test_that("each ROI is processed with its own data when an earlier ROI is invalid", {
  fx <- .err_spec()
  N <- length(fx$ds$dataset$mask)
  for (backend in c("default", "shard")) {
    spec <- fx$mspec
    if (backend == "shard") {
      skip_if_not_installed("shard")
      spec <- use_shard(fx$mspec)
    }
    res <- mvpa_iterate(spec, list(1:10, N + 1:2, 21:30, 41:50),
                        ids = 1:4, batch_size = 4, verbose = FALSE,
                        analysis_type = "regional",
                        processor = .err_index_processor())
    if (backend == "shard") rMVPA:::shard_cleanup(spec$shard_data)
    res <- res[order(res$id), ]
    expect_identical(res$error, c(FALSE, TRUE, FALSE, FALSE), info = backend)
    first <- vapply(res$performance[!res$error], function(p) p[["first"]], numeric(1))
    expect_equal(first, c(1, 21, 41), info = backend)
  }
})

test_that("an invalid last ROI is an error row, not a copy of another ROI", {
  fx <- .err_spec()
  N <- length(fx$ds$dataset$mask)
  res <- mvpa_iterate(fx$mspec, list(1:10, N + 1:2), ids = 1:2, batch_size = 2,
                      verbose = FALSE, analysis_type = "regional",
                      processor = .err_index_processor())
  res <- res[order(res$id), ]
  expect_identical(res$error, c(FALSE, TRUE))
})

test_that("drop_probs keeps ROIs without probabilities and aligns prob_observed", {
  set.seed(8)
  rows <- tibble::tibble(
    result = list(
      NULL,
      structure(list(probs = matrix(c(0.2, 0.8, 0.6, 0.4), 2,
                                    dimnames = list(NULL, c("a", "b"))),
                     observed = factor(c("b", "a"), levels = c("a", "b"))),
                class = c("classification_result", "list"))
    ),
    id = 1:2
  )
  out <- rMVPA:::.drop_roi_probs(rows)
  expect_identical(length(out$prob_observed), 2L)
  expect_null(out$prob_observed[[1]])
  expect_null(out$result[[2]]$probs)

  one <- rMVPA:::.drop_roi_probs(tibble::tibble(result = list(NULL), id = 1L))
  expect_identical(nrow(one), 1L)
  expect_null(one$prob_observed[[1]])
})

test_that("drop_probs = TRUE with return_predictions = FALSE does not fail ROIs", {
  fx <- .err_spec()
  spec <- fx$mspec
  spec$return_predictions <- FALSE
  res <- mvpa_iterate(spec, list(1:20, 21:40), ids = 1:2, batch_size = 2,
                      verbose = FALSE, analysis_type = "regional",
                      drop_probs = TRUE)
  expect_false(any(res$error), info = paste(res$error_message, collapse = "; "))
  expect_true(all(vapply(res$performance, function(p) !is.null(p), logical(1))))
})
