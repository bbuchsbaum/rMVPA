test_that("roi_result constructs success correctly", {
  res <- roi_result(
    metrics = c(accuracy = 0.85, AUC = 0.9),
    indices = 1:10,
    id = 42
  )
  expect_s3_class(res, "roi_result")
  expect_false(res$error)
  expect_equal(res$error_message, "~")
  expect_false(res$warning)
  expect_equal(res$warning_message, "~")
  expect_equal(res$id, 42)
  expect_equal(length(res$indices), 10)
  expect_equal(res$metrics[["accuracy"]], 0.85)
})

test_that("roi_result constructs error correctly", {
  res <- roi_result(
    metrics = NULL,
    indices = 1:5,
    id = 7,
    error = TRUE,
    error_message = "Too few features"
  )
  expect_true(res$error)
  expect_equal(res$error_message, "Too few features")
  expect_null(res$metrics)
})

test_that("roi_result stores warning fields", {
  res <- roi_result(
    metrics = c(r2 = 0.5),
    indices = 1:3,
    id = 1,
    warning = TRUE,
    warning_message = "Near-singular matrix"
  )
  expect_false(res$error)
  expect_true(res$warning)
  expect_equal(res$warning_message, "Near-singular matrix")
})

test_that("roi_result_to_tibble converts success", {
  res <- roi_result(
    metrics = c(accuracy = 0.9),
    indices = 1:4,
    id = 10,
    result = list(observed = 1, predicted = 1)
  )
  tbl <- rMVPA:::roi_result_to_tibble(res)
  expect_s3_class(tbl, "tbl_df")
  expect_equal(nrow(tbl), 1)
  expect_false(tbl$error)
  expect_equal(tbl$error_message, "~")
  expect_equal(tbl$id, 10)
  expect_equal(tbl$performance[[1]], c(accuracy = 0.9))
  expect_false(tbl$warning)
  expect_equal(tbl$warning_message, "~")
})

test_that("roi_result_to_tibble converts error", {
  res <- roi_result(
    metrics = NULL,
    indices = integer(),
    id = 99,
    error = TRUE,
    error_message = "boom",
    warning = TRUE,
    warning_message = "also boom"
  )
  tbl <- rMVPA:::roi_result_to_tibble(res)
  expect_equal(nrow(tbl), 1)
  expect_true(tbl$error)
  expect_equal(tbl$error_message, "boom")
  expect_null(tbl$performance[[1]])
  expect_true(tbl$warning)
  expect_equal(tbl$warning_message, "also boom")
})

test_that("roi_result_to_tibble rejects non-roi_result", {
  expect_error(
    rMVPA:::roi_result_to_tibble(list(a = 1)),
    "inherits"
  )
})

# --- split_results ---

test_that("split_results splits by error column", {
  df <- tibble::tibble(
    id = 1:5,
    error = c(FALSE, TRUE, FALSE, TRUE, FALSE),
    error_message = c("~", "bad1", "~", "bad2", "~"),
    value = 10:14
  )
  sp <- rMVPA:::split_results(df)
  expect_equal(nrow(sp$good), 3)
  expect_equal(nrow(sp$bad), 2)
  expect_equal(sp$good$id, c(1L, 3L, 5L))
  expect_equal(sp$bad$id, c(2L, 4L))
})

test_that("split_results handles all-good", {
  df <- tibble::tibble(id = 1:3, error = c(FALSE, FALSE, FALSE))
  sp <- rMVPA:::split_results(df)
  expect_equal(nrow(sp$good), 3)
  expect_equal(nrow(sp$bad), 0)
})

test_that("split_results handles all-bad", {
  df <- tibble::tibble(id = 1:2, error = c(TRUE, TRUE))
  sp <- rMVPA:::split_results(df)
  expect_equal(nrow(sp$good), 0)
  expect_equal(nrow(sp$bad), 2)
})

test_that("split_results handles missing error column gracefully", {
  df <- tibble::tibble(id = 1:3, value = 4:6)
  sp <- rMVPA:::split_results(df)
  expect_equal(nrow(sp$good), 3)
  expect_equal(nrow(sp$bad), 0)
})

test_that("split_results handles non-data.frame input", {
  sp <- rMVPA:::split_results("not a df")
  expect_equal(sp$good, "not a df")
  expect_s3_class(sp$bad, "tbl_df")
  expect_equal(nrow(sp$bad), 0)
})

# --- summarize_errors ---

test_that("summarize_errors returns 0 for no errors", {
  df <- tibble::tibble(
    id = 1:3,
    error = c(FALSE, FALSE, FALSE),
    error_message = rep("~", 3)
  )
  n <- rMVPA:::summarize_errors(df, "test_context")
  expect_equal(n, 0L)
})

test_that("summarize_errors returns error count", {
  df <- tibble::tibble(
    id = 1:5,
    error = c(FALSE, TRUE, TRUE, FALSE, TRUE),
    error_message = c("~", "err_a", "err_a", "~", "err_b")
  )
  n <- rMVPA:::summarize_errors(df, "test_context")
  expect_equal(n, 3L)
})

test_that("summarize_errors handles non-data.frame", {
  n <- rMVPA:::summarize_errors(NULL, "test")
  expect_equal(n, 0L)
})

test_that("summarize_errors handles missing error_message column", {
  df <- tibble::tibble(id = 1:2, error = c(TRUE, TRUE))
  # Should not error even without error_message column
  n <- rMVPA:::summarize_errors(df, "test")
  expect_equal(n, 2L)
})
