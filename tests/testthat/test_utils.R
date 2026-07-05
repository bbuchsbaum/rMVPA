test_that("matrix utility helpers handle grouping, centering, and variance filters", {
  X <- matrix(c(
    1, 2,
    3, 4,
    5, 6,
    7, 8
  ), nrow = 4, byrow = TRUE)

  row_means <- group_means(X, margin = 1, group = factor(c("a", "a", "b", "b")))
  expect_equal(row_means["a", ], c(2, 3))
  expect_equal(row_means["b", ], c(6, 7))

  col_means <- group_means(X, margin = 2, group = factor(c("left", "right")))
  expect_equal(col_means[, "left"], X[, 1])
  expect_error(group_means(X, margin = 3, group = 1:4), "margin")

  expect_equal(rMVPA:::center_patterns(X, "none"), X)
  centered <- rMVPA:::center_patterns(X, "stimulus_mean")
  expect_equal(colMeans(centered), c(0, 0))

  split_centered <- rMVPA:::center_patterns_train_test(
    X,
    X + 10,
    method = "stimulus_mean"
  )
  expect_equal(colMeans(split_centered$train), c(0, 0))
  expect_equal(split_centered$test, sweep(X + 10, 2, colMeans(X), "-"))
  expect_error(
    rMVPA:::center_patterns_train_test(X, matrix(1, nrow = 2, ncol = 3), "stimulus_mean"),
    "ncol"
  )

  M <- cbind(const = 1, varying = 1:4, has_na = c(1, NA, 1, NA))
  expect_equal(unname(rMVPA:::zeroVarianceColumns(M)), c(1L, 3L))
  expect_equal(unname(rMVPA:::zeroVarianceColumns2(M)), c(TRUE, FALSE, TRUE))
  expect_equal(unname(rMVPA:::na_cols(M)), c(FALSE, FALSE, TRUE))
  expect_equal(unname(rMVPA:::nonzeroVarianceColumns(M)), 2L)
  expect_equal(unname(rMVPA:::nonzeroVarianceColumns2(M)), c(FALSE, TRUE, FALSE))
  expect_equal(rMVPA:::removeZeroVarianceColumns(M), M[, 2, drop = FALSE])
  expect_equal(rMVPA:::removeZeroVarianceColumns(matrix(1, nrow = 3, ncol = 2)),
               matrix(1, nrow = 3, ncol = 2))
})

test_that("coalesce joins and correlation helpers return expected values", {
  x <- data.frame(id = c(1, 2), value = c(NA, 20), keep = c("x1", "x2"))
  y <- data.frame(id = c(1, 2), value = c(10, 30), extra = c("y1", "y2"))

  joined <- rMVPA:::coalesce_join(x, y, by = "id")
  expect_equal(joined$value, c(10, 20))
  expect_equal(joined$extra, c("y1", "y2"))

  joined2 <- rMVPA:::coalesce_join2(x, y, by = "id")
  expect_equal(joined2$value, c(10, 20))
  expect_equal(joined2$keep, c("x1", "x2"))

  no_overlap <- rMVPA:::coalesce_join(
    data.frame(id = 1, a = 2),
    data.frame(id = 1, b = 3),
    by = "id"
  )
  expect_equal(no_overlap, data.frame(id = 1, a = 2, b = 3))

  expect_equal(rMVPA:::spearman_cor(1:5, 1:5), 1)
  expect_equal(rMVPA:::kendall_cor(1:5, 1:5), 1)
})

test_that("logging and sysinfo helpers are callable", {
  expect_false(rMVPA:::.use_crayon_styles())

  old_level <- futile.logger::flog.threshold()
  on.exit(futile.logger::flog.threshold(old_level), add = TRUE)
  expect_equal(set_log_level("WARN"), futile.logger::WARN)
  expect_error(set_log_level("NOPE"), "Unknown log level")

  info <- capture.output(out <- mvpa_sysinfo())
  expect_s3_class(out, "mvpa_sysinfo")
  expect_true(all(c("platform", "dependencies", "parallel_backend") %in% names(out)))
  expect_true(any(grepl("rMVPA System Information", info)))

  minimal <- structure(
    list(
      r_version = "R",
      platform = "",
      os = NULL,
      nodename = character(),
      user = "user",
      locale = "C",
      rmvpa_version = "0",
      parallel_backend = "sequential",
      dependencies = list()
    ),
    class = "mvpa_sysinfo"
  )
  expect_output(print(minimal), "No dependency versions")
})
