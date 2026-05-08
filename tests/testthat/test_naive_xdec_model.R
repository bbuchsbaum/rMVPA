context("naive_xdec_model")

test_that("naive_xdec_model runs and returns performance", {
  skip_on_cran()
  set.seed(7)
  toy <- gen_sample_dataset(D = c(5,5,5), nobs = 60, nlevels = 3, blocks = 3, external_test = TRUE)
  regionMask <- neuroim2::NeuroVol(sample(1:3, size = length(toy$dataset$mask), replace = TRUE),
                                   neuroim2::space(toy$dataset$mask))
  ms <- naive_xdec_model(toy$dataset, toy$design, link_by = NULL)
  res <- run_regional(ms, regionMask)
  expect_s3_class(res, "regional_mvpa_result")
  expect_true(is.data.frame(res$performance_table))
})

test_that("naive_xdec_model honours a custom performance function", {
  skip_on_cran()
  set.seed(11)
  toy <- gen_sample_dataset(D = c(5,5,5), nobs = 60, nlevels = 3, blocks = 3,
                             external_test = TRUE)

  # `seen_test_design` records what the custom function actually receives.
  seen_test_design <- new.env(parent = emptyenv())
  seen_test_design$value <- NULL

  custom_fun <- function(result) {
    seen_test_design$value <- result$test_design
    obs <- as.character(result$observed)
    pred <- as.character(result$predicted)
    c(custom_acc = mean(obs == pred))
  }

  regionMask <- neuroim2::NeuroVol(
    sample(1:2, size = length(toy$dataset$mask), replace = TRUE),
    neuroim2::space(toy$dataset$mask)
  )

  ms <- naive_xdec_model(toy$dataset, toy$design, performance = custom_fun)

  # The custom hook must be annotated as "custom" so the fast-metric kernel
  # cannot bypass it.
  expect_identical(attr(ms$performance, "rmvpa_perf_kind"), "custom")

  res <- run_regional(ms, regionMask)
  expect_s3_class(res, "regional_mvpa_result")
  expect_true("custom_acc" %in% colnames(res$performance_table))
  expect_true(all(is.finite(res$performance_table$custom_acc)))
  # The custom function did receive a test_design (the whole point of the hook).
  expect_true(!is.null(seen_test_design$value))
})

test_that("naive_xdec_model rejects non-function `performance`", {
  skip_on_cran()
  toy <- gen_sample_dataset(D = c(5,5,5), nobs = 30, nlevels = 2, blocks = 3,
                             external_test = TRUE)
  expect_error(
    naive_xdec_model(toy$dataset, toy$design, performance = "not_a_function"),
    "must be a function"
  )
})

