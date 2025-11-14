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

