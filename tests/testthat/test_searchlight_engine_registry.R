library(rMVPA)

test_that("searchlight_engines reports registered engines", {
  tbl <- searchlight_engines()

  expect_true(is.data.frame(tbl))
  expect_true(all(c("engine", "label", "eligible") %in% names(tbl)))
  expect_true(all(c("legacy", "swift", "dual_lda_fast") %in% tbl$engine))
})

test_that("explain_searchlight_engine marks a selected engine", {
  ds <- gen_sample_dataset(c(4, 4, 4), 24, nlevels = 3, blocks = 3)
  mspec <- mvpa_model(
    load_model("corclass"),
    dataset = ds$dataset,
    design = ds$design,
    crossval = blocked_cross_validation(ds$design$block_var)
  )

  tbl <- explain_searchlight_engine(mspec, method = "standard", engine = "auto")

  expect_true(is.data.frame(tbl))
  expect_true(any(tbl$selected))
  expect_identical(tbl$requested[[1]], "auto")
})
