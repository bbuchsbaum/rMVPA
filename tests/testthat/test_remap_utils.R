context("remap_utils helpers")

test_that("summarize_remap_roi returns expected columns from performance_table", {
  skip_on_cran()
  set.seed(11)
  toy <- gen_sample_dataset(D = c(5,5,5), nobs = 60, nlevels = 3, blocks = 3, external_test = TRUE)
  regionMask <- neuroim2::NeuroVol(sample(1:3, size = length(toy$dataset$mask), replace = TRUE),
                                   neuroim2::space(toy$dataset$mask))

  mspec <- remap_rrr_model(
    dataset = toy$dataset,
    design  = toy$design,
    rank    = 0,                # identity to avoid rrpack dependency
    leave_one_key_out = TRUE
  )

  res <- run_regional(mspec, regionMask)
  tab <- summarize_remap_roi(res)
  expect_true(is.data.frame(tab))
  expect_true(all(c("roinum","mean_rank","mean_lambda","mean_roi_improv","mean_delta_frob") %in% names(tab)))
})

test_that("summarize_remap_items works with return_fits", {
  skip_on_cran()
  set.seed(12)
  toy <- gen_sample_dataset(D = c(4,4,4), nobs = 48, nlevels = 3, blocks = 3, external_test = TRUE)
  regionMask <- neuroim2::NeuroVol(sample(1:2, size = length(toy$dataset$mask), replace = TRUE),
                                   neuroim2::space(toy$dataset$mask))

  ms <- remap_rrr_model(dataset = toy$dataset, design = toy$design,
                        rank = 0, leave_one_key_out = TRUE, return_adapter = TRUE)
  res <- run_regional(ms, regionMask, return_fits = TRUE)

  items <- summarize_remap_items(res, roi = 1)
  expect_true(is.data.frame(items))
  expect_true(all(c("item","res_naive","res_remap","res_ratio","n_folds") %in% names(items)))
})

