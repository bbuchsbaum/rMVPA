context("remap_rrr_model predictor fields")

test_that("LOKO predictor carries per-item fields when return_adapter=TRUE", {
  skip_on_cran()
  set.seed(3)
  toy <- gen_sample_dataset(D = c(4,4,4), nobs = 48, nlevels = 3, blocks = 3, external_test = TRUE)
  regionMask <- neuroim2::NeuroVol(sample(1:2, size = length(toy$dataset$mask), replace = TRUE),
                                   neuroim2::space(toy$dataset$mask))

  ms <- remap_rrr_model(
    dataset = toy$dataset,
    design  = toy$design,
    base_classifier = "sda_notune",
    rank = 0,  # identity adapter to keep dependency-light
    leave_one_key_out = TRUE,
    return_adapter = TRUE,
    save_fold_singulars = TRUE,
    min_pairs = 3
  )

  res <- run_regional(ms, regionMask)
  # pull predictor (stored in fits when return_fits=TRUE)
  # For this test, call process_roi directly to avoid aggregation stripping
  vox <- which(toy$dataset$mask > 0)[1:20]
  roi <- extract_roi(vox, toy$dataset)
  row <- process_roi(ms, roi, 1)
  cres <- row$result[[1]]
  pred <- cres$predictor
  expect_true(is.list(pred))
  expect_true(all(c("r2_per_voxel","resid_by_item","rank_by_item","sv1_by_item","keys_scored","skipped_keys","sv_spectra_by_item") %in% names(pred)))
})
