# class_metrics = TRUE adds per-class AUCs; the output schema must declare
# them, or every ROI fails the schema width check.
testthat::skip_if_not_installed("neuroim2")

test_that("class_metrics = TRUE works for multiclass searchlight and regional analyses", {
  set.seed(903)
  ds <- gen_sample_dataset(D = c(4, 4, 4), nobs = 48, nlevels = 3, blocks = 4)
  ms <- mvpa_model(load_model("corclass"), ds$dataset, ds$design, "classification",
                   crossval = blocked_cross_validation(ds$design$block_var),
                   class_metrics = TRUE)
  per_class <- paste0("AUC_", levels(ds$design$y_train))
  expect_identical(names(output_schema(ms)), c("Accuracy", "AUC", per_class))

  sl <- run_searchlight(ms, radius = 2, engine = "legacy", backend = "default")
  expect_true(all(c("Accuracy", "AUC", per_class) %in% names(sl$results)))

  regions <- neuroim2::NeuroVol(rep(1:2, length.out = 64), neuroim2::space(ds$dataset$mask))
  rr <- run_regional(ms, regions)
  expect_true(all(c("Accuracy", "AUC", per_class) %in% names(rr$performance_table)))
  expect_equal(nrow(rr$performance_table), 2L)
})
