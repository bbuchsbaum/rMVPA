build_dual_lda_mspec <- function(D, nobs, nlevels, blocks, gamma = 1e-2, external_test = FALSE) {
  ds <- gen_sample_dataset(
    D = D,
    nobs = nobs,
    nlevels = nlevels,
    blocks = blocks,
    external_test = external_test
  )
  cv <- blocked_cross_validation(ds$design$block_var)
  model <- load_model("dual_lda")
  mspec <- mvpa_model(
    model = model,
    dataset = ds$dataset,
    design = ds$design,
    model_type = "classification",
    crossval = cv,
    tune_grid = data.frame(gamma = gamma)
  )
  list(mspec = mspec, dataset = ds)
}

metric_map_values <- function(res, metric) {
  neuroim2::values(res$results[[metric]])
}
