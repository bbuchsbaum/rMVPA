# Permutation searchlights use the exact engines (prepared once, rescored per
# permutation). The null distribution and p-values must equal those of the
# per-ROI path for the same seed.
testthat::skip_if_not_installed("neuroim2")

perm_pair <- function(model, strategy) {
  set.seed(1301)
  ds <- gen_sample_dataset(c(4, 4, 4), 40, nlevels = 3, blocks = 4)
  ms <- mvpa_model(load_model(model), ds$dataset, ds$design, "classification",
                   crossval = blocked_cross_validation(ds$design$block_var))
  pc <- permutation_control(n_perm = 3, seed = 7, perm_strategy = strategy,
                            subsample = 0.5, diagnose = FALSE)
  obs <- run_searchlight(ms, radius = 2, engine = "legacy", backend = "default")
  list(
    fast = suppressWarnings(run_permutation_searchlight(ms, observed = obs, radius = 2,
                                                        perm_ctrl = pc, metric = "Accuracy")),
    ref = suppressWarnings(run_permutation_searchlight(ms, observed = obs, radius = 2,
                                                       perm_ctrl = pc, metric = "Accuracy",
                                                       engine = "legacy"))
  )
}

for (model in c("corclass", "naive_bayes", "sda_notune")) {
  for (strategy in c("iterate", "searchlight")) {
    local({
      m <- model; st <- strategy
      test_that(sprintf("engine permutations equal the per-ROI path (%s, %s)", m, st), {
        pp <- perm_pair(m, st)
        expect_equal(sort(unlist(pp$fast$adj_null$bin_nulls)),
                     sort(unlist(pp$ref$adj_null$bin_nulls)), tolerance = 1e-12)
        expect_equal(pp$fast$p_values, pp$ref$p_values, tolerance = 1e-12)
        expect_identical(pp$fast$n_perm_used, pp$ref$n_perm_used)
      })
    })
  }
}

test_that("models without an exact engine keep the per-ROI permutation path", {
  set.seed(1302)
  ds <- gen_sample_dataset(c(4, 4, 4), 40, nlevels = 3, blocks = 4)
  ms <- mvpa_model(load_model("dual_lda"), ds$dataset, ds$design, "classification",
                   crossval = blocked_cross_validation(ds$design$block_var))
  expect_null(rMVPA:::.permutation_engine(ms, 2, get_center_ids(ds$dataset), "Accuracy", list()))
})
