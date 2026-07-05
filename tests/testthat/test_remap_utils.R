context("remap_utils helpers")

test_that("summarize_remap_roi uses performance_table and fit diagnostics", {
  perf_res <- list(
    performance_table = data.frame(
      roinum = c(10L, 20L),
      adapter_rank = c(1, 2),
      lambda_mean = c(0.1, 0.2),
      remap_improv = c(0.4, 0.5),
      delta_frob_mean = c(3, 4)
    )
  )
  class(perf_res) <- c("regional_mvpa_result", "list")

  tab <- summarize_remap_roi(perf_res)
  expect_equal(tab$roinum, c(10L, 20L))
  expect_equal(tab$mean_rank, c(1, 2))
  expect_equal(tab$mean_lambda, c(0.1, 0.2))
  expect_equal(tab$mean_roi_improv, c(0.4, 0.5))
  expect_equal(tab$mean_delta_frob, c(3, 4))

  fit_res <- list(
    fits = list(
      list(diag_by_fold = list(
        list(rank_used = 1, lambda_used = 0.1, roi_improv = -0.2, delta_frob = 2),
        list(rank_used = 3, lambda_used = 0.3, roi_improv = 0.6, delta_frob = 4)
      )),
      list(diag_by_fold = list(
        list(rank_used = 2, lambda_used = 0.5, roi_improv = 0.8, delta_frob = 6)
      ))
    )
  )
  class(fit_res) <- c("regional_mvpa_result", "list")

  fit_tab <- summarize_remap_roi(fit_res)
  expect_equal(fit_tab$roinum, c(1L, 2L))
  expect_equal(fit_tab$mean_rank, c(2, 2))
  expect_equal(fit_tab$mean_lambda, c(0.2, 0.5))
  expect_equal(fit_tab$mean_roi_improv, c(0.3, 0.8))
  expect_equal(fit_tab$mean_delta_frob, c(3, 6))
})

test_that("summarize_remap_items aggregates predictor diagnostics and errors clearly", {
  pred <- list(diag_by_fold = list(
    list(
      train_items = c("a", "b"),
      item_res_naive = c(4, 8),
      item_res_remap = c(2, 4)
    ),
    list(
      train_items = c("a", "b"),
      item_res_naive = c(6, 10),
      item_res_remap = c(3, 5)
    )
  ))

  items <- summarize_remap_items(pred)
  expect_equal(items$item, c("a", "b"))
  expect_equal(items$res_naive, c(5, 9))
  expect_equal(items$res_remap, c(2.5, 4.5))
  expect_equal(items$n_folds, c(2L, 2L))
  expect_equal(items$res_ratio, c(0.5, 0.5))

  regional <- list(
    performance_table = data.frame(roinum = 42L),
    fits = list(pred)
  )
  class(regional) <- c("regional_mvpa_result", "list")
  expect_equal(summarize_remap_items(regional, roi = 1)$item, c("a", "b"))
  expect_error(summarize_remap_items(regional), "please supply roi")
  expect_error(summarize_remap_items(list()), "unsupported input")
  expect_error(summarize_remap_items(list(diag_by_fold = list())), "no diag_by_fold")
})

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
