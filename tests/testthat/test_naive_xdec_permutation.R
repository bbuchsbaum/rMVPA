# Permutation inference for naive_xdec_model.
#
# naive_xdec_model caches label-derived keys in `.fast_kernel` at
# construction. A permutation run must rebuild that cache from the permuted
# design; otherwise every permutation scores the observed labels and the null
# reproduces the observed map.

naive_xdec_perm_fixture <- function(seed = 7101L, return_predictions = FALSE) {
  set.seed(seed)
  ds <- gen_sample_dataset(
    D = c(5, 5, 5),
    nobs = 48L,
    nlevels = 2L,
    blocks = 4L,
    external_test = TRUE,
    ntest_obs = 48L
  )
  model <- naive_xdec_model(ds$dataset, ds$design,
                            return_predictions = return_predictions)
  list(ds = ds, model = model)
}

# `.permuted_model_spec` is an unregistered internal generic (like
# `.run_searchlight_engine`): call it from the namespace, as package code does.
permuted_spec <- function(spec, design) {
  eval(quote(.permuted_model_spec(spec, design)),
       list(spec = spec, design = design), asNamespace("rMVPA"))
}

test_that(".permuted_model_spec rebuilds the naive_xdec fast kernel", {
  fix <- naive_xdec_perm_fixture()
  model <- fix$model
  expect_true(is.list(model$.fast_kernel))

  perm_design <- permute_labels(model$design, method = "global", seed = 11L)
  expect_false(identical(as.character(perm_design$y_train),
                         as.character(model$design$y_train)))

  perm_spec <- permuted_spec(model, perm_design)

  expect_identical(perm_spec$design, perm_design)
  expect_identical(as.character(perm_spec$.fast_kernel$key_tr),
                   as.character(perm_design$y_train))
  # The unpermuted spec keeps its own kernel.
  expect_identical(as.character(model$.fast_kernel$key_tr),
                   as.character(model$design$y_train))
})

test_that("default .permuted_model_spec only swaps the design", {
  dset <- gen_sample_dataset(c(4, 4, 4), 24, blocks = 3, nlevels = 2)
  mspec <- mvpa_model(load_model("corclass"), dset$dataset, dset$design,
                      "classification",
                      crossval = blocked_cross_validation(dset$design$block_var))
  perm_design <- permute_labels(mspec$design, method = "global", seed = 3L)
  perm_spec <- permuted_spec(mspec, perm_design)

  expected <- mspec
  expected$design <- perm_design
  expect_identical(perm_spec, expected)
})

test_that("naive_xdec permutation null differs from the observed map", {
  skip_on_cran()
  futile.logger::flog.threshold(futile.logger::ERROR)
  on.exit(futile.logger::flog.threshold(futile.logger::INFO), add = TRUE)

  for (rp in c(FALSE, TRUE)) {
    fix <- naive_xdec_perm_fixture(seed = 7102L, return_predictions = rp)
    pc <- permutation_control(n_perm = 3, subsample = 0.5, seed = 5L,
                              shuffle = "global", null_method = "global",
                              correction = "none", diagnose = FALSE)
    res <- suppressWarnings(
      run_permutation_searchlight(fix$model, radius = 3, perm_ctrl = pc)
    )
    res <- if (inherits(res, "permutation_result_set")) res[[1]] else res
    expect_s3_class(res, "permutation_result")

    null_vals <- unlist(res$adj_null$bin_nulls, use.names = FALSE)
    obs_vals <- as.numeric(res$observed$results[[res$metric]]@.Data)
    obs_vals <- obs_vals[is.finite(obs_vals)]
    null_vals <- null_vals[is.finite(null_vals)]
    expect_gt(length(null_vals), 0)

    # Under the stale-kernel bug every null value was an observed value.
    expect_false(all(null_vals %in% obs_vals),
                 info = sprintf("return_predictions = %s", rp))
  }
})

test_that("naive_xdec with an active link_by refuses permutation up front", {
  set.seed(7103L)
  ds <- gen_sample_dataset(c(4, 4, 4), 24, blocks = 3, nlevels = 2,
                           external_test = TRUE, ntest_obs = 24)
  ds$design$train_design$item <- factor(rep(letters[1:6], 4))
  ds$design$test_design$item <- factor(rep(letters[1:6], 4))
  model <- naive_xdec_model(ds$dataset, ds$design, link_by = "item")

  ran_observed <- FALSE
  testthat::local_mocked_bindings(
    run_searchlight = function(...) { ran_observed <<- TRUE; stop("not reached") },
    .package = "rMVPA"
  )
  expect_error(
    run_permutation_searchlight(model, radius = 3,
                                perm_ctrl = permutation_control(n_perm = 2)),
    "link_by = 'item'"
  )
  expect_false(ran_observed)
})
