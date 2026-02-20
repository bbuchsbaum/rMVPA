testthat::skip_if_not_installed("neuroim2")

build_dual_sampled_mspec <- function(D = c(4, 4, 4), nobs = 54, blocks = 6) {
  ds <- gen_sample_dataset(D = D, nobs = nobs, nlevels = 3, blocks = blocks)

  classes <- letters[1:3]
  stopifnot(nobs %% blocks == 0L)
  per_block <- nobs %/% blocks
  stopifnot(per_block %% length(classes) == 0L)

  y_block <- rep(classes, each = per_block %/% length(classes))
  y <- factor(rep(y_block, times = blocks), levels = classes)
  block <- rep(seq_len(blocks), each = per_block)

  des <- mvpa_design(
    data.frame(y = y, block = block),
    y_train = ~ y,
    block_var = ~ block
  )

  cval <- blocked_cross_validation(des$block_var)
  mvpa_model(
    model = load_model("dual_lda"),
    dataset = ds$dataset,
    design = des,
    model_type = "classification",
    crossval = cval,
    tune_grid = data.frame(gamma = 1e-2)
  )
}

with_dual_sampled_options <- function(code) {
  keys <- c(
    "rMVPA.searchlight_mode",
    "rMVPA.warn_legacy_options",
    "rMVPA.dual_lda_rank_update"
  )
  old <- options()[keys]
  on.exit(options(old), add = TRUE)
  options(
    rMVPA.searchlight_mode = "fast",
    rMVPA.warn_legacy_options = FALSE,
    rMVPA.dual_lda_rank_update = FALSE
  )
  force(code)
}

build_non_dual_sampled_mspec <- function(D = c(4, 4, 4), nobs = 54, blocks = 6) {
  ds <- gen_sample_dataset(D = D, nobs = nobs, nlevels = 3, blocks = blocks)

  classes <- letters[1:3]
  y <- factor(rep(classes, length.out = nobs), levels = classes)
  block <- rep(seq_len(blocks), each = ceiling(nobs / blocks))[seq_len(nobs)]

  des <- mvpa_design(
    data.frame(y = y, block = block),
    y_train = ~ y,
    block_var = ~ block
  )

  cval <- blocked_cross_validation(des$block_var)
  mvpa_model(
    model = load_model("sda_notune"),
    dataset = ds$dataset,
    design = des,
    model_type = "classification",
    crossval = cval
  )
}

test_that("dual_lda_fast engine runs randomized searchlight", {
  with_dual_sampled_options({
    set.seed(4101)
    mspec <- build_dual_sampled_mspec(D = c(4, 4, 4), nobs = 54, blocks = 6)
    res <- run_searchlight(
      mspec,
      radius = 2,
      method = "randomized",
      niter = 1,
      backend = "default",
      engine = "dual_lda_fast"
    )

    expect_s3_class(res, "searchlight_result")
    expect_identical(attr(res, "searchlight_engine"), "dual_lda_fast")
    expect_true("Accuracy" %in% names(res$results))
    expect_true("AUC" %in% names(res$results))
  })
})

test_that("explicit dual_lda_fast request errors when engine is ineligible", {
  with_dual_sampled_options({
    set.seed(4103)
    mspec <- build_non_dual_sampled_mspec(D = c(4, 4, 4), nobs = 54, blocks = 6)
    expect_error(
      run_searchlight(
        mspec,
        radius = 2,
        method = "randomized",
        niter = 1,
        backend = "default",
        engine = "dual_lda_fast"
      ),
      regexp = "not eligible"
    )
  })
})

test_that("dual_lda_fast engine runs resampled searchlight", {
  with_dual_sampled_options({
    set.seed(4102)
    mspec <- build_dual_sampled_mspec(D = c(4, 4, 4), nobs = 54, blocks = 6)
    res <- run_searchlight(
      mspec,
      radius = 2,
      method = "resampled",
      niter = 12,
      backend = "default",
      engine = "dual_lda_fast"
    )

    expect_s3_class(res, "searchlight_result")
    expect_identical(attr(res, "searchlight_engine"), "dual_lda_fast")
    expect_true("Accuracy" %in% names(res$results))
    expect_true("AUC" %in% names(res$results))
  })
})
