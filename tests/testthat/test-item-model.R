library(testthat)
library(rMVPA)

futile.logger::flog.threshold(futile.logger::ERROR)

skip_if_not_installed("fmrilss")

test_that("item_design validates dimensions and stores ITEM metadata", {
  set.seed(9001)

  n_time <- 40
  n_trials <- 12
  train_design <- data.frame(tr = seq_len(n_time))
  X_t <- matrix(rnorm(n_time * n_trials), nrow = n_time, ncol = n_trials)
  run_id <- rep(1:3, each = 4)
  T_target <- rep(c("A", "B"), length.out = n_trials)

  des <- item_design(
    train_design = train_design,
    X_t = X_t,
    T_target = T_target,
    run_id = run_id
  )

  expect_s3_class(des, "item_design")
  expect_s3_class(des, "mvpa_design")
  expect_equal(dim(des$X_t), c(n_time, n_trials))
  expect_equal(length(des$run_id), n_trials)
  expect_true(inherits(des$item_bundle, "item_bundle"))
  expect_equal(nrow(des$T_target), n_trials)

  expect_error(
    item_design(
      train_design = train_design,
      X_t = X_t,
      T_target = T_target,
      run_id = rep(1:3, each = 3)
    ),
    "run_id must have length"
  )
})

test_that("item_model fit_roi path returns scalar metrics", {
  set.seed(9002)

  n_time <- 48
  n_trials <- 12

  ds <- gen_sample_dataset(c(4, 4, 4), nobs = n_time, blocks = 3, nlevels = 2)

  X_t <- matrix(rnorm(n_time * n_trials), nrow = n_time, ncol = n_trials)
  run_id <- rep(1:3, each = 4)
  T_target <- as.numeric(scale(rnorm(n_trials)))

  des <- item_design(
    train_design = ds$design$train_design,
    X_t = X_t,
    T_target = T_target,
    run_id = run_id
  )

  mspec <- item_model(
    dataset = ds$dataset,
    design = des,
    mode = "regression",
    metric = "correlation",
    solver = "svd",
    u_storage = "by_run"
  )

  vox <- sample(which(ds$dataset$mask > 0), 16)
  samp <- data_sample(ds$dataset, vox)
  roi <- as_roi(samp, ds$dataset)

  res <- process_roi(mspec, roi, 1)
  expect_s3_class(res, "tbl_df")
  expect_false(res$error)
  expect_true(is.numeric(res$performance[[1]]))
  expect_equal(
    sort(names(res$performance[[1]])),
    sort(c("item_score_mean", "item_score_sd", "item_n_folds"))
  )
})

test_that("item_model integrates with run_searchlight via schema combiner", {
  set.seed(9003)

  n_time <- 40
  n_trials <- 12

  ds <- gen_sample_dataset(c(3, 3, 3), nobs = n_time, blocks = 3, nlevels = 2)
  X_t <- matrix(rnorm(n_time * n_trials), nrow = n_time, ncol = n_trials)
  run_id <- rep(1:3, each = 4)
  T_target <- as.numeric(scale(rnorm(n_trials)))

  des <- item_design(
    train_design = ds$design$train_design,
    X_t = X_t,
    T_target = T_target,
    run_id = run_id
  )

  mspec <- item_model(
    dataset = ds$dataset,
    design = des,
    mode = "regression",
    metric = "correlation",
    solver = "svd",
    u_storage = "matrix",
    return_predictions = FALSE
  )

  schema <- output_schema(mspec)
  expect_false(is.null(schema))

  out <- run_searchlight(mspec, radius = 1, method = "standard")
  expect_s3_class(out, "searchlight_result")
  expect_equal(sort(names(out$results)), sort(names(schema)))
})
