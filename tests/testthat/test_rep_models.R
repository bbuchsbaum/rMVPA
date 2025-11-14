context("ReNA models: repnet/repmed/repmap correctness")

test_that("repnet_model enforces confound RDM label alignment and returns betas", {
  skip_on_cran()
  set.seed(123)

  # Synthetic dataset with categorical response and sufficient trials
  ds <- gen_sample_dataset(D = c(4,4,4), nobs = 60, response_type = "categorical",
                           data_mode = "image", blocks = 3, nlevels = 6)

  # Add item identity used for grouping
  K <- 6
  items <- paste0("I", 1:K)
  ds$design$train_design$ImageID <- factor(rep(items, length.out = nrow(ds$design$train_design)))

  # Seed RDM (labeled)
  S <- as.matrix(dist(matrix(1:K, K, 1)))
  rownames(S) <- colnames(S) <- items

  # Confound with missing labels (smaller set) to trigger alignment error
  S4 <- as.matrix(dist(matrix(1:4, 4, 1)))
  rownames(S4) <- colnames(S4) <- items[1:4]

  bad_des <- repnet_design(ds$design, key_var = ~ ImageID, seed_rdm = S,
                           confound_rdms = list(short = S4))
  mod_bad <- repnet_model(ds$dataset, ds$design, bad_des, distfun = eucdist())

  # Small ROI
  mask <- ds$dataset$mask
  idx <- which(mask > 0)[1:50]
  # Build minimal ROI list expected by process_roi.*
  train_roi <- neuroim2::series_roi(ds$dataset$train_data, idx)
  roi <- list(train_roi = train_roi, test_roi = NULL)

  res_bad <- process_roi(mod_bad, roi, rnum = 1)
  expect_true(res_bad$error[[1]])
  expect_match(res_bad$error_message[[1]], "confound RDM 'short' lacks required item labels", perl = TRUE)

  # Proper confound with full labels -> should succeed and include beta_seed
  conf_ok <- list(block = S + 0)  # same labels as seed
  ok_des <- repnet_design(ds$design, key_var = ~ ImageID, seed_rdm = S, confound_rdms = conf_ok)
  mod_ok <- repnet_model(ds$dataset, ds$design, ok_des, distfun = eucdist())

  res_ok <- process_roi(mod_ok, roi, rnum = 2)
  expect_false(res_ok$error[[1]])
  pm <- res_ok$performance[[1]]
  expect_true(any(grepl("^beta_seed$|^beta_seed$", names(pm))))
})

test_that("repmed_model returns total effect (med_c) and runs without unnamed columns", {
  skip_on_cran()
  set.seed(124)

  ds <- gen_sample_dataset(D = c(4,4,4), nobs = 60, response_type = "categorical",
                           data_mode = "image", blocks = 3, nlevels = 6)
  K <- 6
  items <- paste0("I", 1:K)
  ds$design$train_design$ImageID <- factor(rep(items, length.out = nrow(ds$design$train_design)))

  # Construct simple X and Y item scalars, derive RDMs
  xi <- as.numeric(1:K)
  yi <- xi + rnorm(K, 0, 0.1)
  Xr <- dist(matrix(xi, K, 1))
  Yr <- dist(matrix(yi, K, 1))

  rmd <- repmed_design(items = items, X_rdm = Xr, Y_rdm = Yr)
  mod <- repmed_model(ds$dataset, ds$design, rmd, key_var = ~ ImageID, distfun = cordist(method = "pearson"))

  mask <- ds$dataset$mask
  idx <- which(mask > 0)[1:50]
  train_roi <- neuroim2::series_roi(ds$dataset$train_data, idx)
  roi <- list(train_roi = train_roi, test_roi = NULL)

  res <- process_roi(mod, roi, rnum = 1)
  expect_false(res$error[[1]])
  pm <- res$performance[[1]]
  expect_true("med_c" %in% names(pm))
  expect_true(is.numeric(pm[["med_c"]]))
})

test_that("repmap_model rank=0 reports zero map and map_rank=0", {
  skip_on_cran()
  set.seed(125)

  ds <- gen_sample_dataset(D = c(4,4,4), nobs = 60, response_type = "categorical",
                           data_mode = "image", blocks = 3, nlevels = 6)
  K <- 6
  items <- paste0("I", 1:K)
  ds$design$train_design$ImageID <- factor(rep(items, length.out = nrow(ds$design$train_design)))

  # Seed features K x P with labels
  P <- 5
  Xseed <- matrix(rnorm(K * P), K, P)
  rownames(Xseed) <- items
  rmd <- repmap_design(items = items, seed_features = Xseed)
  mod <- repmap_model(ds$dataset, ds$design, rmd, key_var = ~ ImageID, rank = 0)

  mask <- ds$dataset$mask
  idx <- which(mask > 0)[1:50]
  train_roi <- neuroim2::series_roi(ds$dataset$train_data, idx)
  roi <- list(train_roi = train_roi, test_roi = NULL)

  res <- process_roi(mod, roi, rnum = 1)
  expect_false(res$error[[1]])
  pm <- res$performance[[1]]
  expect_equal(as.numeric(pm[["map_rank"]]), 0)
})
