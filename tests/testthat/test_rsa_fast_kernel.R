rsa_fast_kernel_fixture <- function(regtype = "lm", distmethod = "pearson", semipartial = FALSE) {
  set.seed(8123)
  n_obs <- 24L
  n_vox <- 40L

  # train_mat passed to train_model.rsa_model (rows = observations, cols = voxels)
  train_mat <- matrix(rnorm(n_obs * n_vox), nrow = n_obs, ncol = n_vox)

  # minimal mvpa_dataset-like object for rsa_model constructor checks
  dset <- list(
    train_data = matrix(rnorm(n_vox * n_obs), nrow = n_vox, ncol = n_obs),
    mask = seq_len(n_vox)
  )
  class(dset) <- "mvpa_dataset"

  block <- rep(seq_len(4L), each = n_obs / 4L)
  D1 <- dist(matrix(rnorm(n_obs * 6L), nrow = n_obs, ncol = 6L))
  D2 <- dist(matrix(rnorm(n_obs * 6L), nrow = n_obs, ncol = 6L))
  rdes <- rsa_design(~ D1 + D2, list(D1 = D1, D2 = D2, block = block), block_var = "block")

  mspec <- rsa_model(
    dataset = dset,
    design = rdes,
    distmethod = distmethod,
    regtype = regtype,
    semipartial = semipartial,
    check_collinearity = FALSE
  )

  list(mspec = mspec, train_mat = train_mat)
}

test_that("rsa fast kernel state is constructed by default", {
  fix <- rsa_fast_kernel_fixture(regtype = "lm", distmethod = "pearson", semipartial = FALSE)
  expect_true(is.list(fix$mspec$.fast_kernel))
  expect_true(is.list(fix$mspec$.fast_kernel$lm))

  semipartial_fix <- rsa_fast_kernel_fixture(regtype = "lm", distmethod = "pearson", semipartial = TRUE)
  expect_true(is.null(semipartial_fix$mspec$.fast_kernel) || is.null(semipartial_fix$mspec$.fast_kernel$lm))
})

test_that("rsa fast kernel preserves lm t-values (differential parity)", {
  fix <- rsa_fast_kernel_fixture(regtype = "lm", distmethod = "pearson", semipartial = FALSE)

  old_opt <- options(
    rMVPA.searchlight_mode = "legacy",
    rMVPA.rsa_fast_kernel = NULL
  )
  on.exit(options(old_opt), add = TRUE)
  base <- rMVPA:::train_model.rsa_model(fix$mspec, fix$train_mat, y = NULL, indices = NULL)

  options(rMVPA.searchlight_mode = "fast", rMVPA.rsa_fast_kernel = NULL)
  fast_fix <- rsa_fast_kernel_fixture(regtype = "lm", distmethod = "pearson", semipartial = FALSE)
  fast <- rMVPA:::train_model.rsa_model(fast_fix$mspec, fast_fix$train_mat, y = NULL, indices = NULL)

  expect_identical(names(fast), names(base))
  expect_equal(unname(fast), unname(base), tolerance = 1e-8)
})

test_that("rsa fast kernel preserves correlation outputs for spearman predictor fit", {
  fix <- rsa_fast_kernel_fixture(regtype = "pearson", distmethod = "spearman", semipartial = FALSE)

  old_opt <- options(
    rMVPA.searchlight_mode = "legacy",
    rMVPA.rsa_fast_kernel = NULL
  )
  on.exit(options(old_opt), add = TRUE)
  base <- rMVPA:::train_model.rsa_model(fix$mspec, fix$train_mat, y = NULL, indices = NULL)

  options(rMVPA.searchlight_mode = "fast", rMVPA.rsa_fast_kernel = NULL)
  fast_fix <- rsa_fast_kernel_fixture(regtype = "pearson", distmethod = "spearman", semipartial = FALSE)
  fast <- rMVPA:::train_model.rsa_model(fast_fix$mspec, fast_fix$train_mat, y = NULL, indices = NULL)

  expect_identical(names(fast), names(base))
  expect_equal(unname(fast), unname(base), tolerance = 1e-10)
})

test_that("rsa fast kernel is invariant to predictor ordering (metamorphic)", {
  set.seed(9211)
  n_obs <- 24L
  n_vox <- 36L
  train_mat <- matrix(rnorm(n_obs * n_vox), nrow = n_obs, ncol = n_vox)

  dset <- list(
    train_data = matrix(rnorm(n_vox * n_obs), nrow = n_vox, ncol = n_obs),
    mask = seq_len(n_vox)
  )
  class(dset) <- "mvpa_dataset"

  block <- rep(seq_len(4L), each = n_obs / 4L)
  D1 <- dist(matrix(rnorm(n_obs * 6L), nrow = n_obs, ncol = 6L))
  D2 <- dist(matrix(rnorm(n_obs * 6L), nrow = n_obs, ncol = 6L))
  data_list <- list(D1 = D1, D2 = D2, block = block)

  old_opt <- options(
    rMVPA.searchlight_mode = "fast",
    rMVPA.rsa_fast_kernel = NULL
  )
  on.exit(options(old_opt), add = TRUE)

  rdes_12 <- rsa_design(~ D1 + D2, data_list, block_var = "block")
  rdes_21 <- rsa_design(~ D2 + D1, data_list, block_var = "block")

  mod_12 <- rsa_model(dset, rdes_12, regtype = "lm", distmethod = "pearson", check_collinearity = FALSE)
  mod_21 <- rsa_model(dset, rdes_21, regtype = "lm", distmethod = "pearson", check_collinearity = FALSE)

  res_12 <- rMVPA:::train_model.rsa_model(mod_12, train_mat, y = NULL, indices = NULL)
  res_21 <- rMVPA:::train_model.rsa_model(mod_21, train_mat, y = NULL, indices = NULL)

  expect_setequal(names(res_12), names(res_21))
  expect_equal(res_12[names(res_21)], res_21, tolerance = 1e-8)
})
