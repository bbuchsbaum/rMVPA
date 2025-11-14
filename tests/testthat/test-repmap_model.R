context("ReNA-Map: repmap_design and repmap_model")

library(testthat)
library(rMVPA)

set.seed(123)

test_that("repmap_design stores items and features", {
  items <- letters[1:6]
  F <- matrix(rnorm(6 * 4), nrow = 6, ncol = 4)
  rownames(F) <- items

  des <- repmap_design(items = items, seed_features = F)
  expect_equal(des$items, items)
  expect_true(is.matrix(des$seed_features))
  expect_equal(rownames(des$seed_features), items)
})


test_that("repmap_model + process_roi.repmap_model return mapping metrics", {
  ds <- gen_sample_dataset(D = c(4,4,4), nobs = 20, blocks = 4, nlevels = 2)
  mvdes <- mvpa_design(
    train_design = ds$design$train_design,
    y_train      = ds$design$y_train,
    block_var    = ds$design$block_var
  )

  items <- as.character(unique(mvdes$train_design$.rownum))
  F <- matrix(rnorm(length(items) * 5), nrow = length(items), ncol = 5)
  rownames(F) <- items

  mp_des <- repmap_design(items = items, seed_features = F)

  mspec <- repmap_model(
    dataset    = ds$dataset,
    design     = mvdes,
    repmap_des = mp_des,
    key_var    = ~ .rownum,
    rank       = 2,
    max_rank   = 3
  )

  vox <- which(ds$dataset$mask > 0)
  samp <- data_sample(ds$dataset, vox)
  roi  <- as_roi(samp, ds$dataset)

  res <- process_roi(mspec, roi, rnum = 1L, center_global_id = NA)

  expect_s3_class(res, "tbl_df")
  expect_false(res$error)
  perf <- res$performance[[1]]
  expect_true(all(c("n_items", "n_seed_feats", "n_vox",
                    "map_rank", "map_r2_mean", "map_frob_norm") %in% names(perf)))
})


## ----------------------------------------------------------------------
## Synthetic mapping tests
## ----------------------------------------------------------------------

simulate_repmap_data <- function(K = 32,
                                 P = 10,
                                 V = 30,
                                 rank_true = 3,
                                 snr = 5) {

  items <- paste0("it", seq_len(K))
  X <- matrix(rnorm(K * P), nrow = K, ncol = P)
  rownames(X) <- items

  U <- matrix(rnorm(P * rank_true), nrow = P, ncol = rank_true)
  Vmat <- matrix(rnorm(V * rank_true), nrow = V, ncol = rank_true)
  B_true <- U %*% t(Vmat)  ## P x V

  signal <- X %*% B_true
  noise <- matrix(rnorm(K * V), nrow = K, ncol = V)
  Y <- signal + noise / snr
  rownames(Y) <- items

  list(
    items   = items,
    X       = X,
    Y       = Y,
    B_true  = B_true,
    rank_true = rank_true
  )
}


test_that("repmap_model: known low-rank mapping is recovered (C1)", {
  skip_on_cran()
  skip_if_not_installed("rrpack")

  sim <- simulate_repmap_data(K = 40, P = 12, V = 24, rank_true = 3, snr = 6)
  X <- sim$X
  Y <- sim$Y

  Xc <- scale(X, center = TRUE, scale = FALSE)
  Yc <- scale(Y, center = TRUE, scale = FALSE)

  fit_full <- rMVPA:::.repmap_fit_rrr(
    Xc,
    Yc,
    rank     = "auto",
    max_rank = 5
  )

  expect_true(is.numeric(fit_full$rank))
  expect_gte(fit_full$rank, 1)
  expect_lte(fit_full$rank, 5)

  Y_hat_full <- Xc %*% fit_full$C
  ss_tot <- colSums(Yc^2)
  ss_res_full <- colSums((Yc - Y_hat_full)^2)
  ss_tot[ss_tot < .Machine$double.eps] <- .Machine$double.eps
  r2_full <- mean(1 - ss_res_full / ss_tot)

  fit_low <- rMVPA:::.repmap_fit_rrr(
    Xc,
    Yc,
    rank     = 1,
    max_rank = 1
  )
  Y_hat_low <- Xc %*% fit_low$C
  ss_res_low <- colSums((Yc - Y_hat_low)^2)
  r2_low <- mean(1 - ss_res_low / ss_tot)

  expect_gt(r2_full, r2_low + 0.05)
  expect_gt(r2_full, 0.2)
})


test_that("repmap_model: pure noise mapping does not hallucinate structure (C2)", {
  skip_on_cran()
  skip_if_not_installed("rrpack")

  K <- 40
  P <- 12
  V <- 24

  items <- paste0("it", seq_len(K))
  X <- matrix(rnorm(K * P), nrow = K, ncol = P)
  Y <- matrix(rnorm(K * V), nrow = K, ncol = V)
  rownames(X) <- rownames(Y) <- items

  Xc <- scale(X, center = TRUE, scale = FALSE)
  Yc <- scale(Y, center = TRUE, scale = FALSE)

  fit_noise <- rMVPA:::.repmap_fit_rrr(
    Xc,
    Yc,
    rank     = "auto",
    max_rank = 5
  )

  Y_hat <- Xc %*% fit_noise$C
  ss_tot <- colSums(Yc^2)
  ss_res <- colSums((Yc - Y_hat)^2)
  ss_tot[ss_tot < .Machine$double.eps] <- .Machine$double.eps
  r2_noise <- mean(1 - ss_res / ss_tot)

  expect_lt(abs(r2_noise), 0.1)
})


test_that("repmap_model: mapping reflects only used feature subspace (C3)", {
  skip_on_cran()
  skip_if_not_installed("rrpack")

  K <- 36
  P1 <- 6
  P2 <- 6
  P <- P1 + P2
  V <- 20

  items <- paste0("it", seq_len(K))
  X1 <- matrix(rnorm(K * P1), nrow = K, ncol = P1)
  X2 <- matrix(rnorm(K * P2), nrow = K, ncol = P2)
  X_all <- cbind(X1, X2)
  rownames(X_all) <- items
  rownames(X2) <- items

  B1 <- matrix(rnorm(P1 * V), nrow = P1, ncol = V)
  signal <- X1 %*% B1
  noise <- matrix(rnorm(K * V), nrow = K, ncol = V)
  Y <- signal + noise / 5
  rownames(Y) <- items

  ## Mapping using all features
  Xc_all <- scale(X_all, center = TRUE, scale = FALSE)
  Yc <- scale(Y, center = TRUE, scale = FALSE)
  fit_all <- rMVPA:::.repmap_fit_rrr(
    Xc_all,
    Yc,
    rank     = "auto",
    max_rank = 6
  )
  Y_hat_all <- Xc_all %*% fit_all$C
  ss_tot <- colSums(Yc^2)
  ss_res_all <- colSums((Yc - Y_hat_all)^2)
  ss_tot[ss_tot < .Machine$double.eps] <- .Machine$double.eps
  r2_all <- mean(1 - ss_res_all / ss_tot)

  ## Mapping using only unused features X2
  Xc_2 <- scale(X2, center = TRUE, scale = FALSE)
  fit_2 <- rMVPA:::.repmap_fit_rrr(
    Xc_2,
    Yc,
    rank     = "auto",
    max_rank = 6
  )
  Y_hat_2 <- Xc_2 %*% fit_2$C
  ss_res_2 <- colSums((Yc - Y_hat_2)^2)
  r2_2 <- mean(1 - ss_res_2 / ss_tot)

  expect_gt(r2_all, 0.2)
  expect_lt(abs(r2_2), 0.1)
  expect_gt(r2_all, r2_2 + 0.1)
})


test_that("repmap_model: orthogonal rotation of seed space preserves mapping quality (C4)", {
  skip_on_cran()
  skip_if_not_installed("rrpack")

  sim <- simulate_repmap_data(K = 40, P = 10, V = 20, rank_true = 3, snr = 6)
  X <- sim$X
  Y <- sim$Y

  Xc <- scale(X, center = TRUE, scale = FALSE)
  Yc <- scale(Y, center = TRUE, scale = FALSE)

  ## Random orthogonal rotation in feature space
  P <- ncol(X)
  Rmat <- qr.Q(qr(matrix(rnorm(P * P), nrow = P, ncol = P)))
  X_rot <- X %*% Rmat
  Xc_rot <- scale(X_rot, center = TRUE, scale = FALSE)

  rank_used <- 3
  fit_orig <- rMVPA:::.repmap_fit_rrr(
    Xc,
    Yc,
    rank     = rank_used,
    max_rank = rank_used
  )
  fit_rot <- rMVPA:::.repmap_fit_rrr(
    Xc_rot,
    Yc,
    rank     = rank_used,
    max_rank = rank_used
  )

  Y_hat_orig <- Xc %*% fit_orig$C
  Y_hat_rot  <- Xc_rot %*% fit_rot$C

  ss_tot <- colSums(Yc^2)
  ss_tot[ss_tot < .Machine$double.eps] <- .Machine$double.eps

  r2_orig <- mean(1 - colSums((Yc - Y_hat_orig)^2) / ss_tot)
  r2_rot  <- mean(1 - colSums((Yc - Y_hat_rot)^2) / ss_tot)

  expect_true(is.finite(r2_orig))
  expect_true(is.finite(r2_rot))
  expect_gt(r2_orig, 0.2)
  expect_gt(r2_rot, 0.2)
  expect_lt(abs(r2_orig - r2_rot), 0.05)
})


test_that("repmap_model: train-test generalization behaves sensibly (C5)", {
  skip_on_cran()
  skip_if_not_installed("rrpack")

  K <- 60
  P <- 10
  V <- 20
  rank_true <- 3

  items <- paste0("it", seq_len(K))
  X <- matrix(rnorm(K * P), nrow = K, ncol = P)
  rownames(X) <- items

  U <- matrix(rnorm(P * rank_true), nrow = P, ncol = rank_true)
  Vmat <- matrix(rnorm(V * rank_true), nrow = V, ncol = rank_true)
  B_true <- U %*% t(Vmat)

  signal <- X %*% B_true
  noise <- matrix(rnorm(K * V), nrow = K, ncol = V)
  Y <- signal + noise / 5
  rownames(Y) <- items

  set.seed(123)
  train_idx <- sort(sample(seq_len(K), size = floor(0.7 * K)))
  test_idx  <- setdiff(seq_len(K), train_idx)

  X_train <- X[train_idx, , drop = FALSE]
  Y_train <- Y[train_idx, , drop = FALSE]
  X_test  <- X[test_idx, , drop = FALSE]
  Y_test  <- Y[test_idx, , drop = FALSE]

  Xc_train <- scale(X_train, center = TRUE, scale = FALSE)
  Yc_train <- scale(Y_train, center = TRUE, scale = FALSE)
  x_mean <- attr(Xc_train, "scaled:center")
  y_mean <- attr(Yc_train, "scaled:center")

  fit <- rMVPA:::.repmap_fit_rrr(
    Xc_train,
    Yc_train,
    rank     = rank_true,
    max_rank = rank_true
  )

  ## In-sample R2
  Y_hat_train <- Xc_train %*% fit$C
  ss_tot_train <- colSums(Yc_train^2)
  ss_tot_train[ss_tot_train < .Machine$double.eps] <- .Machine$double.eps
  ss_res_train <- colSums((Yc_train - Y_hat_train)^2)
  r2_train <- mean(1 - ss_res_train / ss_tot_train)

  ## Out-of-sample R2 (using train centering)
  Xc_test <- sweep(X_test, 2, x_mean, "-")
  Yc_test <- sweep(Y_test, 2, y_mean, "-")
  Y_hat_test <- Xc_test %*% fit$C
  ss_tot_test <- colSums(Yc_test^2)
  ss_tot_test[ss_tot_test < .Machine$double.eps] <- .Machine$double.eps
  ss_res_test <- colSums((Yc_test - Y_hat_test)^2)
  r2_test <- mean(1 - ss_res_test / ss_tot_test)

  expect_gt(r2_train, r2_test)
  expect_gt(r2_test, 0.1)
})
