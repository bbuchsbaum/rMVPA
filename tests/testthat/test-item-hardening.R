library(testthat)
library(rMVPA)

futile.logger::flog.threshold(futile.logger::ERROR)

skip_if_not_installed("fmrilss")

.item_hash <- function(trial_id) {
  getFromNamespace(".item_simple_hash", "fmrilss")(trial_id)
}

.drop_item_diag <- function(x) {
  attr(x, "item_diagnostics") <- NULL
  x
}

test_that("item_compute_u contract supports dense/sparse covariance and precision", {
  set.seed(4101)

  n_time <- 80
  n_trials <- 24
  X_t <- matrix(rnorm(n_time * n_trials), nrow = n_time, ncol = n_trials)

  idx <- seq_len(n_time)
  rho <- 0.20
  V_cov <- rho^abs(outer(idx, idx, "-"))
  V_prec <- solve(V_cov)
  V_prec_sparse <- Matrix::Matrix(V_prec, sparse = TRUE)

  U_cov <- fmrilss::item_compute_u(
    X_t = X_t,
    V = V_cov,
    v_type = "cov",
    method = "svd",
    output = "matrix"
  )
  U_prec <- fmrilss::item_compute_u(
    X_t = X_t,
    V = V_prec,
    v_type = "precision",
    method = "svd",
    output = "matrix"
  )
  U_prec_sparse <- fmrilss::item_compute_u(
    X_t = X_t,
    V = V_prec_sparse,
    v_type = "precision",
    method = "svd",
    output = "matrix"
  )

  expect_equal(dim(U_cov), c(n_trials, n_trials))
  expect_lt(max(abs(U_cov - t(U_cov))), 1e-8)
  evals <- eigen((U_cov + t(U_cov)) / 2, symmetric = TRUE, only.values = TRUE)$values
  expect_gt(min(evals), -1e-7)

  expect_equal(.drop_item_diag(U_cov), .drop_item_diag(U_prec), tolerance = 1e-6)
  expect_equal(.drop_item_diag(U_cov), .drop_item_diag(U_prec_sparse), tolerance = 1e-6)
})

test_that("item_slice_fold submatrices and by-run slicing are consistent", {
  set.seed(4102)

  n_time <- 72
  n_trials <- 24
  run_id <- rep(c(1, 2, 3, 4), each = 6)

  X_t <- matrix(rnorm(n_time * n_trials), nrow = n_time, ncol = n_trials)
  Gamma <- matrix(rnorm(n_trials * 7), nrow = n_trials, ncol = 7)
  T_target <- factor(rep(c("A", "B"), length.out = n_trials), levels = c("A", "B"))

  U_full <- fmrilss::item_compute_u(
    X_t = X_t,
    run_id = run_id,
    method = "svd",
    output = "matrix"
  )
  U_by_run <- fmrilss::item_compute_u(
    X_t = X_t,
    run_id = run_id,
    method = "svd",
    output = "by_run"
  )

  run_levels <- sort(unique(run_id))
  expect_equal(length(U_by_run), length(run_levels))
  for (i in seq_along(run_levels)) {
    expect_equal(dim(U_by_run[[i]]), c(6, 6))
  }

  bundle <- fmrilss::item_build_design(
    X_t = X_t,
    T_target = T_target,
    run_id = run_id,
    validate = TRUE
  )
  bundle$Gamma <- Gamma
  bundle$U <- U_full
  bundle$U_by_run <- NULL

  fold <- fmrilss::item_slice_fold(bundle, test_run = "2", check_hash = FALSE)
  expect_equal(fold$U_test, U_full[fold$test_idx, fold$test_idx], tolerance = 1e-10)
  expect_equal(fold$U_train, U_full[fold$train_idx, fold$train_idx], tolerance = 1e-10)
  expect_equal(length(intersect(fold$train_idx, fold$test_idx)), 0L)
  expect_equal(nrow(fold$Gamma_train) + nrow(fold$Gamma_test), n_trials)

  bundle$U <- NULL
  bundle$U_by_run <- U_by_run
  fold_by_run <- fmrilss::item_slice_fold(bundle, test_run = "2", check_hash = FALSE)
  expect_equal(fold_by_run$U_test, U_by_run[[2]], tolerance = 1e-10)
  run_train <- run_id[fold_by_run$train_idx]
  off_block <- outer(run_train, run_train, "!=")
  expect_lt(max(abs(fold_by_run$U_train[off_block])), 1e-12)
})

test_that("item_cv has deterministic fold order and deterministic tie-breaks", {
  set.seed(4103)

  n_trials <- 12
  run_id <- c(rep(10, 4), rep(2, 4), rep(7, 4))
  U <- diag(n_trials)
  Gamma_tie <- matrix(0, nrow = n_trials, ncol = 5)
  T_target <- factor(rep(c("A", "B"), length.out = n_trials), levels = c("A", "B"))

  res1 <- suppressWarnings(
    fmrilss::item_cv(
      Gamma = Gamma_tie,
      T_target = T_target,
      U = U,
      run_id = run_id,
      mode = "classification",
      method = "svd"
    )
  )
  res2 <- suppressWarnings(
    fmrilss::item_cv(
      Gamma = Gamma_tie,
      T_target = T_target,
      U = U,
      run_id = run_id,
      mode = "classification",
      method = "svd"
    )
  )

  expect_identical(res1$diagnostics$fold_order, as.character(sort(unique(run_id))))
  expect_identical(res1$predictions$predicted_class, res2$predictions$predicted_class)
  expect_true(all(res1$predictions$predicted_class == levels(T_target)[1]))
})

test_that("item_model check_hash enforces trial alignment guard", {
  set.seed(4104)

  n_time <- 48
  n_trials <- 12
  run_id <- rep(1:3, each = 4)
  trial_id <- paste0("trial_", seq_len(n_trials))
  good_hash <- .item_hash(trial_id)

  ds <- gen_sample_dataset(c(4, 4, 4), nobs = n_time, blocks = 3, nlevels = 2)
  X_t <- matrix(rnorm(n_time * n_trials), nrow = n_time, ncol = n_trials)
  T_target <- as.numeric(scale(rnorm(n_trials)))

  des_ok <- item_design(
    train_design = ds$design$train_design,
    X_t = X_t,
    T_target = T_target,
    run_id = run_id,
    trial_id = trial_id,
    trial_hash = good_hash
  )
  mspec_ok <- item_model(
    dataset = ds$dataset,
    design = des_ok,
    mode = "regression",
    metric = "correlation",
    check_hash = TRUE,
    solver = "svd",
    u_storage = "matrix"
  )

  vox <- sample(which(ds$dataset$mask > 0), 18)
  roi <- as_roi(data_sample(ds$dataset, vox), ds$dataset)
  res_ok <- process_roi(mspec_ok, roi, 1)
  expect_false(res_ok$error)

  des_bad <- item_design(
    train_design = ds$design$train_design,
    X_t = X_t,
    T_target = T_target,
    run_id = run_id,
    trial_id = trial_id,
    trial_hash = "definitely_wrong_hash"
  )
  mspec_bad <- item_model(
    dataset = ds$dataset,
    design = des_bad,
    mode = "regression",
    metric = "correlation",
    check_hash = TRUE,
    solver = "svd",
    u_storage = "matrix"
  )
  res_bad <- process_roi(mspec_bad, roi, 2)
  expect_true(res_bad$error)
  expect_match(res_bad$error_message, "Trial hash mismatch")
})

test_that("collinearity triggers solver fallback warning path", {
  set.seed(4105)

  n_trials <- 12
  run_id <- rep(1:3, each = 4)
  U <- diag(n_trials)

  Gamma <- matrix(0, nrow = n_trials, ncol = 4)
  Gamma[, 1] <- 1
  Gamma[, 2] <- Gamma[, 1]
  Gamma[, 3] <- rep(c(1, -1), length.out = n_trials)
  Gamma[, 4] <- Gamma[, 3]

  T_target <- factor(rep(c("A", "B"), length.out = n_trials), levels = c("A", "B"))

  expect_warning(
    out <- fmrilss::item_cv(
      Gamma = Gamma,
      T_target = T_target,
      U = U,
      run_id = run_id,
      mode = "classification",
      ridge = 0,
      method = "chol"
    ),
    "requested 'chol' but used 'svd'"
  )
  expect_s3_class(out, "item_cv_result")
})

test_that("null simulations stay near chance and near-zero regression correlation", {
  set.seed(4106)

  reps <- 20
  n_trials <- 120
  n_vox <- 25
  run_id <- rep(1:4, each = n_trials / 4)
  U <- diag(n_trials)

  acc <- numeric(reps)
  cors <- numeric(reps)

  for (r in seq_len(reps)) {
    Gamma <- matrix(rnorm(n_trials * n_vox), nrow = n_trials, ncol = n_vox)
    y_cls <- factor(sample(c("A", "B"), n_trials, replace = TRUE), levels = c("A", "B"))
    y_reg <- rnorm(n_trials)

    acc[r] <- fmrilss::item_cv(
      Gamma = Gamma,
      T_target = y_cls,
      U = U,
      run_id = run_id,
      mode = "classification",
      method = "svd"
    )$aggregate$mean

    cors[r] <- fmrilss::item_cv(
      Gamma = Gamma,
      T_target = y_reg,
      U = U,
      run_id = run_id,
      mode = "regression",
      method = "svd"
    )$aggregate$mean
  }

  expect_lt(abs(mean(acc) - 0.5), 0.10)
  expect_lt(abs(mean(cors)), 0.12)
})

test_that("signal simulations improve monotonically with SNR", {
  set.seed(4107)

  reps <- 15
  n_trials <- 120
  n_vox <- 20
  run_id <- rep(1:4, each = n_trials / 4)
  U <- diag(n_trials)
  snr_levels <- c(0.00, 0.05, 0.10, 0.20)

  cls_means <- numeric(length(snr_levels))
  reg_means <- numeric(length(snr_levels))

  for (i in seq_along(snr_levels)) {
    snr <- snr_levels[[i]]
    cls_scores <- numeric(reps)
    reg_scores <- numeric(reps)

    for (r in seq_len(reps)) {
      cls <- factor(rep(c("A", "B"), each = n_trials / 2), levels = c("A", "B"))
      cls_code <- ifelse(cls == "A", -1, 1)
      w_cls <- rnorm(n_vox)
      Gamma_cls <- snr * outer(cls_code, w_cls) + matrix(rnorm(n_trials * n_vox), n_trials, n_vox)

      cls_scores[[r]] <- fmrilss::item_cv(
        Gamma = Gamma_cls,
        T_target = cls,
        U = U,
        run_id = run_id,
        mode = "classification",
        method = "svd"
      )$aggregate$mean

      y_reg <- as.numeric(scale(rnorm(n_trials)))
      w_reg <- rnorm(n_vox)
      Gamma_reg <- snr * outer(y_reg, w_reg) + matrix(rnorm(n_trials * n_vox), n_trials, n_vox)

      reg_scores[[r]] <- fmrilss::item_cv(
        Gamma = Gamma_reg,
        T_target = y_reg,
        U = U,
        run_id = run_id,
        mode = "regression",
        method = "svd"
      )$aggregate$mean
    }

    cls_means[[i]] <- mean(cls_scores)
    reg_means[[i]] <- mean(reg_scores)
  }

  expect_true(all(diff(cls_means) >= -0.02))
  expect_true(all(diff(reg_means) >= -0.03))
  expect_gt(tail(cls_means, 1) - cls_means[[1]], 0.15)
  expect_gt(tail(reg_means, 1) - reg_means[[1]], 0.20)
})
