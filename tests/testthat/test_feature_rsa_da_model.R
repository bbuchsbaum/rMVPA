context("feature_rsa_da_model")

.make_frda_dataset <- function(Tenc = 24, Trec = 18, dims = c(3, 3, 2), Dfeat = 8, seed = 101) {
  set.seed(seed)
  toy <- gen_sample_dataset(
    D = dims,
    nobs = Tenc,
    nlevels = 2,
    blocks = 2,
    external_test = TRUE,
    ntest_obs = Trec
  )

  Xenc <- matrix(rnorm(Tenc * Dfeat), Tenc, Dfeat)
  fs_enc <- feature_sets(Xenc, blocks(a = Dfeat / 2, b = Dfeat / 2))

  G <- matrix(runif(Trec * Tenc), Trec, Tenc)
  G <- G / rowSums(G)
  fs_rec <- expected_features(fs_enc, G, drop_null = FALSE, renormalize = TRUE)

  des <- feature_sets_design(fs_enc, fs_rec, block_var_test = rep(1:2, length.out = Trec))
  mask <- neuroim2::NeuroVol(array(1, dims), neuroim2::space(toy$dataset$mask))

  list(dataset = toy$dataset, design = des, mask = mask)
}

.make_frda_linear_dataset <- function(Tenc = 24, Trec = 18, dims = c(3, 3, 2), set_sizes = c(a = 4, b = 4), seed = 301) {
  set.seed(seed)
  Dfeat <- sum(set_sizes)
  V <- prod(dims)

  X_enc <- matrix(rnorm(Tenc * Dfeat), Tenc, Dfeat)
  X_rec <- matrix(rnorm(Trec * Dfeat), Trec, Dfeat)
  B <- matrix(rnorm(Dfeat * V), nrow = Dfeat, ncol = V)

  Y_enc <- X_enc %*% B
  Y_rec <- X_rec %*% B

  arr_enc <- array(0, dim = c(dims, Tenc))
  for (t in seq_len(Tenc)) arr_enc[,,, t] <- array(Y_enc[t, ], dim = dims)
  arr_rec <- array(0, dim = c(dims, Trec))
  for (t in seq_len(Trec)) arr_rec[,,, t] <- array(Y_rec[t, ], dim = dims)

  mask <- neuroim2::NeuroVol(array(1, dim = dims), neuroim2::NeuroSpace(dims))
  train_vec <- neuroim2::NeuroVec(arr_enc, neuroim2::NeuroSpace(c(dims, Tenc)))
  test_vec <- neuroim2::NeuroVec(arr_rec, neuroim2::NeuroSpace(c(dims, Trec)))
  dataset <- mvpa_dataset(train_vec, test_vec, mask)

  spec <- do.call(blocks, as.list(set_sizes))
  fs_enc <- feature_sets(X_enc, spec)
  fs_rec <- feature_sets(X_rec, spec, set_order = fs_enc$set_order)
  des <- feature_sets_design(fs_enc, fs_rec, block_var_test = rep(1:2, length.out = Trec))

  list(dataset = dataset, design = des, mask = mask)
}

.independent_geom_metrics <- function(observed, predicted, rsa_simfun = "spearman") {
  observed <- as.matrix(observed)
  predicted <- as.matrix(predicted)

  out <- c(
    pattern_correlation = NA_real_,
    pattern_discrimination = NA_real_,
    pattern_rank_percentile = NA_real_,
    rdm_correlation = NA_real_,
    voxel_correlation = NA_real_,
    mean_voxelwise_temporal_cor = NA_real_
  )

  if (!all(dim(observed) == dim(predicted))) return(out)

  sd_thresh <- 1e-12
  obs_sd <- apply(observed, 2, stats::sd)
  pred_sd <- apply(predicted, 2, stats::sd)
  valid_col <- which(obs_sd > sd_thresh & pred_sd > sd_thresh)
  if (length(valid_col) < 1L) return(out)

  obs_use <- observed[, valid_col, drop = FALSE]
  pred_use <- predicted[, valid_col, drop = FALSE]

  out["voxel_correlation"] <- tryCatch(
    stats::cor(c(pred_use), c(obs_use)),
    error = function(e) NA_real_
  )

  if (nrow(observed) > 1L) {
    vv <- vapply(seq_len(ncol(obs_use)), function(j) {
      tryCatch(stats::cor(obs_use[, j], pred_use[, j]), error = function(e) NA_real_)
    }, numeric(1))
    out["mean_voxelwise_temporal_cor"] <- mean(vv, na.rm = TRUE)
    if (is.nan(out["mean_voxelwise_temporal_cor"])) out["mean_voxelwise_temporal_cor"] <- NA_real_
  }

  obs_row_sd <- apply(obs_use, 1, stats::sd)
  pred_row_sd <- apply(pred_use, 1, stats::sd)
  valid_row <- which(obs_row_sd > sd_thresh & pred_row_sd > sd_thresh)
  if (length(valid_row) < 2L) return(out)

  obs_r <- obs_use[valid_row, , drop = FALSE]
  pred_r <- pred_use[valid_row, , drop = FALSE]
  cm <- tryCatch(stats::cor(t(pred_r), t(obs_r)), error = function(e) NULL)
  if (is.null(cm)) return(out)

  diag_vals <- diag(cm)
  out["pattern_correlation"] <- mean(diag_vals, na.rm = TRUE)

  n <- nrow(cm)
  if (n > 1L) {
    off <- cm[row(cm) != col(cm)]
    off <- off[is.finite(off)]
    if (length(off) > 0L && is.finite(out["pattern_correlation"])) {
      out["pattern_discrimination"] <- out["pattern_correlation"] - mean(off)
    }
  }

  ranks <- rep(NA_real_, n)
  for (i in seq_len(n)) {
    row_i <- cm[i, ]
    denom <- sum(!is.na(row_i)) - 1L
    if (denom > 0L && is.finite(row_i[i])) {
      ranks[i] <- (sum(row_i <= row_i[i], na.rm = TRUE) - 1L) / denom
    }
  }
  out["pattern_rank_percentile"] <- mean(ranks, na.rm = TRUE)
  if (is.nan(out["pattern_rank_percentile"])) out["pattern_rank_percentile"] <- NA_real_

  if (length(valid_row) >= 3L) {
    pcor <- tryCatch(stats::cor(t(pred_r)), error = function(e) NULL)
    ocor <- tryCatch(stats::cor(t(obs_r)), error = function(e) NULL)
    if (!is.null(pcor) && !is.null(ocor)) {
      prdm <- 1 - pcor
      ordm <- 1 - ocor
      pv <- prdm[upper.tri(prdm)]
      ov <- ordm[upper.tri(ordm)]
      if (length(pv) >= 2L && length(pv) == length(ov)) {
        out["rdm_correlation"] <- tryCatch(
          stats::cor(pv, ov, method = rsa_simfun, use = "complete.obs"),
          error = function(e) NA_real_
        )
      }
    }
  }

  out
}

test_that("feature_rsa_da_model runs and returns target geometry metrics", {
  skip_on_cran()
  dat <- .make_frda_dataset(seed = 201)

  ms <- feature_rsa_da_model(
    dataset = dat$dataset,
    design = dat$design,
    mode = "stacked",
    lambdas = c(a = 1, b = 1),
    return_diagnostics = TRUE
  )

  res <- run_regional(ms, dat$mask)
  expect_s3_class(res, "regional_mvpa_result")
  perf_names <- names(res$performance_table)

  expect_true("target_pattern_correlation" %in% perf_names)
  expect_true("target_pattern_discrimination" %in% perf_names)
  expect_true("target_rdm_correlation" %in% perf_names)
  expect_true("target_voxel_correlation" %in% perf_names)
  expect_true("target_mse_full" %in% perf_names)
  expect_true("target_r2_full" %in% perf_names)
})

test_that("feature_rsa_da geometry metrics match independent implementation", {
  set.seed(401)
  observed <- matrix(rnorm(18 * 7), 18, 7)
  predicted <- observed + matrix(rnorm(18 * 7, sd = 0.2), 18, 7)

  got <- rMVPA:::.frda_geometry_metrics(observed, predicted, rsa_simfun = "spearman")
  ref <- .independent_geom_metrics(observed, predicted, rsa_simfun = "spearman")
  expect_equal(got[names(ref)], ref, tolerance = 1e-10)
})

test_that("feature_rsa_da geometry metrics handle degenerate columns consistently", {
  set.seed(402)
  observed <- matrix(rnorm(12 * 5), 12, 5)
  predicted <- matrix(rnorm(12 * 5), 12, 5)
  observed[, 1] <- 1
  predicted[, 1] <- 1

  got <- rMVPA:::.frda_geometry_metrics(observed, predicted, rsa_simfun = "pearson")
  ref <- .independent_geom_metrics(observed, predicted, rsa_simfun = "pearson")
  expect_equal(got[names(ref)], ref, tolerance = 1e-10)
})

test_that("feature_rsa_da_model recovers exact mapping in stacked mode", {
  dat <- .make_frda_linear_dataset(seed = 302, Tenc = 28, Trec = 20)
  ms <- feature_rsa_da_model(
    dataset = dat$dataset,
    design = dat$design,
    mode = "stacked",
    lambdas = c(a = 0, b = 0),
    alpha_target = 0,
    rsa_simfun = "spearman"
  )

  res <- run_regional(ms, dat$mask)
  perf <- res$performance_table[1, ]
  expect_gt(perf$target_r2_full, 0.99)
  expect_gt(perf$target_rdm_correlation, 0.95)
  expect_gt(perf$target_pattern_correlation, 0.90)
})

test_that("feature_rsa_da_model coupled mode remains stable for extreme rho", {
  dat <- .make_frda_linear_dataset(seed = 305, Tenc = 20, Trec = 14, dims = c(2, 2, 2), set_sizes = c(a = 3, b = 3))
  rhos <- c(0, 1e-4, 1, 1e3, 1e6)

  mses <- vapply(rhos, function(rh) {
    ms <- feature_rsa_da_model(
      dataset = dat$dataset,
      design = dat$design,
      mode = "coupled",
      lambdas = c(a = 0.1, b = 0.1),
      rho = rh
    )
    res <- run_regional(ms, dat$mask)
    as.numeric(res$performance_table$target_mse_full[1])
  }, numeric(1))

  expect_true(all(is.finite(mses)))
})

test_that("feature_rsa_da_model coupled mode handles collinearity and tiny folds", {
  dat <- .make_frda_linear_dataset(seed = 306, Tenc = 12, Trec = 8, dims = c(2, 2, 2), set_sizes = c(a = 3, b = 3))
  dat$design$block_var_test <- rep(1:2, each = 4)

  # Induce severe collinearity in both source and target predictors
  dat$design$X_train$X[, 2] <- dat$design$X_train$X[, 1]
  dat$design$X_train$X[, 6] <- dat$design$X_train$X[, 5]
  dat$design$X_test$X[, 2] <- dat$design$X_test$X[, 1]
  dat$design$X_test$X[, 6] <- dat$design$X_test$X[, 5]

  ms <- feature_rsa_da_model(
    dataset = dat$dataset,
    design = dat$design,
    mode = "coupled",
    lambdas = c(a = 0, b = 0),
    rho = 1e4,
    return_diagnostics = TRUE
  )

  res <- run_regional(ms, dat$mask)
  perf <- res$performance_table[1, ]
  expect_true(is.finite(perf$target_mse_full))
  expect_true("target_rdm_correlation" %in% names(perf))
  expect_equal(length(res$fits[[1]]$folds), 2)
})

test_that("feature_rsa_da_model validates lambda names", {
  dat <- .make_frda_dataset(seed = 303)
  expect_error(
    feature_rsa_da_model(
      dataset = dat$dataset,
      design = dat$design,
      mode = "stacked",
      lambdas = c(a = 1)
    ),
    "missing entries"
  )
})

test_that("feature_rsa_da_model validates non-negative alpha_target", {
  dat <- .make_frda_dataset(seed = 304)
  expect_error(
    feature_rsa_da_model(
      dataset = dat$dataset,
      design = dat$design,
      mode = "stacked",
      lambdas = c(a = 1, b = 1),
      alpha_target = -0.2
    ),
    "alpha_target"
  )
})

test_that("feature_rsa_da_model supports single-run purge gap", {
  skip_on_cran()
  dat <- .make_frda_dataset(seed = 202, Trec = 18)
  dat$design$block_var_test <- rep(1, 18)

  ms <- feature_rsa_da_model(
    dataset = dat$dataset,
    design = dat$design,
    mode = "stacked",
    lambdas = c(a = 1, b = 1),
    recall_nfolds = 3,
    target_gap = 1,
    return_diagnostics = TRUE
  )

  res <- run_regional(ms, dat$mask)
  folds <- res$fits[[1]]$folds
  expect_equal(length(folds), 3)

  for (f in folds) {
    test_idx <- f$test
    purge_lo <- max(1L, min(test_idx) - 1L)
    purge_hi <- min(18L, max(test_idx) + 1L)
    purged <- seq.int(purge_lo, purge_hi)
    expect_false(any(f$train %in% purged))
  }
})

test_that("feature_rsa_da_model computes single-run permutation p-values", {
  skip_on_cran()
  dat <- .make_frda_dataset(seed = 203, Trec = 16)
  dat$design$block_var_test <- rep(1, 16)

  ms <- feature_rsa_da_model(
    dataset = dat$dataset,
    design = dat$design,
    mode = "stacked",
    lambdas = c(a = 1, b = 1),
    recall_nfolds = 4,
    target_nperm = 6,
    target_perm_strategy = "circular_shift",
    return_diagnostics = TRUE
  )

  res <- run_regional(ms, dat$mask)
  perf <- res$performance_table[1, ]
  expect_true("target_perm_p_rdm_correlation" %in% names(perf))
  expect_true("target_perm_p_mse_full" %in% names(perf))
  expect_true(perf$target_perm_p_rdm_correlation >= 0 && perf$target_perm_p_rdm_correlation <= 1)
  expect_true(perf$target_perm_p_mse_full >= 0 && perf$target_perm_p_mse_full <= 1)
  expect_equal(perf$target_perm_n, 6)
})

test_that("feature_rsa_da_model disables permutation null for multi-run targets", {
  dat <- .make_frda_dataset(seed = 204, Trec = 18)
  dat$design$block_var_test <- rep(1:2, each = 9)

  ms <- expect_warning(
    feature_rsa_da_model(
      dataset = dat$dataset,
      design = dat$design,
      mode = "stacked",
      lambdas = c(a = 1, b = 1),
      target_nperm = 5
    ),
    "single-run targets"
  )
  expect_equal(ms$target_nperm, 0L)

  res <- run_regional(ms, dat$mask)
  perf <- res$performance_table[1, ]
  expect_false("target_perm_p_rdm_correlation" %in% names(perf))
})

test_that("feature_rsa_da_model permutation p-values are approximately calibrated under null", {
  skip_on_cran()
  set.seed(501)

  nrep <- 14L
  pvals <- rep(NA_real_, nrep)

  for (i in seq_len(nrep)) {
    dat <- .make_frda_dataset(Tenc = 18, Trec = 14, dims = c(2, 2, 2), Dfeat = 6, seed = 500 + i)
    dat$design$block_var_test <- rep(1, 14)

    # Break any accidental source->target predictor alignment while preserving target autocorrelation structure.
    ord <- sample.int(14)
    dat$design$X_test$X <- dat$design$X_test$X[ord, , drop = FALSE]
    dat$design$X_test$row_weights <- dat$design$X_test$row_weights[ord]

    ms <- feature_rsa_da_model(
      dataset = dat$dataset,
      design = dat$design,
      mode = "stacked",
      lambdas = c(a = 1, b = 1),
      recall_nfolds = 3,
      target_nperm = 20,
      target_perm_strategy = "circular_shift"
    )

    res <- run_regional(ms, dat$mask)
    pvals[i] <- as.numeric(res$performance_table$target_perm_p_rdm_correlation[1])
  }

  pvals <- pvals[is.finite(pvals)]
  expect_gte(length(pvals), 10)

  false_pos_rate <- mean(pvals < 0.05)
  mean_p <- mean(pvals)

  # Wide but meaningful bounds to avoid flakiness while catching severe miscalibration.
  expect_lte(false_pos_rate, 0.25)
  expect_gte(mean_p, 0.30)
  expect_lte(mean_p, 0.70)
})
