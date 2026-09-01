library(rMVPA)
library(neuroim2)

.pred_fixture <- function(seed = 79L,
                          dims = c(3, 3, 3),
                          nobs = 24L,
                          blocks = 3L,
                          n_features = 6L,
                          n_roi = 2L,
                          method = "pca",
                          return_predictions = TRUE,
                          return_rdm_vectors = FALSE,
                          max_retained_mb = 1024,
                          prediction_overflow = "error") {
  set.seed(seed)
  dset <- gen_sample_dataset(dims, nobs, blocks = blocks)
  Fmat <- matrix(rnorm(nobs * n_features), nobs, n_features)
  fdes <- feature_rsa_design(
    F = Fmat,
    labels = paste0("o", seq_len(nobs)),
    max_comps = min(3L, n_features),
    block_var = dset$design$block_var
  )
  model_args <- list(
    dataset = dset$dataset,
    design = fdes,
    method = method,
    crossval = blocked_cross_validation(dset$design$block_var),
    return_predictions = return_predictions,
    return_rdm_vectors = return_rdm_vectors,
    max_retained_mb = max_retained_mb,
    prediction_overflow = prediction_overflow
  )
  if (method %in% c("pls", "pca")) {
    model_args$ncomp_selection <- "max"
  } else if (identical(method, "ridge")) {
    model_args$lambda <- 0.1
    model_args$lambda_selection <- "fixed"
  }
  mspec <- do.call(feature_rsa_model, model_args)
  region_mask <- NeuroVol(
    sample(seq_len(n_roi), length(dset$dataset$mask), replace = TRUE),
    space(dset$dataset$mask)
  )
  list(
    dataset = dset$dataset,
    design = fdes,
    block_var = dset$design$block_var,
    model = mspec,
    region_mask = region_mask,
    nobs = nobs
  )
}

.window_average_by_fold <- function(mat, fold_id, window = 2L) {
  mat <- as.matrix(mat)
  fold_id <- as.vector(fold_id)
  stopifnot(nrow(mat) == length(fold_id), window >= 1L)
  rows <- lapply(split(seq_len(nrow(mat)), fold_id), function(idx) {
    n <- length(idx)
    starts <- seq.int(1L, n, by = window)
    do.call(rbind, lapply(starts, function(start) {
      take <- idx[seq.int(start, min(start + window - 1L, n))]
      colMeans(mat[take, , drop = FALSE])
    }))
  })
  list(
    mat = do.call(rbind, rows),
    fold_id = rep.int(as.integer(names(rows)), vapply(rows, nrow, integer(1)))
  )
}

.pattern_rank <- function(predicted, observed) {
  predicted <- as.matrix(predicted)
  observed <- as.matrix(observed)
  cm <- stats::cor(t(predicted), t(observed), use = "pairwise.complete.obs")
  ranks <- vapply(seq_len(nrow(cm)), function(i) {
    row <- cm[i, ]
    (sum(row <= row[[i]], na.rm = TRUE) - 1L) / (sum(is.finite(row)) - 1L)
  }, numeric(1))
  mean(ranks, na.rm = TRUE)
}

.diagonal_whiten <- function(predicted, observed) {
  resid <- as.matrix(observed) - as.matrix(predicted)
  scale <- apply(resid, 2L, stats::sd, na.rm = TRUE)
  scale[!is.finite(scale) | scale < 1e-8] <- 1
  list(
    predicted = sweep(as.matrix(predicted), 2L, scale, "/"),
    observed = sweep(as.matrix(observed), 2L, scale, "/")
  )
}

.merge_result_set <- function(observed, predicted, folds) {
  tibble::tibble(
    observed = lapply(folds, function(index) observed[index, , drop = FALSE]),
    predicted = lapply(folds, function(index) predicted[index, , drop = FALSE]),
    test_index = folds,
    result = rep(list(list(ncomp = 2L)), length(folds)),
    performance = rep(list(NULL), length(folds)),
    error = FALSE,
    error_message = "~"
  )
}


test_that("feature_rsa_model exposes return_predictions with a safe default", {
  defaults <- formals(feature_rsa_model)
  expect_true("return_predictions" %in% names(defaults))
  expect_false(eval(defaults$return_predictions))
  expect_equal(eval(defaults$max_retained_mb), 1024)
  expect_identical(eval(defaults$prediction_overflow)[[1L]], "error")
  expect_true("feature_rsa_predictions" %in% getNamespaceExports("rMVPA"))
})

test_that("feature_rsa_model validates prediction-retention arguments", {
  dset <- gen_sample_dataset(c(3, 3, 3), 18, blocks = 3)
  Fmat <- matrix(rnorm(18 * 4), 18, 4)
  fdes <- feature_rsa_design(F = Fmat, labels = paste0("o", 1:18), max_comps = 3)

  expect_error(
    feature_rsa_model(
      dset$dataset, fdes, method = "pca",
      crossval = blocked_cross_validation(dset$design$block_var),
      return_predictions = NA
    ),
    "return_predictions must be TRUE or FALSE"
  )
  expect_error(
    feature_rsa_model(
      dset$dataset, fdes, method = "pca",
      crossval = blocked_cross_validation(dset$design$block_var),
      return_predictions = c(TRUE, FALSE)
    ),
    "return_predictions must be TRUE or FALSE"
  )
  expect_error(
    feature_rsa_model(
      dset$dataset, fdes, method = "pca",
      crossval = blocked_cross_validation(dset$design$block_var),
      return_predictions = TRUE,
      max_retained_mb = 0
    ),
    "max_retained_mb must be one positive finite value"
  )
  expect_error(
    feature_rsa_model(
      dset$dataset, fdes, method = "pca",
      crossval = blocked_cross_validation(dset$design$block_var),
      return_predictions = TRUE,
      max_retained_mb = Inf
    ),
    "max_retained_mb must be one positive finite value"
  )
  expect_error(
    feature_rsa_model(
      dset$dataset, fdes, method = "pca",
      crossval = blocked_cross_validation(dset$design$block_var),
      prediction_overflow = "drop"
    ),
    "prediction_overflow"
  )
})

test_that("prediction retention helpers count Yhat and Y and enforce the contract", {
  expect_equal(
    rMVPA:::.feature_rsa_prediction_retention_bytes(10, 5, FALSE),
    0
  )
  expect_equal(
    rMVPA:::.feature_rsa_prediction_retention_bytes(10, 5, TRUE),
    8 * 10 * 5 * 2
  )

  ok <- rMVPA:::.feature_rsa_validate_prediction_retention(
    n_obs = 10, n_voxels = 5, return_predictions = TRUE,
    max_retained_mb = 1, prediction_overflow = "error"
  )
  expect_true(ok$return_predictions)
  expect_equal(ok$retention_estimated_bytes[["total"]], 8 * 10 * 5 * 2)

  expect_error(
    rMVPA:::.feature_rsa_validate_prediction_retention(
      n_obs = 2000, n_voxels = 40000, return_predictions = TRUE,
      max_retained_mb = 1, prediction_overflow = "error"
    ),
    "above max_retained_mb"
  )

  fallback <- rMVPA:::.feature_rsa_validate_prediction_retention(
    n_obs = 2000, n_voxels = 40000, return_predictions = TRUE,
    max_retained_mb = 1, prediction_overflow = "none"
  )
  expect_false(fallback$return_predictions)
  expect_equal(fallback$retention_estimated_bytes[["total"]], 0)
  expect_match(fallback$retention_notice, "disabled")
})

test_that("constructor refuses or disables oversized prediction retention", {
  dset <- gen_sample_dataset(c(3, 3, 3), 24, blocks = 3)
  Fmat <- matrix(rnorm(24 * 5), 24, 5)
  fdes <- feature_rsa_design(F = Fmat, labels = paste0("o", 1:24), max_comps = 3)
  cv <- blocked_cross_validation(dset$design$block_var)

  expect_error(
    feature_rsa_model(
      dset$dataset, fdes, method = "pca", crossval = cv,
      ncomp_selection = "max",
      return_predictions = TRUE,
      max_retained_mb = 1e-6
    ),
    "above max_retained_mb"
  )

  fallback <- feature_rsa_model(
    dset$dataset, fdes, method = "pca", crossval = cv,
    ncomp_selection = "max",
    return_predictions = TRUE,
    max_retained_mb = 1e-6,
    prediction_overflow = "none"
  )
  expect_false(fallback$return_predictions)
  expect_false(isTRUE(fallback$return_fits))
  expect_match(fallback$retention_notice, "disabled")
  expect_equal(fallback$retention_estimated_bytes[["total"]], 0)
})

test_that("merge_results stores aligned OOF predictions and fold ids", {
  set.seed(4422)
  observed <- matrix(rnorm(24 * 8), 24, 8)
  predicted <- observed * 0.6 + matrix(rnorm(24 * 8, sd = 0.4), 24, 8)
  folds <- list(17:24, 1:8, 9:16)
  model <- structure(
    list(
      method = "pls",
      nperm = 0L,
      save_distributions = FALSE,
      return_rdm_vectors = FALSE,
      return_predictions = TRUE
    ),
    class = "feature_rsa_model"
  )
  fold_result <- .merge_result_set(observed, predicted, folds)
  merged <- merge_results(
    model, fold_result, indices = 101:108, id = 4L
  )

  pred <- merged$result[[1L]]$predictor
  expect_false(is.null(pred$predicted))
  expect_equal(pred$predicted, predicted, tolerance = 1e-12)
  expect_equal(pred$observed, observed, tolerance = 1e-12)
  expect_identical(pred$observation_index, seq_len(24))
  expect_identical(pred$fold_id, c(rep(2L, 8), rep(3L, 8), rep(1L, 8)))
  expect_identical(pred$voxel_index, 101:108)
  expect_identical(pred$n_obs, 24L)
  expect_null(pred$predicted_rdm_vec)
  expect_false(is.null(merged$performance[[1L]]))

  optimized <- fold_result
  optimized$observed <- rep(list(NULL), nrow(optimized))
  optimized$predicted <- rep(list(NULL), nrow(optimized))
  attr(optimized, "feature_rsa_oof") <- list(
    observed = observed,
    predicted = predicted,
    test_index = seq_len(24),
    fold_id = c(rep(2L, 8), rep(3L, 8), rep(1L, 8))
  )
  compact <- merge_results(
    model, optimized, indices = 101:108, id = 4L
  )
  expect_equal(compact$result[[1L]]$predictor$predicted, predicted, tolerance = 1e-12)
  expect_identical(
    compact$result[[1L]]$predictor$fold_id,
    merged$result[[1L]]$predictor$fold_id
  )
  expect_equal(
    compact$performance[[1L]],
    merged$performance[[1L]],
    tolerance = 1e-12
  )
})

test_that("merge_results omits prediction matrices unless requested", {
  set.seed(17)
  observed <- matrix(rnorm(12 * 4), 12, 4)
  predicted <- observed + 0.1
  folds <- list(1:4, 5:8, 9:12)

  none <- merge_results(
    structure(
      list(
        method = "pls", nperm = 0L, save_distributions = FALSE,
        return_rdm_vectors = FALSE, return_predictions = FALSE
      ),
      class = "feature_rsa_model"
    ),
    .merge_result_set(observed, predicted, folds),
    indices = 1:4, id = 1L
  )
  expect_null(none$result[[1L]]$predictor)

  rdm_only <- merge_results(
    structure(
      list(
        method = "pls", nperm = 0L, save_distributions = FALSE,
        return_rdm_vectors = TRUE, return_predictions = FALSE
      ),
      class = "feature_rsa_model"
    ),
    .merge_result_set(observed, predicted, folds),
    indices = 1:4, id = 1L
  )
  expect_false(is.null(rdm_only$result[[1L]]$predictor[["predicted_rdm_vec"]]))
  expect_null(rdm_only$result[[1L]]$predictor[["predicted"]])

  both <- merge_results(
    structure(
      list(
        method = "pls", nperm = 0L, save_distributions = FALSE,
        return_rdm_vectors = TRUE, return_predictions = TRUE
      ),
      class = "feature_rsa_model"
    ),
    .merge_result_set(observed, predicted, folds),
    indices = 1:4, id = 1L
  )
  expect_false(is.null(both$result[[1L]]$predictor[["predicted_rdm_vec"]]))
  expect_equal(both$result[[1L]]$predictor[["predicted"]], predicted, tolerance = 1e-12)
})

test_that("default feature RSA regional results still omit prediction payloads", {
  fix <- .pred_fixture(return_predictions = FALSE)
  res <- run_regional(fix$model, fix$region_mask)

  expect_s3_class(res, "regional_mvpa_result")
  expect_null(res$fits)
  expect_null(res$prediction_table)
  expect_error(feature_rsa_predictions(res), "return_predictions=TRUE")
  expect_true("pattern_correlation" %in% names(res$performance_table))
})

test_that("run_regional retains per-ROI OOF predictions with fold geometry", {
  fix <- .pred_fixture(method = "pca", n_roi = 3L)
  res <- run_regional(fix$model, fix$region_mask)

  expect_null(res$prediction_table)
  expect_true(!is.null(res$fits))
  preds <- feature_rsa_predictions(res)
  expect_equal(nrow(preds), nrow(res$performance_table))
  expect_setequal(preds$roinum, res$performance_table$roinum)

  for (i in seq_len(nrow(preds))) {
    yhat <- preds$predicted[[i]]
    yobs <- preds$observed[[i]]
    expect_true(is.matrix(yhat))
    expect_identical(dim(yhat), dim(yobs))
    expect_equal(nrow(yhat), fix$nobs)
    expect_equal(preds$n_obs[[i]], fix$nobs)
    expect_identical(preds$observation_index[[i]], seq_len(fix$nobs))
    expect_length(preds$fold_id[[i]], fix$nobs)
    expect_false(anyNA(preds$fold_id[[i]]))
    expect_equal(length(preds$voxel_index[[i]]), ncol(yhat))
    expect_true(all(preds$voxel_index[[i]] > 0))
    expect_equal(length(unique(preds$voxel_index[[i]])), ncol(yhat))
  }

  expect_true(all(vapply(
    preds$fold_id,
    function(fid) identical(fid, preds$fold_id[[1L]]),
    logical(1)
  )))
})

test_that("returned predictions reconstruct the stored performance metrics", {
  for (method in c("pca", "pls", "ridge")) {
    fix <- .pred_fixture(seed = 800 + nchar(method), method = method, n_roi = 2L)
    res <- run_regional(fix$model, fix$region_mask)
    preds <- feature_rsa_predictions(res)
    perf <- res$performance_table

    for (i in seq_len(nrow(preds))) {
      row <- perf[perf$roinum == preds$roinum[[i]], , drop = FALSE]
      rebuilt <- evaluate_model.feature_rsa_model(
        object = fix$model,
        predicted = preds$predicted[[i]],
        observed = preds$observed[[i]],
        fold_id = preds$fold_id[[i]]
      )
      expect_equal(
        rebuilt$pattern_correlation,
        row$pattern_correlation[[1L]],
        tolerance = 1e-12,
        info = method
      )
      expect_equal(
        rebuilt$pattern_discrimination,
        row$pattern_discrimination[[1L]],
        tolerance = 1e-12,
        info = method
      )
      expect_equal(
        rebuilt$pattern_rank_percentile,
        row$pattern_rank_percentile[[1L]],
        tolerance = 1e-12,
        info = method
      )
      expect_equal(
        rebuilt$rdm_correlation,
        row$rdm_correlation[[1L]],
        tolerance = 1e-12,
        info = method
      )
      expect_equal(rebuilt$mse, row$mse[[1L]], tolerance = 1e-12, info = method)
    }
  }
})

test_that("blocked CV fold ids stay aligned with acquisition blocks", {
  fix <- .pred_fixture(nobs = 30L, blocks = 3L)
  res <- run_regional(fix$model, fix$region_mask)
  preds <- feature_rsa_predictions(res)
  obs_index <- preds$observation_index[[1L]]
  fold_id <- preds$fold_id[[1L]]
  block <- fix$block_var[obs_index]

  expect_equal(length(unique(fold_id)), length(unique(block)))
  fold_by_block <- tapply(fold_id, block, function(x) length(unique(x)))
  expect_true(all(fold_by_block == 1L))
  block_by_fold <- tapply(block, fold_id, function(x) length(unique(x)))
  expect_true(all(block_by_fold == 1L))
})

test_that("fold ids prevent cross-fold identification scoring", {
  fix <- .pred_fixture(nobs = 24L, blocks = 3L, n_roi = 1L)
  res <- run_regional(fix$model, fix$region_mask)
  preds <- feature_rsa_predictions(res)
  yhat <- preds$predicted[[1L]]
  yobs <- preds$observed[[1L]]
  fold_id <- preds$fold_id[[1L]]

  within_fold <- evaluate_model.feature_rsa_model(
    object = fix$model,
    predicted = yhat,
    observed = yobs,
    fold_id = fold_id
  )
  pooled <- evaluate_model.feature_rsa_model(
    object = fix$model,
    predicted = yhat,
    observed = yobs,
    fold_id = rep.int(1L, nrow(yhat))
  )

  expect_equal(
    within_fold$pattern_rank_percentile,
    res$performance_table$pattern_rank_percentile[[1L]],
    tolerance = 1e-12
  )
  expect_false(isTRUE(all.equal(
    within_fold$pattern_rank_percentile,
    pooled$pattern_rank_percentile,
    tolerance = 1e-12
  )))
  expect_false(isTRUE(all.equal(
    within_fold$rdm_correlation,
    pooled$rdm_correlation,
    tolerance = 1e-12
  )))
})

test_that("returned predictions support fold-aware temporal window scoring", {
  fix <- .pred_fixture(nobs = 24L, blocks = 3L, n_roi = 1L)
  res <- run_regional(fix$model, fix$region_mask)
  preds <- feature_rsa_predictions(res)
  yhat <- preds$predicted[[1L]]
  yobs <- preds$observed[[1L]]
  fold_id <- preds$fold_id[[1L]]
  obs_index <- preds$observation_index[[1L]]

  fold_windows <- .window_average_by_fold(yhat, fold_id, window = 2L)
  obs_windows <- .window_average_by_fold(yobs, fold_id, window = 2L)
  expect_identical(fold_windows$fold_id, obs_windows$fold_id)
  expect_true(is.finite(.pattern_rank(fold_windows$mat, obs_windows$mat)))

  source_rows <- split(seq_len(nrow(yhat)), fold_id)
  expect_true(all(vapply(source_rows, function(idx) {
    length(unique(fold_id[idx])) == 1L
  }, logical(1))))

  # Adjacent stored rows can belong to different folds after sorting by
  # observation_index. A naive fixed window then mixes training and test
  # geometry; fold ids make that detectable.
  fold_changes <- which(fold_id[-1L] != fold_id[-length(fold_id)])
  expect_true(length(fold_changes) >= 1L)
  boundary <- fold_changes[[1L]]
  naive_rows <- seq.int(max(1L, boundary - 1L), min(length(fold_id), boundary + 2L))
  expect_gt(length(unique(fold_id[naive_rows])), 1L)
  expect_equal(length(unique(fix$block_var[obs_index[naive_rows]])), 2L)
})

test_that("returned predictions support residual-whitened pattern scoring", {
  fix <- .pred_fixture(n_roi = 1L, method = "pca")
  res <- run_regional(fix$model, fix$region_mask)
  preds <- feature_rsa_predictions(res)
  raw <- evaluate_model.feature_rsa_model(
    object = fix$model,
    predicted = preds$predicted[[1L]],
    observed = preds$observed[[1L]],
    fold_id = preds$fold_id[[1L]]
  )
  whitened <- .diagonal_whiten(preds$predicted[[1L]], preds$observed[[1L]])
  white <- evaluate_model.feature_rsa_model(
    object = fix$model,
    predicted = whitened$predicted,
    observed = whitened$observed,
    fold_id = preds$fold_id[[1L]]
  )

  expect_equal(
    raw$pattern_correlation,
    res$performance_table$pattern_correlation[[1L]],
    tolerance = 1e-12
  )
  expect_true(is.finite(white$pattern_correlation))
  expect_true(is.finite(white$pattern_rank_percentile))
  expect_false(isTRUE(all.equal(
    raw$pattern_correlation,
    white$pattern_correlation,
    tolerance = 1e-12
  )))
})

test_that("return_predictions can be combined with return_rdm_vectors", {
  fix <- .pred_fixture(return_predictions = TRUE, return_rdm_vectors = TRUE)
  res <- run_regional(fix$model, fix$region_mask)

  preds <- feature_rsa_predictions(res)
  vecs <- feature_rsa_rdm_vectors(res)
  expect_equal(preds$roinum, vecs$roinum)
  expect_equal(preds$observation_index, vecs$observation_index)
  expect_equal(preds$fold_id, vecs$fold_id)
  expect_equal(length(vecs$rdm_vec[[1L]]), fix$nobs * (fix$nobs - 1L) / 2L)
})

test_that("feature_rsa_predictions rejects RDM-only results and accepts tibbles", {
  rdm_fix <- .pred_fixture(return_predictions = FALSE, return_rdm_vectors = TRUE)
  rdm_res <- run_regional(rdm_fix$model, rdm_fix$region_mask)
  expect_error(feature_rsa_predictions(rdm_res), "return_predictions=TRUE")
  expect_error(feature_rsa_predictions(list()), "regional_mvpa_result")

  ready <- tibble::tibble(
    roinum = 7L,
    predicted = list(matrix(1, 2, 2))
  )
  expect_equal(feature_rsa_predictions(ready), ready)
})

test_that("run_regional re-checks the allocation contract after construction", {
  fix <- .pred_fixture(return_predictions = TRUE)
  expect_true(fix$model$return_predictions)

  too_small <- fix$model
  too_small$max_retained_mb <- 1e-6
  expect_error(
    run_regional(too_small, fix$region_mask),
    "above max_retained_mb"
  )

  fallback <- fix$model
  fallback$max_retained_mb <- 1e-6
  fallback$prediction_overflow <- "none"
  res <- run_regional(fallback, fix$region_mask)
  expect_s3_class(res, "regional_mvpa_result")
  expect_error(feature_rsa_predictions(res), "return_predictions=TRUE")
})

test_that("searchlight refuses feature RSA prediction retention", {
  fix <- .pred_fixture(return_predictions = TRUE)
  expect_error(
    run_searchlight(fix$model, radius = 2),
    "return_predictions = TRUE"
  )
  expect_false(isTRUE(.pred_fixture(return_predictions = FALSE)$model$return_predictions))
})

test_that("print.feature_rsa_model reports prediction retention", {
  on <- .pred_fixture(return_predictions = TRUE)$model
  off <- .pred_fixture(return_predictions = FALSE)$model
  expect_match(paste(capture.output(print(on)), collapse = "\n"), "Return predictions:\\s+Yes")
  expect_match(paste(capture.output(print(off)), collapse = "\n"), "Return predictions:\\s+No")
})
