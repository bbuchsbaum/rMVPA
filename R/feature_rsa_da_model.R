#' Feature-RSA Domain Adaptation Model
#'
#' Learns a source->target mapping with grouped ridge domain adaptation and evaluates
#' target-domain representational geometry (feature-RSA style) on held-out target folds.
#' This model is intended for settings where source rows have higher SNR/coverage and
#' target rows are the analysis domain of interest.
#'
#' @param dataset An `mvpa_dataset` with external `test_data`.
#' @param design A `feature_sets_design` carrying `X_train` and `X_test`.
#' @param mode `"stacked"` or `"coupled"`.
#' @param lambdas Named ridge penalties per feature set.
#' @param alpha_recall Non-negative scalar weighting target-domain rows.
#' @param alpha_target Optional alias for `alpha_recall`; if provided, overrides.
#' @param rho Coupling strength for `mode = "coupled"`.
#' @param recall_folds Optional explicit test-domain folds (list of `train`/`test` indices).
#' @param target_folds Optional alias for `recall_folds`; if provided, overrides.
#' @param recall_nfolds Number of contiguous folds for single-run targets.
#' @param target_nfolds Optional alias for `recall_nfolds`; if provided, overrides.
#' @param recall_gap Non-negative purge gap (TRs) for single-run contiguous folds.
#' @param target_gap Optional alias for `recall_gap`; if provided, overrides.
#' @param target_nperm Number of single-run target permutations for a null model.
#' @param target_perm_strategy `"circular_shift"` or `"block_shuffle"` for target permutations.
#' @param target_perm_block Optional block size (TRs) for `"block_shuffle"`.
#' @param rsa_simfun Similarity for target RDM correlation, `"spearman"` or `"pearson"`.
#' @param return_diagnostics If TRUE, store fold/parameter diagnostics in fits.
#' @param ... Additional fields stored on the model spec.
#'
#' @return A model spec of class `feature_rsa_da_model`.
#' @export
feature_rsa_da_model <- function(dataset,
                                 design,
                                 mode = c("stacked", "coupled"),
                                 lambdas,
                                 alpha_recall = 0.2,
                                 alpha_target = NULL,
                                 rho = 5,
                                 recall_folds = NULL,
                                 target_folds = NULL,
                                 recall_nfolds = 5L,
                                 target_nfolds = NULL,
                                 recall_gap = 0L,
                                 target_gap = NULL,
                                 target_nperm = 0L,
                                 target_perm_strategy = c("circular_shift", "block_shuffle"),
                                 target_perm_block = NULL,
                                 rsa_simfun = c("spearman", "pearson"),
                                 return_diagnostics = FALSE,
                                 ...) {
  mode <- match.arg(mode)
  rsa_simfun <- match.arg(rsa_simfun)
  target_perm_strategy <- match.arg(target_perm_strategy)

  assertthat::assert_that(inherits(dataset, "mvpa_dataset"))
  assertthat::assert_that(inherits(design, "feature_sets_design"))

  if (!has_test_set(dataset) || is.null(dataset$test_data) || is.null(design$X_test)) {
    stop("feature_rsa_da_model requires an external test set (dataset$test_data + design$X_test).", call. = FALSE)
  }

  if (!is.numeric(alpha_recall) || length(alpha_recall) != 1L || !is.finite(alpha_recall) || alpha_recall < 0) {
    stop("feature_rsa_da_model: alpha_recall must be a finite non-negative scalar.", call. = FALSE)
  }
  if (!is.null(alpha_target)) {
    if (!is.numeric(alpha_target) || length(alpha_target) != 1L || !is.finite(alpha_target) || alpha_target < 0) {
      stop("feature_rsa_da_model: alpha_target must be a finite non-negative scalar.", call. = FALSE)
    }
    if (!missing(alpha_recall) && !isTRUE(all.equal(alpha_recall, alpha_target))) {
      rlang::warn("feature_rsa_da_model: both `alpha_recall` and `alpha_target` supplied; using `alpha_target`.")
    }
  }
  if (mode == "coupled") {
    if (!is.numeric(rho) || length(rho) != 1L || !is.finite(rho) || rho < 0) {
      stop("feature_rsa_da_model: rho must be a finite non-negative scalar.", call. = FALSE)
    }
  }
  alpha <- if (!is.null(alpha_target)) alpha_target else alpha_recall

  sets <- names(design$X_train$indices)
  if (is.null(names(lambdas)) || !is.numeric(lambdas)) {
    stop("feature_rsa_da_model: lambdas must be a *named* numeric vector.", call. = FALSE)
  }
  if (!all(sets %in% names(lambdas))) {
    missing <- setdiff(sets, names(lambdas))
    stop("feature_rsa_da_model: lambdas missing entries for sets: ", paste(missing, collapse = ", "), call. = FALSE)
  }

  folds_input <- recall_folds
  if (!is.null(target_folds)) {
    if (!is.null(recall_folds)) {
      rlang::warn("feature_rsa_da_model: both `recall_folds` and `target_folds` supplied; using `target_folds`.")
    }
    folds_input <- target_folds
  }

  nfolds_input <- recall_nfolds
  if (!is.null(target_nfolds)) {
    if (!missing(recall_nfolds) && !isTRUE(all.equal(recall_nfolds, target_nfolds))) {
      rlang::warn("feature_rsa_da_model: both `recall_nfolds` and `target_nfolds` supplied; using `target_nfolds`.")
    }
    nfolds_input <- target_nfolds
  }

  gap_input <- recall_gap
  if (!is.numeric(gap_input) || length(gap_input) != 1L || !is.finite(gap_input) || gap_input < 0) {
    stop("feature_rsa_da_model: recall_gap must be a finite non-negative scalar.", call. = FALSE)
  }
  if (!is.null(target_gap)) {
    if (!is.numeric(target_gap) || length(target_gap) != 1L || !is.finite(target_gap) || target_gap < 0) {
      stop("feature_rsa_da_model: target_gap must be a finite non-negative scalar.", call. = FALSE)
    }
    if (!missing(recall_gap) && !isTRUE(all.equal(recall_gap, target_gap))) {
      rlang::warn("feature_rsa_da_model: both `recall_gap` and `target_gap` supplied; using `target_gap`.")
    }
    gap_input <- target_gap
  }

  if (!is.numeric(target_nperm) || length(target_nperm) != 1L || !is.finite(target_nperm) || target_nperm < 0) {
    stop("feature_rsa_da_model: target_nperm must be a finite non-negative scalar.", call. = FALSE)
  }
  target_nperm <- as.integer(target_nperm)

  if (!is.null(target_perm_block)) {
    if (!is.numeric(target_perm_block) || length(target_perm_block) != 1L || !is.finite(target_perm_block) || target_perm_block < 1) {
      stop("feature_rsa_da_model: target_perm_block must be a finite positive scalar when provided.", call. = FALSE)
    }
    target_perm_block <- as.integer(target_perm_block)
  }

  n_blocks_target <- if (is.null(design$block_var_test)) 1L else length(unique(design$block_var_test))
  if (target_nperm > 0L && n_blocks_target >= 2L) {
    rlang::warn("feature_rsa_da_model: target_nperm > 0 requested with multi-run target data; permutation null is only used for single-run targets. Disabling target permutations.")
    target_nperm <- 0L
  }

  folds <- if (!is.null(folds_input)) {
    folds_input
  } else {
    .br_make_recall_folds(design$block_var_test, nrow(design$X_test$X), nfolds_input, gap = gap_input)
  }

  create_model_spec(
    "feature_rsa_da_model",
    dataset = dataset,
    design = design,
    mode = mode,
    lambdas = lambdas,
    alpha_recall = alpha,
    alpha_target = alpha,
    rho = rho,
    recall_folds = folds,
    target_folds = folds,
    recall_gap = as.integer(gap_input),
    target_gap = as.integer(gap_input),
    target_nperm = target_nperm,
    target_perm_strategy = target_perm_strategy,
    target_perm_block = target_perm_block,
    rsa_simfun = rsa_simfun,
    compute_performance = TRUE,
    performance = get_regression_perf(design$split_groups),
    return_predictions = isTRUE(return_diagnostics),
    return_fits = isTRUE(return_diagnostics),
    return_diagnostics = return_diagnostics,
    ...
  )
}

.frda_geometry_metrics <- function(observed, predicted, rsa_simfun = "spearman") {
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

  if (nrow(observed) != nrow(predicted) || ncol(observed) != ncol(predicted)) {
    return(out)
  }

  sd_thresh <- 1e-12
  obs_sd <- apply(observed, 2, stats::sd)
  pred_sd <- apply(predicted, 2, stats::sd)
  valid_col <- which(obs_sd > sd_thresh & pred_sd > sd_thresh)
  if (length(valid_col) < 1L) return(out)

  out["voxel_correlation"] <- tryCatch(
    stats::cor(as.vector(predicted[, valid_col, drop = FALSE]),
               as.vector(observed[, valid_col, drop = FALSE])),
    error = function(e) NA_real_
  )

  if (nrow(observed) > 1L) {
    vt <- vapply(valid_col, function(j) {
      tryCatch(stats::cor(observed[, j], predicted[, j]), error = function(e) NA_real_)
    }, numeric(1))
    out["mean_voxelwise_temporal_cor"] <- mean(vt, na.rm = TRUE)
    if (is.nan(out["mean_voxelwise_temporal_cor"])) {
      out["mean_voxelwise_temporal_cor"] <- NA_real_
    }
  }

  obs_row_sd <- apply(observed[, valid_col, drop = FALSE], 1, stats::sd)
  pred_row_sd <- apply(predicted[, valid_col, drop = FALSE], 1, stats::sd)
  valid_row <- which(obs_row_sd > sd_thresh & pred_row_sd > sd_thresh)

  if (length(valid_row) >= 2L) {
    pmat <- predicted[valid_row, valid_col, drop = FALSE]
    omat <- observed[valid_row, valid_col, drop = FALSE]
    cormat_cond <- tryCatch(stats::cor(t(pmat), t(omat)), error = function(e) NULL)
    if (!is.null(cormat_cond)) {
      diag_cors <- diag(cormat_cond)
      out["pattern_correlation"] <- mean(diag_cors, na.rm = TRUE)

      nc <- nrow(cormat_cond)
      if (nc > 1L) {
        off_vals <- cormat_cond[row(cormat_cond) != col(cormat_cond)]
        off_vals <- off_vals[!is.na(off_vals)]
        off_diag <- if (length(off_vals) > 0L) mean(off_vals) else NA_real_
        if (is.finite(out["pattern_correlation"]) && is.finite(off_diag)) {
          out["pattern_discrimination"] <- out["pattern_correlation"] - off_diag
        }
      }

      ranks <- numeric(nc)
      for (i in seq_len(nc)) {
        row_cors <- cormat_cond[i, ]
        denom <- sum(!is.na(row_cors)) - 1L
        ranks[i] <- if (denom > 0L && is.finite(row_cors[i])) {
          (sum(row_cors <= row_cors[i], na.rm = TRUE) - 1L) / denom
        } else {
          NA_real_
        }
      }
      out["pattern_rank_percentile"] <- mean(ranks, na.rm = TRUE)
      if (is.nan(out["pattern_rank_percentile"])) out["pattern_rank_percentile"] <- NA_real_
    }
  }

  if (length(valid_row) >= 3L) {
    pmat <- predicted[valid_row, valid_col, drop = FALSE]
    omat <- observed[valid_row, valid_col, drop = FALSE]
    pc <- tryCatch(stats::cor(t(pmat)), error = function(e) NULL)
    oc <- tryCatch(stats::cor(t(omat)), error = function(e) NULL)
    if (!is.null(pc) && !is.null(oc)) {
      prdm <- 1 - pc
      ordm <- 1 - oc
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

#' @method output_schema feature_rsa_da_model
#' @export
output_schema.feature_rsa_da_model <- function(model) {
  nms <- c(
    "target_pattern_correlation",
    "target_pattern_discrimination",
    "target_pattern_rank_percentile",
    "target_rdm_correlation",
    "target_voxel_correlation",
    "target_mean_voxelwise_temporal_cor",
    "target_mse_full",
    "target_r2_full"
  )

  if (as.integer(model$target_nperm %||% 0L) > 0L) {
    base_suffixes <- c(
      "pattern_correlation", "pattern_discrimination", "pattern_rank_percentile",
      "rdm_correlation", "voxel_correlation", "mean_voxelwise_temporal_cor",
      "mse_full", "r2_full"
    )
    perm_nms <- c(
      paste0("target_perm_p_", base_suffixes),
      paste0("target_perm_z_", base_suffixes),
      "target_perm_n"
    )
    nms <- c(nms, perm_nms)
  }

  setNames(rep("scalar", length(nms)), nms)
}

#' @export
compute_performance.feature_rsa_da_model <- function(obj, result) {
  if (is.numeric(result) && !is.null(names(result))) {
    return(result)
  }
  if (is.list(result) && !is.null(result$performance) &&
      is.numeric(result$performance) && !is.null(names(result$performance))) {
    return(result$performance)
  }
  if (is.list(result) && !is.null(result$predictor) && is.null(result$observed)) {
    return(NULL)
  }
  obj$performance(result)
}

#' @rdname fit_roi
#' @export
fit_roi.feature_rsa_da_model <- function(model, roi_data, context, ...) {
  ind <- roi_data$indices
  id <- context$id

  if (is.null(roi_data$test_data)) {
    return(roi_result(
      metrics = NULL, indices = ind, id = id,
      error = TRUE,
      error_message = "feature_rsa_da_model requires external test set"
    ))
  }

  Xtrain_fs <- model$design$X_train
  Xtest_fs <- model$design$X_test

  Ye <- roi_data$train_data
  Yr <- roi_data$test_data

  if (nrow(Ye) < 2L || nrow(Yr) < 2L || ncol(Ye) < 1L) {
    return(roi_result(
      metrics = NULL, indices = ind, id = id,
      error = TRUE,
      error_message = "feature_rsa_da_model: insufficient observations or voxels"
    ))
  }

  Xe <- Xtrain_fs$X
  Xr <- Xtest_fs$X
  if (nrow(Xe) != nrow(Ye)) {
    return(roi_result(
      metrics = NULL, indices = ind, id = id,
      error = TRUE,
      error_message = "feature_rsa_da_model: nrow(X_train) must match nrow(Y_train)"
    ))
  }
  if (nrow(Xr) != nrow(Yr)) {
    return(roi_result(
      metrics = NULL, indices = ind, id = id,
      error = TRUE,
      error_message = "feature_rsa_da_model: nrow(X_test) must match nrow(Y_test)"
    ))
  }

  pen_full <- .br_penalty_vec(Xtrain_fs, model$lambdas)
  if (any(!is.finite(pen_full)) || any(pen_full < 0)) {
    return(roi_result(
      metrics = NULL, indices = ind, id = id,
      error = TRUE,
      error_message = "feature_rsa_da_model: invalid lambdas (must be finite and non-negative)"
    ))
  }

  folds <- model$target_folds
  w_rec <- Xtest_fs$row_weights
  if (is.null(w_rec)) {
    w_rec <- rep(1, nrow(Xr))
  } else if (length(w_rec) != nrow(Xr)) {
    stop(sprintf("feature_rsa_da_model: row_weights length (%d) != nrow(X_test) (%d)",
                 length(w_rec), nrow(Xr)), call. = FALSE)
  }

  mode <- model$mode
  alpha <- model$alpha_target
  rho <- model$rho

  fit_fold <- function(keep_cols, lambdas_pen, Yr_use = Yr) {
    Xe_k <- Xe[, keep_cols, drop = FALSE]
    Xr_k <- Xr[, keep_cols, drop = FALSE]
    pen_k <- lambdas_pen[keep_cols]

    fold_scores <- lapply(seq_along(folds), function(fi) {
      tr <- folds[[fi]]$train
      te <- folds[[fi]]$test
      if (length(te) < 2L || length(tr) < 2L) {
        return(c(
          target_pattern_correlation = NA_real_,
          target_pattern_discrimination = NA_real_,
          target_pattern_rank_percentile = NA_real_,
          target_rdm_correlation = NA_real_,
          target_voxel_correlation = NA_real_,
          target_mean_voxelwise_temporal_cor = NA_real_,
          target_mse_full = NA_real_,
          target_r2_full = NA_real_
        ))
      }

      Xr_tr <- Xr_k[tr, , drop = FALSE]
      Yr_tr <- Yr_use[tr, , drop = FALSE]
      Xr_te <- Xr_k[te, , drop = FALSE]
      Yr_te <- Yr_use[te, , drop = FALSE]
      w_tr <- w_rec[tr]

      fit <- if (mode == "stacked") {
        .br_fit_stacked(Xe_k, Ye, Xr_tr, Yr_tr, pen_k, alpha, w_tr)
      } else {
        .br_fit_coupled(Xe_k, Ye, Xr_tr, Yr_tr, pen_k, alpha, rho, w_tr)
      }

      Yhat <- .br_predict(Xr_te, fit, mode = mode)
      geom <- .frda_geometry_metrics(Yr_te, Yhat, rsa_simfun = model$rsa_simfun)
      c(
        target_pattern_correlation = unname(geom["pattern_correlation"]),
        target_pattern_discrimination = unname(geom["pattern_discrimination"]),
        target_pattern_rank_percentile = unname(geom["pattern_rank_percentile"]),
        target_rdm_correlation = unname(geom["rdm_correlation"]),
        target_voxel_correlation = unname(geom["voxel_correlation"]),
        target_mean_voxelwise_temporal_cor = unname(geom["mean_voxelwise_temporal_cor"]),
        target_mse_full = .br_mse(Yr_te, Yhat),
        target_r2_full = .br_r2_mean(Yr_te, Yhat)
      )
    })

    fold_mat <- do.call(rbind, fold_scores)
    colnames(fold_mat) <- names(fold_scores[[1]])
    fold_mean <- colMeans(fold_mat, na.rm = TRUE)
    fold_mean[is.nan(fold_mean)] <- NA_real_
    list(mean = fold_mean, by_fold = fold_mat)
  }

  full_cols <- seq_len(ncol(Xe))
  full_fit <- fit_fold(full_cols, pen_full)
  perf <- full_fit$mean

  perm_n <- as.integer(model$target_nperm %||% 0L)
  perm_diag <- NULL
  if (perm_n > 0L) {
    T_rec <- nrow(Yr)
    perm_strategy <- model$target_perm_strategy %||% "circular_shift"
    perm_block <- model$target_perm_block %||% NULL

    null_mat <- matrix(NA_real_, nrow = perm_n, ncol = length(perf))
    colnames(null_mat) <- names(perf)
    for (pi in seq_len(perm_n)) {
      perm_idx <- .br_perm_indices(T_rec, strategy = perm_strategy, block_size = perm_block)
      Yr_perm <- Yr[perm_idx, , drop = FALSE]
      null_mat[pi, ] <- fit_fold(full_cols, pen_full, Yr_use = Yr_perm)$mean
    }

    for (nm in names(perf)) {
      obs <- perf[[nm]]
      null <- null_mat[, nm]
      ok <- is.finite(null)
      n_ok <- sum(ok)
      suffix <- sub("^target_", "", nm)
      pnm <- paste0("target_perm_p_", suffix)
      znm <- paste0("target_perm_z_", suffix)

      if (n_ok < 1L || !is.finite(obs)) {
        perf[pnm] <- NA_real_
        perf[znm] <- NA_real_
        next
      }

      if (identical(nm, "target_mse_full")) {
        perf[pnm] <- (sum(null[ok] <= obs) + 1) / (n_ok + 1)
      } else {
        perf[pnm] <- (sum(null[ok] >= obs) + 1) / (n_ok + 1)
      }

      if (n_ok > 1L) {
        mu <- mean(null[ok])
        sdv <- stats::sd(null[ok])
        if (identical(nm, "target_mse_full")) {
          perf[znm] <- (mu - obs) / max(sdv, .Machine$double.eps)
        } else {
          perf[znm] <- (obs - mu) / max(sdv, .Machine$double.eps)
        }
      } else {
        perf[znm] <- NA_real_
      }
    }
    perf["target_perm_n"] <- perm_n
    perm_diag <- list(null_metrics = null_mat)
  }

  diag_obj <- NULL
  if (isTRUE(model$return_diagnostics)) {
    diag_obj <- list(
      folds = folds,
      full_metrics_by_fold = full_fit$by_fold,
      lambdas = model$lambdas,
      alpha_target = alpha,
      target_gap = model$target_gap %||% 0L,
      target_nperm = perm_n,
      target_perm_strategy = model$target_perm_strategy %||% "circular_shift",
      target_perm_block = model$target_perm_block %||% NA_integer_,
      rsa_simfun = model$rsa_simfun,
      rho = if (mode == "coupled") rho else NA_real_,
      mode = mode,
      perm = perm_diag
    )
  }

  roi_result(
    metrics = perf,
    indices = ind,
    id = id,
    result = list(predictor = diag_obj)
  )
}
