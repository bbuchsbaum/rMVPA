#' REMAP-RRR: Residual low-rank, domain-adaptive cross-decoding
#'
#' Domain-adaptive cross-decoding that models memory as "perception + low-rank correction"
#' in a jointly whitened space. Given paired item/class prototypes from a source
#' domain (e.g., perception) and a target domain (e.g., memory), REMAP learns a
#' reduced-rank residual map \eqn{\Delta} by regressing the residuals \eqn{R_w = Y_w - X_w}
#' on the whitened source prototypes \eqn{X_w}. Predicted target templates are then
#' \deqn{\hat Y_w = X_w + \lambda X_w \Delta,}
#' where \eqn{\lambda \in [0,1]} is chosen on the training items (prototype self-retrieval).
#' Classification of held-out target trials is performed by correlating their whitened
#' patterns with \eqn{\hat Y_w}. When the data do not support a reliable residual, the
#' method automatically shrinks to the naive cross-decoder (\eqn{\lambda = 0}).
#'
#' @param dataset mvpa_dataset with `train_data`, optional `test_data`, and `mask`.
#' @param design mvpa_design with `y_train`, optional `y_test`, and design tables.
#' @param base_classifier Name in the rMVPA model registry (kept for compatibility; ignored by REMAP which uses correlation to templates internally).
#' @param link_by Optional column present in both `train_design` and `test_design` used to pair
#'   source and target trials (e.g., "ImageID"). If `NULL`, pairs by class labels (`y`).
#' @param center logical; center prototypes before whitening (default TRUE).
#' @param shrink_whiten logical; use `corpcor::cov.shrink` for joint whitening (default TRUE).
#' @param rank "auto" for `rrpack::cv.rrr` rank selection, integer for fixed rank, or `0`
#'   for identity/no adaptation fallback.
#' @param max_rank Upper bound on rank search (default 20).
#' @param ridge_rrr_lambda Optional numeric lambda for `rrpack::rrs.fit` (ridge RRR) on the residuals.
#' @param leave_one_key_out logical; if TRUE, use leave-one-key-out (LOKO) over items for adaptor learning (default FALSE).
#' @param min_pairs Minimum number of paired prototypes required for adaptor fitting (default 5).
#' @param save_fold_singulars logical; if TRUE, save singular values from each fold's adaptor (default FALSE).
#' @param return_adapter logical; if TRUE, store diagnostics (e.g., per-voxel R2, per-fold stats) in `result$predictor`.
#' @param lambda_grid Numeric vector of candidate \eqn{\lambda} values for shrinkage toward naive; defaults to `c(0, .25, .5, .75, 1)`.
#' @param compute_naive_baseline logical; if TRUE, also compute performance with \eqn{\lambda = 0} (no adaptor) to provide baseline metrics (default FALSE).
#' @param calibration_rows Optional distinct row indices, or a logical selector, in
#'   the target dataset. Fit on these labeled rows and score only their complement.
#'   Select whole independent runs when evaluating run generalization. With
#'   this split, `compute_naive_baseline` adds strictly source-only metrics
#'   prefixed `source_only_`, and `return_adapter` retains a reusable fit.
#'   Incompatible with leave_one_key_out.
#' @param ... Additional arguments passed to `create_model_spec`.
#'
#' @return A model spec of class `remap_rrr_model` compatible with `run_regional()` / `run_searchlight()`.
#'
#' @section Performance metrics (added to base metrics):
#' \itemize{
#'   \item Base metrics: whatever `compute_performance()` returns for the response type (e.g., accuracy/AUC for classification or regression metrics).
#'   \item Naive baseline metrics (if `compute_naive_baseline = TRUE`): base metrics recomputed with \eqn{\lambda = 0} (no adaptor), prefixed `naive_` (e.g., `naive_Accuracy`, `naive_AUC`).
#'   \item `adapter_rank`: Mean adaptor rank used (LOKO folds or single fit).
#'   \item `adapter_sv1`: Mean leading singular value of the adaptor across folds.
#'   \item `adapter_mean_r2`: Mean per-voxel cross-validated R\eqn{^2} (LOKO only).
#'   \item `remap_improv`: Mean relative reduction in residual sum of squares vs. naive (\code{Y_w - X_w}).
#'   \item `delta_frob_mean`: Mean Frobenius norm of \eqn{\lambda \Delta} (magnitude of the learned correction).
#'   \item `lambda_mean`: Mean selected \eqn{\lambda} across folds or the single fit.
#'   \item `n_pairs_used`: Number of paired prototypes used.
#'   \item `n_skipped_keys`: Keys skipped during LOKO due to insufficient pairs or fit issues.
#' }
#' If `return_adapter = TRUE`, the result also includes ROI-level diagnostics (e.g., per-voxel R\eqn{^2}, per-item residuals, ranks, singular spectra) in `result$predictor`.
#'
#' @section Key ideas:
#' \itemize{
#'   \item Joint whitening (shared covariance) puts perception and memory in the same coordinates.
#'   \item Residual RRR learns a low-rank correction \eqn{\Delta} on \eqn{R_w = Y_w - X_w}.
#'   \item Automatic shrinkage to naive via \eqn{\lambda} selection on training items.
#'   \item Prediction/classification by correlation to \eqn{\hat Y_w} in the whitened space.
#'   \item Rank zero applies no residual correction. Positive and automatic ranks require rrpack.
#' }
#'
#' @section LOKO vs. single-fit adapters:
#' When `leave_one_key_out = TRUE`, REMAP performs a leave-one-key-out (LOKO)
#' loop over items to estimate adaptor diagnostics (e.g. per-voxel R2 and
#' improvement over naive) in an item-wise cross-validated fashion. This is
#' computationally heavy but yields less optimistic adaptor-level summaries.
#' Without a calibration split, the single-fit mode also scores the target
#' rows used to estimate its prototypes. In key-wise mode the observed key
#' selects the adapter, so its trial scores are conditional on that key.
#' For independent prediction, use calibration_rows or \code{\link{fit_remap}}.
#' @examples
#' ds <- gen_sample_dataset(c(5, 5, 5), 20, external_test = TRUE)
#' model <- remap_rrr_model(ds$dataset, ds$design, rank = 5)
#' class(model)
#' @export
remap_rrr_model <- function(dataset,
                            design,
                            base_classifier = "corclass",
                            link_by = NULL,
                            center = TRUE,
                            shrink_whiten = TRUE,
                            rank = "auto",
                            max_rank = 20,
                            ridge_rrr_lambda = NULL,
                            leave_one_key_out = FALSE,
                            min_pairs = 5,
                            save_fold_singulars = FALSE,
                            return_adapter = FALSE,
                            lambda_grid = c(0, 0.25, 0.5, 0.75, 1),
                            compute_naive_baseline = FALSE,
                            calibration_rows = NULL,
                            ...) {

  if (!is.null(calibration_rows)) {
    n <- nrow(design$test_design)
    if (is.null(n) || isTRUE(leave_one_key_out))
      stop("calibration_rows requires a target design and leave_one_key_out = FALSE.")
    if (is.logical(calibration_rows)) {
      if (length(calibration_rows) != n || anyNA(calibration_rows))
        stop("Logical calibration_rows must have one nonmissing entry per target row.")
      calibration_rows <- which(calibration_rows)
    }
    if (!is.numeric(calibration_rows) || !length(calibration_rows) ||
        any(!is.finite(calibration_rows)) || any(calibration_rows != trunc(calibration_rows)) ||
        any(calibration_rows < 1 | calibration_rows > n) || anyDuplicated(calibration_rows) ||
        length(calibration_rows) >= n)
      stop("calibration_rows must identify a nonempty proper subset of distinct target rows.")
  }

  # Decide performance function based on response type
  perf_fun <- if (is.numeric(design$y_train)) {
    get_regression_perf(design$split_groups)
  } else if (length(levels(design$y_train)) > 2) {
    get_multiclass_perf(design$split_groups, class_metrics = FALSE)
  } else {
    get_binary_perf(design$split_groups)
  }

  create_model_spec(
    "remap_rrr_model",
    dataset = dataset,
    design = design,
    base_classifier = base_classifier,
    link_by = link_by,
    calibration_rows = calibration_rows,
    .preserve_roi_features = !is.null(calibration_rows),
    center = center,
    shrink_whiten = shrink_whiten,
    rank = rank,
    max_rank = max_rank,
    ridge_rrr_lambda = ridge_rrr_lambda,
    leave_one_key_out = leave_one_key_out,
    min_pairs = min_pairs,
    save_fold_singulars = save_fold_singulars,
    return_adapter = return_adapter,
    lambda_grid = lambda_grid,
    compute_naive_baseline = compute_naive_baseline,
    performance = perf_fun,
    compute_performance = TRUE,
    return_predictions = TRUE,
    ...
  )
}


# --- helpers ---------------------------------------------------------------

# ----------------------------------------------------------------------------
# Residual REMAP helpers (joint whitening + low-rank correction)
# ----------------------------------------------------------------------------

#' Joint whitening on stacked X/Y prototypes to place them in common coords
#' Returns whitened fits for X and Y plus the whitening transform.
#' @return A list with whitened data matrices.
#' @keywords internal
.joint_whiten <- function(X_fit, Y_fit, do_shrink = TRUE, center = TRUE, eps = 1e-7) {
  Z <- rbind(X_fit, Y_fit)
  if (center) {
    mu <- colMeans(Z)
    Zc <- sweep(Z, 2, mu, "-")
  } else {
    mu <- rep(0, ncol(Z))
    Zc <- Z
  }
  S <- if (isTRUE(do_shrink) && nrow(Zc) == 2L) {
    # Two prototypes supply one covariance degree of freedom. Use the
    # spherical target directly; corpcor cannot estimate shrinkage here.
    diag(mean(matrixStats::colVars(Zc)), nrow = ncol(Zc))
  } else if (isTRUE(do_shrink) && ncol(Zc) > 1L) {
    require_package("corpcor", "for shrinkage whitening")
    corpcor::cov.shrink(Zc, verbose = FALSE)
  } else {
    stats::cov(Zc)
  }
  ee <- eigen(S, symmetric = TRUE)
  W  <- ee$vectors %*% diag(1 / sqrt(pmax(ee$values, eps)), nrow = ncol(Zc)) %*% t(ee$vectors)
  Zw <- Zc %*% W
  nX <- nrow(X_fit)
  list(
    Xw_fit = Zw[seq_len(nX), , drop = FALSE],
    Yw_fit = Zw[nX + seq_len(nrow(Y_fit)), , drop = FALSE],
    W = W,
    mu = mu
  )
}

#' Residual reduced-rank regression: (Yw - Xw) ~ Xw
#' @return A list with fitted residual RRR components.
#' @keywords internal
.fit_residual_rrr <- function(Xw, Yw, rank = "auto", max_rank = 20, ridge_lambda = NULL) {
  whole <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x) && x == trunc(x)
  auto <- identical(rank, "auto")
  if (!auto && (!whole(rank) || rank < 0)) stop("rank must be a nonnegative integer or 'auto'.")
  if (!whole(max_rank) || max_rank < 1L) stop("max_rank must be a positive integer.")
  if (!is.null(ridge_lambda) && (!is.numeric(ridge_lambda) || length(ridge_lambda) != 1L ||
      !is.finite(ridge_lambda) || ridge_lambda < 0)) stop("Invalid ridge penalty.")
  p <- ncol(Xw)
  if (!auto && rank == 0) {
    return(list(Delta = matrix(0, p, p), rank = 0L, Ad = numeric(0)))
  }
  require_package("rrpack", "for positive or automatic REMAP rank")
  stopifnot(nrow(Xw) == nrow(Yw))
  Rw <- Yw - Xw
  n <- nrow(Xw)
  feasible <- min(n, p, ncol(Yw))
  max_rank <- min(max_rank, feasible, n - 1L)
  if (!auto && rank > feasible) stop("Requested rank exceeds the training dimensions.")
  if (!is.null(ridge_lambda)) {
    fit <- rrpack::rrs.fit(Y = Rw, X = Xw, nrank = if (identical(rank, "auto")) max_rank else as.integer(rank),
                           lambda = ridge_lambda)
  } else if (identical(rank, "auto")) {
    # rrpack::cv.rrr drops singleton folds to vectors. For tiny designs,
    # use explicit leave-one-key-out fits with dimensions preserved.
    if (n < 4L || p == 1L) {
      errors <- numeric(max_rank + 1L)
      errors[1L] <- sum(Rw^2)
      for (r in seq_len(max_rank)) {
        for (i in seq_len(n)) {
          inner <- rrpack::rrr.fit(Rw[-i, , drop = FALSE], Xw[-i, , drop = FALSE], nrank = r)
          residual <- Rw[i, , drop = FALSE] - Xw[i, , drop = FALSE] %*% inner$coef
          errors[r + 1L] <- errors[r + 1L] + sum(residual^2)
        }
      }
      return(.fit_residual_rrr(Xw, Yw, which.min(errors) - 1L, max_rank))
    }
    folds <- min(5L, max(2L, floor(n / 3)))
    largest_fold <- max(round(n / folds), n - (folds - 1L) * round(n / folds))
    fit <- rrpack::cv.rrr(Y = Rw, X = Xw, nfold = folds, maxrank = min(max_rank, n - largest_fold))
  } else {
    fit <- rrpack::rrr.fit(Y = Rw, X = Xw, nrank = as.integer(rank))
  }

  Delta <- fit$coef
  if (!is.matrix(Delta) || !identical(dim(Delta), c(p, ncol(Yw))) || any(!is.finite(Delta)))
    stop("The REMAP solver returned invalid coefficient estimates.")
  list(Delta = Delta, rank = fit$rank %||% fit$nrank, Ad = as.numeric(fit$Ad))
}

#' Select lambda by prototype self-retrieval on training items
#' @return A numeric scalar: the selected lambda value.
#' @keywords internal
.select_lambda <- function(Xw_fit, Yw_fit, Delta, lambda_grid) {
  if (length(lambda_grid) <= 1L) return(list(lambda = lambda_grid, accs = 1))
  n_train <- nrow(Xw_fit)
  accs <- numeric(length(lambda_grid))
  for (i in seq_along(lambda_grid)) {
    lam <- lambda_grid[i]
    Yhat_train <- Xw_fit + lam * (Xw_fit %*% Delta)
    scores <- .remap_scores(Yw_fit, Yhat_train)
    accs[i] <- mean(max.col(scores, ties.method = "first") == seq_len(n_train))
  }
  best <- which(accs == max(accs))
  lam_opt <- min(lambda_grid[best])
  list(lambda = lam_opt, accs = accs)
}

#' @keywords internal
.remap_paired <- function(train_mat, test_mat, train_des, test_des, link_by, y_train, y_test) {
  if (!is.null(link_by) && link_by %in% colnames(train_des) && link_by %in% colnames(test_des)) {
    key_tr <- factor(train_des[[link_by]])
    key_te <- factor(test_des[[link_by]])
  } else {
    key_tr <- factor(y_train)
    key_te <- factor(y_test)
  }
  # source (train) prototypes
  Xp <- .nx_rowsum_mean(train_mat, key_tr)
  # target (test) prototypes (average any repeats per key)
  Ym <- .nx_rowsum_mean(test_mat, key_te)
  # Align common keys only
  common <- intersect(rownames(Xp), rownames(Ym))
  if (length(common) == 0L) {
    return(list(Xp = matrix(0, 0, ncol(train_mat)), Ym = matrix(0, 0, ncol(test_mat)), keys = character(0)))
  }
  Xp <- Xp[common, , drop = FALSE]
  Ym <- Ym[common, , drop = FALSE]
  list(Xp = Xp, Ym = Ym, keys = common)
}


# --- S3 methods ------------------------------------------------------------

#' @export
compute_performance.remap_rrr_model <- function(obj, result) {
  # Delegate to performance function created in constructor
  obj$performance(result)
}


#' @export
output_schema.remap_rrr_model <- function(model) {
  NULL
}


#' Per-ROI REMAP fit
#'
#' Implements the end-to-end REMAP flow for one ROI/searchlight: build paired
#' prototypes, fit the residual transformation, and classify target patterns
#' by correlation with the transformed source templates.
#' Emits a standard classification_result plus ROI-level performance.
#'
#' @return An \code{roi_result} object.
#' @keywords internal
#' @export
fit_roi.remap_rrr_model <- function(model, roi_data, context, ...) {
  Xtrain <- roi_data$train_data  # n_train x p
  ind    <- roi_data$indices
  id     <- context$id

  # We need test data for forward transfer
  if (is.null(roi_data$test_data)) {
    return(roi_result(
      metrics = NULL, indices = ind, id = id,
      error = TRUE,
      error_message = "remap_rrr_model requires an external test set (dataset$test_data + design$y_test)"
    ))
  }

  Xtest <- roi_data$test_data   # n_test x q (q==p typically)
  des   <- model$design
  ytr   <- y_train(des)
  yte   <- y_test(des)

  # Guard: ensure non-degenerate ROI
  if (ncol(Xtrain) < 1L || ncol(Xtest) < 1L) {
    return(roi_result(
      metrics = NULL, indices = ind, id = id,
      error = TRUE,
      error_message = "ROI has fewer than 1 feature"
    ))
  }

  .log_debug("REMAP ROI %s: start (train_dim=%s x %s, test_dim=%s x %s)",
                            id, nrow(Xtrain), ncol(Xtrain), nrow(Xtest), ncol(Xtest))

  if (!is.null(model$calibration_rows)) {
    if (!is.null(model$link_by)) {
      if (!model$link_by %in% names(des$train_design) || !model$link_by %in% names(des$test_design))
        stop("link_by must identify a column in both designs.")
      ytr <- des$train_design[[model$link_by]]
      yte <- des$test_design[[model$link_by]]
    }
    return(.remap_assess(model, Xtrain, ytr, Xtest, yte, ind, id))
  }

  # Build paired prototypes and key factors
  pairs <- .remap_paired(Xtrain, Xtest, des$train_design, des$test_design,
                         model$link_by, ytr, yte)
  Xp <- pairs$Xp; Ym <- pairs$Ym
  n_pairs <- nrow(Xp)
  .log_debug("REMAP ROI %s: paired prototypes = %s", id, n_pairs)
  # Hard error if fewer than 2 paired items
  if (n_pairs < 2L || nrow(Ym) < 2L) {
    return(roi_result(
      metrics = NULL, indices = ind, id = id,
      error = TRUE,
      error_message = "REMAP: Need at least 2 paired items (prototypes) to fit adapter"
    ))
  }
  # Soft warning if below recommended minimum
  soft_warn <- n_pairs < (model$min_pairs %||% 5L)

  # Determine keys for prototypes and observed labels
  if (!is.null(model$link_by) && model$link_by %in% colnames(des$train_design) && model$link_by %in% colnames(des$test_design)) {
    key_tr_full <- factor(des$train_design[[model$link_by]])
    key_te_full <- factor(des$test_design[[model$link_by]])
  } else {
    key_tr_full <- factor(ytr)
    key_te_full <- factor(yte)
  }

  # We'll predict over the common item set; set levels accordingly for probs
  unique_keys <- intersect(rownames(Xp), rownames(Ym))
  levs <- unique_keys
  prob_all <- matrix(0, nrow(Xtest), length(levs), dimnames = list(NULL, levs))
  pred_all <- rep(NA_character_, nrow(Xtest))
  use_naive <- isTRUE(model$compute_naive_baseline)
  if (use_naive) {
    prob_naive <- matrix(0, nrow(Xtest), length(levs), dimnames = list(NULL, levs))
    pred_naive <- rep(NA_character_, nrow(Xtest))
  }

  # CV R^2 accumulator per voxel
  sse_vox <- rep(0, ncol(Xtrain))
  ss_vox  <- rep(0, ncol(Xtrain))

  use_loko <- isTRUE(model$leave_one_key_out) && length(unique_keys) >= 3L
  .log_debug("REMAP ROI %s: use_loko=%s, unique_keys=%s", id, use_loko, length(unique_keys))

  adapter_svals <- c()
  adapter_ranks <- c()
  skipped_keys <- character(0)
  keys_scored  <- character(0)
  # Residual diagnostics accumulators (LOKO folds)
  roi_improv_vals <- c()
  delta_frob_vals <- c()
  lambda_vals     <- c()

  if (use_loko) {
    # Collect per-item diagnostics
    resid_by_item <- setNames(numeric(0), character(0))
    rank_by_item  <- setNames(numeric(0), character(0))
    sv1_by_item   <- setNames(numeric(0), character(0))
    sv_spectra_by_item <- if (isTRUE(model$save_fold_singulars)) vector("list", length(unique_keys)) else NULL
    jj_map <- list()

    for (k in unique_keys) {
      keep <- setdiff(unique_keys, k)
      if (length(keep) < 2L) { skipped_keys <- c(skipped_keys, k); next }
      Xp_k <- Xp[keep, , drop = FALSE]
      Ym_k <- Ym[keep, , drop = FALSE]

      fit <- .remap_model_templates(model, Xp, Ym_k)
      JW <- list(mu = fit$mu, W = fit$W,
                 Xw_fit = sweep(Xp_k, 2L, fit$mu, "-") %*% fit$W,
                 Yw_fit = sweep(Ym_k, 2L, fit$mu, "-") %*% fit$W)
      Delta <- fit$Delta
      lam_opt <- fit$lambda

      # ---- diagnostics in whitened space ----
      R_naive <- JW$Yw_fit - JW$Xw_fit
      Yhat_fit <- JW$Xw_fit + lam_opt * (JW$Xw_fit %*% Delta)
      R_remap  <- JW$Yw_fit - Yhat_fit
      ss_naive <- sum(R_naive^2)
      ss_remap <- sum(R_remap^2)
      roi_improv <- if (ss_naive > .Machine$double.eps) 1 - (ss_remap / ss_naive) else NA_real_
      item_res_naive <- sqrt(rowSums(R_naive^2))
      item_res_remap <- sqrt(rowSums(R_remap^2))
      delta_frob <- sqrt(sum((lam_opt * Delta)^2))

      # accumulate fold-level stats
      roi_improv_vals <- c(roi_improv_vals, roi_improv)
      delta_frob_vals <- c(delta_frob_vals, delta_frob)
      lambda_vals     <- c(lambda_vals, lam_opt)

      rows_k <- which(as.character(key_te_full) == k)
      if (length(rows_k)) {
        test <- Xtest[rows_k, , drop = FALSE]
        prob_all[rows_k, ] <- predict(fit, test, type = "prob")
        pred_all[rows_k] <- unique_keys[max.col(prob_all[rows_k, , drop = FALSE], ties.method = "first")]
        if (use_naive) {
          raw <- sweep(Xp, 2L, fit$mu, "-") %*% fit$W
          yw <- sweep(test, 2L, fit$mu, "-") %*% fit$W
          prob_naive[rows_k, ] <- .nx_softmax(.remap_scores(yw, raw))
          pred_naive[rows_k] <- unique_keys[max.col(prob_naive[rows_k, , drop = FALSE], ties.method = "first")]
        }
        keys_scored <- c(keys_scored, k)
      }

      # CV R^2 for left-out item in voxel space
      if (k %in% rownames(Xp) && k %in% rownames(Ym)) {
        # Unwhiten predicted template for held-out key
        Winv <- try(solve(JW$W), silent = TRUE); if (inherits(Winv, "try-error")) Winv <- MASS::ginv(JW$W)
        # Recompute Xw for single row k using training JW
        Xw_k <- (Xp[k, , drop = FALSE] - matrix(JW$mu, nrow = 1)) %*% JW$W
        Yhat_k_vox <- Xw_k %*% (diag(ncol(Delta)) + lam_opt * Delta) %*% Winv + matrix(JW$mu, nrow = 1)
        resid  <- Ym[k, , drop = FALSE] - Yhat_k_vox
        sse_vox <- sse_vox + as.numeric(resid^2)
        Ym_center <- Ym[k, , drop = FALSE] - matrix(colMeans(Ym_k), nrow = 1)
        ss_vox  <- ss_vox + as.numeric(Ym_center^2)
        resid_by_item[k] <- mean(as.numeric(resid)^2)
      }

      # Diagnostics
      rk <- fit$rank
      sv <- fit$singular_values
      adapter_ranks <- c(adapter_ranks, rk)
      if (!is.null(sv) && length(sv)) adapter_svals <- c(adapter_svals, sv[1])
      rank_by_item[k] <- rk
      sv1_by_item[k]  <- if (!is.null(sv) && length(sv)) sv[1] else NA_real_
      if (isTRUE(model$save_fold_singulars)) {
        if (is.null(jj_map[[k]])) jj_map[[k]] <- length(jj_map) + 1L
        sv_spectra_by_item[[jj_map[[k]]]] <- sv
        names(sv_spectra_by_item)[jj_map[[k]]] <- k
      }

      # Optional per-fold diagnostics package (only if return_adapter)
      if (isTRUE(model$return_adapter)) {
        # create list-once
        if (!exists("remap_diag_by_fold", inherits = FALSE)) remap_diag_by_fold <- list()
        remap_diag_by_fold[[length(remap_diag_by_fold) + 1L]] <- list(
          heldout_item   = k,
          rank_used      = rk,
          lambda_used    = lam_opt,
          roi_improv     = roi_improv,
          delta_frob     = delta_frob,
          train_items    = keep,
          item_res_naive = item_res_naive,
          item_res_remap = item_res_remap
        )
      }
    }

    # Compute per-voxel CV R^2
    ss_vox[ss_vox < .Machine$double.eps] <- .Machine$double.eps
    r2_cv_vox <- 1 - (sse_vox / ss_vox)

    prob <- prob_all
    pred <- factor(pred_all, levels = levs)
    if (any(is.na(pred))) {
      mc <- max.col(prob, ties.method = "first")
      pred <- factor(levs[mc], levels = levs)
    }
    if (use_naive) {
      prob_naive_full <- prob_naive
      pred_naive_fac <- factor(pred_naive, levels = levs)
      if (any(is.na(pred_naive_fac))) {
        mc_naive <- max.col(prob_naive_full, ties.method = "first")
        pred_naive_fac <- factor(levs[mc_naive], levels = levs)
      }
    }

    # Observed labels are the key factor
    obs <- factor(as.character(key_te_full), levels = levs)

    predictor_obj <- NULL
    if (isTRUE(model$return_adapter)) {
      predictor_obj <- list(r2_per_voxel = r2_cv_vox,
                            resid_by_item = resid_by_item,
                            rank_by_item  = rank_by_item,
                            sv1_by_item   = sv1_by_item,
                            keys_scored   = keys_scored,
                            skipped_keys  = skipped_keys)
      if (isTRUE(model$save_fold_singulars)) predictor_obj$sv_spectra_by_item <- sv_spectra_by_item
      if (exists("remap_diag_by_fold", inherits = FALSE)) predictor_obj$diag_by_fold <- remap_diag_by_fold
    }
    # filter out test rows whose observed key not in training/common set
    keep_rows <- which(!is.na(obs))
    cres <- classification_result(obs[keep_rows], pred[keep_rows], prob[keep_rows, , drop = FALSE],
                                  testind = keep_rows,
                                  test_design = des$test_design[keep_rows, , drop = FALSE],
                                  predictor = predictor_obj)
    base_perf <- compute_performance(model, cres)
    base_perf_naive <- NULL
    if (use_naive) {
      cres_naive <- classification_result(obs[keep_rows], pred_naive_fac[keep_rows],
                                          prob_naive_full[keep_rows, , drop = FALSE],
                                          testind = keep_rows,
                                          test_design = des$test_design[keep_rows, , drop = FALSE],
                                          predictor = NULL)
      base_perf_naive <- compute_performance(model, cres_naive)
      base_perf_naive <- setNames(base_perf_naive, paste0("naive_", names(base_perf_naive)))
    }
    extra <- c(adapter_rank    = suppressWarnings(mean(adapter_ranks, na.rm = TRUE)),
               adapter_sv1     = suppressWarnings(mean(adapter_svals, na.rm = TRUE)),
               adapter_mean_r2 = mean(r2_cv_vox, na.rm = TRUE),
               remap_improv    = suppressWarnings(mean(roi_improv_vals, na.rm = TRUE)),
               delta_frob_mean = suppressWarnings(mean(delta_frob_vals, na.rm = TRUE)),
               lambda_mean     = suppressWarnings(mean(lambda_vals, na.rm = TRUE)))
    perf_vec <- if (use_naive) c(base_perf, base_perf_naive, extra) else c(base_perf, extra)

    return(roi_result(
      metrics = c(perf_vec, n_pairs_used = n_pairs, n_skipped_keys = length(skipped_keys)),
      indices = ind,
      id = id,
      result = cres
    ))
  }

  # Non-LOKO single-fit path using joint whitening + residual RRR
  fit <- .remap_model_templates(model, Xp, Ym)
  JW <- list(mu = fit$mu, W = fit$W,
             Xw_fit = sweep(Xp, 2L, fit$mu, "-") %*% fit$W,
             Yw_fit = sweep(Ym, 2L, fit$mu, "-") %*% fit$W)
  Delta <- fit$Delta
  lam_opt <- fit$lambda

  # Residual diagnostics (global fit)
  R_naive <- JW$Yw_fit - JW$Xw_fit
  Yhat_fit <- JW$Xw_fit + lam_opt * (JW$Xw_fit %*% Delta)
  R_remap  <- JW$Yw_fit - Yhat_fit
  ss_naive <- sum(R_naive^2)
  ss_remap <- sum(R_remap^2)
  roi_improv <- if (ss_naive > .Machine$double.eps) 1 - (ss_remap / ss_naive) else NA_real_
  delta_frob <- sqrt(sum((lam_opt * Delta)^2))
  .log_debug("REMAP ROI %s: roi_improv=%.3f, delta_frob=%.3f", id, roi_improv, delta_frob)

  # Templates for all keys and correlation-based classification
  Xw_all <- JW$Xw_fit
  prob <- predict(fit, Xtest, type = "prob")
  pred <- factor(unique_keys[max.col(prob, ties.method = "first")], levels = unique_keys)
  if (use_naive) {
    yw <- sweep(Xtest, 2L, fit$mu, "-") %*% fit$W
    prob_naive <- .nx_softmax(.remap_scores(yw, Xw_all))
    pred_naive <- factor(unique_keys[max.col(prob_naive, ties.method = "first")], levels = unique_keys)
  }
  obs <- factor(as.character(key_te_full), levels = unique_keys)

  # Diagnostics on prototypes in voxel space
  Winv <- try(solve(JW$W), silent = TRUE); if (inherits(Winv, "try-error")) Winv <- MASS::ginv(JW$W)
  Yhat_vox <- (Xw_all %*% (diag(ncol(Delta)) + lam_opt * Delta)) %*% Winv + matrix(JW$mu, nrow = nrow(Xw_all), ncol = ncol(Xw_all), byrow = TRUE)
  Ym_centered <- scale(Ym[unique_keys, , drop = FALSE], center = TRUE, scale = FALSE)
  denom <- colSums(Ym_centered^2); denom[denom < .Machine$double.eps] <- .Machine$double.eps
  r2 <- 1 - colSums((Ym[unique_keys, , drop = FALSE] - Yhat_vox)^2) / denom
  rnk <- fit$rank
  svals <- fit$singular_values

  keep_rows <- which(!is.na(obs))
  cres <- classification_result(obs[keep_rows], pred[keep_rows], prob[keep_rows, , drop = FALSE],
                                testind = keep_rows,
                                test_design = des$test_design[keep_rows, , drop = FALSE],
                                predictor = if (isTRUE(model$return_adapter)) list(rank = rnk, singvals = svals, r2_per_voxel = r2) else NULL)
  base_perf <- compute_performance(model, cres)
  base_perf_naive <- NULL
  if (use_naive) {
    cres_naive <- classification_result(obs[keep_rows], pred_naive[keep_rows], prob_naive[keep_rows, , drop = FALSE],
                                        testind = keep_rows,
                                        test_design = des$test_design[keep_rows, , drop = FALSE],
                                        predictor = NULL)
    base_perf_naive <- compute_performance(model, cres_naive)
    base_perf_naive <- setNames(base_perf_naive, paste0("naive_", names(base_perf_naive)))
  }
  extra <- c(adapter_rank    = as.numeric(rnk %||% NA_real_),
             adapter_sv1     = if (!is.null(svals) && length(svals) > 0) svals[1] else NA_real_,
             adapter_mean_r2 = mean(r2, na.rm = TRUE),
             remap_improv    = roi_improv,
             delta_frob_mean = delta_frob,
             lambda_mean     = lam_opt)
  perf_vec <- if (use_naive) c(base_perf, base_perf_naive, extra) else c(base_perf, extra)

  roi_result(
    metrics = c(perf_vec, n_pairs_used = n_pairs, n_skipped_keys = 0L),
    indices = ind,
    id = id,
    result = cres
  )
}



.remap_model_templates <- function(model, source, target) {
  .remap_fit_templates(source, target, model$rank, model$max_rank,
                       model$ridge_rrr_lambda, model$lambda_grid,
                       model$center, model$shrink_whiten)
}

# Assess one user-specified target split. Only the final performance calculation
# sees assessment labels. Keep original target indices for prediction ledgers.
.remap_assess <- function(model, source, source_key, target, target_key, indices, id) {
  cal <- model$calibration_rows
  test <- setdiff(seq_len(nrow(target)), cal)
  fit <- fit_remap(source, source_key, target[cal, , drop = FALSE], target_key[cal],
                   rank = model$rank, max_rank = model$max_rank,
                   ridge_rrr_lambda = model$ridge_rrr_lambda, lambda_grid = model$lambda_grid,
                   center = model$center, shrink_whiten = model$shrink_whiten)
  assessment <- target[test, , drop = FALSE]
  measured <- function(fitted) {
    prob <- predict(fitted, assessment, type = "prob")
    keys <- colnames(prob)
    predicted <- factor(keys[max.col(prob, ties.method = "first")], levels = keys)
    observed <- factor(as.character(target_key[test]), levels = keys)
    if (anyNA(observed)) stop("Assessment labels must refer to known source keys.")
    classification_result(observed, predicted, prob, testind = test,
                          test_design = model$design$test_design[test, , drop = FALSE],
                          predictor = if (isTRUE(model$return_adapter)) fitted else NULL)
  }
  result <- measured(fit)
  metrics <- compute_performance(model, result)
  if (isTRUE(model$compute_naive_baseline)) {
    reference <- fit_remap(source, source_key, center = model$center,
                           shrink_whiten = model$shrink_whiten)
    baseline <- compute_performance(model, measured(reference))
    metrics <- c(metrics, setNames(baseline, paste0("source_only_", names(baseline))))
  }
  roi_result(metrics = c(metrics, adapter_rank = fit$rank, lambda_mean = fit$lambda,
                         n_pairs_used = length(fit$calibration_keys)),
             indices = indices, id = id, result = result)
}
