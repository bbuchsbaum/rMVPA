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
#' method automatically shrinks to the naïve cross-decoder (\eqn{\lambda = 0}).
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
#' @param forward_adapt Logical; kept for compatibility, currently ignored.
#' @param ridge_rrr_lambda Optional numeric lambda for `rrpack::rrs.fit` (ridge RRR) on the residuals.
#' @param leave_one_key_out logical; if TRUE, use leave-one-key-out (LOKO) over items for adaptor learning (default TRUE).
#' @param min_pairs Minimum number of paired prototypes required for adaptor fitting (default 5).
#' @param save_fold_singulars logical; if TRUE, save singular values from each fold's adaptor (default FALSE).
#' @param return_adapter logical; if TRUE, store diagnostics (e.g., per-voxel R2, per-fold stats) in `result$predictor`.
#' @param lambda_grid Numeric vector of candidate \eqn{\lambda} values for shrinkage toward naïve; defaults to `c(0, .25, .5, .75, 1)`.
#' @param return_diag Logical; reserved for future detailed diagnostics (currently unused).
#' @param ... Additional arguments passed to `create_model_spec`.
#'
#' @return A model spec of class `remap_rrr_model` compatible with `run_regional()` / `run_searchlight()`.
#'
#' @section Key ideas:
#' \itemize{
#'   \item Joint whitening (shared covariance) puts perception and memory in the same coordinates.
#'   \item Residual RRR learns a low-rank correction \eqn{\Delta} on \eqn{R_w = Y_w - X_w}.
#'   \item Automatic shrinkage to naïve via \eqn{\lambda} selection on training items.
#'   \item Prediction/classification by correlation to \eqn{\hat Y_w} in the whitened space.
#'   \item Fallback to naïve cross-decoding when `rank = 0` or `rrpack` is unavailable.
#' }
#' @export
remap_rrr_model <- function(dataset,
                            design,
                            base_classifier = "corclass",
                            link_by = NULL,
                            center = TRUE,
                            shrink_whiten = TRUE,
                            rank = "auto",
                            max_rank = 20,
                            forward_adapt = TRUE,
                            ridge_rrr_lambda = NULL,
                            leave_one_key_out = TRUE,
                            min_pairs = 5,
                            save_fold_singulars = FALSE,
                            return_adapter = TRUE,
                            lambda_grid = c(0, 0.25, 0.5, 0.75, 1),
                            return_diag = FALSE,
                            ...) {

  # Decide performance function based on response type
  perf_fun <- if (is.numeric(design$y_train)) {
    get_regression_perf(design$split_groups)
  } else if (length(levels(design$y_train)) > 2) {
    get_multiclass_perf(design$split_groups, class_metrics = TRUE)
  } else {
    get_binary_perf(design$split_groups)
  }

  create_model_spec(
    "remap_rrr_model",
    dataset = dataset,
    design = design,
    base_classifier = base_classifier,
    link_by = link_by,
    center = center,
    shrink_whiten = shrink_whiten,
    rank = rank,
    max_rank = max_rank,
    forward_adapt = forward_adapt,
    ridge_rrr_lambda = ridge_rrr_lambda,
    leave_one_key_out = leave_one_key_out,
    min_pairs = min_pairs,
    save_fold_singulars = save_fold_singulars,
    return_adapter = return_adapter,
    lambda_grid = lambda_grid,
    return_diag = return_diag,
    performance = perf_fun,
    compute_performance = TRUE,
    return_predictions = TRUE,
    ...
  )
}


# --- helpers ---------------------------------------------------------------

#' @keywords internal
#' @importFrom stats cov
.remap_whiten <- function(X, do_shrink = TRUE, center = TRUE) {
  if (center) X <- scale(X, center = TRUE, scale = FALSE)
  S <- try({
    if (!isTRUE(do_shrink)) {
      cov(X)
    } else {
      # Suppress verbose console output from corpcor during shrinkage
      corpcor::cov.shrink(X, verbose = FALSE)
    }
  }, silent = TRUE)
  if (inherits(S, "try-error")) S <- cov(X)
  E <- eigen(S, symmetric = TRUE)
  # ZCA whitening matrix
  W <- E$vectors %*% diag(1 / sqrt(pmax(E$values, .Machine$double.eps))) %*% t(E$vectors)
  list(Xw = X %*% W, W = W, mean = attr(X, "scaled:center"))
}

#' @keywords internal
.remap_fit_rrr <- function(Xw, Yw, rank = "auto", max_rank = 20, ridge_lambda = NULL) {
  stopifnot(nrow(Xw) == nrow(Yw))
  n <- nrow(Xw); p <- ncol(Xw); q <- ncol(Yw)
  max_rank <- max(0L, min(max_rank, p, q, max(1L, n - 1L)))

  if (is.numeric(rank) && as.integer(rank) == 0L) {
    return(list(method = "identity", Cw = matrix(0, nrow = p, ncol = q), rank = 0L, Ad = numeric(0)))
  }

  if (!requireNamespace("rrpack", quietly = TRUE)) {
    return(list(method = "unavailable", Cw = matrix(0, nrow = p, ncol = q), rank = 0L, Ad = numeric(0)))
  }

  if (!is.null(ridge_lambda)) {
    fit <- rrpack::rrs.fit(Yw, Xw, nrank = if (identical(rank, "auto")) max_rank else as.integer(rank),
                           lambda = ridge_lambda)
  } else if (identical(rank, "auto")) {
    fit <- rrpack::cv.rrr(Y = Yw, X = Xw, nfold = min(5, max(2, floor(n/3))), maxrank = max_rank)
  } else {
    fit <- rrpack::rrr.fit(Y = Yw, X = Xw, nrank = as.integer(rank))
  }
  fit
}

# ----------------------------------------------------------------------------
# Residual REMAP helpers (joint whitening + low-rank correction)
# ----------------------------------------------------------------------------

#' @keywords internal
.row_center <- function(X) sweep(X, 2, colMeans(X), "-")

#' Joint whitening on stacked X/Y prototypes to place them in common coords
#' Returns whitened fits for X and Y plus the whitening transform.
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
  S <- try({
    if (!isTRUE(do_shrink)) stats::cov(Zc) else corpcor::cov.shrink(Zc)
  }, silent = TRUE)
  if (inherits(S, "try-error")) S <- stats::cov(Zc)
  ee <- eigen(S, symmetric = TRUE)
  W  <- ee$vectors %*% diag(1 / sqrt(pmax(ee$values, eps))) %*% t(ee$vectors)
  Zw <- Zc %*% W
  nX <- nrow(X_fit)
  list(
    Xw_fit = Zw[seq_len(nX), , drop = FALSE],
    Yw_fit = Zw[(nX + 1):nrow(Zw), , drop = FALSE],
    W = W,
    mu = mu
  )
}

#' Row-wise correlation of matrix rows with a vector
#' @keywords internal
.row_cor <- function(A, v) {
  A0 <- .row_center(A)
  v0 <- v - mean(v)
  den <- sqrt(rowSums(A0^2)) * sqrt(sum(v0^2))
  as.numeric((A0 %*% v0) / pmax(den, .Machine$double.eps))
}

#' Residual reduced-rank regression: (Yw - Xw) ~ Xw
#' @keywords internal
.fit_residual_rrr <- function(Xw, Yw, rank = "auto", max_rank = 20, ridge_lambda = NULL) {
  stopifnot(nrow(Xw) == nrow(Yw))
  Rw <- Yw - Xw
  n <- nrow(Xw); p <- ncol(Xw); q <- ncol(Yw)
  max_rank <- max(1L, min(max_rank, p, q, max(1L, n - 1L)))

  if (is.numeric(rank) && as.integer(rank) <= 0L) {
    return(list(Delta = matrix(0, nrow = p, ncol = q), rank = 0L, Ad = numeric(0)))
  }
  if (!requireNamespace("rrpack", quietly = TRUE)) {
    return(list(Delta = matrix(0, nrow = p, ncol = q), rank = 0L, Ad = numeric(0)))
  }

  if (!is.null(ridge_lambda)) {
    fit <- rrpack::rrs.fit(Y = Rw, X = Xw, nrank = if (identical(rank, "auto")) max_rank else as.integer(rank),
                           lambda = ridge_lambda)
  } else if (identical(rank, "auto")) {
    fit <- rrpack::cv.rrr(Y = Rw, X = Xw, nfold = min(5, max(2, floor(n/3))), maxrank = max_rank)
  } else {
    fit <- rrpack::rrr.fit(Y = Rw, X = Xw, nrank = as.integer(rank))
  }

  # Try to extract coefficients and leading singulars consistently
  Delta <- try({
    if (!is.null(fit$coef)) fit$coef else stats::coef(fit)
  }, silent = TRUE)
  if (inherits(Delta, "try-error")) Delta <- matrix(0, nrow = p, ncol = q)
  Ad <- tryCatch(as.numeric(fit$Ad), error = function(...) numeric(0))
  r_used <- tryCatch({ if (!is.null(fit$rank)) fit$rank else if (length(Ad)) length(Ad) else NA_integer_ }, error = function(...) NA_integer_)
  list(Delta = Delta, rank = r_used, Ad = Ad)
}

#' Select lambda by prototype self-retrieval on training items
#' @keywords internal
.select_lambda <- function(Xw_fit, Yw_fit, Delta, lambda_grid) {
  if (length(lambda_grid) <= 1L) return(list(lambda = lambda_grid, accs = 1))
  n_train <- nrow(Xw_fit)
  accs <- numeric(length(lambda_grid))
  for (i in seq_along(lambda_grid)) {
    lam <- lambda_grid[i]
    Yhat_train <- Xw_fit + lam * (Xw_fit %*% Delta)
    correct <- logical(n_train)
    for (j in seq_len(n_train)) {
      scores <- .row_cor(Yhat_train, Yw_fit[j, ])
      pred   <- which.max(scores)
      correct[j] <- (pred == j)
    }
    accs[i] <- mean(correct)
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
  Xp <- rowsum(train_mat, key_tr, reorder = TRUE) / as.vector(table(key_tr))
  # target (test) prototypes (average any repeats per key)
  Ym <- rowsum(test_mat, key_te, reorder = TRUE) / pmax(1, as.vector(table(key_te)))
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


#' Per-ROI REMAP processing
#'
#' Implements the end-to-end REMAP flow for one ROI/searchlight: build paired
#' prototypes, whiten, fit reduced-rank mapping, transform source trials, train
#' base classifier in the target domain, and classify target-domain trials.
#' Emits a standard classification_result plus ROI-level performance.
#'
#' @keywords internal
#' @export
process_roi.remap_rrr_model <- function(mod_spec, roi, rnum, ...) {
  # Pull ROI matrices directly from the ROI container (dataset is stripped in workers)
  Xtrain <- as.matrix(neuroim2::values(roi$train_roi))  # n_train x p
  ind    <- neuroim2::indices(roi$train_roi)

  # We need test data for forward transfer
  if (!has_test_set(mod_spec)) {
    return(tibble::tibble(
      result = list(NULL), indices = list(ind), performance = list(NULL), id = rnum,
      error = TRUE, error_message = "remap_rrr_model requires an external test set (dataset$test_data + design$y_test)"
    ))
  }

  Xtest <- as.matrix(neuroim2::values(roi$test_roi))   # n_test  x q (q==p typically)
  des   <- mod_spec$design
  ytr   <- y_train(des)
  yte   <- y_test(des)

  # Guard: ensure non-degenerate ROI
  if (ncol(Xtrain) < 2L || ncol(Xtest) < 2L) {
    return(tibble::tibble(
      result = list(NULL), indices = list(ind), performance = list(NULL), id = rnum,
      error = TRUE, error_message = "ROI has fewer than 2 features"
    ))
  }

  # Build paired prototypes and key factors
  pairs <- .remap_paired(Xtrain, Xtest, des$train_design, des$test_design,
                         mod_spec$link_by, ytr, yte)
  Xp <- pairs$Xp; Ym <- pairs$Ym
  n_pairs <- nrow(Xp)
  # Hard error if fewer than 2 paired items
  if (n_pairs < 2L || nrow(Ym) < 2L) {
    return(tibble::tibble(
      result = list(NULL), indices = list(ind), performance = list(NULL), id = rnum,
      warning = FALSE, warning_message = "~",
      error = TRUE, error_message = "REMAP: Need at least 2 paired items (prototypes) to fit adapter"
    ))
  }
  # Soft warning if below recommended minimum
  soft_warn <- n_pairs < (mod_spec$min_pairs %||% 5L)

  # Determine keys for prototypes and observed labels
  if (!is.null(mod_spec$link_by) && mod_spec$link_by %in% colnames(des$train_design) && mod_spec$link_by %in% colnames(des$test_design)) {
    key_tr_full <- factor(des$train_design[[mod_spec$link_by]])
    key_te_full <- factor(des$test_design[[mod_spec$link_by]])
  } else {
    key_tr_full <- factor(ytr)
    key_te_full <- factor(yte)
  }

  # We'll predict over the common item set; set levels accordingly for probs
  unique_keys <- intersect(rownames(Xp), rownames(Ym))
  levs <- unique_keys
  prob_all <- matrix(0, nrow(Xtest), length(levs), dimnames = list(NULL, levs))
  pred_all <- rep(NA_character_, nrow(Xtest))

  # CV R^2 accumulator per voxel
  sse_vox <- rep(0, ncol(Xtrain))
  ss_vox  <- rep(0, ncol(Xtrain))

  use_loko <- isTRUE(mod_spec$leave_one_key_out) && length(unique_keys) >= 3L

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
    sv_spectra_by_item <- if (isTRUE(mod_spec$save_fold_singulars)) vector("list", length(unique_keys)) else NULL
    jj_map <- list()

    for (k in unique_keys) {
      keep <- setdiff(unique_keys, k)
      if (length(keep) < 2L) { skipped_keys <- c(skipped_keys, k); next }
      Xp_k <- Xp[keep, , drop = FALSE]
      Ym_k <- Ym[keep, , drop = FALSE]

      # Joint whitening and residual RRR on training keys
      JW <- .joint_whiten(Xp_k, Ym_k, do_shrink = mod_spec$shrink_whiten, center = mod_spec$center)
      rf <- .fit_residual_rrr(JW$Xw_fit, JW$Yw_fit, rank = mod_spec$rank,
                              max_rank = mod_spec$max_rank,
                              ridge_lambda = mod_spec$ridge_rrr_lambda)
      Delta <- rf$Delta

      # Lambda selection on training keys
      lam_res <- .select_lambda(JW$Xw_fit, JW$Yw_fit, Delta, mod_spec$lambda_grid)
      lam_opt <- lam_res$lambda

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

      # Build predicted templates for all items in whitened space
      Xw_all <- (sweep(Xp[unique_keys, , drop = FALSE], 2, JW$mu, "-")) %*% JW$W
      Yhat_all <- Xw_all + lam_opt * (Xw_all %*% Delta)

      # Predict rows in test belonging to key k via correlation in whitened space
      rows_k <- which(as.character(key_te_full) == k)
      if (length(rows_k) > 0L) {
        for (rr in rows_k) {
          y_w <- (Xtest[rr, , drop = FALSE] - matrix(JW$mu, nrow = 1)) %*% JW$W
          scores <- .row_cor(Yhat_all, y_w[1, ])
          # stable softmax over scores
          sshift <- scores - max(scores)
          w <- exp(sshift)
          probs <- w / sum(w)
          prob_all[rr, ] <- probs
          pred_all[rr] <- unique_keys[which.max(scores)]
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
      rk <- tryCatch({ if (!is.null(rf$rank)) rf$rank else if (!is.null(rf$Ad)) length(rf$Ad) else NA_integer_ }, error = function(...) NA_integer_)
      sv <- tryCatch(as.numeric(rf$Ad), error = function(...) NULL)
      adapter_ranks <- c(adapter_ranks, rk)
      if (!is.null(sv) && length(sv)) adapter_svals <- c(adapter_svals, sv[1])
      rank_by_item[k] <- rk
      sv1_by_item[k]  <- if (!is.null(sv) && length(sv)) sv[1] else NA_real_
      if (isTRUE(mod_spec$save_fold_singulars)) {
        if (is.null(jj_map[[k]])) jj_map[[k]] <- length(jj_map) + 1L
        sv_spectra_by_item[[jj_map[[k]]]] <- sv
        names(sv_spectra_by_item)[jj_map[[k]]] <- k
      }

      # Optional per-fold diagnostics package (only if return_adapter)
      if (isTRUE(mod_spec$return_adapter)) {
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

    # Handle any test rows whose key not in unique_keys: fit on all keys
    leftover <- which(!(as.character(key_te_full) %in% unique_keys))
    if (length(leftover) > 0L) {
      JW <- .joint_whiten(Xp, Ym, do_shrink = mod_spec$shrink_whiten, center = mod_spec$center)
      rf <- .fit_residual_rrr(JW$Xw_fit, JW$Yw_fit, rank = mod_spec$rank,
                              max_rank = mod_spec$max_rank,
                              ridge_lambda = mod_spec$ridge_rrr_lambda)
      Delta <- rf$Delta
      lam_opt <- .select_lambda(JW$Xw_fit, JW$Yw_fit, Delta, mod_spec$lambda_grid)$lambda
      Xw_all <- (sweep(Xp[unique_keys, , drop = FALSE], 2, JW$mu, "-")) %*% JW$W
      Yhat_all <- Xw_all + lam_opt * (Xw_all %*% Delta)
      for (rr in leftover) {
        y_w <- (Xtest[rr, , drop = FALSE] - matrix(JW$mu, nrow = 1)) %*% JW$W
        scores <- .row_cor(Yhat_all, y_w[1, ])
        sshift <- scores - max(scores)
        w <- exp(sshift)
        probs <- w / sum(w)
        prob_all[rr, ] <- probs
        pred_all[rr] <- unique_keys[which.max(scores)]
      }
    }

    # Compute per-voxel CV R^2
    ss_vox[ss_vox < .Machine$double.eps] <- .Machine$double.eps
    r2_cv_vox <- 1 - (sse_vox / ss_vox)

    prob <- prob_all
    pred <- factor(pred_all, levels = levs)
    if (any(is.na(pred))) {
      mc <- max.col(prob)
      pred <- factor(levs[mc], levels = levs)
    }

    # Observed labels are the key factor
    obs <- factor(as.character(key_te_full), levels = levs)

    predictor_obj <- NULL
    if (isTRUE(mod_spec$return_adapter)) {
      predictor_obj <- list(r2_per_voxel = r2_cv_vox,
                            resid_by_item = resid_by_item,
                            rank_by_item  = rank_by_item,
                            sv1_by_item   = sv1_by_item,
                            keys_scored   = keys_scored,
                            skipped_keys  = skipped_keys)
      if (isTRUE(mod_spec$save_fold_singulars)) predictor_obj$sv_spectra_by_item <- sv_spectra_by_item
      if (exists("remap_diag_by_fold", inherits = FALSE)) predictor_obj$diag_by_fold <- remap_diag_by_fold
    }
    # filter out test rows whose observed key not in training/common set
    keep_rows <- which(!is.na(obs))
    cres <- classification_result(obs[keep_rows], pred[keep_rows], prob[keep_rows, , drop = FALSE],
                                  testind = keep_rows,
                                  test_design = des$test_design[keep_rows, , drop = FALSE],
                                  predictor = predictor_obj)

    base_perf <- compute_performance(mod_spec, cres)
    extra <- c(adapter_rank    = suppressWarnings(mean(adapter_ranks, na.rm = TRUE)),
               adapter_sv1     = suppressWarnings(mean(adapter_svals, na.rm = TRUE)),
               adapter_mean_r2 = mean(r2_cv_vox, na.rm = TRUE),
               remap_improv    = suppressWarnings(mean(roi_improv_vals, na.rm = TRUE)),
               delta_frob_mean = suppressWarnings(mean(delta_frob_vals, na.rm = TRUE)),
               lambda_mean     = suppressWarnings(mean(lambda_vals, na.rm = TRUE)))
    perf_vec <- c(base_perf, extra)

    return(tibble::tibble(
      result = list(cres),
      indices = list(ind),
      performance = list(c(perf_vec,
                           n_pairs_used = n_pairs,
                           n_skipped_keys = length(skipped_keys))),
      id = rnum,
      error = FALSE,
      error_message = "~",
      warning = soft_warn || length(skipped_keys) > 0,
      warning_message = if (soft_warn) sprintf("REMAP: only %d paired items (<%d).", n_pairs, mod_spec$min_pairs) else if (length(skipped_keys) > 0) sprintf("REMAP LOKO skipped %d keys due to insufficient training pairs or fit issues.", length(skipped_keys)) else "~"
    ))
  }

  # Non-LOKO single-fit path using joint whitening + residual RRR
  JW <- .joint_whiten(Xp, Ym, do_shrink = mod_spec$shrink_whiten, center = mod_spec$center)
  rf <- .fit_residual_rrr(JW$Xw_fit, JW$Yw_fit, rank = mod_spec$rank,
                          max_rank = mod_spec$max_rank,
                          ridge_lambda = mod_spec$ridge_rrr_lambda)
  Delta <- rf$Delta
  lam_opt <- .select_lambda(JW$Xw_fit, JW$Yw_fit, Delta, mod_spec$lambda_grid)$lambda

  # Residual diagnostics (global fit)
  R_naive <- JW$Yw_fit - JW$Xw_fit
  Yhat_fit <- JW$Xw_fit + lam_opt * (JW$Xw_fit %*% Delta)
  R_remap  <- JW$Yw_fit - Yhat_fit
  ss_naive <- sum(R_naive^2)
  ss_remap <- sum(R_remap^2)
  roi_improv <- if (ss_naive > .Machine$double.eps) 1 - (ss_remap / ss_naive) else NA_real_
  delta_frob <- sqrt(sum((lam_opt * Delta)^2))

  # Templates for all keys and correlation-based classification
  Xw_all <- (sweep(Xp[unique_keys, , drop = FALSE], 2, JW$mu, "-")) %*% JW$W
  Yhat_all <- Xw_all + lam_opt * (Xw_all %*% Delta)
  prob <- matrix(0, nrow(Xtest), length(unique_keys), dimnames = list(NULL, unique_keys))
  pred_chr <- character(nrow(Xtest))
  for (i in seq_len(nrow(Xtest))) {
    y_w <- (Xtest[i, , drop = FALSE] - matrix(JW$mu, nrow = 1)) %*% JW$W
    scores <- .row_cor(Yhat_all, y_w[1, ])
    sshift <- scores - max(scores)
    w <- exp(sshift); probs <- w / sum(w)
    prob[i, ] <- probs
    pred_chr[i] <- unique_keys[which.max(scores)]
  }
  pred <- factor(pred_chr, levels = unique_keys)
  obs <- factor(as.character(key_te_full), levels = unique_keys)

  # Diagnostics on prototypes in voxel space
  Winv <- try(solve(JW$W), silent = TRUE); if (inherits(Winv, "try-error")) Winv <- MASS::ginv(JW$W)
  Yhat_vox <- (Xw_all %*% (diag(ncol(Delta)) + lam_opt * Delta)) %*% Winv + matrix(JW$mu, nrow = nrow(Xw_all), ncol = ncol(Xw_all), byrow = TRUE)
  Ym_centered <- scale(Ym[unique_keys, , drop = FALSE], center = TRUE, scale = FALSE)
  denom <- colSums(Ym_centered^2); denom[denom < .Machine$double.eps] <- .Machine$double.eps
  r2 <- 1 - colSums((Ym[unique_keys, , drop = FALSE] - Yhat_vox)^2) / denom
  rnk <- tryCatch({ if (!is.null(rf$rank)) rf$rank else if (!is.null(rf$Ad)) length(rf$Ad) else NA_integer_ }, error = function(...) NA_integer_)
  svals <- tryCatch(as.numeric(rf$Ad), error = function(...) NULL)

  keep_rows <- which(!is.na(obs))
  cres <- classification_result(obs[keep_rows], pred[keep_rows], prob[keep_rows, , drop = FALSE],
                                testind = keep_rows,
                                test_design = des$test_design[keep_rows, , drop = FALSE],
                                predictor = if (isTRUE(mod_spec$return_adapter)) list(rank = rnk, singvals = svals, r2_per_voxel = r2) else NULL)
  base_perf <- compute_performance(mod_spec, cres)
  extra <- c(adapter_rank    = as.numeric(rnk %||% NA_real_),
             adapter_sv1     = if (!is.null(svals) && length(svals) > 0) svals[1] else NA_real_,
             adapter_mean_r2 = mean(r2, na.rm = TRUE),
             remap_improv    = roi_improv,
             delta_frob_mean = delta_frob,
             lambda_mean     = lam_opt)
  perf_vec <- c(base_perf, extra)

  tibble::tibble(
    result = list(cres),
    indices = list(ind),
    performance = list(c(perf_vec, n_pairs_used = n_pairs, n_skipped_keys = 0L)),
    id = rnum,
    error = FALSE,
    error_message = "~",
    warning = soft_warn,
    warning_message = if (soft_warn) sprintf("REMAP: only %d paired items (<%d).", n_pairs, mod_spec$min_pairs) else "~"
  )
}
