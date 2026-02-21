#' @keywords internal
#' @noRd
.swift_searchlight_enabled <- function() {
  TRUE
}

#' @keywords internal
#' @noRd
.swift_whitening_mode <- function(mode = NULL) {
  allowed <- c("none", "zscore", "diag_shrink")

  if (!is.null(mode)) {
    return(match.arg(as.character(mode)[1], allowed))
  }
  "zscore"
}

#' @keywords internal
#' @noRd
.is_swift_fast_path <- function(model_spec, method) {
  if (!inherits(model_spec, "mvpa_model")) {
    return(FALSE)
  }
  if (!(method %in% c("standard", "randomized", "resampled"))) {
    return(FALSE)
  }
  if (!isTRUE(has_crossval(model_spec))) {
    return(FALSE)
  }
  if (isTRUE(has_test_set(model_spec))) {
    return(FALSE)
  }

  ds <- model_spec$dataset
  if (!inherits(ds, "mvpa_image_dataset") || inherits(ds, "mvpa_multibasis_image_dataset")) {
    return(FALSE)
  }

  y <- y_train(model_spec)
  if (!is.factor(y) || nlevels(y) < 3L) {
    return(FALSE)
  }

  perf_fun <- model_spec$performance
  perf_kind <- attr(perf_fun, "rmvpa_perf_kind", exact = TRUE)
  identical(perf_kind, "multiclass")
}

#' @keywords internal
#' @noRd
.swift_extract_fold_idx <- function(split_obj) {
  if (is.null(split_obj)) {
    return(integer(0))
  }

  idx <- NULL
  if (is.list(split_obj) && !is.null(split_obj$idx)) {
    idx <- split_obj$idx
  } else if (inherits(split_obj, "resample") && !is.null(split_obj$idx)) {
    idx <- split_obj$idx
  } else if (is.integer(split_obj) || is.numeric(split_obj)) {
    idx <- split_obj
  }

  as.integer(idx %||% integer(0))
}

#' @keywords internal
#' @noRd
.swift_whiten_train_test <- function(x_train, x_test, mode = c("none", "zscore", "diag_shrink")) {
  mode <- match.arg(mode)

  if (identical(mode, "none")) {
    return(list(train = x_train, test = x_test))
  }

  mu <- colMeans(x_train)
  x_train_centered <- sweep(x_train, 2, mu, "-")
  x_test_centered <- sweep(x_test, 2, mu, "-")

  sdev <- matrixStats::colSds(x_train_centered)
  sdev[!is.finite(sdev)] <- 0

  if (identical(mode, "diag_shrink")) {
    med_sd <- stats::median(sdev[sdev > 0], na.rm = TRUE)
    if (!is.finite(med_sd) || med_sd <= 0) {
      med_sd <- 1
    }
    alpha <- 0.1
    sdev <- sqrt((1 - alpha) * (sdev^2) + alpha * (med_sd^2))
  }

  sdev[sdev <= .Machine$double.eps] <- 1

  list(
    train = sweep(x_train_centered, 2, sdev, "/"),
    test = sweep(x_test_centered, 2, sdev, "/")
  )
}

#' @keywords internal
#' @noRd
.swift_class_means <- function(x, y, classes) {
  y <- factor(y, levels = classes)
  design <- stats::model.matrix(~ y - 1)
  counts <- colSums(design)
  means <- crossprod(design, x)
  means <- sweep(means, 1, counts, "/")
  rownames(means) <- classes
  means
}

#' @keywords internal
#' @noRd
.swift_multiclass_evidence <- function(x_train, x_test, y_train, y_test, classes) {
  stats <- .swift_multiclass_fold_stats(
    x_train = x_train,
    x_test = x_test,
    y_train = y_train,
    y_test = y_test,
    classes = classes
  )
  if (is.null(stats)) {
    return(NULL)
  }
  stats$evidence
}

#' @keywords internal
#' @noRd
.swift_multiclass_fold_stats <- function(x_train, x_test, y_train, y_test, classes) {
  y_train <- factor(y_train, levels = classes)
  y_test <- factor(y_test, levels = classes)

  n_train <- tabulate(as.integer(y_train), nbins = length(classes))
  n_test <- tabulate(as.integer(y_test), nbins = length(classes))
  if (any(is.na(n_train)) || any(is.na(n_test)) || any(n_train <= 0) || any(n_test <= 0)) {
    return(NULL)
  }

  mu_train <- .swift_class_means(x_train, y_train, classes)
  mu_test <- .swift_class_means(x_test, y_test, classes)

  mu_train_global <- colMeans(x_train)
  mu_test_global <- colMeans(x_test)

  d_train <- sweep(mu_train, 2, mu_train_global, "-")
  d_test <- sweep(mu_test, 2, mu_test_global, "-")

  d_train <- sweep(d_train, 1, sqrt(n_train), "*")
  d_test <- sweep(d_test, 1, sqrt(n_test), "*")

  list(
    mu_train = mu_train,
    mu_test = mu_test,
    n_train = n_train,
    n_test = n_test,
    evidence = colSums(d_train * d_test)
  )
}

#' @keywords internal
#' @noRd
.swift_precompute_center_cols <- function(centers, center_coords, offsets,
                                          dims, mask_active, col_lookup) {
  lapply(seq_along(centers), function(i) {
    ids <- .offset_to_indices(
      center_coord = center_coords[i, ],
      offsets = offsets,
      dims = dims,
      mask_active = mask_active
    )
    cols <- col_lookup[ids]
    unique(cols[cols > 0L])
  })
}

#' @keywords internal
#' @noRd
.swift_score_centers <- function(evidence, center_cols) {
  vapply(center_cols, function(cols) {
    if (length(cols) == 0L) {
      NA_real_
    } else {
      sum(evidence[cols])
    }
  }, numeric(1))
}

#' @keywords internal
#' @noRd
.swift_center_sparse_matrix <- function(center_cols, n_voxels) {
  nnz <- sum(lengths(center_cols))
  if (nnz == 0L) {
    return(Matrix::sparseMatrix(i = integer(0), j = integer(0), x = numeric(0), dims = c(length(center_cols), n_voxels)))
  }

  i <- rep.int(seq_along(center_cols), lengths(center_cols))
  j <- unlist(center_cols, use.names = FALSE)
  Matrix::sparseMatrix(
    i = as.integer(i),
    j = as.integer(j),
    x = rep.int(1, length(j)),
    dims = c(length(center_cols), n_voxels)
  )
}

#' @keywords internal
#' @noRd
.swift_metric_mode <- function(model_spec, classes) {
  perf_fun <- model_spec$performance
  perf_kind <- attr(perf_fun, "rmvpa_perf_kind", exact = TRUE)
  if (!identical(perf_kind, "multiclass") || length(classes) < 3L) {
    return(list(enabled = FALSE))
  }

  list(
    enabled = TRUE,
    split_groups = model_spec$design$split_groups,
    class_metrics = isTRUE(attr(perf_fun, "rmvpa_class_metrics", exact = TRUE))
  )
}

#' @keywords internal
#' @noRd
.swift_metric_block_size <- function() {
  256L
}

#' @keywords internal
#' @noRd
.swift_auc_matrix_from_scores <- function(score_mat, positive) {
  positive <- as.logical(positive)
  n_cols <- ncol(score_mat)
  auc <- rep(NA_real_, n_cols)
  n_pos <- sum(positive)
  n_neg <- length(positive) - n_pos
  if (n_pos == 0L || n_neg == 0L || n_cols == 0L) {
    return(auc)
  }

  finite_cols <- colSums(!is.finite(score_mat)) == 0L
  if (!any(finite_cols)) {
    return(auc)
  }

  ranks <- matrixStats::colRanks(
    score_mat[, finite_cols, drop = FALSE],
    ties.method = "average",
    preserveShape = TRUE
  )
  sum_pos <- colSums(ranks[positive, , drop = FALSE])
  auc_vals <- (sum_pos - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
  auc[finite_cols] <- auc_vals
  auc
}

#' @keywords internal
#' @noRd
.swift_metric_core_matrix <- function(observed, probs, classes, class_metrics = FALSE) {
  n_obs <- length(observed)
  n_centers <- dim(probs)[2]
  k <- length(classes)

  obs_idx <- as.integer(factor(observed, levels = classes))
  prob_flat <- matrix(probs, nrow = n_obs * n_centers, ncol = k)

  pred_flat <- max.col(prob_flat)
  pred_idx <- matrix(pred_flat, nrow = n_obs, ncol = n_centers)
  acc <- colMeans(pred_idx == obs_idx)

  auc_centered <- matrix(NA_real_, nrow = n_centers, ncol = k)
  for (kk in seq_len(k)) {
    auc_k <- .swift_auc_matrix_from_scores(probs[, , kk, drop = TRUE], obs_idx == kk)
    auc_centered[, kk] <- 2 * auc_k - 1
  }
  mean_auc <- rowMeans(auc_centered, na.rm = TRUE)
  mean_auc[is.nan(mean_auc)] <- NA_real_

  out <- cbind(Accuracy = acc, AUC = mean_auc)
  if (isTRUE(class_metrics)) {
    colnames(auc_centered) <- paste0("AUC_", classes)
    out <- cbind(out, auc_centered)
  }
  out
}

#' @keywords internal
#' @noRd
.swift_metric_with_splits_matrix <- function(observed, probs, test_idx, classes,
                                             split_groups = NULL, class_metrics = FALSE) {
  base <- .swift_metric_core_matrix(
    observed = observed,
    probs = probs,
    classes = classes,
    class_metrics = class_metrics
  )
  base_names <- colnames(base)

  if (is.null(split_groups) || length(split_groups) == 0L) {
    return(base)
  }

  tags <- names(split_groups)
  if (is.null(tags) || length(tags) == 0L) {
    return(base)
  }

  split_blocks <- lapply(tags, function(tag) {
    res_idx <- which(test_idx %in% split_groups[[tag]])
    if (length(res_idx) == 0L) {
      vals <- matrix(NA_real_, nrow = nrow(base), ncol = ncol(base))
      colnames(vals) <- paste0(base_names, "_", tag)
      return(vals)
    }

    vals <- .swift_metric_core_matrix(
      observed = observed[res_idx],
      probs = probs[res_idx, , , drop = FALSE],
      classes = classes,
      class_metrics = class_metrics
    )
    colnames(vals) <- paste0(base_names, "_", tag)
    vals
  })

  do.call(cbind, c(list(base), split_blocks))
}

#' @keywords internal
#' @noRd
.swift_compute_metric_matrix <- function(fold_cache,
                                         center_t,
                                         center_lookup_full,
                                         observed_union,
                                         test_idx_union,
                                         classes,
                                         split_groups = NULL,
                                         class_metrics = FALSE,
                                         block_size = 256L,
                                         verbose = FALSE) {
  n_centers <- ncol(center_t)
  if (n_centers == 0L) {
    return(matrix(numeric(0), nrow = 0L, ncol = 0L))
  }

  n_union <- length(observed_union)
  k <- length(classes)
  block_starts <- seq.int(1L, n_centers, by = max(1L, as.integer(block_size)))
  out_blocks <- vector("list", length(block_starts))
  metric_names <- NULL
  fold_row_idx <- lapply(fold_cache, function(fold) match(fold$test_idx, test_idx_union))

  for (b in seq_along(block_starts)) {
    b_start <- block_starts[b]
    b_end <- min(n_centers, b_start + block_size - 1L)
    b_cols <- b_start:b_end
    b_size <- length(b_cols)

    center_t_block <- center_t[, b_cols, drop = FALSE]
    if (!inherits(center_t_block, "dgCMatrix")) {
      center_t_block <- methods::as(center_t_block, "dgCMatrix")
    }
    center_t_rows <- center_t_block@i + 1L
    center_t_x <- center_t_block@x
    center_full <- center_lookup_full[b_cols]
    prob_sum <- array(0, dim = c(n_union, b_size, k))
    prob_n <- matrix(0L, nrow = n_union, ncol = b_size)

    for (f in seq_along(fold_cache)) {
      fold <- fold_cache[[f]]
      row_idx <- fold_row_idx[[f]]
      keep_rows <- which(!is.na(row_idx))
      if (length(keep_rows) == 0L) {
        next
      }

      row_idx <- row_idx[keep_rows]
      if (length(keep_rows) == nrow(fold$x_test)) {
        x_test <- fold$x_test
      } else {
        x_test <- fold$x_test[keep_rows, , drop = FALSE]
      }
      n_test <- nrow(x_test)
      if (n_test == 0L) {
        next
      }

      score_flat <- matrix(0, nrow = n_test * b_size, ncol = k)
      mu_sq_center_block <- fold$mu_sq_center[center_full, , drop = FALSE]
      scaled_block <- center_t_block
      for (kk in seq_len(k)) {
        scaled_block@x <- center_t_x * fold$mu_train[kk, center_t_rows]
        dot <- as.matrix(x_test %*% scaled_block)
        bias <- fold$log_priors[kk] - 0.5 * mu_sq_center_block[, kk]
        dot <- sweep(dot, 2, bias, "+")
        score_flat[, kk] <- dot
      }

      max_scores <- matrixStats::rowMaxs(score_flat)
      probs_flat <- exp(score_flat - max_scores)
      probs_flat <- probs_flat / rowSums(probs_flat)

      for (kk in seq_len(k)) {
        prob_k <- matrix(probs_flat[, kk], nrow = n_test, ncol = b_size)
        prob_sum[row_idx, , kk] <- prob_sum[row_idx, , kk] + prob_k
      }
      prob_n[row_idx, ] <- prob_n[row_idx, ] + 1L
    }

    probs <- array(1 / k, dim = c(n_union, b_size, k))
    nz <- prob_n > 0L
    for (kk in seq_len(k)) {
      val <- prob_sum[, , kk]
      val[nz] <- val[nz] / prob_n[nz]
      probs[, , kk] <- val
    }

    block_perf <- .swift_metric_with_splits_matrix(
      observed = observed_union,
      probs = probs,
      test_idx = test_idx_union,
      classes = classes,
      split_groups = split_groups,
      class_metrics = class_metrics
    )

    if (is.null(metric_names)) {
      metric_names <- colnames(block_perf)
    } else if (!identical(metric_names, colnames(block_perf))) {
      stop("SWIFT fast path: inconsistent performance metric names across center blocks.")
    }
    out_blocks[[b]] <- block_perf

    if (isTRUE(verbose) && (b %% 4L == 0L || b == length(block_starts))) {
      futile.logger::flog.info(
        "SWIFT fast path: scored metric blocks %d/%d",
        b,
        length(block_starts)
      )
    }
  }

  out <- do.call(rbind, out_blocks)
  colnames(out) <- metric_names
  out
}

#' @keywords internal
#' @noRd
run_searchlight_swift_fast <- function(model_spec, radius, whitening = NULL, verbose = FALSE, ...) {
  ds <- model_spec$dataset
  if (!inherits(ds, "mvpa_image_dataset")) {
    stop("SWIFT fast path currently supports mvpa_image_dataset only.")
  }
  if (inherits(ds, "mvpa_multibasis_image_dataset")) {
    stop("SWIFT fast path does not currently support multibasis datasets.")
  }

  y_all <- y_train(model_spec)
  if (!is.factor(y_all)) {
    stop("SWIFT fast path requires factor responses.")
  }
  classes <- levels(y_all)
  if (length(classes) < 3L) {
    stop("SWIFT fast path currently targets multiclass (>= 3 classes).")
  }

  whitening_mode <- .swift_whitening_mode(whitening)

  space_obj <- resolve_volume_space(ds)
  dims <- spatial_dim_shape(space_obj)
  spacing <- neuroim2::spacing(space_obj)[1:3]

  mask_indices <- ds$mask_indices
  if (is.null(mask_indices)) {
    mask_indices <- compute_mask_indices(ds$mask)
  }
  if (length(mask_indices) == 0L) {
    return(empty_searchlight_result(ds))
  }

  centers <- get_center_ids(ds)
  centers <- intersect(centers, mask_indices)
  if (length(centers) == 0L) {
    return(empty_searchlight_result(ds))
  }

  offsets <- .searchlight_offsets(radius, spacing = spacing)
  if (nrow(offsets) == 0L) {
    return(empty_searchlight_result(ds))
  }

  n_space <- prod(dims)
  mask_active <- logical(n_space)
  mask_active[mask_indices] <- TRUE

  col_lookup <- integer(n_space)
  col_lookup[mask_indices] <- seq_along(mask_indices)

  center_coords <- neuroim2::index_to_grid(space_obj, centers)
  center_cols <- .swift_precompute_center_cols(
    centers = centers,
    center_coords = center_coords,
    offsets = offsets,
    dims = dims,
    mask_active = mask_active,
    col_lookup = col_lookup
  )
  center_sparse <- .swift_center_sparse_matrix(center_cols, n_voxels = length(mask_indices))
  center_t <- methods::as(Matrix::t(center_sparse), "dgCMatrix")

  x_all <- as.matrix(neuroim2::series(ds$train_data, mask_indices))
  if (nrow(x_all) != length(y_all)) {
    stop("SWIFT fast path: mismatch between training rows and y_train length.")
  }

  folds <- generate_folds(
    model_spec$crossval,
    tibble::tibble(.row = seq_len(nrow(x_all))),
    y_all
  )
  if (nrow(folds) == 0L) {
    stop("SWIFT fast path: cross-validation produced zero folds.")
  }

  score_sum <- numeric(length(centers))
  score_n <- integer(length(centers))
  metric_mode <- .swift_metric_mode(model_spec, classes)
  fold_cache <- vector("list", nrow(folds))
  n_valid_folds <- 0L
  bad_rows <- list()

  for (i in seq_len(nrow(folds))) {
    tr_idx <- .swift_extract_fold_idx(folds$train[[i]])
    te_idx <- .swift_extract_fold_idx(folds$test[[i]])

    if (length(tr_idx) == 0L || length(te_idx) == 0L) {
      bad_rows[[length(bad_rows) + 1L]] <- tibble::tibble(
        id = NA_integer_,
        error = TRUE,
        error_message = sprintf("SWIFT fold %d has empty train/test indices.", i)
      )
      next
    }

    y_train_fold <- factor(y_all[tr_idx], levels = classes)
    y_test_fold <- factor(y_all[te_idx], levels = classes)

    if (any(table(y_train_fold) == 0L) || any(table(y_test_fold) == 0L)) {
      bad_rows[[length(bad_rows) + 1L]] <- tibble::tibble(
        id = NA_integer_,
        error = TRUE,
        error_message = sprintf(
          "SWIFT fold %d missing at least one class in train/test; skipping fold.",
          i
        )
      )
      next
    }

    x_train <- x_all[tr_idx, , drop = FALSE]
    x_test <- x_all[te_idx, , drop = FALSE]
    xw <- .swift_whiten_train_test(
      x_train = x_train,
      x_test = x_test,
      mode = whitening_mode
    )

    fold_stats <- .swift_multiclass_fold_stats(
      x_train = xw$train,
      x_test = xw$test,
      y_train = y_train_fold,
      y_test = y_test_fold,
      classes = classes
    )
    evidence <- if (is.null(fold_stats)) NULL else fold_stats$evidence

    if (is.null(evidence) || any(!is.finite(evidence))) {
      bad_rows[[length(bad_rows) + 1L]] <- tibble::tibble(
        id = NA_integer_,
        error = TRUE,
        error_message = sprintf("SWIFT fold %d produced invalid evidence.", i)
      )
      next
    }

    fold_scores <- as.numeric(center_sparse %*% evidence)
    keep <- is.finite(fold_scores)
    if (!any(keep)) {
      bad_rows[[length(bad_rows) + 1L]] <- tibble::tibble(
        id = NA_integer_,
        error = TRUE,
        error_message = sprintf("SWIFT fold %d produced no finite center scores.", i)
      )
      next
    }

    score_sum[keep] <- score_sum[keep] + fold_scores[keep]
    score_n[keep] <- score_n[keep] + 1L

    n_valid_folds <- n_valid_folds + 1L
    if (isTRUE(metric_mode$enabled)) {
      priors <- fold_stats$n_train / sum(fold_stats$n_train)
      priors[!is.finite(priors) | priors <= 0] <- .Machine$double.eps
      priors <- priors / sum(priors)
      mu_sq_center <- as.matrix(center_sparse %*% Matrix::t(fold_stats$mu_train^2))

      fold_cache[[n_valid_folds]] <- list(
        test_idx = te_idx,
        x_test = xw$test,
        mu_train = fold_stats$mu_train,
        log_priors = log(priors),
        mu_sq_center = mu_sq_center
      )
    }

    if (isTRUE(verbose) && (i %% 2L == 0L || i == nrow(folds))) {
      futile.logger::flog.info(
        "SWIFT fast path: processed %d/%d folds",
        i,
        nrow(folds)
      )
    }
  }

  if (n_valid_folds == 0L) {
    out <- empty_searchlight_result(ds)
    attr(out, "bad_results") <- if (length(bad_rows) > 0L) dplyr::bind_rows(bad_rows) else tibble::tibble()
    return(out)
  }

  avg_scores <- rep(NA_real_, length(centers))
  keep <- score_n > 0L
  avg_scores[keep] <- score_sum[keep] / score_n[keep]

  good_idx <- which(is.finite(avg_scores))
  if (length(good_idx) == 0L) {
    out <- empty_searchlight_result(ds)
    attr(out, "bad_results") <- if (length(bad_rows) > 0L) dplyr::bind_rows(bad_rows) else tibble::tibble()
    return(out)
  }

  perf_mat <- NULL
  if (isTRUE(metric_mode$enabled)) {
    fold_cache <- fold_cache[seq_len(n_valid_folds)]
    test_idx_union <- sort(unique(unlist(lapply(fold_cache, function(x) x$test_idx))))
    observed_union <- factor(y_all[test_idx_union], levels = classes)

    metric_mat <- .swift_compute_metric_matrix(
      fold_cache = fold_cache,
      center_t = center_t[, good_idx, drop = FALSE],
      center_lookup_full = good_idx,
      observed_union = observed_union,
      test_idx_union = test_idx_union,
      classes = classes,
      split_groups = metric_mode$split_groups,
      class_metrics = metric_mode$class_metrics,
      block_size = .swift_metric_block_size(),
      verbose = verbose
    )
    perf_mat <- cbind(metric_mat, SWIFT_Info = avg_scores[good_idx])
  } else {
    perf_mat <- matrix(avg_scores[good_idx], ncol = 1)
    colnames(perf_mat) <- "SWIFT_Info"
  }

  out <- wrap_out(perf_mat, ds, ids = centers[good_idx])
  attr(out, "bad_results") <- if (length(bad_rows) > 0L) dplyr::bind_rows(bad_rows) else tibble::tibble()
  attr(out, "searchlight_engine") <- "swift"
  attr(out, "swift_whitening") <- whitening_mode
  out
}

#' @keywords internal
#' @noRd
.swift_sampled_combiner <- function(combiner) {
  .engine_sampled_combiner(combiner, "SWIFT")
}

#' @keywords internal
#' @noRd
.swift_extract_roi_indices <- .engine_extract_roi_indices

#' @keywords internal
#' @noRd
.swift_extract_roi_ids <- .engine_extract_roi_ids

#' @keywords internal
#' @noRd
.swift_randomized_iteration <- function(slight,
                                        fold_base,
                                        metric_mode,
                                        observed_union,
                                        test_idx_union,
                                        classes,
                                        mask_indices,
                                        col_lookup,
                                        verbose = FALSE) {
  if (length(slight) == 0L) {
    return(NULL)
  }

  roi_indices <- .swift_extract_roi_indices(slight)
  roi_ids <- .swift_extract_roi_ids(slight)
  roi_cols <- lapply(roi_indices, function(ids) {
    cols <- col_lookup[ids]
    unique(cols[cols > 0L])
  })

  roi_sparse <- .swift_center_sparse_matrix(roi_cols, n_voxels = length(mask_indices))
  n_roi <- nrow(roi_sparse)
  if (n_roi == 0L) {
    return(NULL)
  }

  score_sum <- numeric(n_roi)
  score_n <- integer(n_roi)
  fold_cache <- vector("list", length(fold_base))

  for (f in seq_along(fold_base)) {
    fold <- fold_base[[f]]
    fold_scores <- as.numeric(roi_sparse %*% fold$evidence)
    keep <- is.finite(fold_scores)
    if (any(keep)) {
      score_sum[keep] <- score_sum[keep] + fold_scores[keep]
      score_n[keep] <- score_n[keep] + 1L
    }

    if (isTRUE(metric_mode$enabled)) {
      fold_cache[[f]] <- list(
        test_idx = fold$test_idx,
        x_test = fold$x_test,
        mu_train = fold$mu_train,
        log_priors = fold$log_priors,
        mu_sq_center = as.matrix(roi_sparse %*% Matrix::t(fold$mu_train^2))
      )
    }
  }

  avg_scores <- rep(NA_real_, n_roi)
  keep <- score_n > 0L
  avg_scores[keep] <- score_sum[keep] / score_n[keep]
  good_idx <- which(is.finite(avg_scores))
  if (length(good_idx) == 0L) {
    return(NULL)
  }

  perf_mat <- NULL
  if (isTRUE(metric_mode$enabled)) {
    center_t <- methods::as(Matrix::t(roi_sparse[good_idx, , drop = FALSE]), "dgCMatrix")
    metric_mat <- .swift_compute_metric_matrix(
      fold_cache = fold_cache,
      center_t = center_t,
      center_lookup_full = good_idx,
      observed_union = observed_union,
      test_idx_union = test_idx_union,
      classes = classes,
      split_groups = metric_mode$split_groups,
      class_metrics = metric_mode$class_metrics,
      block_size = .swift_metric_block_size(),
      verbose = verbose
    )
    perf_mat <- cbind(metric_mat, SWIFT_Info = avg_scores[good_idx])
  } else {
    perf_mat <- matrix(avg_scores[good_idx], ncol = 1)
    colnames(perf_mat) <- "SWIFT_Info"
  }

  perf_list <- lapply(seq_len(nrow(perf_mat)), function(i) {
    vals <- as.numeric(perf_mat[i, ])
    names(vals) <- colnames(perf_mat)
    vals
  })

  tibble::tibble(
    result = rep(list(NULL), nrow(perf_mat)),
    indices = roi_indices[good_idx],
    performance = perf_list,
    id = roi_ids[good_idx],
    error = FALSE,
    error_message = "~"
  )
}

#' @keywords internal
#' @noRd
run_searchlight_swift_sampled_fast <- function(model_spec,
                                               radius,
                                               method = c("randomized", "resampled"),
                                               niter = 4L,
                                               combiner = "average",
                                               whitening = NULL,
                                               verbose = FALSE,
                                               ...) {
  method <- match.arg(method)
  combiner_fun <- .swift_sampled_combiner(combiner)
  niter <- as.integer(niter)[1]
  if (!is.finite(niter) || is.na(niter) || niter < 1L) {
    stop("SWIFT sampled fast path requires niter >= 1.", call. = FALSE)
  }

  ds <- model_spec$dataset
  if (!inherits(ds, "mvpa_image_dataset")) {
    stop("SWIFT fast path currently supports mvpa_image_dataset only.")
  }
  if (inherits(ds, "mvpa_multibasis_image_dataset")) {
    stop("SWIFT fast path does not currently support multibasis datasets.")
  }

  y_all <- y_train(model_spec)
  if (!is.factor(y_all)) {
    stop("SWIFT fast path requires factor responses.")
  }
  classes <- levels(y_all)
  if (length(classes) < 3L) {
    stop("SWIFT fast path currently targets multiclass (>= 3 classes).")
  }

  whitening_mode <- .swift_whitening_mode(whitening)

  mask_indices <- ds$mask_indices
  if (is.null(mask_indices)) {
    mask_indices <- compute_mask_indices(ds$mask)
  }
  if (length(mask_indices) == 0L) {
    return(empty_searchlight_result(ds))
  }

  space_obj <- resolve_volume_space(ds)
  dims <- spatial_dim_shape(space_obj)
  n_space <- prod(dims)
  col_lookup <- integer(n_space)
  col_lookup[mask_indices] <- seq_along(mask_indices)

  x_all <- as.matrix(neuroim2::series(ds$train_data, mask_indices))
  if (nrow(x_all) != length(y_all)) {
    stop("SWIFT sampled fast path: mismatch between training rows and y_train length.")
  }

  folds <- generate_folds(
    model_spec$crossval,
    tibble::tibble(.row = seq_len(nrow(x_all))),
    y_all
  )
  if (nrow(folds) == 0L) {
    stop("SWIFT sampled fast path: cross-validation produced zero folds.")
  }

  fold_base <- vector("list", nrow(folds))
  n_valid_folds <- 0L
  bad_rows <- list()

  for (i in seq_len(nrow(folds))) {
    tr_idx <- .swift_extract_fold_idx(folds$train[[i]])
    te_idx <- .swift_extract_fold_idx(folds$test[[i]])

    if (length(tr_idx) == 0L || length(te_idx) == 0L) {
      bad_rows[[length(bad_rows) + 1L]] <- tibble::tibble(
        id = NA_integer_,
        error = TRUE,
        error_message = sprintf("SWIFT sampled fold %d has empty train/test indices.", i)
      )
      next
    }

    y_train_fold <- factor(y_all[tr_idx], levels = classes)
    y_test_fold <- factor(y_all[te_idx], levels = classes)
    if (any(table(y_train_fold) == 0L) || any(table(y_test_fold) == 0L)) {
      bad_rows[[length(bad_rows) + 1L]] <- tibble::tibble(
        id = NA_integer_,
        error = TRUE,
        error_message = sprintf(
          "SWIFT sampled fold %d missing at least one class in train/test; skipping fold.",
          i
        )
      )
      next
    }

    x_train <- x_all[tr_idx, , drop = FALSE]
    x_test <- x_all[te_idx, , drop = FALSE]
    xw <- .swift_whiten_train_test(
      x_train = x_train,
      x_test = x_test,
      mode = whitening_mode
    )

    fold_stats <- .swift_multiclass_fold_stats(
      x_train = xw$train,
      x_test = xw$test,
      y_train = y_train_fold,
      y_test = y_test_fold,
      classes = classes
    )
    if (is.null(fold_stats) || any(!is.finite(fold_stats$evidence))) {
      bad_rows[[length(bad_rows) + 1L]] <- tibble::tibble(
        id = NA_integer_,
        error = TRUE,
        error_message = sprintf("SWIFT sampled fold %d produced invalid evidence.", i)
      )
      next
    }

    priors <- fold_stats$n_train / sum(fold_stats$n_train)
    priors[!is.finite(priors) | priors <= 0] <- .Machine$double.eps
    priors <- priors / sum(priors)

    n_valid_folds <- n_valid_folds + 1L
    fold_base[[n_valid_folds]] <- list(
      test_idx = te_idx,
      x_test = xw$test,
      mu_train = fold_stats$mu_train,
      evidence = fold_stats$evidence,
      log_priors = log(priors)
    )
  }

  if (n_valid_folds == 0L) {
    out <- empty_searchlight_result(ds)
    attr(out, "bad_results") <- if (length(bad_rows) > 0L) dplyr::bind_rows(bad_rows) else tibble::tibble()
    return(out)
  }
  fold_base <- fold_base[seq_len(n_valid_folds)]

  metric_mode <- .swift_metric_mode(model_spec, classes)
  test_idx_union <- NULL
  observed_union <- NULL
  if (isTRUE(metric_mode$enabled)) {
    test_idx_union <- sort(unique(unlist(lapply(fold_base, `[[`, "test_idx"))))
    observed_union <- factor(y_all[test_idx_union], levels = classes)
  }

  collect_good <- list()
  iter_counter <- 0L
  iter_seq <- if (identical(method, "randomized")) seq_len(niter) else 1L

  for (iter in iter_seq) {
    slight <- if (identical(method, "randomized")) {
      get_searchlight(ds, type = "randomized", radius = radius)
    } else {
      get_searchlight(ds, type = "resampled", radius = radius, iter = niter)
    }

    iter_tbl <- .swift_randomized_iteration(
      slight = slight,
      fold_base = fold_base,
      metric_mode = metric_mode,
      observed_union = observed_union,
      test_idx_union = test_idx_union,
      classes = classes,
      mask_indices = mask_indices,
      col_lookup = col_lookup,
      verbose = verbose
    )
    if (!is.null(iter_tbl) && nrow(iter_tbl) > 0L) {
      iter_counter <- iter_counter + 1L
      collect_good[[iter_counter]] <- iter_tbl
    }

    if (identical(method, "resampled")) {
      break
    }
  }

  if (iter_counter == 0L) {
    out <- empty_searchlight_result(ds)
    attr(out, "bad_results") <- if (length(bad_rows) > 0L) dplyr::bind_rows(bad_rows) else tibble::tibble()
    return(out)
  }

  good_results <- dplyr::bind_rows(collect_good)
  bad_results <- if (length(bad_rows) > 0L) dplyr::bind_rows(bad_rows) else tibble::tibble()

  out <- combiner_fun(model_spec, good_results, bad_results)
  attr(out, "bad_results") <- bad_results
  attr(out, "searchlight_engine") <- "swift"
  attr(out, "swift_whitening") <- whitening_mode
  out
}
