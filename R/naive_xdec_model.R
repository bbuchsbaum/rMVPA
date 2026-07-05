#' Naive Cross-Decoding (correlation to source prototypes)
#'
#' Baseline analysis that classifies target-domain trials by correlating them
#' directly with source-domain prototypes (means over repeats) without any
#' adaptation. Serves as a comparator for REMAP-RRR.
#'
#' The term "naive" is a technical designation meaning "without domain adaptation"
#' (i.e., direct prototype matching). It serves as the principled baseline for
#' evaluating domain-adaptive methods like REMAP-RRR. Also known as "zero-shot
#' transfer" or "direct transfer" in machine learning literature.
#'
#' @param dataset mvpa_dataset with train_data (source) and test_data (target)
#' @param design  mvpa_design with y_train/y_test and optional `link_by` column
#'   present in both train/test designs. If `link_by` is NULL, prototypes are
#'   formed by the training labels `y_train` and evaluated against `y_test`.
#' @param link_by Optional character; if provided, prototypes and observed labels
#'   are keyed by this column instead of y.
#' @param return_predictions logical; keep per-ROI predictions.
#' @param performance Optional user-supplied performance function. When
#'   non-`NULL`, must be a function that accepts the classification result
#'   object (with fields `observed`, `predicted`, `probs`, `testind`,
#'   `test_design`) and returns a named numeric vector of metrics. Routed
#'   through [get_custom_perf()], which annotates it as `"custom"`; this
#'   disables the optimised fast-metric kernel and ensures the user's
#'   function is always called on the full result object. When `NULL`, the
#'   default binary / multiclass / regression performance helper is used.
#' @param ... Additional arguments (currently unused).
#'
#' @return model spec object of class `naive_xdec_model` for use with
#'   `run_regional()` or `run_searchlight()`.
#' @examples
#' \dontrun{
#'   # Requires dataset with train_data and test_data
#'   ds <- gen_sample_dataset(c(5,5,5), 20, external_test=TRUE)
#'   model <- naive_xdec_model(ds$dataset, ds$design)
#'
#'   # Custom metric that uses a column from test_design
#'   custom_fun <- function(result) {
#'     vivid <- result$test_design$RateVivid
#'     probs <- as.matrix(result$probs)
#'     obs   <- as.character(result$observed)
#'     true_p <- probs[cbind(seq_along(obs), match(obs, colnames(probs)))]
#'     c(vivid_spearman = stats::cor(vivid, true_p, method = "spearman"))
#'   }
#'   model <- naive_xdec_model(ds$dataset, ds$design, performance = custom_fun)
#' }
#' @export
naive_xdec_model <- function(dataset, design, link_by = NULL,
                              return_predictions = TRUE,
                              performance = NULL, ...) {
  # Choose performance function: user-supplied custom takes precedence; else
  # pick the default helper based on the training response type.
  perf_fun <- if (!is.null(performance)) {
    if (!is.function(performance)) {
      stop("`performance` must be a function or NULL", call. = FALSE)
    }
    get_custom_perf(performance, design$split_groups)
  } else if (is.numeric(design$y_train)) {
    get_regression_perf(design$split_groups)
  } else if (length(levels(design$y_train)) > 2) {
    get_multiclass_perf(design$split_groups, class_metrics = FALSE)
  } else {
    get_binary_perf(design$split_groups)
  }

  fast_kernel <- NULL
  if (.naive_xdec_fast_kernel_enabled()) {
    fast_kernel <- .naive_xdec_prepare_fast_kernel(design = design, link_by = link_by)
  }

  create_model_spec(
    "naive_xdec_model",
    dataset = dataset,
    design  = design,
    link_by = link_by,
    .fast_kernel = fast_kernel,
    performance = perf_fun,
    compute_performance = TRUE,
    return_predictions = return_predictions,
    ...
  )
}

#' @export
print.naive_xdec_model <- function(x, ...) {
  cat("Naive Cross-Decoding Model\n")
  cat("  link_by:            ", x$link_by %||% "(class labels)", "\n")
  cat("  training levels:    ", paste(levels(x$design$y_train), collapse = ", "), "\n")
  cat("  n_train:            ", length(x$design$y_train), "\n")
  if (!is.null(x$design$y_test)) {
    cat("  n_test:             ", length(x$design$y_test), "\n")
  }
  cat("  return_predictions: ", x$return_predictions, "\n")
  invisible(x)
}

#' @keywords internal
.nx_rowsum_mean <- function(X, f) {
  f <- factor(f)
  rowsum(X, f, reorder = TRUE) / as.vector(table(f))
}

#' @keywords internal
.naive_xdec_fast_kernel_enabled <- function() {
  TRUE
}

#' @keywords internal
.naive_xdec_fast_metric_mode <- function(model, classes) {
  perf_fun <- model$performance
  perf_kind <- attr(perf_fun, "rmvpa_perf_kind", exact = TRUE)

  if (is.null(perf_kind)) {
    return(list(enabled = FALSE))
  }

  # User-supplied custom performance functions must always run on the full
  # classification result; never short-circuit through the fast-metric kernel.
  if (identical(perf_kind, "custom")) {
    return(list(enabled = FALSE))
  }

  if (identical(perf_kind, "binary") && length(classes) == 2L) {
    return(list(
      enabled = TRUE,
      kind = "binary",
      class_metrics = FALSE,
      split_groups = model$design$split_groups
    ))
  }

  if (identical(perf_kind, "multiclass") && length(classes) > 2L) {
    return(list(
      enabled = TRUE,
      kind = "multiclass",
      class_metrics = isTRUE(attr(perf_fun, "rmvpa_class_metrics", exact = TRUE)),
      split_groups = model$design$split_groups
    ))
  }

  list(enabled = FALSE)
}

#' @keywords internal
.naive_xdec_fast_metrics <- function(model, core, test_idx) {
  if (anyNA(core$obs) || anyNA(core$probs) || any(!is.finite(core$probs))) {
    return(NULL)
  }

  classes <- levels(core$obs)
  metric_mode <- .naive_xdec_fast_metric_mode(model, classes = classes)
  if (!isTRUE(metric_mode$enabled)) {
    return(NULL)
  }

  if (!exists(".dual_lda_metric_with_splits", mode = "function", inherits = TRUE)) {
    return(NULL)
  }
  metric_fun <- get(".dual_lda_metric_with_splits", mode = "function", inherits = TRUE)

  metric_fun(
    observed = core$obs,
    probs = core$probs,
    test_idx = test_idx,
    classes = classes,
    kind = metric_mode$kind,
    split_groups = metric_mode$split_groups,
    class_metrics = metric_mode$class_metrics
  )
}

#' @keywords internal
.nx_resolve_keys <- function(design, link_by = NULL) {
  if (!is.null(link_by) &&
      link_by %in% colnames(design$train_design) &&
      link_by %in% colnames(design$test_design)) {
    key_tr <- factor(design$train_design[[link_by]])
    key_te <- factor(design$test_design[[link_by]])
  } else {
    key_tr <- factor(design$y_train)
    key_te <- factor(design$y_test)
  }
  list(key_tr = key_tr, key_te = key_te)
}

#' @keywords internal
.naive_xdec_prepare_fast_kernel <- function(design, link_by = NULL) {
  keys <- .nx_resolve_keys(design = design, link_by = link_by)
  key_tr <- keys$key_tr
  key_te <- keys$key_te
  levs <- levels(key_tr)

  if (length(levs) < 1L) {
    return(NULL)
  }

  group_mat <- model.matrix(~ key_tr - 1)
  if (!is.matrix(group_mat) || nrow(group_mat) != length(key_tr) || ncol(group_mat) != length(levs)) {
    return(NULL)
  }
  storage.mode(group_mat) <- "double"
  colnames(group_mat) <- levs

  counts <- as.numeric(table(key_tr))
  if (length(counts) != length(levs) || any(counts <= 0)) {
    return(NULL)
  }

  list(
    levels = levs,
    key_tr = key_tr,
    obs = factor(as.character(key_te), levels = levs),
    group_mat = group_mat,
    counts = counts,
    n_train = length(key_tr),
    n_test = length(key_te)
  )
}

#' @keywords internal
.nx_row_cor <- function(Xa, Xb) {
  Xa <- as.matrix(Xa)
  Xb <- as.matrix(Xb)

  if (ncol(Xa) != ncol(Xb)) {
    stop("Row correlation requires matching feature dimensions.")
  }

  p <- ncol(Xa)
  out <- matrix(NA_real_, nrow = nrow(Xa), ncol = nrow(Xb))
  if (p < 2L || nrow(Xa) < 1L || nrow(Xb) < 1L) {
    return(out)
  }

  out <- suppressWarnings(cor(t(Xa), t(Xb)))
  out[!is.finite(out)] <- NA_real_
  out
}

#' @keywords internal
.naive_xdec_fit_core <- function(Xtr, Xte, kernel = NULL, key_tr = NULL, key_te = NULL,
                                 return_probs = TRUE, return_scores = FALSE) {
  if (!is.null(kernel)) {
    if (nrow(Xtr) == kernel$n_train && nrow(Xte) == kernel$n_test) {
      Xp <- crossprod(kernel$group_mat, Xtr)
      Xp <- sweep(Xp, 1, kernel$counts, "/")
      storage.mode(Xp) <- "double"
      levs <- kernel$levels
      rownames(Xp) <- levs
      obs <- kernel$obs
      scores <- .nx_row_cor(Xte, Xp)
    } else {
      # Defensive fallback if design/test rows differ from cached metadata.
      key_tr <- kernel$key_tr
      key_te <- kernel$obs
      Xp <- .nx_rowsum_mean(Xtr, key_tr)
      levs <- rownames(Xp)
      obs <- factor(as.character(key_te), levels = levs)
      scores <- cor(t(Xte), t(Xp))
    }
  } else {
    Xp <- .nx_rowsum_mean(Xtr, key_tr)
    levs <- rownames(Xp)
    obs <- factor(as.character(key_te), levels = levs)
    scores <- cor(t(Xte), t(Xp))
  }

  if (is.null(colnames(scores))) {
    colnames(scores) <- levs
  }
  pred <- factor(levs[max.col(scores)], levels = levs)

  probs <- NULL
  if (isTRUE(return_probs)) {
    probs <- .nx_softmax(scores)
    colnames(probs) <- levs
  }

  ret <- list(obs = obs, pred = pred, probs = probs)
  if (isTRUE(return_scores)) {
    ret$scores <- scores
  }
  ret
}

#' @keywords internal
.nx_softmax <- function(scores) {
  # numerically-stable row-wise softmax
  m <- matrixStats::rowMaxs(scores)
  z <- exp(sweep(scores, 1, m, "-"))
  sweep(z, 1, rowSums(z), "/")
}

#' @keywords internal
.nx_row_logsumexp <- function(scores) {
  m <- matrixStats::rowMaxs(scores)
  m + log(rowSums(exp(sweep(scores, 1, m, "-"))))
}

#' @keywords internal
.naive_xdec_metric_from_scores <- function(observed, pred, scores, classes, kind,
                                           class_metrics = FALSE) {
  observed <- factor(observed, levels = classes)
  pred <- factor(pred, levels = classes)
  acc <- mean(pred == observed)

  if (identical(kind, "binary")) {
    score <- if (ncol(scores) >= 2L) {
      scores[, 2L] - scores[, 1L]
    } else {
      scores[, 1L]
    }
    auc <- .dual_lda_auc_from_scores(score, observed == classes[2L])
    auc_centered <- if (is.na(auc)) NA_real_ else 2 * auc - 1
    return(c(Accuracy = acc, AUC = auc_centered))
  }

  k <- length(classes)
  log_denom <- .nx_row_logsumexp(scores)
  aucres <- rep(NA_real_, k)
  for (i in seq_len(k)) {
    score <- scores[, i] - log_denom
    aucres[i] <- .dual_lda_auc_from_scores(score, observed == classes[i])
  }

  auc_centered <- 2 * aucres - 1
  mean_auc <- mean(auc_centered, na.rm = TRUE)
  if (is.nan(mean_auc)) {
    mean_auc <- NA_real_
  }

  metrics <- c(Accuracy = acc, AUC = mean_auc)
  if (isTRUE(class_metrics)) {
    names(auc_centered) <- paste0("AUC_", classes)
    c(metrics, auc_centered)
  } else {
    metrics
  }
}

#' @keywords internal
.naive_xdec_metric_from_scores_with_splits <- function(observed, pred, scores,
                                                       test_idx, classes, kind,
                                                       split_groups = NULL,
                                                       class_metrics = FALSE) {
  base_vals <- .naive_xdec_metric_from_scores(
    observed = observed,
    pred = pred,
    scores = scores,
    classes = classes,
    kind = kind,
    class_metrics = class_metrics
  )

  if (is.null(split_groups) || length(split_groups) == 0L) {
    return(base_vals)
  }

  tags <- names(split_groups)
  if (is.null(tags) || length(tags) == 0L) {
    return(base_vals)
  }

  split_vecs <- lapply(tags, function(tag) {
    design_idx <- split_groups[[tag]]
    res_idx <- which(test_idx %in% design_idx)

    if (!length(res_idx)) {
      vals <- rep(NA_real_, length(base_vals))
      names(vals) <- paste0(names(base_vals), "_", tag)
      return(vals)
    }

    vals <- .naive_xdec_metric_from_scores(
      observed = observed[res_idx],
      pred = pred[res_idx],
      scores = scores[res_idx, , drop = FALSE],
      classes = classes,
      kind = kind,
      class_metrics = class_metrics
    )
    names(vals) <- paste0(names(vals), "_", tag)
    vals
  })

  c(base_vals, unlist(split_vecs))
}

#' @keywords internal
.naive_xdec_fast_metrics_from_scores <- function(model, core, test_idx) {
  if (is.null(core$scores) ||
      anyNA(core$obs) ||
      anyNA(core$scores) ||
      any(!is.finite(core$scores))) {
    return(NULL)
  }

  classes <- levels(core$obs)
  metric_mode <- .naive_xdec_fast_metric_mode(model, classes = classes)
  if (!isTRUE(metric_mode$enabled)) {
    return(NULL)
  }

  if (!exists(".dual_lda_auc_from_scores", mode = "function", inherits = TRUE)) {
    return(NULL)
  }

  .naive_xdec_metric_from_scores_with_splits(
    observed = core$obs,
    pred = core$pred,
    scores = core$scores,
    test_idx = test_idx,
    classes = classes,
    kind = metric_mode$kind,
    split_groups = metric_mode$split_groups,
    class_metrics = metric_mode$class_metrics
  )
}

#' @keywords internal
.is_naive_xdec_fast_path <- function(model_spec, method) {
  if (!inherits(model_spec, "naive_xdec_model") ||
      !identical(method, "standard") ||
      !isTRUE(has_test_set(model_spec)) ||
      isTRUE(model_spec$return_predictions) ||
      is.null(model_spec$.fast_kernel) ||
      inherits(model_spec$dataset, "mvpa_multibasis_image_dataset") ||
      !inherits(model_spec$dataset, "mvpa_image_dataset") ||
      is.null(model_spec$dataset$test_data)) {
    return(FALSE)
  }

  classes <- model_spec$.fast_kernel$levels
  metric_mode <- .naive_xdec_fast_metric_mode(model_spec, classes = classes)
  isTRUE(metric_mode$enabled)
}

#' @keywords internal
.naive_xdec_expand_searchlight_ids <- function(roi) {
  ids <- tryCatch({
    coords <- tryCatch(roi@coords, error = function(...) NULL)
    sp <- tryCatch(roi@space, error = function(...) NULL)
    if (!is.null(coords) && !is.null(sp)) {
      neuroim2::grid_to_index(sp, coords)
    } else {
      roi
    }
  }, error = function(...) roi)

  ids <- suppressWarnings(as.integer(ids))
  ids[is.finite(ids) & ids > 0L]
}

#' @keywords internal
.naive_xdec_filter_feature_mask <- function(Xtr, feature_ids, center_id = NULL,
                                            min_voxels = 1L) {
  stats <- .filter_roi_stats(Xtr)
  keep <- !stats$nas & stats$sdnonzero

  if (!is.null(center_id) && length(center_id) == 1L && !is.na(center_id)) {
    kp <- match(center_id, feature_ids)
    if (!is.na(kp) && !keep[kp]) {
      keep[kp] <- TRUE
    }
  }

  if (sum(keep) < min_voxels) {
    return(NULL)
  }
  keep
}

#' @keywords internal
.naive_xdec_image_feature_matrices <- function(dataset) {
  mask_idx <- which(dataset$mask > 0)
  if (!length(mask_idx)) {
    stop("naive_xdec fast searchlight requires at least one active mask voxel.", call. = FALSE)
  }

  list(
    train = as.matrix(neuroim2::series(dataset$train_data, mask_idx)),
    test = as.matrix(neuroim2::series(dataset$test_data, mask_idx)),
    mask_idx = as.integer(mask_idx)
  )
}

#' @keywords internal
.naive_xdec_bad_searchlight_row <- function(id, message) {
  tibble::tibble(
    result = list(NULL),
    indices = list(NULL),
    performance = list(NULL),
    id = id,
    error = TRUE,
    error_message = message,
    warning = TRUE,
    warning_message = message
  )
}

#' @keywords internal
.naive_xdec_standard_searchlight_fast <- function(model_spec, radius, k = NULL,
                                                  verbose = FALSE) {
  dataset <- model_spec$dataset
  slight <- get_searchlight(dataset, "standard", radius, k = k)
  center_ids <- get_center_ids(dataset)
  if (length(center_ids) != length(slight)) {
    stop("naive_xdec fast searchlight: center id count does not match searchlight count.",
         call. = FALSE)
  }

  mats <- .naive_xdec_image_feature_matrices(dataset)
  col_map <- integer(max(mats$mask_idx))
  col_map[mats$mask_idx] <- seq_along(mats$mask_idx)

  perf_rows <- vector("list", length(slight))
  good_ids <- integer(length(slight))
  bad_rows <- vector("list", length(slight))
  n_good <- 0L
  n_bad <- 0L
  test_idx <- seq_len(nrow(mats$test))

  for (i in seq_along(slight)) {
    center_id <- as.integer(center_ids[[i]])
    roi_ids <- .naive_xdec_expand_searchlight_ids(slight[[i]])
    if (!length(roi_ids)) {
      n_bad <- n_bad + 1L
      bad_rows[[n_bad]] <- .naive_xdec_bad_searchlight_row(center_id, "empty ROI")
      next
    }

    valid <- roi_ids <= length(col_map) & col_map[roi_ids] > 0L
    roi_ids <- roi_ids[valid]
    cols <- col_map[roi_ids]
    if (!length(cols)) {
      n_bad <- n_bad + 1L
      bad_rows[[n_bad]] <- .naive_xdec_bad_searchlight_row(center_id, "ROI outside mask")
      next
    }

    Xtr <- mats$train[, cols, drop = FALSE]
    keep <- .naive_xdec_filter_feature_mask(
      Xtr = Xtr,
      feature_ids = roi_ids,
      center_id = center_id,
      min_voxels = 1L
    )
    if (is.null(keep)) {
      n_bad <- n_bad + 1L
      bad_rows[[n_bad]] <- .naive_xdec_bad_searchlight_row(center_id, "ROI filtered out")
      next
    }

    Xtr <- Xtr[, keep, drop = FALSE]
    Xte <- mats$test[, cols[keep], drop = FALSE]

    perf_error <- NULL
    perf <- tryCatch({
      core <- .naive_xdec_fit_core(
        Xtr = Xtr,
        Xte = Xte,
        kernel = model_spec$.fast_kernel,
        return_probs = FALSE,
        return_scores = TRUE
      )
      out <- .naive_xdec_fast_metrics_from_scores(model_spec, core, test_idx = test_idx)
      if (is.null(out)) {
        core <- .naive_xdec_fit_core(Xtr = Xtr, Xte = Xte, kernel = model_spec$.fast_kernel)
        out <- .naive_xdec_fast_metrics(model_spec, core, test_idx = test_idx)
      }
      out
    }, error = function(e) {
      perf_error <<- conditionMessage(e)
      NULL
    })

    if (is.null(perf)) {
      n_bad <- n_bad + 1L
      bad_rows[[n_bad]] <- .naive_xdec_bad_searchlight_row(
        center_id,
        perf_error %||% "metric computation failed"
      )
      next
    }

    n_good <- n_good + 1L
    perf_rows[[n_good]] <- perf
    good_ids[[n_good]] <- center_id
  }

  if (n_good == 0L) {
    stop("No valid results for naive_xdec fast searchlight: all ROIs failed to process.",
         call. = FALSE)
  }

  perf_mat <- do.call(rbind, perf_rows[seq_len(n_good)])
  good_ids <- good_ids[seq_len(n_good)]
  out <- wrap_out(perf_mat, dataset, good_ids)

  bad_results <- if (n_bad > 0L) {
    dplyr::bind_rows(bad_rows[seq_len(n_bad)])
  } else {
    tibble::tibble(
      result = list(),
      indices = list(),
      performance = list(),
      id = integer(),
      error = logical(),
      error_message = character(),
      warning = logical(),
      warning_message = character()
    )
  }
  attr(out, "bad_results") <- bad_results
  attr(out, "searchlight_engine") <- "naive_xdec_fast"
  out
}

#' @keywords internal
#' @noRd
.run_searchlight_engine.naive_xdec_model <- function(model_spec, radius, method,
                                                     engine = "auto",
                                                     niter = 4L,
                                                     combiner = "average",
                                                     drop_probs = FALSE,
                                                     fail_fast = FALSE,
                                                     backend = c("default", "shard", "auto"),
                                                     incremental = TRUE,
                                                     gamma = NULL,
                                                     verbose = FALSE,
                                                     k = NULL,
                                                     ...) {
  requested <- .match_searchlight_engine(engine)
  backend <- match.arg(backend)

  if (identical(requested, "legacy")) {
    return(list(handled = FALSE, result = NULL, engine = "legacy"))
  }

  eligible <- .is_naive_xdec_fast_path(model_spec, method) &&
    requested %in% c("auto", "naive_xdec_fast") &&
    backend %in% c("default", "auto") &&
    !isTRUE(drop_probs) &&
    !isTRUE(fail_fast) &&
    is.character(combiner) &&
    combiner[1] %in% c("average", "standard")

  if (!isTRUE(eligible)) {
    if (identical(requested, "naive_xdec_fast")) {
      stop("Requested searchlight engine 'naive_xdec_fast' is not eligible for this analysis.",
           call. = FALSE)
    }
    return(list(handled = FALSE, result = NULL, engine = "legacy"))
  }

  list(
    handled = TRUE,
    result = .naive_xdec_standard_searchlight_fast(
      model_spec = model_spec,
      radius = radius,
      k = k,
      verbose = verbose
    ),
    engine = "naive_xdec_fast"
  )
}

#' @export
compute_performance.naive_xdec_model <- function(obj, result) {
  obj$performance(result)
}

#' @rdname fit_roi
#' @method fit_roi naive_xdec_model
#' @export
fit_roi.naive_xdec_model <- function(model, roi_data, context, ...) {
  if (!has_test_set(model)) {
    return(roi_result(
      metrics = NULL,
      indices = roi_data$indices,
      id = context$id,
      error = TRUE,
      error_message = "naive_xdec_model requires external test set"
    ))
  }

  Xtr <- roi_data$train_data
  Xte <- roi_data$test_data
  des <- model$design

  if (is.null(Xte) || nrow(Xtr) < 2 || ncol(Xtr) < 1 || nrow(Xte) < 1) {
    return(roi_result(
      metrics = NULL,
      indices = roi_data$indices,
      id = context$id,
      error = TRUE,
      error_message = "naive_xdec_model: insufficient data"
    ))
  }

  key_tr <- NULL
  key_te <- NULL
  if (is.null(model$.fast_kernel)) {
    keys <- .nx_resolve_keys(design = des, link_by = model$link_by)
    key_tr <- keys$key_tr
    key_te <- keys$key_te
  }

  classes <- if (!is.null(model$.fast_kernel)) {
    model$.fast_kernel$levels
  } else {
    levels(key_tr)
  }
  metric_mode <- .naive_xdec_fast_metric_mode(model, classes = classes)
  score_fast_metrics <- !isTRUE(model$return_predictions) && isTRUE(metric_mode$enabled)

  core <- .naive_xdec_fit_core(
    Xtr = Xtr,
    Xte = Xte,
    kernel = model$.fast_kernel,
    key_tr = key_tr,
    key_te = key_te,
    return_probs = !score_fast_metrics,
    return_scores = score_fast_metrics
  )

  test_ind <- seq_len(nrow(Xte))
  fast_perf <- if (score_fast_metrics) {
    .naive_xdec_fast_metrics_from_scores(model, core, test_idx = test_ind)
  } else {
    .naive_xdec_fast_metrics(model, core, test_idx = test_ind)
  }

  if (score_fast_metrics && is.null(fast_perf)) {
    core <- .naive_xdec_fit_core(
      Xtr = Xtr,
      Xte = Xte,
      kernel = model$.fast_kernel,
      key_tr = key_tr,
      key_te = key_te
    )
    fast_perf <- .naive_xdec_fast_metrics(model, core, test_idx = test_ind)
  }

  need_result_object <- isTRUE(model$return_predictions) || is.null(fast_perf)
  cres <- NULL
  if (need_result_object) {
    cres <- classification_result(core$obs, core$pred, core$probs,
                                  testind = test_ind,
                                  test_design = des$test_design,
                                  predictor = NULL)
  }

  perf <- if (is.null(fast_perf)) {
    compute_performance(model, cres)
  } else {
    fast_perf
  }

  roi_result(
    metrics = perf,
    indices = roi_data$indices,
    id = context$id,
    result = if (isTRUE(model$return_predictions)) cres else NULL
  )
}

#' @rdname output_schema
#' @method output_schema naive_xdec_model
#' @export
output_schema.naive_xdec_model <- function(model) {
  NULL
}
