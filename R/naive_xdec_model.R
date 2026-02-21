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
#' @param ... Additional arguments (currently unused).
#'
#' @return model spec object of class `naive_xdec_model` for use with
#'   `run_regional()` or `run_searchlight()`.
#' @examples
#' \dontrun{
#'   # Requires dataset with train_data and test_data
#'   ds <- gen_sample_dataset(c(5,5,5), 20, external_test=TRUE)
#'   model <- naive_xdec_model(ds$dataset, ds$design)
#' }
#' @export
naive_xdec_model <- function(dataset, design, link_by = NULL, return_predictions = TRUE, ...) {
  # Choose performance function based on training response
  perf_fun <- if (is.numeric(design$y_train)) {
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

  if (anyNA(Xa) || anyNA(Xb) || any(!is.finite(Xa)) || any(!is.finite(Xb))) {
    return(cor(t(Xa), t(Xb)))
  }

  ma <- rowMeans(Xa)
  mb <- rowMeans(Xb)
  Ac <- sweep(Xa, 1, ma, "-")
  Bc <- sweep(Xb, 1, mb, "-")

  sa <- sqrt(rowSums(Ac * Ac) / (p - 1))
  sb <- sqrt(rowSums(Bc * Bc) / (p - 1))
  cov <- tcrossprod(Ac, Bc) / (p - 1)
  denom <- tcrossprod(sa, sb)

  out <- cov / denom
  out[!is.finite(out)] <- NA_real_
  out
}

#' @keywords internal
.naive_xdec_fit_core <- function(Xtr, Xte, kernel = NULL, key_tr = NULL, key_te = NULL) {
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
  probs <- .nx_softmax(scores)
  colnames(probs) <- levs
  pred <- factor(levs[max.col(scores)], levels = levs)

  list(obs = obs, pred = pred, probs = probs)
}

#' @keywords internal
.nx_softmax <- function(scores) {
  # numerically-stable row-wise softmax
  m <- apply(scores, 1, max)
  z <- exp(sweep(scores, 1, m, "-"))
  sweep(z, 1, rowSums(z), "/")
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

  core <- .naive_xdec_fit_core(
    Xtr = Xtr,
    Xte = Xte,
    kernel = model$.fast_kernel,
    key_tr = key_tr,
    key_te = key_te
  )

  test_ind <- seq_len(nrow(Xte))
  fast_perf <- .naive_xdec_fast_metrics(model, core, test_idx = test_ind)

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
