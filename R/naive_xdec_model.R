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

  create_model_spec(
    "naive_xdec_model",
    dataset = dataset,
    design  = design,
    link_by = link_by,
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

  # Keys for prototypes and observed labels
  if (!is.null(model$link_by) && model$link_by %in% colnames(des$train_design) && model$link_by %in% colnames(des$test_design)) {
    key_tr <- factor(des$train_design[[model$link_by]])
    key_te <- factor(des$test_design[[model$link_by]])
  } else {
    key_tr <- factor(des$y_train)
    key_te <- factor(des$y_test)
  }

  # Align on common keys and build prototypes
  Xp <- .nx_rowsum_mean(Xtr, key_tr)
  levs <- rownames(Xp)
  obs <- factor(as.character(key_te), levels = levs)

  # Correlate each test row to prototypes
  scores <- cor(t(Xte), t(Xp))
  if (is.null(colnames(scores))) colnames(scores) <- levs
  probs <- .nx_softmax(scores)
  colnames(probs) <- levs
  pred <- factor(levs[max.col(scores)], levels = levs)

  cres <- classification_result(obs, pred, probs,
                                testind = seq_len(nrow(Xte)),
                                test_design = des$test_design,
                                predictor = NULL)

  perf <- compute_performance(model, cres)

  roi_result(
    metrics = perf,
    indices = roi_data$indices,
    id = context$id,
    result = cres
  )
}

#' @rdname output_schema
#' @method output_schema naive_xdec_model
#' @export
output_schema.naive_xdec_model <- function(model) {
  NULL
}

