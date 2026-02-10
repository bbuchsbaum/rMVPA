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

#' Naive cross-decoding per ROI
#' @return A tibble row with columns \code{result}, \code{indices}, \code{performance}, and \code{id}.
#' @examples
#' \dontrun{
#'   # Internal method called by run_searchlight/run_regional
#'   # See naive_xdec_model examples for usage
#' }
#' @keywords internal
#' @export
process_roi.naive_xdec_model <- function(mod_spec, roi, rnum, ...) {
  if (!has_test_set(mod_spec)) {
    return(tibble::tibble(
      result = list(NULL), indices = list(neuroim2::indices(roi$train_roi)),
      performance = list(NULL), id = rnum,
      error = TRUE, error_message = "naive_xdec_model requires external test set"
    ))
  }

  Xtr <- as.matrix(neuroim2::values(roi$train_roi))
  Xte <- as.matrix(neuroim2::values(roi$test_roi))
  des <- mod_spec$design

  # Keys for prototypes and observed labels
  if (!is.null(mod_spec$link_by) && mod_spec$link_by %in% colnames(des$train_design) && mod_spec$link_by %in% colnames(des$test_design)) {
    key_tr <- factor(des$train_design[[mod_spec$link_by]])
    key_te <- factor(des$test_design[[mod_spec$link_by]])
  } else {
    key_tr <- factor(des$y_train)
    key_te <- factor(des$y_test)
  }

  # Align on common keys and build prototypes
  Xp <- .nx_rowsum_mean(Xtr, key_tr)   # prototypes per key
  levs <- rownames(Xp)
  # Ensure observed uses same level set
  obs <- factor(as.character(key_te), levels = levs)

  # Correlate each test row to prototypes (columns = levs)
  # Use base cor for correctness; speed is adequate at ROI scale
  scores <- cor(t(Xte), t(Xp))
  if (is.null(colnames(scores))) colnames(scores) <- levs
  probs <- .nx_softmax(scores)
  colnames(probs) <- levs
  pred <- factor(levs[max.col(scores)], levels = levs)

  cres <- classification_result(obs, pred, probs,
                                testind = seq_len(nrow(Xte)),
                                test_design = des$test_design,
                                predictor = NULL)

  perf <- compute_performance(mod_spec, cres)
  tibble::tibble(
    result = list(cres),
    indices = list(neuroim2::indices(roi$train_roi)),
    performance = list(perf),
    id = rnum,
    error = FALSE,
    error_message = "~"
  )
}
