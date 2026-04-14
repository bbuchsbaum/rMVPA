#' Validate an MVPA Analysis for Common Methodological Issues
#'
#' Performs static checks on a model specification (or design + CV pair)
#' to detect potential data leakage, class imbalance, inadequate sample
#' sizes, and other common pitfalls \emph{before} running a potentially
#' expensive analysis.
#'
#' @param x An object to validate. Typically an \code{mvpa_model} or
#'   an \code{mvpa_design} (in which case \code{crossval} must also be
#'   supplied).
#' @param crossval A cross-validation specification (only needed when
#'   \code{x} is an \code{mvpa_design}).
#' @param dataset An optional \code{mvpa_dataset}, used for feature-count
#'   checks.
#' @param verbose Logical; if \code{TRUE} (the default), print the
#'   validation report to the console.
#' @param ... Additional arguments (currently unused).
#'
#' @return A \code{validation_result} object (invisibly) containing a list
#'   of individual check results, each with \code{name}, \code{status}
#'   (\code{"pass"}, \code{"warn"}, or \code{"fail"}), and \code{message}.
#'
#' @details
#' The following checks are performed:
#'
#' \describe{
#'   \item{cv_block_alignment}{Verifies that the CV scheme respects the
#'     blocking structure (e.g. run-level blocking for fMRI).}
#'   \item{block_structure}{Checks that blocks are well-formed — enough
#'     blocks, roughly equal sizes, and all classes represented.}
#'   \item{class_balance}{Checks class balance across folds and flags
#'     severe imbalance.}
#'   \item{single_class_folds}{Detects folds where train or test sets
#'     contain only one class (guaranteed failure for classification).}
#'   \item{fold_sizes}{Warns about very small test sets or extreme
#'     train/test ratios.}
#'   \item{sample_size}{Warns when the total number of observations per
#'     class is very low.}
#'   \item{temporal_leakage}{For non-blocked CV, detects when adjacent
#'     observations (likely temporally autocorrelated) are split across
#'     train and test.}
#' }
#'
#' @examples
#' # Create a simple design
#' des_df <- data.frame(
#'   condition = factor(rep(c("A","B"), each = 50)),
#'   run = rep(1:5, each = 20)
#' )
#' des <- mvpa_design(des_df, y_train = ~ condition, block_var = ~ run)
#' cv  <- blocked_cross_validation(des$block_var)
#'
#' # Validate — should pass all checks
#' validate_analysis(des, crossval = cv)
#'
#' # Risky: kfold CV ignoring run structure
#' cv_bad <- kfold_cross_validation(len = 100, nfolds = 5)
#' validate_analysis(des, crossval = cv_bad)
#'
#' @export
validate_analysis <- function(x, ...) {
  UseMethod("validate_analysis")
}

#' @rdname validate_analysis
#' @export
validate_analysis.mvpa_model <- function(x, verbose = TRUE, ...) {
  .do_validate(
    design  = x$design,
    crossval = x$crossval,
    dataset = x$dataset,
    verbose = verbose
  )
}

#' @rdname validate_analysis
#' @export
validate_analysis.model_spec <- function(x, verbose = TRUE, ...) {
  .do_validate(
    design = x$design,
    crossval = x$crossval,
    dataset = x$dataset,
    verbose = verbose
  )
}

#' @rdname validate_analysis
#' @export
validate_analysis.mvpa_design <- function(x, crossval = NULL,
                                          dataset = NULL,
                                          verbose = TRUE, ...) {
  if (is.null(crossval)) {
    stop("'crossval' must be provided when validating an mvpa_design.",
         call. = FALSE)
  }
  .do_validate(
    design   = x,
    crossval = crossval,
    dataset  = dataset,
    verbose  = verbose
  )
}

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

#' @keywords internal
#' @noRd
.make_check <- function(name, status, message) {
  list(name = name, status = status, message = message)
}

#' @keywords internal
#' @noRd
.do_validate <- function(design, crossval, dataset = NULL, verbose = TRUE) {
  checks <- list()
  y <- y_train(design)
  block_var <- design$block_var

  # ---- Check 1: CV-block alignment ----
  checks <- c(checks, .check_cv_block_alignment(crossval, block_var))

  # ---- Check 2: Block structure ----
  if (!is.null(block_var)) {
    checks <- c(checks, .check_block_structure(block_var, y))
  }

  # ---- Check 3: Minimum observations per class ----
  checks <- c(checks, .check_min_observations(y))

  # ---- Checks 4-7: Fold-level checks (require generating folds) ----
  if (!is.null(y)) {
    n <- if (is.matrix(y)) nrow(y) else length(y)
    dummy_data <- data.frame(dummy = seq_len(n))
    folds <- tryCatch(
      crossval_samples(crossval, dummy_data, y),
      error = function(e) NULL
    )

    if (!is.null(folds)) {
      checks <- c(checks, .check_class_balance(folds, y))
      checks <- c(checks, .check_single_class_folds(folds, y))
      checks <- c(checks, .check_fold_sizes(folds, y))
      checks <- c(checks, .check_temporal_leakage(folds, block_var))
    }
  }

  n_pass <- sum(vapply(checks, function(ch) ch$status == "pass", logical(1)))
  n_warn <- sum(vapply(checks, function(ch) ch$status == "warn", logical(1)))
  n_fail <- sum(vapply(checks, function(ch) ch$status == "fail", logical(1)))

  result <- structure(
    list(checks = checks, n_pass = n_pass, n_warn = n_warn, n_fail = n_fail),
    class = "validation_result"
  )

  if (verbose) print(result)
  invisible(result)
}

#' @keywords internal
#' @noRd
.validation_result_payload <- function(x) {
  if (is.null(x) || !inherits(x, "validation_result")) {
    return(NULL)
  }

  list(
    n_pass = x$n_pass,
    n_warn = x$n_warn,
    n_fail = x$n_fail,
    checks = lapply(x$checks, function(ch) {
      list(
        name = ch$name,
        status = ch$status,
        message = ch$message
      )
    })
  )
}

#' @keywords internal
#' @noRd
.format_preflight_message <- function(x, context = "analysis preflight") {
  payload <- .validation_result_payload(x)
  if (is.null(payload)) {
    return(sprintf("%s failed for an unknown reason.", context))
  }

  flagged <- Filter(function(ch) ch$status %in% c("fail", "warn"), payload$checks)
  head_flagged <- flagged[seq_len(min(3L, length(flagged)))]

  detail <- if (length(head_flagged) == 0L) {
    "No failing or warning checks were recorded."
  } else {
    paste(
      vapply(
        head_flagged,
        function(ch) sprintf("%s: %s", ch$name, ch$message),
        character(1)
      ),
      collapse = "\n"
    )
  }

  sprintf(
    "%s reported %d failure(s) and %d warning(s).\n%s",
    context,
    payload$n_fail,
    payload$n_warn,
    detail
  )
}

#' @keywords internal
#' @noRd
.preflight_supported <- function(model_spec) {
  inherits(model_spec, "model_spec") &&
    !is.null(model_spec$design) &&
    !is.null(model_spec$crossval)
}

#' @keywords internal
#' @noRd
.apply_analysis_preflight <- function(model_spec,
                                      preflight = c("warn", "error", "off"),
                                      context = "analysis") {
  preflight <- match.arg(preflight)
  if (identical(preflight, "off") || !.preflight_supported(model_spec)) {
    return(NULL)
  }

  vres <- tryCatch(
    validate_analysis(model_spec, verbose = FALSE),
    error = function(e) e
  )

  if (inherits(vres, "error")) {
    msg <- sprintf("%s preflight could not be completed: %s", context, conditionMessage(vres))
    if (identical(preflight, "error")) {
      stop(msg, call. = FALSE)
    }
    warning(msg, call. = FALSE)
    return(NULL)
  }

  if (identical(preflight, "error") && vres$n_fail > 0L) {
    stop(.format_preflight_message(vres, context = paste(context, "preflight")), call. = FALSE)
  }

  if (identical(preflight, "warn") && (vres$n_fail > 0L || vres$n_warn > 0L)) {
    warning(.format_preflight_message(vres, context = paste(context, "preflight")), call. = FALSE)
  }

  vres
}

# ---------------------------------------------------------------------------
# Individual check functions
# ---------------------------------------------------------------------------

#' @keywords internal
#' @noRd
.check_cv_block_alignment <- function(crossval, block_var) {
  checks <- list()
  cv_class <- class(crossval)[1]

  # Design has block_var but CV ignores it

  if (!is.null(block_var) && cv_class == "kfold_cross_validation") {
    checks <- c(checks, list(.make_check(
      "cv_block_alignment", "fail",
      paste0(
        "Design has a blocking variable (e.g. fMRI run) but CV scheme is ",
        "kfold_cross_validation, which ignores block structure. ",
        "Observations from the same run may appear in both train and test, ",
        "causing temporal autocorrelation leakage. ",
        "Use blocked_cross_validation(block_var) instead."
      )
    )))
  } else if (is.null(block_var) && !cv_class %in% c("custom_cross_validation",
                                                      "kfold_cross_validation")) {
    checks <- c(checks, list(.make_check(
      "no_block_var", "warn",
      paste0(
        "No blocking variable in the design. For fMRI data, blocking by ",
        "run is strongly recommended to prevent temporal leakage."
      )
    )))
  } else {
    checks <- c(checks, list(.make_check(
      "cv_block_alignment", "pass",
      "CV scheme is consistent with the blocking structure."
    )))
  }

  # Check for mismatched block_vars between design and CV spec
  if (cv_class == "blocked_cross_validation" && !is.null(block_var)) {
    cv_bv <- crossval$block_var
    if (!is.null(cv_bv) && length(cv_bv) == length(block_var)) {
      if (!identical(as.integer(cv_bv), as.integer(block_var))) {
        checks <- c(checks, list(.make_check(
          "cv_block_mismatch", "warn",
          paste0(
            "The block variable in the CV spec differs from the design's ",
            "block_var. Verify this is intentional."
          )
        )))
      }
    }
  }

  checks
}


#' @keywords internal
#' @noRd
.check_block_structure <- function(block_var, y) {
  checks <- list()
  block_tab <- table(block_var)
  n_blocks <- length(block_tab)

  if (n_blocks < 2) {
    checks <- c(checks, list(.make_check(
      "block_count", "fail",
      "Only 1 block found. Cross-validation requires at least 2 blocks."
    )))
  } else if (n_blocks == 2) {
    checks <- c(checks, list(.make_check(
      "block_count", "warn",
      paste0(
        "Only 2 blocks found. Leave-one-block-out CV will have only 2 folds, ",
        "providing limited evaluation with high variance estimates."
      )
    )))
  } else {
    checks <- c(checks, list(.make_check(
      "block_count", "pass",
      sprintf("%d blocks found.", n_blocks)
    )))
  }

  # Block size imbalance
  block_sizes <- as.integer(block_tab)
  if (n_blocks >= 2) {
    size_ratio <- max(block_sizes) / max(min(block_sizes), 1L)
    if (size_ratio > 3) {
      checks <- c(checks, list(.make_check(
        "block_size_imbalance", "warn",
        sprintf(
          "Block sizes are highly unequal (range: %d to %d, ratio: %.1fx). ",
          min(block_sizes), max(block_sizes), size_ratio
        )
      )))
    }
  }

  # Class representation within blocks
  if (is.factor(y) && n_blocks >= 2) {
    block_class_tab <- table(block_var, y)
    blocks_missing_class <- apply(block_class_tab, 1, function(row) any(row == 0))
    if (any(blocks_missing_class)) {
      n_bad <- sum(blocks_missing_class)
      checks <- c(checks, list(.make_check(
        "block_class_coverage", "warn",
        sprintf(
          "%d of %d block(s) are missing at least one class. ",
          n_bad, n_blocks
        )
      )))
    }
  }

  checks
}


#' @keywords internal
#' @noRd
.check_min_observations <- function(y) {
  checks <- list()

  if (is.factor(y)) {
    class_counts <- table(y)
    min_count <- min(class_counts)
    n_classes <- length(class_counts)

    if (min_count < 2) {
      checks <- c(checks, list(.make_check(
        "min_observations", "fail",
        sprintf(
          "Class '%s' has only %d observation(s). Need at least 2 per class for CV.",
          names(which.min(class_counts)), min_count
        )
      )))
    } else if (min_count < 5) {
      checks <- c(checks, list(.make_check(
        "min_observations", "warn",
        sprintf(
          "Smallest class '%s' has only %d observations. Results may be unreliable.",
          names(which.min(class_counts)), min_count
        )
      )))
    } else {
      checks <- c(checks, list(.make_check(
        "min_observations", "pass",
        sprintf(
          "%d classes with %d to %d observations each.",
          n_classes, min_count, max(class_counts)
        )
      )))
    }
  } else {
    n <- if (is.matrix(y)) nrow(y) else length(y)
    if (n < 10) {
      checks <- c(checks, list(.make_check(
        "min_observations", "warn",
        sprintf("Only %d observations for regression. Results may be unreliable.", n)
      )))
    } else {
      checks <- c(checks, list(.make_check(
        "min_observations", "pass",
        sprintf("%d observations.", n)
      )))
    }
  }

  checks
}


#' @keywords internal
#' @noRd
.check_class_balance <- function(folds, y) {
  checks <- list()
  if (!is.factor(y)) return(checks)

  n_folds <- nrow(folds)
  worst_ratio <- 1
  worst_fold <- 1L

  for (i in seq_len(n_folds)) {
    train_idx <- folds$train[[i]]$idx
    train_y <- y[train_idx]
    tab <- table(train_y)
    if (length(tab) > 1) {
      ratio <- max(tab) / max(min(tab), 1L)
      if (ratio > worst_ratio) {
        worst_ratio <- ratio
        worst_fold <- i
      }
    }
  }

  if (worst_ratio > 5) {
    checks <- c(checks, list(.make_check(
      "class_balance", "warn",
      sprintf(
        "Severe class imbalance in training fold %d (%.1f:1 ratio). Consider using balance_partitions().",
        worst_fold, worst_ratio
      )
    )))
  } else {
    checks <- c(checks, list(.make_check(
      "class_balance", "pass",
      sprintf("Class balance OK across folds (worst ratio: %.1f:1).", worst_ratio)
    )))
  }

  checks
}


#' @keywords internal
#' @noRd
.check_single_class_folds <- function(folds, y) {
  checks <- list()
  if (!is.factor(y)) return(checks)

  n_folds <- nrow(folds)
  bad_train <- integer(0)
  bad_test  <- integer(0)

  for (i in seq_len(n_folds)) {
    train_idx <- folds$train[[i]]$idx
    test_idx  <- folds$test[[i]]$idx
    if (length(unique(y[train_idx])) < 2) bad_train <- c(bad_train, i)
    if (length(unique(y[test_idx]))  < 2) bad_test  <- c(bad_test, i)
  }

  if (length(bad_train) > 0) {
    checks <- c(checks, list(.make_check(
      "single_class_train", "fail",
      sprintf(
        "%d fold(s) have only one class in the training set (folds: %s). Classification will fail.",
        length(bad_train), paste(bad_train, collapse = ", ")
      )
    )))
  }

  if (length(bad_test) > 0) {
    checks <- c(checks, list(.make_check(
      "single_class_test", "warn",
      sprintf(
        "%d fold(s) have only one class in the test set (folds: %s). Metrics will be unreliable.",
        length(bad_test), paste(bad_test, collapse = ", ")
      )
    )))
  }

  if (length(bad_train) == 0 && length(bad_test) == 0) {
    checks <- c(checks, list(.make_check(
      "single_class_folds", "pass",
      "All folds have multiple classes in both train and test sets."
    )))
  }

  checks
}


#' @keywords internal
#' @noRd
.check_fold_sizes <- function(folds, y) {
  checks <- list()
  n_folds <- nrow(folds)

  test_sizes <- vapply(seq_len(n_folds), function(i) {
    length(folds$test[[i]]$idx)
  }, integer(1))

  train_sizes <- vapply(seq_len(n_folds), function(i) {
    length(folds$train[[i]]$idx)
  }, integer(1))

  min_test  <- min(test_sizes)
  min_train <- min(train_sizes)

  if (min_test < 2) {
    checks <- c(checks, list(.make_check(
      "fold_test_size", "fail",
      sprintf(
        "At least one fold has only %d test observation(s). Metrics cannot be reliably computed.",
        min_test
      )
    )))
  } else if (min_test < 5) {
    checks <- c(checks, list(.make_check(
      "fold_test_size", "warn",
      sprintf(
        "Smallest test fold has only %d observations. Consider fewer folds or more data.",
        min_test
      )
    )))
  } else {
    checks <- c(checks, list(.make_check(
      "fold_sizes", "pass",
      sprintf(
        "Fold sizes OK (train: %d-%d, test: %d-%d).",
        min_train, max(train_sizes), min_test, max(test_sizes)
      )
    )))
  }

  checks
}


#' @keywords internal
#' @noRd
.check_temporal_leakage <- function(folds, block_var) {
  checks <- list()
  # Only relevant when there is no block_var (i.e., CV doesn't respect blocks)
  if (!is.null(block_var)) return(checks)

  n_folds <- nrow(folds)
  n_adjacent_leaks <- 0L

  # Check a few folds to avoid expensive full scan
  for (i in seq_len(min(n_folds, 3L))) {
    train_idx <- sort(folds$train[[i]]$idx)
    test_idx  <- sort(folds$test[[i]]$idx)

    all_idx  <- sort(c(train_idx, test_idx))
    in_train <- all_idx %in% train_idx

    # Count transitions where consecutive indices switch between train/test
    if (length(all_idx) > 1) {
      consecutive <- diff(all_idx) == 1L
      switches    <- diff(as.integer(in_train)) != 0L
      n_adjacent_leaks <- n_adjacent_leaks + sum(consecutive & switches)
    }
  }

  if (n_adjacent_leaks > 0) {
    checks <- c(checks, list(.make_check(
      "temporal_leakage", "warn",
      sprintf(
        paste0(
          "%d adjacent observation pairs split across train/test ",
          "(checked %d fold(s)). If these are fMRI time-series, ",
          "temporal autocorrelation may inflate accuracy."
        ),
        n_adjacent_leaks, min(n_folds, 3L)
      )
    )))
  }

  checks
}


# ---------------------------------------------------------------------------
# Print method
# ---------------------------------------------------------------------------

#' @export
#' @method print validation_result
print.validation_result <- function(x, ...) {
  has_crayon <- requireNamespace("crayon", quietly = TRUE)
  header_fn  <- if (has_crayon) crayon::bold$cyan else identity
  pass_fn    <- if (has_crayon) crayon::green     else identity
  warn_fn    <- if (has_crayon) crayon::yellow    else identity
  fail_fn    <- if (has_crayon) crayon::red       else identity

  status_icon <- function(s) {
    switch(s,
      pass = pass_fn("[PASS]"),
      warn = warn_fn("[WARN]"),
      fail = fail_fn("[FAIL]"),
      paste0("[", toupper(s), "]")
    )
  }

  cat("\n", header_fn("MVPA Analysis Validation"), "\n\n")

  for (ch in x$checks) {
    cat("  ", status_icon(ch$status), " ", ch$name, "\n")
    msg_lines <- strwrap(ch$message, width = 70, prefix = "    ")
    cat(paste(msg_lines, collapse = "\n"), "\n\n")
  }

  parts <- character(0)
  if (x$n_pass > 0) parts <- c(parts, pass_fn(sprintf("%d passed", x$n_pass)))
  if (x$n_warn > 0) parts <- c(parts, warn_fn(sprintf("%d warnings", x$n_warn)))
  if (x$n_fail > 0) parts <- c(parts, fail_fn(sprintf("%d failures", x$n_fail)))

  cat("  ", header_fn("Summary:"), paste(parts, collapse = ", "), "\n\n")
  invisible(x)
}
