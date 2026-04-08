#' Feature-sets design (mvpa_design extension for continuous regression)
#'
#' `feature_sets_design()` defines a small `mvpa_design` extension that attaches
#' grouped continuous predictors for stimulus->brain regression.
#'
#' This mirrors the pattern in `hrfdecoder_design()`: rMVPA's core infrastructure
#' expects `mvpa_design` fields like `train_design`, `y_train`, etc. For continuous
#' regression with external test domains, those fields are largely bookkeeping,
#' while the \emph{actual} predictors are carried explicitly as `feature_sets` objects.
#'
#' @details
#' \strong{Where the data live.}
#' \itemize{
#'   \item Predictors live on the design: `design$X_train` and `design$X_test` (both `feature_sets`).
#'   \item Responses live on the dataset: `dataset$train_data` and `dataset$test_data` (TRxvoxel matrices per ROI).
#' }
#'
#' \strong{cv_labels semantics.}
#' The returned object passes `cv_labels = 1:T_enc` to `mvpa_design()`. These integer
#' indices are used only for length validation and fold-construction bookkeeping in the
#' generic rMVPA iterators; they are not meaningful training targets. The actual
#' predictors are carried as `design$X_train` / `design$X_test`, and models built for
#' this design (e.g. `grouped_ridge_da_model()` / `banded_ridge_da_model()`) retrieve
#' them from those fields rather than from `cv_labels`.
#'
#' \strong{Recall blocks.}
#' `block_var_test` is stored for convenience as `design$block_var_test`. Models can use it
#' to define blocked cross-validation on test/target time (e.g. leave-one-run-out), and fall back
#' to contiguous folds when only a single run is present.
#'
#' \strong{Fold-aware target builders.}
#' For unbiased domain-adaptation workflows, `feature_sets_design()` can store a
#' `target_builder` callback instead of a single fixed `X_test`. The callback is
#' invoked separately for each outer target fold and receives the source
#' predictors plus the target train/test row indices. It must return target-side
#' predictors in the original target row order, either as a `feature_sets`
#' object, a numeric matrix, or a list containing `gamma`, `X`, or `X_test`.
#' This allows upstream matching/alignment to be re-fit on target-train rows
#' only, while held-out target rows remain untouched until evaluation.
#'
#' @param X_train A `feature_sets` object (train/source predictors).
#' @param X_test Optional `feature_sets` object (test/target predictors).
#' @param block_var_test Optional integer/factor vector of length T_rec defining
#'   test/target blocks (typically runs). Used by models that do test-time CV.
#' @param target_builder Optional function that rebuilds target-domain predictors
#'   per outer target fold. It is called with named arguments including
#'   `X_train`, `X_test`, `train_idx`, `test_idx`, `fold_id`, `block_var_test`,
#'   `n_test`, `builder_data`, and `design`; callbacks may accept only the
#'   subset they need. The return value must resolve to target predictors with
#'   the same feature/set layout as `X_train`.
#' @param target_builder_data Optional object passed through to
#'   `target_builder` as `builder_data`.
#' @param n_test Optional target-domain row count used when `X_test` is not
#'   supplied. If omitted, it is inferred from `X_test` or `block_var_test`.
#' @param split_by Optional split definition passed to `mvpa_design`.
#'
#' @return An object inheriting from `mvpa_design` with class `feature_sets_design`.
#' @export
#' @seealso
#'   \code{\link{feature_sets}}, \code{\link{expected_features}},
#'   \code{\link{banded_ridge_da_model}}, \code{\link{grouped_ridge_da_model}},
#'   \code{\link{banded_ridge_da}}, \code{\link{grouped_ridge_da}}
#' @examples
#' # Train predictors (TR x features), split into named sets:
#' X_enc <- matrix(rnorm(20 * 8), 20, 8)
#' fs_enc <- feature_sets(X_enc, blocks(low = 3, sem = 5))
#'
#' # Test predictors (TR x features), for example from a soft alignment:
#' gamma <- matrix(runif(10 * 20), 10, 20)
#' gamma <- gamma / rowSums(gamma)
#' fs_rec <- expected_features(fs_enc, gamma, drop_null = FALSE, renormalize = TRUE)
#'
#' des <- feature_sets_design(fs_enc, fs_rec, block_var_test = rep(1:2, each = 5))
feature_sets_design <- function(X_train,
                                X_test = NULL,
                                block_var_test = NULL,
                                target_builder = NULL,
                                target_builder_data = NULL,
                                n_test = NULL,
                                split_by = NULL) {
  assertthat::assert_that(inherits(X_train, "feature_sets"))
  if (!is.null(target_builder) && !is.function(target_builder)) {
    stop("feature_sets_design: target_builder must be a function when provided.", call. = FALSE)
  }

  if (!is.null(X_test)) {
    X_test <- .feature_sets_validate_target_fs(
      X_train,
      X_test,
      n_test = NULL,
      context = "feature_sets_design: X_test"
    )
  }

  T_enc <- nrow(X_train$X)
  train_df <- data.frame(.dummy = rep(1L, T_enc))

  T_rec <- NULL
  if (!is.null(X_test)) {
    T_rec <- nrow(X_test$X)
  } else if (!is.null(n_test)) {
    if (!is.numeric(n_test) || length(n_test) != 1L || !is.finite(n_test) || n_test < 1) {
      stop("feature_sets_design: n_test must be a finite positive scalar when provided.", call. = FALSE)
    }
    T_rec <- as.integer(n_test)
  } else if (!is.null(block_var_test)) {
    T_rec <- length(block_var_test)
  }

  if (!is.null(block_var_test) && !is.null(T_rec) && length(block_var_test) != T_rec) {
    stop("feature_sets_design: block_var_test must have length matching target rows.", call. = FALSE)
  }
  if (!is.null(target_builder) && is.null(T_rec)) {
    stop("feature_sets_design: supply X_test, n_test, or block_var_test when using target_builder.", call. = FALSE)
  }

  test_df <- NULL
  y_test <- NULL
  if (!is.null(T_rec)) {
    test_df <- data.frame(.dummy = rep(1L, T_rec))
    y_test <- seq_len(T_rec)
    if (!is.null(block_var_test)) {
      test_df$.block <- block_var_test
    }
  }

  mvdes <- mvpa_design(
    train_design = train_df,
    test_design = test_df,
    cv_labels = seq_len(T_enc),
    y_test = y_test,
    targets = list(X_train = X_train, X_test = X_test),
    split_by = split_by
  )

  mvdes$X_train <- X_train
  mvdes$X_test <- X_test
  mvdes$block_var_test <- block_var_test
  mvdes$target_builder <- target_builder
  mvdes$target_builder_data <- target_builder_data
  mvdes$n_test <- T_rec

  class(mvdes) <- c("feature_sets_design", class(mvdes))
  mvdes
}

.feature_sets_validate_target_fs <- function(train_fs,
                                             target_fs,
                                             n_test = NULL,
                                             context = "target feature set") {
  if (!inherits(target_fs, "feature_sets")) {
    stop(context, " must be a feature_sets object.", call. = FALSE)
  }

  if (ncol(train_fs$X) != ncol(target_fs$X)) {
    stop(context, " must have the same number of columns (features) as X_train.", call. = FALSE)
  }
  if (!identical(levels(train_fs$set), levels(target_fs$set))) {
    stop(context, " must have identical set names/order as X_train.", call. = FALSE)
  }
  if (!is.null(n_test) && nrow(target_fs$X) != n_test) {
    stop(context, " must have nrow equal to the target row count.", call. = FALSE)
  }

  if (is.null(target_fs$row_weights)) {
    target_fs$row_weights <- rep(1, nrow(target_fs$X))
  }
  if (!is.numeric(target_fs$row_weights) || length(target_fs$row_weights) != nrow(target_fs$X)) {
    stop(context, " row_weights must be numeric with length nrow(X_test).", call. = FALSE)
  }
  if (any(!is.finite(target_fs$row_weights)) || any(target_fs$row_weights < 0)) {
    stop(context, " row_weights must be finite and non-negative.", call. = FALSE)
  }

  target_fs
}

.feature_sets_target_template <- function(train_fs, X, row_weights = NULL) {
  Xt <- as.matrix(X)
  if (!is.numeric(Xt)) {
    stop("target feature matrices must be numeric.", call. = FALSE)
  }
  if (ncol(Xt) != ncol(train_fs$X)) {
    stop("target feature matrices must have the same number of columns as X_train.", call. = FALSE)
  }

  if (is.null(row_weights)) {
    row_weights <- rep(1, nrow(Xt))
  }
  if (!is.numeric(row_weights) || length(row_weights) != nrow(Xt)) {
    stop("target feature row_weights must be numeric with length nrow(X).", call. = FALSE)
  }
  if (any(!is.finite(row_weights)) || any(row_weights < 0)) {
    stop("target feature row_weights must be finite and non-negative.", call. = FALSE)
  }

  out <- train_fs
  out$X <- Xt
  out$row_weights <- as.numeric(row_weights)
  out
}

.feature_sets_materialize_target <- function(train_fs,
                                             out,
                                             n_test,
                                             context = "target_builder") {
  override_row_weights <- NULL

  if (inherits(out, "feature_sets")) {
    fs <- out
  } else if (is.matrix(out)) {
    fs <- .feature_sets_target_template(train_fs, out)
  } else if (is.list(out)) {
    if (!is.null(out$row_weights)) {
      override_row_weights <- out$row_weights
    }

    if (!is.null(out$X_test)) {
      fs <- .feature_sets_materialize_target(train_fs, out$X_test, n_test, context)
    } else if (!is.null(out$X)) {
      fs <- .feature_sets_materialize_target(train_fs, out$X, n_test, context)
    } else if (!is.null(out$gamma)) {
      drop_null <- if (!is.null(out$drop_null)) out$drop_null else TRUE
      renormalize <- if (!is.null(out$renormalize)) out$renormalize else FALSE
      fs <- expected_features(
        train_fs,
        out$gamma,
        drop_null = drop_null,
        renormalize = renormalize
      )
    } else {
      stop(
        context,
        " must return a feature_sets object, a numeric matrix, or a list with `X_test`, `X`, or `gamma`.",
        call. = FALSE
      )
    }
  } else {
    stop(
      context,
      " must return a feature_sets object, a numeric matrix, or a compatible list.",
      call. = FALSE
    )
  }

  if (!is.null(override_row_weights)) {
    fs$row_weights <- override_row_weights
  }

  .feature_sets_validate_target_fs(train_fs, fs, n_test = n_test, context = context)
}

.feature_sets_invoke_builder <- function(builder, args) {
  fmls <- tryCatch(names(formals(builder)), error = function(e) NULL)
  if (is.null(fmls) || "..." %in% fmls) {
    return(do.call(builder, args))
  }
  do.call(builder, args[intersect(names(args), fmls)])
}

.feature_sets_build_target_fold <- function(design, train_idx, test_idx, fold_id = NULL) {
  if (is.null(design$target_builder)) {
    stop(".feature_sets_build_target_fold: design does not have a target_builder.", call. = FALSE)
  }
  if (is.null(design$n_test)) {
    stop(".feature_sets_build_target_fold: design is missing n_test.", call. = FALSE)
  }

  builder_args <- list(
    X_train = design$X_train,
    X_test = design$X_test,
    train_idx = train_idx,
    test_idx = test_idx,
    fold_id = fold_id,
    block_var_test = design$block_var_test,
    n_test = design$n_test,
    builder_data = design$target_builder_data,
    design = design
  )

  out <- .feature_sets_invoke_builder(design$target_builder, builder_args)
  .feature_sets_materialize_target(
    design$X_train,
    out,
    n_test = design$n_test,
    context = "feature_sets_design: target_builder"
  )
}

.feature_sets_precompute_fold_targets <- function(design, folds) {
  if (is.null(design$target_builder)) {
    return(NULL)
  }

  lapply(seq_along(folds), function(i) {
    .feature_sets_build_target_fold(
      design = design,
      train_idx = folds[[i]]$train,
      test_idx = folds[[i]]$test,
      fold_id = i
    )
  })
}

#' Print method for feature_sets_design
#'
#' @param x A feature_sets_design object
#' @param ... ignored
#' @return Invisibly returns the input object \code{x} (called for side effects).
#' @examples
#' \dontrun{
#'   # print method called on feature_sets_design object
#' }
#' @export
print.feature_sets_design <- function(x, ...) {
  cat("feature_sets_design\n")
  cat("=================\n\n")
  if (!is.null(x$X_train) && inherits(x$X_train, "feature_sets")) {
    cat(sprintf("Train (encoding): %d x %d (sets: %s)\n",
                nrow(x$X_train$X), ncol(x$X_train$X),
                paste(names(x$X_train$indices), collapse = ", ")))
  }
  if (!is.null(x$X_test) && inherits(x$X_test, "feature_sets")) {
    cat(sprintf("Test (recall):    %d x %d\n",
                nrow(x$X_test$X), ncol(x$X_test$X)))
    if (!is.null(x$block_var_test)) {
      cat(sprintf("Recall blocks:    %d\n", length(unique(x$block_var_test))))
    }
  } else if (!is.null(x$n_test)) {
    cat(sprintf("Test (recall):    %d x %d <built per fold>\n",
                x$n_test, ncol(x$X_train$X)))
    if (!is.null(x$block_var_test)) {
      cat(sprintf("Recall blocks:    %d\n", length(unique(x$block_var_test))))
    }
  } else {
    cat("Test (recall):    <none>\n")
  }
  if (!is.null(x$target_builder)) {
    cat("Target builder:   yes\n")
  }
  invisible(x)
}
