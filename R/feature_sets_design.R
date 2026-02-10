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
#' \strong{Why `y_train` is a dummy.}
#' The returned object includes a dummy `y_train = 1:T_enc` (and `y_test = 1:T_rec` when
#' `X_test` is present). This is only to satisfy length checks and indexing in the
#' generic rMVPA iterators. Models built for this design (e.g. `grouped_ridge_da_model()` /
#' `banded_ridge_da_model()`)
#' should ignore `y_train` and use `design$X_train` / `design$X_test` instead.
#'
#' \strong{Recall blocks.}
#' `block_var_test` is stored for convenience as `design$block_var_test`. Models can use it
#' to define blocked cross-validation on test/target time (e.g. leave-one-run-out), and fall back
#' to contiguous folds when only a single run is present.
#'
#' @param X_train A `feature_sets` object (train/source predictors).
#' @param X_test Optional `feature_sets` object (test/target predictors).
#' @param block_var_test Optional integer/factor vector of length T_rec defining
#'   test/target blocks (typically runs). Used by models that do test-time CV.
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
                                split_by = NULL) {
  assertthat::assert_that(inherits(X_train, "feature_sets"))

  if (!is.null(X_test)) {
    assertthat::assert_that(inherits(X_test, "feature_sets"))
    # Must have the same set layout / feature dimension for paired-domain analyses
    if (ncol(X_train$X) != ncol(X_test$X)) {
      stop("feature_sets_design: X_train and X_test must have the same number of columns (features).", call. = FALSE)
    }
    if (!identical(levels(X_train$set), levels(X_test$set))) {
      stop("feature_sets_design: X_train and X_test must have identical set names/order.", call. = FALSE)
    }
  }

  T_enc <- nrow(X_train$X)
  train_df <- data.frame(.y = seq_len(T_enc))

  test_df <- NULL
  y_test <- NULL
  if (!is.null(X_test)) {
    T_rec <- nrow(X_test$X)
    test_df <- data.frame(.y = seq_len(T_rec))
    y_test <- seq_len(T_rec)
    if (!is.null(block_var_test)) {
      if (length(block_var_test) != T_rec) {
        stop("feature_sets_design: block_var_test must have length nrow(X_test$X).", call. = FALSE)
      }
      test_df$.block <- block_var_test
    }
  }

  mvdes <- mvpa_design(
    train_design = train_df,
    test_design = test_df,
    y_train = seq_len(T_enc),
    y_test = y_test,
    split_by = split_by
  )

  mvdes$X_train <- X_train
  mvdes$X_test <- X_test
  mvdes$block_var_test <- block_var_test

  class(mvdes) <- c("feature_sets_design", class(mvdes))
  mvdes
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
  } else {
    cat("Test (recall):    <none>\n")
  }
  invisible(x)
}
