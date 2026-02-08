#' Grouped (banded) ridge domain-adaptation model (continuous predictors -> brain)
#'
#' Fit multivariate linear models from continuous predictors (TR x features)
#' to ROI/searchlight activity (TR x voxels), with predictors organized into named
#' \emph{feature sets}. This supports both:
#' \itemize{
#'   \item \strong{domain adaptation}: allow or enforce similarity between train/source and test/target mappings, and
#'   \item \strong{feature-set competition}: quantify what each feature set adds beyond the others.
#' }
#'
#' @details
#' \strong{Data model.}
#' For each ROI/searchlight, rMVPA provides:
#' \itemize{
#'   \item Train/source responses \eqn{Y_{train}} from `dataset$train_data` (rows = TRs, columns = voxels).
#'   \item Test/target responses \eqn{Y_{test}} from `dataset$test_data` (rows = TRs, columns = voxels).
#' }
#' Predictors are stored on a `feature_sets_design`:
#' \itemize{
#'   \item Train/source predictors \eqn{X_{train}} in `design$X_train$X`.
#'   \item Test/target predictors \eqn{X_{test}} in `design$X_test$X`.
#' }
#'
#' \strong{Grouped (banded) ridge.}
#' Each predictor column belongs to a named feature set (e.g. low/mid/high/sem). The
#' ridge penalty is applied per-column based on that set, i.e. a diagonal penalty
#' \eqn{P} where entries for columns in set \eqn{g} take value \eqn{\lambda_g}.
#'
#' \strong{Two training modes.}
#' \itemize{
#'   \item \emph{stacked}: fit a single \eqn{\beta} by stacking train and test rows
#'         (with test down-weighted by `alpha_recall`/`alpha_target` and optional test TR weights from `design$X_test$row_weights`).
#'   \item \emph{coupled}: fit \eqn{\beta_{train}} and \eqn{\beta_{test}} with coupling strength `rho`,
#'         allowing a controlled train->test shift while still borrowing strength across domains.
#' }
#'
#' \strong{Test/target-time cross-validation (default).}
#' If test blocks/runs are provided (`design$block_var_test`) and there are at least two
#' unique blocks, evaluation uses leave-one-block-out on test time. If only one block
#' is present, evaluation falls back to contiguous folds over test TRs (`recall_nfolds`/`target_nfolds`).
#' In all cases, \emph{all train TRs are always included in training} for each fold.
#' This avoids transductive evaluation when test data are used in training.
#'
#' \strong{Feature-set attribution via} \eqn{\Delta R^2}.
#' When `compute_delta_r2 = TRUE`, the model computes leave-one-set-out unique contribution
#' on held-out test:
#' \deqn{\Delta R^2_g = R^2_{\mathrm{full}} - R^2_{-g},}
#' where \eqn{R^2_{-g}} is obtained by refitting the model without feature set \eqn{g}.
#'
#' @param dataset An `mvpa_dataset` with `train_data` (train/source) and `test_data` (test/target).
#' @param design A `feature_sets_design` created by `feature_sets_design()`.
#' @param mode `"stacked"` or `"coupled"`.
#' @param lambdas Named numeric vector of ridge penalties per feature set.
#'   Names must match `names(design$X_train$indices)`.
#' @param alpha_recall Non-negative scalar weighting the test/target loss (default 0.2).
#'   This argument is kept for backward compatibility with the original encoding/recall use case.
#' @param alpha_target Optional non-negative scalar weighting the test/target loss.
#'   If provided, overrides `alpha_recall`.
#' @param rho Non-negative scalar coupling strength (only used when `mode="coupled"`).
#' @param recall_folds Optional explicit list of folds, each a list with `train` and `test`
#'   integer indices over test rows. If NULL, folds are derived from `design$block_var_test`
#'   when available; otherwise, contiguous K-fold splits are used.
#' @param target_folds Optional explicit folds over test rows. If provided, overrides `recall_folds`.
#' @param recall_nfolds Integer number of contiguous folds to use when only one test block is present.
#' @param target_nfolds Integer number of contiguous folds to use when only one test block is present.
#'   If provided, overrides `recall_nfolds`.
#' @param compute_delta_r2 Logical; if TRUE, compute leave-one-set-out unique contribution
#'   \eqn{\Delta R^2} on held-out test (default TRUE).
#' @param delta_sets Optional character vector of set names to compute \eqn{\Delta R^2} for
#'   (default: all sets).
#' @param return_diagnostics Logical; if TRUE, store fold-level diagnostics (fold definitions,
#'   per-fold R-squared/MSE, hyperparameters) in `regional_mvpa_result$fits` when running
#'   `run_regional()` (default FALSE).
#' @param ... Additional arguments (currently unused).
#'
#' @return A model spec of class `banded_ridge_da_model` compatible with `run_regional()` and `run_searchlight()`.
#' @export
#' @seealso
#'   \code{\link{feature_sets}}, \code{\link{expected_features}}, \code{\link{feature_sets_design}},
#'   \code{\link{banded_ridge_da}}, \code{\link{grouped_ridge_da}}
#' @examples
#' \dontrun{
#' # Build encoding predictors and declare feature sets
#' fs_enc <- feature_sets(X_enc, blocks(low = 100, mid = 100, high = 100, sem = 100))
#'
#' # Build recall predictors from a soft alignment posterior gamma
#' fs_rec <- expected_features(fs_enc, gamma, drop_null = TRUE, renormalize = FALSE)
#'
#' # Create design and model spec
#' des <- feature_sets_design(fs_enc, fs_rec, block_var_test = recall_runs)
#' ms <- grouped_ridge_da_model(
#'   dataset = dset,
#'   design = des,
#'   mode = "coupled",
#'   lambdas = c(low = 10, mid = 10, high = 10, sem = 10),
#'   alpha_recall = 0.2,
#'   rho = 5,
#'   compute_delta_r2 = TRUE,
#'   return_diagnostics = TRUE
#' )
#'
#' res <- run_regional(ms, region_mask)
#' }
banded_ridge_da_model <- function(dataset,
                                  design,
                                  mode = c("stacked", "coupled"),
                                  lambdas,
                                  alpha_recall = 0.2,
                                  alpha_target = NULL,
                                  rho = 5,
                                  recall_folds = NULL,
                                  target_folds = NULL,
                                  recall_nfolds = 5L,
                                  target_nfolds = NULL,
                                  compute_delta_r2 = TRUE,
                                  delta_sets = NULL,
                                  return_diagnostics = FALSE,
                                  ...) {

  mode <- match.arg(mode)
  assertthat::assert_that(inherits(dataset, "mvpa_dataset"))
  assertthat::assert_that(inherits(design, "feature_sets_design"))

  if (!has_test_set(dataset) || is.null(dataset$test_data) || is.null(design$X_test)) {
    stop("banded_ridge_da_model requires an external test set (dataset$test_data + design$X_test).", call. = FALSE)
  }

  if (!is.numeric(alpha_recall) || length(alpha_recall) != 1L || !is.finite(alpha_recall) || alpha_recall < 0) {
    stop("banded_ridge_da_model: alpha_recall must be a finite non-negative scalar.", call. = FALSE)
  }
  if (!is.null(alpha_target)) {
    if (!is.numeric(alpha_target) || length(alpha_target) != 1L || !is.finite(alpha_target) || alpha_target < 0) {
      stop("banded_ridge_da_model: alpha_target must be a finite non-negative scalar.", call. = FALSE)
    }
    if (!missing(alpha_recall) && !isTRUE(all.equal(alpha_recall, alpha_target))) {
      rlang::warn("banded_ridge_da_model: both `alpha_recall` and `alpha_target` supplied; using `alpha_target`.")
    }
  }
  if (mode == "coupled") {
    if (!is.numeric(rho) || length(rho) != 1L || !is.finite(rho) || rho < 0) {
      stop("banded_ridge_da_model: rho must be a finite non-negative scalar.", call. = FALSE)
    }
  }

  alpha <- if (!is.null(alpha_target)) alpha_target else alpha_recall

  sets <- names(design$X_train$indices)
  if (is.null(delta_sets)) delta_sets <- sets

  if (is.null(names(lambdas)) || !is.numeric(lambdas)) {
    stop("banded_ridge_da_model: lambdas must be a *named* numeric vector.", call. = FALSE)
  }
  if (!all(sets %in% names(lambdas))) {
    missing <- setdiff(sets, names(lambdas))
    stop("banded_ridge_da_model: lambdas missing entries for sets: ", paste(missing, collapse = ", "), call. = FALSE)
  }

  folds_input <- recall_folds
  if (!is.null(target_folds)) {
    if (!is.null(recall_folds)) {
      rlang::warn("banded_ridge_da_model: both `recall_folds` and `target_folds` supplied; using `target_folds`.")
    }
    folds_input <- target_folds
  }

  nfolds_input <- recall_nfolds
  if (!is.null(target_nfolds)) {
    if (!missing(recall_nfolds) && !isTRUE(all.equal(recall_nfolds, target_nfolds))) {
      rlang::warn("banded_ridge_da_model: both `recall_nfolds` and `target_nfolds` supplied; using `target_nfolds`.")
    }
    nfolds_input <- target_nfolds
  }

  # Precompute test/target folds once at spec construction time
  folds <- if (!is.null(folds_input)) {
    folds_input
  } else {
    .br_make_recall_folds(design$block_var_test, nrow(design$X_test$X), nfolds_input)
  }

  # Default regression-style "performance" fields for rMVPA maps
  perf_fun <- get_regression_perf(design$split_groups)

  create_model_spec(
    "banded_ridge_da_model",
    dataset = dataset,
    design = design,
    mode = mode,
    lambdas = lambdas,
    alpha_recall = alpha,
    alpha_target = alpha,
    rho = rho,
    recall_folds = folds,
    target_folds = folds,
    compute_delta_r2 = compute_delta_r2,
    delta_sets = delta_sets,
    return_diagnostics = return_diagnostics,
    performance = perf_fun,
    compute_performance = TRUE,
    # NOTE: rMVPA drops per-ROI `result` objects unless model_spec$return_predictions is TRUE.
    # We keep results only when diagnostics are requested, and we override run_regional
    # to avoid building a prediction table from these results.
    return_predictions = isTRUE(return_diagnostics),
    return_fits = isTRUE(return_diagnostics),
    ...
  )
}

#' Grouped (banded) ridge domain-adaptation model (alias)
#'
#' Preferred name for `banded_ridge_da_model()`. See that function for full details.
#'
#' @return A model spec of class `banded_ridge_da_model` compatible with `run_regional()` and `run_searchlight()`.
#' @rdname banded_ridge_da_model
#' @export
grouped_ridge_da_model <- function(dataset,
                                   design,
                                   mode = c("stacked", "coupled"),
                                   lambdas,
                                   alpha_recall = 0.2,
                                   alpha_target = NULL,
                                   rho = 5,
                                   recall_folds = NULL,
                                   target_folds = NULL,
                                   recall_nfolds = 5L,
                                   target_nfolds = NULL,
                                   compute_delta_r2 = TRUE,
                                   delta_sets = NULL,
                                   return_diagnostics = FALSE,
                                   ...) {
  banded_ridge_da_model(
    dataset = dataset,
    design = design,
    mode = mode,
    lambdas = lambdas,
    alpha_recall = alpha_recall,
    alpha_target = alpha_target,
    rho = rho,
    recall_folds = recall_folds,
    target_folds = target_folds,
    recall_nfolds = recall_nfolds,
    target_nfolds = target_nfolds,
    compute_delta_r2 = compute_delta_r2,
    delta_sets = delta_sets,
    return_diagnostics = return_diagnostics,
    ...
  )
}

#' Convenience wrapper: build a grouped-ridge domain-adaptation model from matrices
#'
#' `banded_ridge_da()` is a convenience wrapper that builds:
#' \itemize{
#'   \item a `feature_sets` object for train predictors,
#'   \item a test `feature_sets` object (from `X_test` or from `gamma` via `expected_features()`),
#'   \item a `feature_sets_design`, and
#'   \item the final `banded_ridge_da_model` spec.
#' }
#'
#' Use this when you already have \code{X_train} (TR x features) as a single matrix
#' and you want to declare sets via `blocks()` or `by_set()`.
#'
#' @param dataset mvpa_dataset with train_data/test_data.
#' @param X_train Train predictor matrix (T_train x D) or a `feature_sets` object.
#' @param spec Feature-set spec for matrix inputs, created by `blocks()` or `by_set()`.
#'   Ignored if `X_train` is already a `feature_sets`.
#' @param X_test Optional test predictor matrix (T_test x D) or a `feature_sets` object.
#' @param gamma Optional alignment matrix used when `X_test` is NULL. See `expected_features()`.
#' @param drop_null,renormalize Passed to `expected_features()` when using `gamma`.
#' @param block_var_test Optional test run/block vector (length T_test).
#' @param ... Passed through to `banded_ridge_da_model()`.
#'
#' @return A model spec of class `banded_ridge_da_model`.
#' @export
#' @seealso \code{\link{banded_ridge_da_model}}, \code{\link{grouped_ridge_da_model}},
#'   \code{\link{blocks}}, \code{\link{by_set}}
#' @examples
#' \dontrun{
#' ms <- grouped_ridge_da(
#'   dataset = dset,
#'   X_train = X_enc,
#'   spec = blocks(low = 100, mid = 100, high = 100, sem = 100),
#'   gamma = gamma,
#'   block_var_test = recall_runs,
#'   mode = "stacked",
#'   lambdas = c(low = 10, mid = 10, high = 10, sem = 10),
#'   alpha_recall = 0.2
#' )
#' }
banded_ridge_da <- function(dataset,
                            X_train,
                            spec = NULL,
                            X_test = NULL,
                            gamma = NULL,
                            drop_null = TRUE,
                            renormalize = FALSE,
                            block_var_test = NULL,
                            ...) {
  fs_train <- if (inherits(X_train, "feature_sets")) {
    X_train
  } else {
    feature_sets(X_train, spec = spec)
  }

  fs_test <- NULL
  if (!is.null(X_test)) {
    fs_test <- if (inherits(X_test, "feature_sets")) {
      X_test
    } else if (inherits(fs_train, "feature_sets") && is.matrix(X_test)) {
      Xt <- as.matrix(X_test)
      if (ncol(Xt) != ncol(fs_train$X)) {
        stop("banded_ridge_da: ncol(X_test) must match ncol(X_train) when X_train is a feature_sets.", call. = FALSE)
      }
      tmp <- fs_train
      tmp$X <- Xt
      tmp$row_weights <- rep(1, nrow(Xt))
      tmp
    } else {
      feature_sets(X_test, spec = spec, set_order = fs_train$set_order)
    }
  } else if (!is.null(gamma)) {
    fs_test <- expected_features(fs_train, gamma, drop_null = drop_null, renormalize = renormalize)
  }

  des <- feature_sets_design(
    X_train = fs_train,
    X_test = fs_test,
    block_var_test = block_var_test
  )

  banded_ridge_da_model(dataset = dataset, design = des, ...)
}

#' Grouped (banded) ridge DA convenience wrapper (alias)
#'
#' Preferred name for `banded_ridge_da()`. See that function for full details.
#'
#' @return A model spec of class `banded_ridge_da_model`.
#' @rdname banded_ridge_da
#' @export
grouped_ridge_da <- function(dataset,
                             X_train,
                             spec = NULL,
                             X_test = NULL,
                             gamma = NULL,
                             drop_null = TRUE,
                             renormalize = FALSE,
                             block_var_test = NULL,
                             ...) {
  banded_ridge_da(
    dataset = dataset,
    X_train = X_train,
    spec = spec,
    X_test = X_test,
    gamma = gamma,
    drop_null = drop_null,
    renormalize = renormalize,
    block_var_test = block_var_test,
    ...
  )
}
# ---- internal helpers ------------------------------------------------------

.br_make_recall_folds <- function(block_var_test, T_rec, nfolds) {
  if (!is.numeric(T_rec) || length(T_rec) != 1L || T_rec < 2L) {
    stop(".br_make_recall_folds: T_rec must be >= 2.", call. = FALSE)
  }

  if (!is.null(block_var_test)) {
    if (length(block_var_test) != T_rec) {
      stop(".br_make_recall_folds: block_var_test length mismatch with T_rec.", call. = FALSE)
    }
    blocks <- unique(block_var_test)
    if (length(blocks) >= 2L) {
      folds <- lapply(blocks, function(b) {
        test <- which(block_var_test == b)
        train <- setdiff(seq_len(T_rec), test)
        list(train = train, test = test)
      })
      return(folds)
    }
  }

  # Single-block fallback: contiguous folds over time
  nfolds <- as.integer(nfolds)
  if (!is.finite(nfolds) || nfolds < 2L) {
    stop(".br_make_recall_folds: nfolds must be >= 2.", call. = FALSE)
  }
  # ensure each test fold has at least 2 TRs when possible
  nfolds <- min(nfolds, max(2L, floor(T_rec / 2L)))
  groups <- split(seq_len(T_rec), cut(seq_len(T_rec), nfolds, labels = FALSE))
  folds <- lapply(groups, function(test) {
    train <- setdiff(seq_len(T_rec), test)
    list(train = train, test = test)
  })
  folds
}

.br_center_scale_fit <- function(X) {
  mu <- colMeans(X)
  sd <- apply(X, 2, stats::sd)
  sd[!is.finite(sd) | sd == 0] <- 1
  list(mu = mu, sd = sd)
}

.br_center_scale_apply <- function(X, mu, sd) {
  sweep(sweep(X, 2, mu, "-"), 2, sd, "/")
}

.br_center_Y <- function(Y) {
  mu <- colMeans(Y)
  list(mu = mu, Yc = sweep(Y, 2, mu, "-"))
}

.br_penalty_vec <- function(fs, lambdas) {
  set_names <- as.character(fs$set)
  as.numeric(lambdas[set_names])
}

.br_r2_mean <- function(Y_true, Y_pred, eps = 1e-12) {
  y_mu <- colMeans(Y_true)
  ss_tot <- colSums((sweep(Y_true, 2, y_mu, "-"))^2) + eps
  ss_res <- colSums((Y_true - Y_pred)^2)
  r2 <- 1 - (ss_res / ss_tot)
  mean(r2, na.rm = TRUE)
}

.br_mse <- function(Y_true, Y_pred) {
  mean((Y_true - Y_pred)^2, na.rm = TRUE)
}

.br_fit_stacked <- function(Xe, Ye, Xr_tr, Yr_tr, pen, alpha_recall, w_rec_tr) {
  # Weight recall rows
  w <- sqrt(alpha_recall) * sqrt(pmax(w_rec_tr, 0))
  if (all(w == 0)) {
    X_stack <- Xe
    Y_stack <- Ye
  } else {
    X_stack <- rbind(Xe, sweep(Xr_tr, 1, w, "*"))
    Y_stack <- rbind(Ye, sweep(Yr_tr, 1, w, "*"))
  }

  xs <- .br_center_scale_fit(X_stack)
  Xs <- .br_center_scale_apply(X_stack, xs$mu, xs$sd)

  yc <- .br_center_Y(Y_stack)
  Ys <- yc$Yc

  A <- crossprod(Xs) + diag(pen, ncol(Xs))
  B <- crossprod(Xs, Ys)
  beta <- solve(A, B)

  list(beta = beta, x_mu = xs$mu, x_sd = xs$sd, y_mu = yc$mu)
}

.br_fit_coupled <- function(Xe, Ye, Xr_tr, Yr_tr, pen, alpha_recall, rho, w_rec_tr) {
  D <- ncol(Xe)
  I <- diag(1, D)
  P <- diag(pen, D)

  # Standardize X jointly (unweighted)
  xs <- .br_center_scale_fit(rbind(Xe, Xr_tr))
  Xe_s <- .br_center_scale_apply(Xe, xs$mu, xs$sd)
  Xr_s <- .br_center_scale_apply(Xr_tr, xs$mu, xs$sd)

  # Center Y jointly (unweighted)
  yc <- .br_center_Y(rbind(Ye, Yr_tr))
  Ye_c <- sweep(Ye, 2, yc$mu, "-")
  Yr_c <- sweep(Yr_tr, 2, yc$mu, "-")

  # Apply recall weights into the recall normal equations
  w <- sqrt(alpha_recall) * sqrt(pmax(w_rec_tr, 0))
  Xr_w <- sweep(Xr_s, 1, w, "*")
  Yr_w <- sweep(Yr_c, 1, w, "*")

  A11 <- crossprod(Xe_s) + P + rho * I
  A22 <- crossprod(Xr_w) + P + rho * I
  A12 <- -rho * I
  A21 <- -rho * I

  A <- rbind(
    cbind(A11, A12),
    cbind(A21, A22)
  )

  B1 <- crossprod(Xe_s, Ye_c)
  B2 <- crossprod(Xr_w, Yr_w)
  B <- rbind(B1, B2)

  sol <- solve(A, B)
  beta_enc <- sol[seq_len(D), , drop = FALSE]
  beta_rec <- sol[(D + 1L):(2L * D), , drop = FALSE]

  list(beta_enc = beta_enc, beta_rec = beta_rec, x_mu = xs$mu, x_sd = xs$sd, y_mu = yc$mu)
}

.br_predict <- function(X, fit, mode = c("stacked", "coupled")) {
  mode <- match.arg(mode)
  Xs <- .br_center_scale_apply(X, fit$x_mu, fit$x_sd)
  beta <- if (mode == "coupled") fit$beta_rec else fit$beta
  sweep(Xs %*% beta, 2, fit$y_mu, "+")
}

# ---- per-ROI processing ----------------------------------------------------

#' @export
compute_performance.banded_ridge_da_model <- function(obj, result) {
  # This model typically computes performance inside process_roi and returns it
  # directly. Provide a safe fallback for callers that invoke compute_performance().
  if (is.numeric(result) && !is.null(names(result))) {
    return(result)
  }
  if (is.list(result) && !is.null(result$performance) &&
      is.numeric(result$performance) && !is.null(names(result$performance))) {
    return(result$performance)
  }
  if (is.list(result) && !is.null(result$predictor) && is.null(result$observed)) {
    return(NULL)
  }
  obj$performance(result)
}

#' @keywords internal
#' @export
process_roi.banded_ridge_da_model <- function(mod_spec,
                                              roi,
                                              rnum,
                                              center_global_id = NA,
                                              ...) {
  if (!has_test_set(mod_spec)) {
    return(tibble::tibble(
      result = list(NULL),
      indices = list(neuroim2::indices(roi$train_roi)),
      performance = list(NULL),
      id = rnum,
      error = TRUE,
      error_message = "banded_ridge_da_model requires external test set"
    ))
  }

  Xtrain_fs <- mod_spec$design$X_train
  Xtest_fs <- mod_spec$design$X_test

  Ye <- as.matrix(neuroim2::values(roi$train_roi))
  Yr <- as.matrix(neuroim2::values(roi$test_roi))

  if (nrow(Ye) < 2L || nrow(Yr) < 2L || ncol(Ye) < 1L) {
    return(tibble::tibble(
      result = list(NULL),
      indices = list(neuroim2::indices(roi$train_roi)),
      performance = list(NULL),
      id = rnum,
      error = TRUE,
      error_message = "banded_ridge_da_model: insufficient observations or voxels"
    ))
  }

  Xe <- Xtrain_fs$X
  Xr <- Xtest_fs$X
  if (nrow(Xe) != nrow(Ye)) {
    return(tibble::tibble(
      result = list(NULL),
      indices = list(neuroim2::indices(roi$train_roi)),
      performance = list(NULL),
      id = rnum,
      error = TRUE,
      error_message = "banded_ridge_da_model: nrow(X_train) must match nrow(Y_enc) for this subject"
    ))
  }
  if (nrow(Xr) != nrow(Yr)) {
    return(tibble::tibble(
      result = list(NULL),
      indices = list(neuroim2::indices(roi$train_roi)),
      performance = list(NULL),
      id = rnum,
      error = TRUE,
      error_message = "banded_ridge_da_model: nrow(X_test) must match nrow(Y_rec) for this subject"
    ))
  }

  pen_full <- .br_penalty_vec(Xtrain_fs, mod_spec$lambdas)
  if (any(!is.finite(pen_full)) || any(pen_full < 0)) {
    return(tibble::tibble(
      result = list(NULL),
      indices = list(neuroim2::indices(roi$train_roi)),
      performance = list(NULL),
      id = rnum,
      error = TRUE,
      error_message = "banded_ridge_da_model: invalid lambdas (must be finite and non-negative)"
    ))
  }

  folds <- mod_spec$recall_folds
  w_rec <- Xtest_fs$row_weights
  if (is.null(w_rec) || length(w_rec) != nrow(Xr)) w_rec <- rep(1, nrow(Xr))

  mode <- mod_spec$mode
  alpha <- mod_spec$alpha_recall
  rho <- mod_spec$rho

  # Full model fold scores
  r2_full <- numeric(0)
  mse_full <- numeric(0)
  per_fold_diag <- list()

  fit_fold <- function(keep_cols, lambdas_pen) {
    Xe_k <- Xe[, keep_cols, drop = FALSE]
    Xr_k <- Xr[, keep_cols, drop = FALSE]
    pen_k <- lambdas_pen[keep_cols]

    fold_scores <- lapply(seq_along(folds), function(fi) {
      tr <- folds[[fi]]$train
      te <- folds[[fi]]$test
      if (length(te) < 2L || length(tr) < 2L) {
        return(list(r2 = NA_real_, mse = NA_real_))
      }

      Xr_tr <- Xr_k[tr, , drop = FALSE]
      Yr_tr <- Yr[tr, , drop = FALSE]
      Xr_te <- Xr_k[te, , drop = FALSE]
      Yr_te <- Yr[te, , drop = FALSE]
      w_tr <- w_rec[tr]

      fit <- if (mode == "stacked") {
        .br_fit_stacked(Xe_k, Ye, Xr_tr, Yr_tr, pen_k, alpha, w_tr)
      } else {
        .br_fit_coupled(Xe_k, Ye, Xr_tr, Yr_tr, pen_k, alpha, rho, w_tr)
      }

      Yhat <- .br_predict(Xr_te, fit, mode = mode)
      list(
        r2 = .br_r2_mean(Yr_te, Yhat),
        mse = .br_mse(Yr_te, Yhat)
      )
    })

    r2s <- vapply(fold_scores, `[[`, numeric(1), "r2")
    mses <- vapply(fold_scores, `[[`, numeric(1), "mse")
    list(r2_mean = mean(r2s, na.rm = TRUE), mse_mean = mean(mses, na.rm = TRUE),
         r2_folds = r2s, mse_folds = mses)
  }

  full_cols <- seq_len(ncol(Xe))
  full_fit_res <- fit_fold(full_cols, pen_full)
  r2_full_mean <- full_fit_res$r2_mean
  mse_full_mean <- full_fit_res$mse_mean
  r2_full_folds <- full_fit_res$r2_folds

  perf <- c(recall_r2_full = r2_full_mean,
            recall_mse_full = mse_full_mean,
            target_r2_full = r2_full_mean,
            target_mse_full = mse_full_mean)

  # Leave-one-set-out unique contribution on recall
  if (isTRUE(mod_spec$compute_delta_r2)) {
    delta_sets <- mod_spec$delta_sets
    delta_sets <- intersect(delta_sets, names(Xtrain_fs$indices))
    for (set_name in delta_sets) {
      drop_cols <- Xtrain_fs$indices[[set_name]]
      keep_cols <- setdiff(full_cols, drop_cols)
      if (length(keep_cols) < 1L) {
        perf[paste0("delta_r2_", set_name)] <- NA_real_
        next
      }
      drop_fit_res <- fit_fold(keep_cols, pen_full)
      # Compute delta per-fold to keep fold pairing consistent
      delta_fold <- r2_full_folds - drop_fit_res$r2_folds
      perf[paste0("delta_r2_", set_name)] <- mean(delta_fold, na.rm = TRUE)
    }
  }

  diag_obj <- NULL
  if (isTRUE(mod_spec$return_diagnostics)) {
    diag_obj <- list(
      folds = folds,
      r2_full_folds = r2_full_folds,
      mse_full_folds = full_fit_res$mse_folds,
      lambdas = mod_spec$lambdas,
      alpha_recall = alpha,
      alpha_target = alpha,
      rho = if (mode == "coupled") rho else NA_real_,
      mode = mode
    )
  }

  tibble::tibble(
    result = list(list(predictor = diag_obj)),
    indices = list(neuroim2::indices(roi$train_roi)),
    performance = list(perf),
    id = rnum,
    error = FALSE,
    error_message = "~"
  )
}

#' @rdname run_regional-methods
#' @export
run_regional.banded_ridge_da_model <- function(model_spec,
                                               region_mask,
                                               return_fits = model_spec$return_fits,
                                               compute_performance = TRUE,
                                               coalesce_design_vars = FALSE,
                                               processor = NULL,
                                               verbose = FALSE,
                                               ...) {
  # Keep per-ROI results for optional fits, but do not attempt to build a prediction table.
  run_regional_base(
    model_spec,
    region_mask,
    coalesce_design_vars = coalesce_design_vars,
    processor = processor,
    verbose = verbose,
    compute_performance = compute_performance,
    return_predictions = FALSE,
    return_fits = return_fits,
    ...
  )
}
