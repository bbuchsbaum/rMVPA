# pattern_model: pattern-first spatial reduced-rank MVPA
#
# One forward model, x = A C' y_w + eps, fitted on a whole domain (global) or
# on each region (regional / searchlight), yields classification, multivariate
# decoding, encoding, and forward patterns from the same estimate. This file
# is the rMVPA integration layer: constructor, fit_roi(), output_schema(),
# run_global(), and the domain evaluator they share. Numerics live in
# pattern_core.R / pattern_noise.R / pattern_predict.R.

#' Pattern-First Spatial Reduced-Rank MVPA Model
#'
#' Creates a \code{pattern_model} specification. The model fits a low-rank
#' forward model \eqn{x = A C^\top y + \epsilon} with a structured residual
#' covariance and derives condition classification, multivariate target
#' decoding, and encoding predictions from that single fit. It works with
#' categorical targets (an \code{\link{mvpa_design}} built from
#' \code{y_train}) and with vector-valued continuous targets
#' (\code{targets = } a numeric matrix, a \code{\link{feature_sets_design}},
#' or a \code{\link{feature_rsa_design}}), read through
#' \code{\link{model_targets}}.
#'
#' @section What is fitted:
#' Targets are centred and whitened on the training rows; \eqn{C} has
#' orthonormal columns so all scale lives in the spatial patterns \eqn{A}.
#' The residual covariance \eqn{\Psi = D + UU^\top} is estimated once from a
#' training-only pilot fit and held fixed. Rank is chosen by nested,
#' block-aware cross-validation on held-out decoding loss (log loss for
#' categorical targets, normalized squared error for continuous targets).
#' Spatial penalties are declared through \code{penalty} but are not yet
#' available in this version: any non-\code{NULL} penalty is rejected.
#'
#' @section Outputs:
#' Regional and searchlight runs report scalar metrics per ROI (see
#' \code{\link{output_schema}}): \code{Accuracy}, \code{AUC}, \code{logloss},
#' and \code{rank_mean} for categorical targets; \code{R2}, \code{RMSE},
#' \code{cor}, and \code{rank_mean} for continuous targets. \code{rank_mean}
#' is the mean rank selected across folds and need not be a whole number. With
#' \code{return_predictions = TRUE} a regional run also returns the usual
#' out-of-fold prediction table (categorical targets), and with either
#' \code{return_predictions = TRUE} or \code{return_fits = TRUE} each ROI keeps
#' its prediction ledger (plus the fold fits when \code{return_fits = TRUE})
#' under \code{$fits} of the regional result. \code{\link{run_global}} fits the
#' whole domain once and returns a \code{pattern_global_result}.
#'
#' Two ledgers are kept. The fold-resolved one records every prediction with
#' the fold that produced it. The pooled one holds a single record per tested
#' observation, in sorted order, with repeats averaged; it is what the metrics
#' and the prediction table are built from, so a cross-validation scheme that
#' tests a row several times (\code{\link{twofold_blocked_cross_validation}},
#' sequential, or bootstrap) does not give that row extra weight. This matches
#' how every other rMVPA model aggregates repeated predictions.
#'
#' @param dataset An \code{\link{mvpa_dataset}} (image, multibasis, surface,
#'   or clustered).
#' @param design A design with a \code{\link{model_targets}} method.
#' @param crossval A cross-validation specification. Defaults to
#'   \code{\link{blocked_cross_validation}} on the design's block variable,
#'   or 5-fold cross-validation when no blocks are available.
#' @param rank \code{"auto"} (selected by nested cross-validation) or a fixed
#'   positive integer, capped at the eligible rank.
#' @param max_rank Largest rank considered when \code{rank = "auto"}.
#' @param penalty Spatial penalty specification, a list with elements
#'   \code{sparse}, \code{signed_smooth}, and \code{support_smooth}. Reserved:
#'   this version accepts only \code{NULL}.
#' @param noise Residual covariance specification passed to
#'   \code{\link{pattern_control}}.
#' @param control A \code{\link{pattern_control}} object. When supplied it
#'   overrides \code{max_rank} and \code{noise}.
#' @param return_predictions Retain out-of-fold prediction ledgers per ROI.
#' @param return_fits Retain the fold fits per ROI (implies ledgers).
#' @param refit Also fit the model on all training rows (a descriptive
#'   deployment fit, distinct from the cross-validated evidence).
#' @param ... Additional fields stored on the specification.
#' @return A \code{pattern_model} specification (class
#'   \code{c("pattern_model", "model_spec")}).
#' @examples
#' ds <- gen_sample_dataset(c(6, 6, 4), 60, nlevels = 3, blocks = 3)
#' spec <- pattern_model(ds$dataset, ds$design, max_rank = 2)
#' res <- run_global(spec)
#' res$performance_table
#' @seealso \code{\link{run_global}}, \code{\link{run_regional}},
#'   \code{\link{predict.pattern_fit}}, \code{\link{pattern_control}}
#' @export
pattern_model <- function(dataset, design, crossval = NULL, rank = "auto", max_rank = 8L,
                          penalty = NULL,
                          noise = list(type = "diag_lowrank", rank = "auto"),
                          control = NULL,
                          return_predictions = FALSE, return_fits = FALSE, refit = FALSE, ...) {
  assertthat::assert_that(inherits(dataset, "mvpa_dataset"))
  targets_train <- model_targets(design, "train")
  if (is.null(targets_train)) {
    stop("pattern_model: the design provides no training targets.", call. = FALSE)
  }
  if (!targets_train$type %in% c("categorical", "continuous", "matrix")) {
    stop("pattern_model: unsupported target type.", call. = FALSE)
  }

  if (!is.null(penalty)) {
    penalty <- as.list(penalty)
    allowed <- c("sparse", "signed_smooth", "support_smooth")
    bad <- setdiff(names(penalty), allowed)
    if (length(bad)) {
      stop(sprintf("pattern_model: unknown penalty element(s): %s.", paste(bad, collapse = ", ")),
           call. = FALSE)
    }
    active <- vapply(penalty, function(v) !is.null(v) && !(is.numeric(v) && all(v == 0)), logical(1))
    if (any(active)) {
      stop("pattern_model: spatial penalties (sparse, signed_smooth, support_smooth) are not ",
           "implemented in this version; pass penalty = NULL.", call. = FALSE)
    }
  }

  if (is.null(control)) {
    control <- pattern_control(max_rank = max_rank, noise = noise)
  } else if (!inherits(control, "pattern_control")) {
    stop("pattern_model: 'control' must come from pattern_control().", call. = FALSE)
  }

  if (!identical(rank, "auto")) {
    if (!is.numeric(rank) || length(rank) != 1L || rank < 1) {
      stop("pattern_model: 'rank' must be \"auto\" or a positive integer.", call. = FALSE)
    }
    rank <- as.integer(rank)
  }

  if (is.null(crossval)) {
    # feature_sets_design stores its training blocks in block_var_train, the
    # same convention banded_ridge_model follows.
    bv <- design$block_var %||% design$block_var_train
    if (!is.null(bv)) {
      crossval <- blocked_cross_validation(bv)
    } else {
      warning("pattern_model: the design has no block variable; falling back to ",
              "5-fold cross-validation, which does not respect run structure. ",
              "Supply 'crossval' explicitly for time-series or run-structured data.",
              call. = FALSE)
      crossval <- kfold_cross_validation(length(targets_train$observation_ids), 5)
    }
  }

  if (!is.null(targets_train$row_weights) && any(targets_train$row_weights != 1)) {
    warning("pattern_model: observation row weights are not used by this version ",
            "of the estimator; every observation is weighted equally.", call. = FALSE)
  }

  has_test <- !is.null(dataset$test_data) && !is.null(model_targets(design, "test"))
  targets_test <- if (has_test) model_targets(design, "test") else NULL
  target_type <- if (targets_train$type == "categorical") "categorical" else "continuous"

  # Regional prediction tables (rMVPA's classification_result path) exist only
  # for categorical targets. For continuous targets the out-of-fold ledger is
  # retained under $fits instead, which the iteration engine keeps whenever
  # return_fits is set on the specification.
  keep_predictions <- isTRUE(return_predictions) && target_type == "categorical"
  retain_results <- isTRUE(return_predictions) || isTRUE(return_fits)

  create_model_spec(
    "pattern_model", dataset, design,
    return_predictions = keep_predictions,
    compute_performance = TRUE,
    crossval = crossval,
    rank = rank,
    penalty = penalty,
    control = control,
    return_fits = retain_results,
    keep_fold_fits = isTRUE(return_fits),
    refit = isTRUE(refit),
    has_test_set = has_test,
    targets_train = targets_train,
    targets_test = targets_test,
    target_type = target_type,
    ...
  )
}

#' @export
print.pattern_model <- function(x, ...) {
  tt <- x$targets_train
  q <- if (is.matrix(tt$values)) ncol(tt$values) else if (is.factor(tt$values)) nlevels(tt$values) else 1L
  cat("pattern_model specification\n")
  cat(sprintf("  targets: %s (%d observations x %d responses)\n", x$target_type,
              length(tt$observation_ids), q))
  cat(sprintf("  rank: %s (max %d), noise: %s\n",
              if (identical(x$rank, "auto")) "auto" else as.character(x$rank),
              x$control$max_rank, x$control$noise$type))
  cat(sprintf("  crossval: %s, external test set: %s\n",
              class(x$crossval)[1], if (isTRUE(x$has_test_set)) "yes" else "no"))
  cat(sprintf("  return_predictions: %s, return_fits: %s, refit: %s\n",
              x$return_predictions, x$return_fits, x$refit))
  invisible(x)
}

# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

.pattern_subset_rows <- function(t, idx) {
  if (is.matrix(t)) t[idx, , drop = FALSE] else t[idx]
}

# Outer fold train/test index pairs from the model's crossval spec.
.pattern_outer_folds <- function(crossval, n) {
  samples <- crossval_samples(crossval, data.frame(.row = seq_len(n)), seq_len(n))
  lapply(seq_len(nrow(samples)), function(k) {
    list(train = .extract_sample_indices(samples$train[[k]]),
         test = .extract_sample_indices(samples$test[[k]]))
  })
}

# Inner folds on a subset of training rows: leave-one-block-out (blocks
# merged into at most `max_folds` groups) when blocks are available, random
# k-fold otherwise. Returns index vectors relative to the subset.
.pattern_inner_folds <- function(n, blocks = NULL, max_folds = 5L) {
  if (!is.null(blocks) && length(unique(blocks)) >= 2L) {
    ids <- as.integer(factor(blocks))
    nb <- max(ids)
    grp <- if (nb > max_folds) ((ids - 1L) %% max_folds) + 1L else ids
  } else {
    k <- min(max_folds, n)
    grp <- sample(rep(seq_len(k), length.out = n))
  }
  lapply(sort(unique(grp)), function(g) {
    te <- which(grp == g)
    list(train = setdiff(seq_len(n), te), test = te)
  })
}

# Widen a fold's probability matrix to the full set of classes. A fold whose
# training rows omit a class cannot predict it, so that column is exactly zero:
# the row counts as an error for accuracy and takes the clamped log loss.
.pattern_pad_probs <- function(P, all_levels) {
  if (identical(colnames(P), all_levels)) return(P)
  out <- matrix(0, nrow(P), length(all_levels), dimnames = list(NULL, all_levels))
  shared <- intersect(colnames(P), all_levels)
  out[, shared] <- P[, shared, drop = FALSE]
  out
}

# Held-out decoding loss (lower is better). Returns NA when the fold carries no
# scorable rows so the caller can drop it rather than propagate the NA.
.pattern_loss <- function(fit, X_test, truth) {
  if (fit$y_transform$type == "categorical") {
    P <- predict(fit, X_test, type = "prob")
    truth <- factor(as.character(truth), levels = fit$y_transform$levels)
    # Rows whose class never appeared in this fold's training data carry no
    # information about rank; scoring them would make every rank equally bad.
    ok <- !is.na(truth)
    if (!any(ok)) return(NA_real_)
    p_true <- P[cbind(which(ok), as.integer(truth[ok]))]
    mean(-log(pmax(p_true, 1e-12)))
  } else {
    Yhat <- predict(fit, X_test, type = "decode")
    Y <- as.matrix(truth)
    mu <- fit$y_transform$mu           # training-mean baseline (see .pattern_score_ledger)
    sst <- colSums(sweep(Y, 2L, mu, "-")^2)
    sse <- colSums((Y - Yhat)^2)
    mean(sse / pmax(sst, .Machine$double.eps))
  }
}

# Nested rank selection: exact reduced-rank path per inner fold, summed loss.
.pattern_select_rank <- function(X, targets, blocks, control) {
  n <- nrow(X)
  folds <- .pattern_inner_folds(n, blocks)
  losses <- NULL
  n_used <- 0L
  for (f in folds) {
    if (length(f$train) < 3L || length(f$test) < 1L) next
    tr_targets <- .pattern_subset_rows(targets, f$train)
    if (is.factor(tr_targets) && nlevels(droplevels(tr_targets)) < 2L) next
    path <- tryCatch(.pattern_fit(X[f$train, , drop = FALSE], tr_targets, rank = "path", control = control),
                     error = function(e) NULL)
    if (is.null(path)) next
    l <- vapply(path, function(fit) .pattern_loss(fit, X[f$test, , drop = FALSE],
                                                 .pattern_subset_rows(targets, f$test)), numeric(1))
    if (!all(is.finite(l))) next
    if (is.null(losses)) {
      losses <- l
    } else {
      m <- min(length(losses), length(l))
      losses <- losses[seq_len(m)] + l[seq_len(m)]
    }
    n_used <- n_used + 1L
  }
  if (is.null(losses) || !length(losses) || !any(is.finite(losses))) {
    return(list(rank = 1L, losses = NULL))
  }
  # Ties go to the smaller rank: when an extra dimension buys no held-out
  # accuracy, prefer the more parsimonious model.
  list(rank = as.integer(unname(which.min(losses))), losses = unname(losses) / max(n_used, 1L))
}

# ---------------------------------------------------------------------------
# Domain evaluator shared by fit_roi() and run_global()
# ---------------------------------------------------------------------------

# Fits and evaluates the model on one domain (an ROI or the whole brain).
# Cross-validated when no external test set is present; otherwise a single
# train -> test evaluation. Returns metrics, a prediction ledger, and
# optionally the fold fits and a full-data refit.
.pattern_evaluate_domain <- function(model, X, X_test = NULL, keep_fits = FALSE, refit = FALSE) {
  X <- as.matrix(X)
  n <- nrow(X)
  tt <- model$targets_train
  targets <- tt$values
  n_targets <- if (is.matrix(targets)) nrow(targets) else length(targets)
  if (n_targets != n) {
    stop(sprintf("pattern_model: %d target rows but %d observations.", n_targets, n), call. = FALSE)
  }
  blocks <- model$crossval$block_var
  if (!is.null(blocks) && length(blocks) != n) blocks <- NULL
  control <- model$control
  categorical <- identical(model$target_type, "categorical")
  all_levels <- if (categorical) {
    lv <- model$targets_train$response_ids
    if (is.null(lv)) levels(as.factor(targets)) else lv
  } else {
    NULL
  }

  if (!is.null(X_test)) {
    if (is.null(model$targets_test)) stop("pattern_model: test data without test targets.", call. = FALSE)
    folds <- list(list(train = seq_len(n), test = seq_len(nrow(X_test))))
    external <- TRUE
  } else {
    folds <- .pattern_outer_folds(model$crossval, n)
    external <- FALSE
  }

  fold_fits <- vector("list", length(folds))
  ranks <- integer(length(folds))
  ledger_fold <- integer(0); ledger_obs <- integer(0)
  ledger_truth <- NULL; ledger_pred <- NULL; ledger_base <- NULL

  for (k in seq_along(folds)) {
    tr <- folds[[k]]$train; te <- folds[[k]]$test
    tr_targets <- .pattern_subset_rows(targets, tr)
    if (identical(model$rank, "auto")) {
      sel <- .pattern_select_rank(X[tr, , drop = FALSE], tr_targets, blocks[tr], control)
      r_k <- sel$rank
    } else {
      r_k <- model$rank
    }
    fit_k <- .pattern_fit(X[tr, , drop = FALSE], tr_targets, rank = r_k, control = control, cap_rank = TRUE)
    ranks[k] <- fit_k$rank
    if (keep_fits) fold_fits[[k]] <- fit_k

    X_te <- if (external) as.matrix(X_test)[te, , drop = FALSE] else X[te, , drop = FALSE]
    truth <- if (external) .pattern_subset_rows(model$targets_test$values, te) else .pattern_subset_rows(targets, te)
    pred <- if (categorical) {
      .pattern_pad_probs(predict(fit_k, X_te, type = "prob"), all_levels)
    } else {
      predict(fit_k, X_te, type = "decode")
    }

    ledger_fold <- c(ledger_fold, rep(k, length(te)))
    ledger_obs <- c(ledger_obs, te)
    ledger_truth <- if (is.null(ledger_truth)) truth else {
      if (is.matrix(truth)) rbind(ledger_truth, truth) else
        factor(c(as.character(ledger_truth), as.character(truth)), levels = levels(ledger_truth))
    }
    ledger_pred <- if (is.null(ledger_pred)) pred else rbind(ledger_pred, pred)
    if (!categorical) {
      # this fold's training-target mean, the baseline its predictions are scored against
      b <- matrix(fit_k$y_transform$mu, nrow = length(te), ncol = length(fit_k$y_transform$mu),
                  byrow = TRUE)
      ledger_base <- if (is.null(ledger_base)) b else rbind(ledger_base, b)
    }
  }

  if (categorical) {
    ledger_truth <- factor(as.character(ledger_truth), levels = all_levels)
  }
  ledger <- structure(
    list(fold = ledger_fold, observation = ledger_obs, truth = ledger_truth,
         prediction = ledger_pred, baseline = ledger_base,
         partition = if (external) "external" else "cv",
         type = model$target_type),
    class = c("pattern_ledger", "list")
  )
  pooled <- .pattern_pool_ledger(ledger)
  metrics <- .pattern_score_ledger(pooled)
  metrics <- c(metrics, rank_mean = mean(ranks))

  refit_obj <- NULL
  if (isTRUE(refit)) {
    r_all <- if (identical(model$rank, "auto")) {
      .pattern_select_rank(X, targets, blocks, control)$rank
    } else {
      model$rank
    }
    refit_obj <- .pattern_fit(X, targets, rank = r_all, control = control, cap_rank = TRUE)
  }

  list(metrics = metrics, ledger = ledger, pooled = pooled, ranks = ranks,
       fold_fits = if (keep_fits) fold_fits else NULL, refit = refit_obj)
}

# ---------------------------------------------------------------------------
# Metrics from a ledger
# ---------------------------------------------------------------------------

.pattern_score_ledger <- function(ledger) {
  if (identical(ledger$type, "categorical")) {
    obs <- ledger$truth
    P <- ledger$prediction
    P <- P[, levels(obs), drop = FALSE]
    pred <- factor(colnames(P)[max.col(P, ties.method = "first")], levels = levels(obs))
    base <- if (nlevels(obs) == 2L) binary_perf(obs, pred, P) else multiclass_perf(obs, pred, P)
    p_true <- P[cbind(seq_along(obs), as.integer(obs))]
    c(Accuracy = unname(base["Accuracy"]), AUC = unname(base["AUC"]),
      logloss = mean(-log(pmax(p_true, 1e-12))))
  } else {
    Y <- as.matrix(ledger$truth); Yhat <- as.matrix(ledger$prediction)
    sse <- colSums((Y - Yhat)^2)
    # Predictive R^2 against the training-mean baseline (the same baseline the
    # rank-selection loss uses), not the test mean.
    base <- ledger$baseline
    sst <- if (is.null(base)) {
      colSums(sweep(Y, 2L, colMeans(Y), "-")^2)
    } else {
      colSums((Y - base)^2)
    }
    r2 <- 1 - sse / pmax(sst, .Machine$double.eps)
    cors <- vapply(seq_len(ncol(Y)), function(j) {
      suppressWarnings(stats::cor(Y[, j], Yhat[, j]))
    }, numeric(1))
    c(R2 = mean(r2), RMSE = mean(sqrt(sse / nrow(Y))), cor = mean(cors, na.rm = TRUE))
  }
}

#' Pool a fold-resolved ledger into one record per observation.
#'
#' Cross-validation schemes that test a row more than once (twofold with
#' \code{nreps}, sequential, bootstrap) produce several predictions for the
#' same observation. \code{wrap_result()} -- the path every other rMVPA model
#' takes -- averages those repeats and returns one row per tested observation
#' in sorted order. This does the same for a \code{pattern_ledger}: class
#' probabilities are averaged (equivalently summed and renormalized) and
#' continuous predictions are divided by their repeat count, so a repeatedly
#' tested observation does not get extra weight in the reported metrics.
#'
#' Pooling is a no-op for schemes that test each row once, which includes
#' blocked and k-fold cross-validation and the external-test path.
#' @keywords internal
#' @noRd
.pattern_pool_ledger <- function(ledger) {
  obs <- ledger$observation
  uo <- sort(unique(obs))
  if (identical(obs, uo)) {
    # already one sorted record per observation: pooling is the identity
    ledger$pooled <- TRUE
    ledger$n_repeats <- rep(1L, length(obs))
    return(ledger)
  }
  idx <- match(obs, uo)
  counts <- as.integer(tabulate(idx, nbins = length(uo)))

  pool <- function(M) {
    M <- as.matrix(M)
    out <- rowsum(M, group = idx, reorder = TRUE) / counts
    dimnames(out) <- list(NULL, colnames(M))
    out
  }
  pred <- pool(ledger$prediction)
  if (identical(ledger$type, "categorical")) {
    pred <- pred / rowSums(pred)          # rows already sum to 1; guards drift
  }
  base <- if (is.null(ledger$baseline)) NULL else pool(ledger$baseline)

  first <- match(uo, obs)
  truth <- if (is.matrix(ledger$truth)) {
    ledger$truth[first, , drop = FALSE]
  } else {
    ledger$truth[first]
  }
  # repeats of one observation must agree about the truth
  chk <- if (is.matrix(ledger$truth)) {
    all(abs(ledger$truth - truth[idx, , drop = FALSE]) < 1e-9)
  } else {
    all(as.character(ledger$truth) == as.character(truth)[idx])
  }
  if (!isTRUE(chk)) {
    stop("pattern_model: repeated predictions for one observation disagree about its target.",
         call. = FALSE)
  }

  # `fold` is dropped: a pooled record can come from several folds. The
  # fold-resolved ledger keeps that information.
  structure(
    list(fold = rep(NA_integer_, length(uo)), observation = uo, truth = truth,
         prediction = pred, baseline = base, partition = ledger$partition,
         type = ledger$type, pooled = TRUE, n_repeats = counts),
    class = c("pattern_ledger", "list")
  )
}

#' @export
print.pattern_ledger <- function(x, ...) {
  if (isTRUE(x$pooled) && !is.null(x$n_repeats) && any(x$n_repeats > 1L)) {
    cat(sprintf("pattern_ledger [%s, pooled]: %d observations, %.2f predictions each on average (%s targets)\n",
                x$partition, length(x$observation), mean(x$n_repeats), x$type))
  } else if (isTRUE(x$pooled)) {
    cat(sprintf("pattern_ledger [%s, pooled]: %d observations, tested once each (%s targets)\n",
                x$partition, length(x$observation), x$type))
  } else {
    cat(sprintf("pattern_ledger [%s]: %d predictions over %d folds (%s targets)\n",
                x$partition, length(x$observation), length(unique(x$fold)), x$type))
  }
  invisible(x)
}

# ---------------------------------------------------------------------------
# Plugin contract
# ---------------------------------------------------------------------------

#' @rdname output_schema
#' @export
output_schema.pattern_model <- function(model) {
  if (identical(model$target_type, "categorical")) {
    list(Accuracy = "scalar", AUC = "scalar", logloss = "scalar", rank_mean = "scalar")
  } else {
    list(R2 = "scalar", RMSE = "scalar", cor = "scalar", rank_mean = "scalar")
  }
}

#' @rdname fit_roi
#' @export
fit_roi.pattern_model <- function(model, roi_data, context, ...) {
  X <- roi_data$train_data
  X_test <- if (isTRUE(model$has_test_set) && !is.null(roi_data$test_data)) roi_data$test_data else NULL
  keep <- isTRUE(model$keep_fold_fits)

  out <- tryCatch(
    .pattern_evaluate_domain(model, X, X_test = X_test, keep_fits = keep, refit = isTRUE(model$refit)),
    error = function(e) e
  )
  if (inherits(out, "error")) {
    return(roi_result(metrics = NULL, indices = roi_data$indices, id = context$id,
                      error = TRUE, error_message = conditionMessage(out)))
  }

  schema_names <- names(output_schema(model))
  metrics <- out$metrics[schema_names]
  names(metrics) <- schema_names
  metrics[!is.finite(metrics)] <- NA_real_

  predictor <- structure(
    list(ledger = out$ledger, pooled_ledger = out$pooled, ranks = out$ranks,
         fold_fits = out$fold_fits, refit = out$refit,
         indices = roi_data$indices, id = context$id),
    class = c("pattern_roi_fits", "list")
  )

  # Categorical targets ride in a classification_result so the regional
  # prediction table, pooling, and $fits all work through the standard path.
  # Continuous multi-response targets keep the ledger under $predictor only.
  # The pooled ledger has one row per tested observation in sorted order, the
  # same shape wrap_result() gives every other model, so the regional
  # prediction table has unique .rownum values.
  led <- out$pooled
  result <- if (identical(model$target_type, "categorical")) {
    P <- led$prediction[, levels(led$truth), drop = FALSE]
    pred <- factor(colnames(P)[max.col(P, ties.method = "first")], levels = levels(led$truth))
    test_design <- if (identical(led$partition, "external")) model$design$test_design else model$design$train_design
    classification_result(led$truth, pred, P, testind = led$observation,
                          test_design = test_design, predictor = predictor)
  } else {
    list(predictor = predictor)
  }
  roi_result(metrics = metrics, indices = roi_data$indices, id = context$id, result = result)
}

# ---------------------------------------------------------------------------
# Global (whole-domain) analysis
# ---------------------------------------------------------------------------

# Feature matrix of the test data, mirroring get_feature_matrix() for train.
.pattern_test_matrix <- function(dataset) {
  if (is.null(dataset$test_data)) return(NULL)
  if (inherits(dataset, "mvpa_clustered_dataset")) return(as.matrix(dataset$test_data@ts))
  idx <- which(dataset$mask > 0)
  if (inherits(dataset, "mvpa_multibasis_image_dataset")) {
    return(do.call(cbind, lapply(dataset$test_data, function(v) neuroim2::series(v, idx))))
  }
  neuroim2::series(dataset$test_data, idx)
}

#' @rdname run_global
#' @param return_fits Retain the per-fold \code{pattern_fit} objects.
#' @param refit Also fit the model on all training rows (descriptive fit).
#' @param preflight Preflight validation level (\code{"warn"}, \code{"error"},
#'   or \code{"off"}); the pattern model runs a lightweight specification
#'   check.
#' @export
run_global.pattern_model <- function(model_spec, return_fits = isTRUE(model_spec$keep_fold_fits),
                                     refit = model_spec$refit,
                                     preflight = c("warn", "error", "off"), ...) {
  preflight <- match.arg(preflight)
  preflight_result <- .apply_analysis_preflight(model_spec, preflight, context = "run_global")
  dataset <- model_spec$dataset
  X <- get_feature_matrix(dataset)
  X_test <- if (isTRUE(model_spec$has_test_set)) .pattern_test_matrix(dataset) else NULL
  feature_ids <- feature_ids_for_dataset(dataset, ncol(X))

  out <- .pattern_evaluate_domain(model_spec, X, X_test = X_test,
                                  keep_fits = isTRUE(return_fits), refit = isTRUE(refit))

  structure(
    list(
      performance_table = tibble::as_tibble(as.list(out$metrics)),
      ledger = out$pooled,
      fold_ledger = out$ledger,
      ranks = out$ranks,
      fold_fits = out$fold_fits,
      refit = out$refit,
      feature_ids = feature_ids,
      n_features = ncol(X),
      model_spec = model_spec,
      preflight = preflight_result
    ),
    class = c("pattern_global_result", "list")
  )
}

#' @export
print.pattern_global_result <- function(x, ...) {
  cat("pattern_global_result\n")
  cat(sprintf("  domain: %d features, %d observations, %s targets\n",
              x$n_features, length(x$model_spec$targets_train$observation_ids),
              x$model_spec$target_type))
  cat(sprintf("  evaluation: %s over %d fold(s); rank selected per fold: %s\n",
              x$ledger$partition, length(x$ranks), paste(x$ranks, collapse = " ")))
  if (!is.null(x$ledger$n_repeats) && any(x$ledger$n_repeats > 1L)) {
    cat(sprintf("  repeated testing: %.2f predictions per observation, averaged\n",
                mean(x$ledger$n_repeats)))
  }
  cat("  performance:\n")
  print(x$performance_table)
  cat(sprintf("  fold fits retained: %s; refit: %s\n",
              if (is.null(x$fold_fits)) "no" else "yes",
              if (is.null(x$refit)) "no" else "yes"))
  invisible(x)
}

#' @rdname performance-methods
#' @export
performance.pattern_global_result <- function(x, ...) {
  x$performance_table
}
