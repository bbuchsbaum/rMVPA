# ---------- get_feature_matrix methods ----------

#' @rdname get_feature_matrix
#' @export
get_feature_matrix.mvpa_clustered_dataset <- function(dataset, ...) {
  as.matrix(dataset$train_data@ts)
}

#' @rdname get_feature_matrix
#' @export
get_feature_matrix.mvpa_image_dataset <- function(dataset, ...) {
  idx <- which(dataset$mask > 0)
  neuroim2::series(dataset$train_data, idx)
}

#' @rdname get_feature_matrix
#' @export
get_feature_matrix.mvpa_surface_dataset <- function(dataset, ...) {
  idx <- which(dataset$mask > 0)
  neuroim2::series(dataset$train_data, idx)
}

#' @rdname get_feature_matrix
#' @export
get_feature_matrix.matrix <- function(dataset, ...) {
  dataset
}


# ---------- feature_ids_for_dataset ----------

#' Get feature IDs for building output maps
#' @keywords internal
#' @noRd
feature_ids_for_dataset <- function(dataset, P) {
  if (inherits(dataset, "mvpa_clustered_dataset")) {
    get_center_ids(dataset)
  } else if (inherits(dataset, "mvpa_image_dataset")) {
    which(dataset$mask > 0)
  } else if (inherits(dataset, "mvpa_surface_dataset")) {
    which(dataset$mask > 0)
  } else {
    seq_len(P)
  }
}


# ---------- cv_run_global ----------

#' Run Cross-Validated Global Analysis
#'
#' Shared CV loop used by \code{run_global} and \code{region_importance}.
#' Trains a model on the full feature matrix in each fold and collects
#' predictions, optional weight matrices, optional per-fold fits, and
#' optional per-fold importance vectors.
#'
#' @param mspec An \code{mvpa_model} specification (may have dataset stripped).
#' @param X A T x P feature matrix.
#' @param y Response vector (factor or numeric).
#' @param crossval A cross-validation specification.
#' @param feature_ids Integer vector of feature IDs for \code{train_model}.
#' @param return_fits Logical; if TRUE, store per-fold model fits.
#' @param extract_weights_fn Logical; if TRUE, extract weight matrices per fold.
#' @param compute_importance Logical; if TRUE, compute per-fold importance via
#'   \code{model_importance}. Each fold uses its own training data.
#' @param summary_fun Optional summary function for Haufe importance (passed
#'   to \code{model_importance} methods for linear models).
#' @return A list with components \code{result_table} (tibble), \code{fold_fits}
#'   (list, possibly NULL entries), \code{fold_weights} (list, possibly NULL
#'   entries), and \code{fold_importance} (list, possibly NULL entries).
#' @keywords internal
#' @noRd
cv_run_global <- function(mspec, X, y, crossval, feature_ids,
                          return_fits = FALSE, extract_weights_fn = TRUE,
                          compute_importance = FALSE, summary_fun = NULL) {
  P <- ncol(X)
  X_df <- as.data.frame(X)
  cv_samples <- crossval_samples(crossval, X_df, y)
  nfolds <- nrow(cv_samples)

  fold_fits       <- vector("list", nfolds)
  fold_weights    <- vector("list", nfolds)
  fold_importance <- vector("list", nfolds)
  result_rows     <- vector("list", nfolds)

  for (k in seq_len(nfolds)) {
    train_idx <- as.integer(cv_samples$train[[k]])
    test_idx  <- as.integer(cv_samples$test[[k]])
    y_train_k <- cv_samples$ytrain[[k]]
    y_test_k  <- cv_samples$ytest[[k]]

    X_train <- X[train_idx, , drop = FALSE]
    X_test  <- X[test_idx, , drop = FALSE]

    mfit <- tryCatch(
      train_model(mspec,
                  as.data.frame(X_train),
                  y_train_k,
                  indices = feature_ids),
      error = function(e) {
        futile.logger::flog.warn("cv_run_global fold %d: training error: %s", k, e$message)
        NULL
      }
    )

    if (is.null(mfit)) {
      if (is.factor(y)) {
        result_rows[[k]] <- tibble::tibble(
          class = list(NULL), probs = list(NULL),
          y_true = list(y_test_k), test_ind = list(test_idx),
          fit = list(NULL), error = TRUE, error_message = "training error"
        )
      } else {
        result_rows[[k]] <- tibble::tibble(
          preds = list(NULL),
          y_true = list(y_test_k), test_ind = list(test_idx),
          fit = list(NULL), error = TRUE, error_message = "training error"
        )
      }
      next
    }

    if (return_fits) {
      fold_fits[[k]] <- mfit
    }

    pred <- predict(mfit, X_test)

    if (extract_weights_fn) {
      raw_fit <- mfit$fit
      W_fold <- tryCatch(extract_weights(raw_fit), error = function(e) NULL)

      if (!is.null(W_fold)) {
        W_full <- matrix(0, nrow = P, ncol = ncol(W_fold))
        fm <- mfit$feature_mask
        if (!is.null(fm) && is.logical(fm)) {
          W_full[fm, ] <- W_fold
        } else {
          W_full <- W_fold
        }
        fold_weights[[k]] <- W_full
      }
    }

    # Per-fold importance via model_importance generic
    if (compute_importance) {
      raw_fit <- mfit$fit
      fm <- mfit$feature_mask

      # Build the training data as seen by the model (after feature masking)
      X_train_masked <- if (!is.null(fm) && is.logical(fm)) {
        X_train[, fm, drop = FALSE]
      } else {
        X_train
      }

      imp_args <- list(object = raw_fit, X_train = X_train_masked)
      if (!is.null(summary_fun)) {
        imp_args$summary_fun <- summary_fun
      }

      imp_k <- tryCatch(
        do.call(model_importance, imp_args),
        error = function(e) {
          futile.logger::flog.warn(
            "cv_run_global fold %d: model_importance failed: %s", k, e$message)
          NULL
        }
      )

      # Zero-pad to full P-length vector if feature mask was applied
      if (!is.null(imp_k) && !is.null(fm) && is.logical(fm)) {
        imp_full <- numeric(P)
        imp_full[fm] <- imp_k
        fold_importance[[k]] <- imp_full
      } else {
        fold_importance[[k]] <- imp_k
      }
    }

    if (is.factor(y)) {
      plist <- lapply(pred, list)
      plist$y_true <- list(y_test_k)
      plist$test_ind <- list(test_idx)
      plist$fit <- if (return_fits) list(mfit) else list(NULL)
      plist$error <- FALSE
      plist$error_message <- "~"
      result_rows[[k]] <- tibble::as_tibble(plist, .name_repair = .name_repair)
    } else {
      result_rows[[k]] <- tibble::tibble(
        preds = list(pred$preds),
        y_true = list(y_test_k), test_ind = list(test_idx),
        fit = if (return_fits) list(mfit) else list(NULL),
        error = FALSE, error_message = "~"
      )
    }
  }

  result_table <- dplyr::bind_rows(result_rows)

  list(
    result_table    = result_table,
    fold_fits       = fold_fits,
    fold_weights    = fold_weights,
    fold_importance = fold_importance
  )
}


# ---------- run_global ----------

#' @rdname run_global
#' @param model_spec An \code{mvpa_model} specification.
#' @param X Optional pre-computed T x P feature matrix. If NULL, extracted
#'   from \code{model_spec$dataset} via \code{get_feature_matrix}.
#' @param summary_fun Function to summarize activation pattern matrix rows
#'   into a scalar importance per feature. Default: L2 norm.
#' @param return_fits Logical; if TRUE, store per-fold model fits.
#' @param ... Additional arguments (currently unused).
#'
#' @return A \code{global_mvpa_result} object.
#' @export
run_global.mvpa_model <- function(model_spec, X = NULL, summary_fun = NULL,
                                   return_fits = FALSE, ...) {

  dataset <- model_spec$dataset
  design  <- model_spec$design

  # --- 1. Feature matrix ---
  if (is.null(X)) {
    X <- get_feature_matrix(dataset)
  }
  stopifnot(is.matrix(X))

  P <- ncol(X)
  N <- nrow(X)
  y <- y_train(design)
  stopifnot(length(y) == N)

  feature_ids <- feature_ids_for_dataset(dataset, P)

  # --- 2. Run CV via shared helper (with per-fold importance) ---
  cv_out <- cv_run_global(model_spec, X, y, model_spec$crossval, feature_ids,
                          return_fits = return_fits, extract_weights_fn = TRUE,
                          compute_importance = TRUE, summary_fun = summary_fun)

  # --- 3. Merge results ---
  cres <- tryCatch(
    wrap_result(cv_out$result_table, design),
    error = function(e) {
      futile.logger::flog.warn("run_global: wrap_result failed: %s", e$message)
      NULL
    }
  )

  # --- 4. Compute performance ---
  perf_table <- if (!is.null(cres)) {
    perf <- compute_performance(model_spec, cres)
    tibble::as_tibble(as.list(perf))
  } else {
    tibble::tibble()
  }

  # --- 5. Per-fold importance (Haufe-then-average for linear, native for others) ---
  valid_imp <- Filter(Negate(is.null), cv_out$fold_importance)
  importance_vec <- NULL
  importance_map <- NULL

  if (length(valid_imp) > 0) {
    importance_vec <- Reduce("+", valid_imp) / length(valid_imp)

    importance_map <- tryCatch(
      build_output_map(dataset, importance_vec, feature_ids),
      error = function(e) {
        futile.logger::flog.warn("run_global: build_output_map failed: %s", e$message)
        NULL
      }
    )
  }

  # --- 5b. Backward-compat: averaged weights + activation patterns (linear only) ---
  valid_weights <- Filter(Negate(is.null), cv_out$fold_weights)
  W_avg <- NULL
  A_avg <- NULL

  if (length(valid_weights) > 0) {
    max_D <- max(sapply(valid_weights, ncol))
    valid_weights <- lapply(valid_weights, function(w) {
      if (ncol(w) < max_D) {
        cbind(w, matrix(0, nrow = nrow(w), ncol = max_D - ncol(w)))
      } else {
        w
      }
    })
    W_avg <- Reduce("+", valid_weights) / length(valid_weights)

    # Compute averaged activation pattern matrix for linear models
    Sigma_x <- cov(X)
    haufe_args <- list(W = W_avg, Sigma_x = Sigma_x)
    if (!is.null(summary_fun)) {
      haufe_args$summary_fun <- summary_fun
    }
    haufe_result <- tryCatch(
      do.call(haufe_importance, haufe_args),
      error = function(e) NULL
    )
    if (!is.null(haufe_result)) {
      A_avg <- haufe_result$A
    }
  }

  # --- 6. Build result ---
  global_mvpa_result(
    performance_table   = perf_table,
    result              = cres,
    importance_map      = importance_map,
    importance_vector   = importance_vec,
    activation_patterns = A_avg,
    raw_weights         = W_avg,
    fold_fits           = if (return_fits) cv_out$fold_fits else NULL,
    model_spec          = model_spec
  )
}


# ---------- global_mvpa_result constructor ----------

#' Construct a Global MVPA Result
#'
#' @param performance_table Tibble of cross-validated performance metrics.
#' @param result The merged classification/regression result object.
#' @param importance_map Spatial object (NeuroVol/NeuroSurface) of per-feature importance.
#' @param importance_vector Numeric vector of per-feature importance.
#' @param activation_patterns The P x D activation pattern matrix A.
#' @param raw_weights The averaged P x D weight matrix W.
#' @param fold_fits Optional list of per-fold model_fit objects.
#' @param model_spec The input model specification.
#'
#' @return An S3 object of class \code{global_mvpa_result}.
#' @export
global_mvpa_result <- function(performance_table, result,
                                importance_map, importance_vector,
                                activation_patterns, raw_weights,
                                fold_fits, model_spec) {
  structure(
    list(
      performance_table   = performance_table,
      result              = result,
      importance_map      = importance_map,
      importance_vector   = importance_vector,
      activation_patterns = activation_patterns,
      raw_weights         = raw_weights,
      fold_fits           = fold_fits,
      model_spec          = model_spec
    ),
    class = "global_mvpa_result"
  )
}

#' @export
#' @method print global_mvpa_result
print.global_mvpa_result <- function(x, ...) {
  cat("\n  Global MVPA Result\n\n")

  if (nrow(x$performance_table) > 0) {
    cat("  Performance:\n")
    print(x$performance_table, n = 20)
    cat("\n")
  }

  if (!is.null(x$importance_vector)) {
    cat("  Importance vector: length", length(x$importance_vector), "\n")
    cat("    range: [", round(min(x$importance_vector), 4), ",",
        round(max(x$importance_vector), 4), "]\n")
  }

  if (!is.null(x$importance_map)) {
    cat("  Importance map: ", class(x$importance_map)[1], "\n")
  }

  if (!is.null(x$raw_weights)) {
    cat("  Weight matrix: ", nrow(x$raw_weights), "x", ncol(x$raw_weights), "\n")
  }

  if (!is.null(x$fold_fits)) {
    cat("  Fold fits: ", length(x$fold_fits), "\n")
  }

  cat("\n")
}
