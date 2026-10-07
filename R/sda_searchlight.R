# Exact searchlight engine for sda_notune (the native sda fit, R/sda_native.R).
#
# The sda fit splits into per-column statistics (class means, centring,
# moments, standardisation, lambda.var ingredients) and cross-column steps
# (shrinkage intensities, the Gram matrix and its Woodbury solve). Per fold,
# the per-column statistics are computed once for every voxel in the mask.
# Each sphere then runs only the cross-column steps, on its columns in the
# general path's voxel order. Per-column arithmetic does not depend on the
# other columns, so every sphere's fit and posteriors are bit-identical to
# the general per-sphere path, without its per-sphere bookkeeping or repeated
# column work.
#
# Data outside the regime (missing values, identical voxel columns, a class
# absent from a training fold) raise "rmvpa_engine_ineligible", and the run
# falls back to the general path. A sphere whose fit would be delegated to
# sda::sda() is fitted that way, as the general path does.

#' @keywords internal
#' @noRd
.is_sda_fast_path <- function(model_spec, method) {
  if (!inherits(model_spec, "mvpa_model")) return(FALSE)
  if (!identical(method, "standard")) return(FALSE)
  if (!identical(model_spec$model$label, "sda_notune")) return(FALSE)
  if (!isTRUE(has_crossval(model_spec)) || isTRUE(has_test_set(model_spec))) return(FALSE)
  ds <- model_spec$dataset
  if (!inherits(ds, "mvpa_image_dataset") ||
      inherits(ds, "mvpa_multibasis_image_dataset")) return(FALSE)
  if (!is.null(model_spec$feature_selector)) return(FALSE)
  y <- y_train(model_spec)
  if (!is.factor(y) || nlevels(y) < 2L) return(FALSE)
  grid <- model_spec$tune_grid
  if (!is.null(grid) && (!is.data.frame(grid) || nrow(grid) != 1L)) return(FALSE)
  if (!isTRUE(model_spec$compute_performance)) return(FALSE)
  perf <- model_spec$performance
  kind <- attr(perf, "rmvpa_perf_kind", exact = TRUE)
  if (is.null(kind) || !(kind %in% c("multiclass", "binary"))) return(FALSE)
  if (!is.null(attr(perf, "rmvpa_split_list", exact = TRUE))) return(FALSE)
  if (identical(kind, "binary") != (nlevels(y) == 2L)) return(FALSE)
  TRUE
}

#' @keywords internal
#' @noRd
run_searchlight_sda_fast <- function(model_spec, radius, verbose = FALSE, ...) {
  ds <- model_spec$dataset
  y_all <- y_train(model_spec)
  classes <- levels(y_all)
  K <- length(classes)
  kind <- attr(model_spec$performance, "rmvpa_perf_kind", exact = TRUE)
  class_metrics <- isTRUE(attr(model_spec$performance, "rmvpa_class_metrics", exact = TRUE))

  nb <- .aggregate_neighbourhoods(ds, radius)
  n_centres <- length(nb$centers)
  if (n_centres == 0L) return(empty_searchlight_result(ds))

  x_all <- as.matrix(neuroim2::series(ds$train_data, nb$mask_indices))
  if (nrow(x_all) != length(y_all)) {
    stop("sda_fast: mismatch between train rows and y_train length.")
  }
  if (!all(is.finite(x_all))) {
    stop(.aggregate_ineligible("data contain missing or non-finite values"))
  }

  folds <- generate_folds(model_spec$crossval, tibble::tibble(.row = seq_len(nrow(x_all))), y_all)
  fold_list <- lapply(seq_len(nrow(folds)), function(i) {
    tr <- as.integer(.extract_sample_indices(folds$train[[i]]))
    te <- as.integer(.extract_sample_indices(folds$test[[i]]))
    ytr <- factor(y_all[tr], levels = classes)
    if (any(table(ytr) == 0L)) {
      stop(.aggregate_ineligible("a class is absent from a training fold"))
    }
    if (length(tr) < 3L) {
      stop(.aggregate_ineligible("fewer than three training observations"))
    }
    valid <- nonzeroVarianceColumns2(x_all[tr, , drop = FALSE])
    if (anyDuplicated(t(x_all[tr, valid, drop = FALSE])) > 0L) {
      stop(.aggregate_ineligible("identical voxel columns in a training fold"))
    }
    list(train = tr, test = te, valid = valid, ytr = ytr,
         stats = .sda_column_stats(x_all[tr, , drop = FALSE], ytr))
  })
  testind <- sort(unique(unlist(lapply(fold_list, `[[`, "test"))))
  observed <- y_all[testind]
  n_obs <- length(testind)

  # Each sphere's columns in the general path's voxel order.
  sphere_cols <- .aggregate_generic_order(ds, radius, nb$centers, nb$mask_indices)

  pooled <- matrix(0, n_obs * n_centres, K)
  ok <- rep(TRUE, n_centres)
  for (f in fold_list) {
    x_test <- x_all[f$test, , drop = FALSE]
    dest_rows <- match(f$test, testind)
    for (c in seq_len(n_centres)) {
      if (!ok[c]) next
      cols <- sphere_cols[[c]]
      cols <- cols[f$valid[cols]]
      if (length(cols) < 2L) { ok[c] <- FALSE; next }
      fit <- .sda_fit_columns(f$stats, cols)
      pr <- if (is.null(fit)) {
        nbm <- MVPAModels$sda_notune
        dfit <- nbm$fit(x_all[f$train, cols, drop = FALSE], f$ytr, NULL, NULL, classes, NULL, NULL, TRUE)
        nbm$prob(dfit, x_test[, cols, drop = FALSE])
      } else {
        .sda_native_posterior(fit, x_test[, cols, drop = FALSE])
      }
      dest <- (c - 1L) * n_obs + dest_rows
      pooled[dest, ] <- pooled[dest, , drop = FALSE] + pr
    }
  }

  perf <- .engine_pooled_metrics(pooled, observed, classes, kind, class_metrics, ok)
  good <- !is.na(perf[, "Accuracy"])
  if (!any(good)) return(empty_searchlight_result(ds))
  out <- wrap_out(perf[good, , drop = FALSE], ds, ids = nb$centers[good])
  attr(out, "bad_results") <- tibble::tibble()
  out
}
