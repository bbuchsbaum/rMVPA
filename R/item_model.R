#' ITEM decoding model for ROI/searchlight analysis
#'
#' Integrates ITEM-style trial-wise decoding into the `fit_roi` architecture by
#' delegating trial-level estimation and covariance-aware decoding to `fmrilss`.
#'
#' Conceptually, ITEM differs from `hrfdecoder_model()` in where the decoder is
#' fit:
#' - ITEM (`item_model`) first estimates trial-wise responses (`Gamma`) via LS-A
#'   and then performs covariance-aware decoding on those trial estimates.
#' - `hrfdecoder_model` fits a continuous-time decoder directly on TR-level data
#'   and aggregates TR predictions to events afterward.
#'
#' Use ITEM when you want an explicit trial-estimation stage and direct control
#' over trial covariance handling (`U`), especially for trial-level diagnostics.
#'
#' @param dataset An `mvpa_dataset`.
#' @param design An `item_design` object.
#' @param mode Decoding mode: `"classification"` or `"regression"`.
#' @param metric Optional ITEM metric.
#'   Classification: `"accuracy"`, `"balanced_accuracy"`.
#'   Regression: `"correlation"`, `"rmse"`.
#' @param ridge_u Non-negative ridge used when computing `U`.
#' @param ridge_w Non-negative ridge used when fitting ITEM weights per fold.
#' @param lsa_method LS-A backend for `fmrilss::lsa()` (`"r"` or `"cpp"`).
#' @param solver Solver preference for ITEM linear solves (`"chol"`, `"svd"`, `"pinv"`).
#' @param u_storage Store trial covariance as full matrix (`"matrix"`) or run blocks (`"by_run"`).
#' @param class_levels Optional class order for classification.
#' @param check_hash Validate trial hash before CV when available.
#' @param return_predictions Keep trial-level prediction tables in ROI results.
#' @param compute_performance Reserved for API compatibility. ITEM always
#'   computes scalar ROI metrics for mapping.
#' @param ... Additional fields stored on the model spec.
#'
#' @return A model spec of class `item_model` compatible with
#'   `run_regional()` and `run_searchlight()`.
#' @export
item_model <- function(dataset,
                       design,
                       mode = c("classification", "regression"),
                       metric = NULL,
                       ridge_u = 0,
                       ridge_w = 1e-4,
                       lsa_method = c("r", "cpp"),
                       solver = c("chol", "svd", "pinv"),
                       u_storage = c("matrix", "by_run"),
                       class_levels = NULL,
                       check_hash = FALSE,
                       return_predictions = TRUE,
                       compute_performance = TRUE,
                       ...) {
  .item_require_fmrilss("item_model")

  mode <- match.arg(mode)
  lsa_method <- match.arg(lsa_method)
  solver <- match.arg(solver)
  u_storage <- match.arg(u_storage)

  assertthat::assert_that(inherits(dataset, "mvpa_dataset"))
  assertthat::assert_that(inherits(design, "item_design"), msg = "design must inherit class 'item_design'.")

  if (!isTRUE(compute_performance)) {
    warning("item_model always computes scalar mapping metrics; forcing compute_performance=TRUE.", call. = FALSE)
  }
  compute_performance <- TRUE

  if (!is.numeric(ridge_u) || length(ridge_u) != 1L || is.na(ridge_u) || ridge_u < 0) {
    stop("ridge_u must be a single non-negative number.", call. = FALSE)
  }
  if (!is.numeric(ridge_w) || length(ridge_w) != 1L || is.na(ridge_w) || ridge_w < 0) {
    stop("ridge_w must be a single non-negative number.", call. = FALSE)
  }

  n_time <- nobs(dataset)
  if (nrow(design$X_t) != n_time) {
    stop(
      sprintf(
        "nrow(design$X_t) (%d) must match dataset observations (%d).",
        nrow(design$X_t),
        n_time
      ),
      call. = FALSE
    )
  }

  bundle <- design$item_bundle

  if (identical(mode, "classification")) {
    if (!is.null(class_levels)) {
      class_levels <- as.character(class_levels)
      if (length(class_levels) != ncol(bundle$T_target)) {
        stop(
          sprintf("class_levels must have length %d.", ncol(bundle$T_target)),
          call. = FALSE
        )
      }
      colnames(bundle$T_target) <- class_levels
      bundle$meta$class_levels <- class_levels
    } else if (!is.null(bundle$meta$class_levels)) {
      class_levels <- as.character(bundle$meta$class_levels)
    }
  } else {
    class_levels <- NULL
  }

  u_obj <- fmrilss::item_compute_u(
    X_t = bundle$X_t,
    V = design$V,
    v_type = design$v_type,
    ridge = ridge_u,
    method = solver,
    run_id = bundle$run_id,
    output = u_storage
  )

  if (identical(u_storage, "matrix")) {
    bundle$U <- u_obj
    bundle$U_by_run <- NULL
  } else {
    bundle$U <- NULL
    bundle$U_by_run <- u_obj
  }

  bundle$diagnostics <- utils::modifyList(
    bundle$diagnostics,
    list(
      u = attr(u_obj, "item_diagnostics")
    )
  )

  create_model_spec(
    "item_model",
    dataset = dataset,
    design = design,
    mode = mode,
    metric = metric,
    ridge_u = ridge_u,
    ridge_w = ridge_w,
    lsa_method = lsa_method,
    solver = solver,
    u_storage = u_storage,
    class_levels = class_levels,
    check_hash = isTRUE(check_hash),
    return_predictions = isTRUE(return_predictions),
    compute_performance = compute_performance,
    has_test_set = FALSE,
    item_bundle = bundle,
    Z = design$Z,
    Nuisance = design$Nuisance,
    ...
  )
}

#' @rdname fit_roi
#' @method fit_roi item_model
#' @export
fit_roi.item_model <- function(model, roi_data, context, ...) {
  Y <- roi_data$train_data
  if (is.null(Y) && !is.null(roi_data$train_roi)) {
    Y <- as.matrix(neuroim2::values(roi_data$train_roi))
  }
  Y <- as.matrix(Y)

  if (!is.numeric(Y)) {
    return(roi_result(
      metrics = NULL,
      indices = roi_data$indices,
      id = context$id,
      error = TRUE,
      error_message = "item_model: ROI data must be numeric."
    ))
  }

  bundle <- model$item_bundle
  if (nrow(Y) != nrow(bundle$X_t)) {
    return(roi_result(
      metrics = NULL,
      indices = roi_data$indices,
      id = context$id,
      error = TRUE,
      error_message = sprintf(
        "item_model: ROI rows (%d) must match nrow(X_t) (%d).",
        nrow(Y),
        nrow(bundle$X_t)
      )
    ))
  }

  gamma_hat <- tryCatch(
    fmrilss::lsa(
      Y = Y,
      X = bundle$X_t,
      Z = model$Z,
      Nuisance = model$Nuisance,
      method = model$lsa_method
    ),
    error = function(e) e
  )

  if (inherits(gamma_hat, "error")) {
    return(roi_result(
      metrics = NULL,
      indices = roi_data$indices,
      id = context$id,
      error = TRUE,
      error_message = paste("item_model: lsa failed:", conditionMessage(gamma_hat))
    ))
  }

  u_obj <- if (!is.null(bundle$U_by_run)) bundle$U_by_run else bundle$U

  cv_res <- tryCatch(
    fmrilss::item_cv(
      Gamma = gamma_hat,
      T_target = bundle$T_target,
      U = u_obj,
      run_id = bundle$run_id,
      mode = model$mode,
      metric = model$metric,
      ridge = model$ridge_w,
      method = model$solver,
      class_levels = model$class_levels,
      trial_id = bundle$trial_id,
      trial_hash = bundle$trial_hash,
      check_hash = isTRUE(model$check_hash)
    ),
    error = function(e) e
  )

  if (inherits(cv_res, "error")) {
    return(roi_result(
      metrics = NULL,
      indices = roi_data$indices,
      id = context$id,
      error = TRUE,
      error_message = paste("item_model: item_cv failed:", conditionMessage(cv_res))
    ))
  }

  metric_mean <- as.numeric(cv_res$aggregate$mean)
  metric_sd <- as.numeric(cv_res$aggregate$sd)
  n_folds <- as.numeric(cv_res$aggregate$n_folds)

  metrics <- c(
    item_score_mean = metric_mean,
    item_score_sd = metric_sd,
    item_n_folds = n_folds
  )

  if (!isTRUE(model$return_predictions) && !is.null(cv_res$predictions)) {
    cv_res$predictions <- NULL
  }

  roi_result(
    metrics = metrics,
    indices = roi_data$indices,
    id = context$id,
    result = cv_res
  )
}

#' @rdname output_schema
#' @method output_schema item_model
#' @export
output_schema.item_model <- function(model) {
  list(
    item_score_mean = "scalar",
    item_score_sd = "scalar",
    item_n_folds = "scalar"
  )
}

#' @keywords internal
#' @export
process_roi.item_model <- function(mod_spec, roi, rnum, center_global_id = NA, ...) {
  NextMethod()
}
