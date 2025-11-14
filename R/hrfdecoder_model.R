#' hrfdecoder Model (continuous-time MVPA)
#'
#' Adapter model that trains a continuous-time decoder on TR-level data within each
#' ROI/fold and aggregates to event-level predictions. This model integrates with
#' existing rMVPA searchlight/cross-validation iterators and delegates the core
#' solver to the external 'hrfdecoder' package.
#'
#' The hrfdecoder approach differs from traditional MVPA by operating directly on
#' continuous TR-level data without requiring trial-level beta estimation. It jointly
#' estimates decoder weights (W), TR-level soft labels (P), and HRF parameters via
#' an alternating least squares optimization that incorporates temporal smoothness
#' and HRF conformity penalties.
#'
#' Usage notes:
#' - Requires the 'hrfdecoder' package at runtime; we check via requireNamespace.
#' - Use the explicit `hrfdecoder_design()` constructor (subclasses `mvpa_design`).
#'   This keeps event metadata explicit without adding ad-hoc members to a plain design.
#' - Cross-validation operates at the TR level (typically blocked by run), but
#'   performance metrics are computed on event-level predictions after aggregation.
#'
#' @param dataset An `mvpa_dataset` with continuous TR x V data and mask.
#' @param design An `hrfdecoder_design` (subclass of `mvpa_design`) created by
#'   `hrfdecoder_design(event_model, events, block_var, ...)`.
#' @param lambda_W Decoder ridge penalty (default: 10). Controls regularization of
#'   decoder weights to prevent overfitting.
#' @param lambda_HRF Weight for HRF-conformity of soft labels (default: 1). Higher
#'   values force soft labels closer to the HRF-convolved event design.
#' @param lambda_smooth Temporal smoothness penalty for soft labels (default: 5).
#'   Penalizes rapid changes in soft labels between adjacent TRs.
#' @param basis Optional HRF basis (e.g., fmrihrf::spmg1()). If NULL, uses the hrfdecoder default.
#' @param window Numeric length-2 (start,end) seconds after onset for event aggregation
#'   (default: c(4, 8)). Defines the time window for aggregating TR-level predictions
#'   to event-level.
#' @param nonneg Logical; whether to project soft labels to non-negative (default: TRUE).
#' @param max_iter Integer; maximum ALS iterations (default: 10).
#' @param tol Convergence tolerance (default: 1e-4).
#' @param performance Optional custom performance function that takes a classification_result.
#' @param metrics Character vector of fast metrics to compute on event-level probs
#'   when available (via hrfdecoder::hrf_metrics). Examples: c("auc_ovo","accuracy").
#'   These are added non-invasively to the result row and do not affect the
#'   legacy rMVPA performance vector.
#' @param primary_metric Optional character naming which metric to surface as
#'   the primary summary (stored as `primary_metric` and `primary_value` in ROI
#'   outputs). Defaults to "auc_ovo" if available.
#' @param crossval Optional cross-validation spec. If NULL and `design$block_var` exists, uses blocked CV.
#' @param return_predictions Logical; keep per-fold aggregated event predictions (default: TRUE).
#' @param return_fits Logical; keep per-fold fit objects (default: FALSE).
#'
#' @return An S3 model spec object of class `hrfdecoder_model` compatible with run_searchlight.
#'
#' @examples
#' \dontrun{
#' library(fmridesign)
#' library(fmrihrf)
#' library(hrfdecoder)
#'
#' # 1. Create event table
#' events_df <- data.frame(
#'   onset = seq(10, 290, by = 20),  # Events every 20 seconds
#'   condition = rep(c("A", "B", "C"), length.out = 15),
#'   run = rep(1:3, each = 5)
#' )
#'
#' # 2. Define temporal structure
#' sframe <- sampling_frame(blocklens = c(100, 100, 100), TR = 2)
#'
#' # 3. Build event model
#' evmod <- event_model(
#'   onset ~ hrf(condition, basis = "spmg3"),
#'   data = events_df,
#'   block = ~run,
#'   sampling_frame = sframe
#' )
#'
#' # 4. Create rMVPA dataset
#' mask <- neuroim2::NeuroVol(...)
#' fmri_data <- neuroim2::NeuroVec(...)  # 300 TRs
#' dset <- mvpa_dataset(train_data = fmri_data, mask = mask)
#'
#' # 5. Create hrfdecoder design
#' block_var <- rep(1:3, each = 100)
#' design <- hrfdecoder_design(
#'   event_model = evmod,
#'   events = events_df,
#'   block_var = block_var
#' )
#'
#' # 6. Specify model with custom hyperparameters
#' mspec <- hrfdecoder_model(
#'   dataset = dset,
#'   design = design,
#'   lambda_W = 10,
#'   lambda_HRF = 1,
#'   lambda_smooth = 5,
#'   basis = fmrihrf::spmg1(),
#'   window = c(4, 8),
#'   max_iter = 15
#' )
#'
#' # 7. Run searchlight analysis
#' results <- run_searchlight(mspec, radius = 8, method = "randomized", niter = 4)
#' }
#' @export
hrfdecoder_model <- function(dataset, design,
                             lambda_W = 10,
                             lambda_HRF = 1,
                             lambda_smooth = 5,
                             basis = NULL,                 # use hrfdecoder default if NULL
                             window = c(4, 8),
                             nonneg = TRUE,
                             max_iter = 10,
                             tol = 1e-4,
                             performance = NULL,           # classification_result -> named metrics
                             crossval = NULL,              # default: blocked by run
                             return_predictions = TRUE,
                             return_fits = FALSE,
                             metrics = c("auc_ovo", "accuracy"),
                             primary_metric = "auc_ovo") {

  assertthat::assert_that(inherits(dataset, "mvpa_dataset"))
  assertthat::assert_that(inherits(design,  "hrfdecoder_design"),
                          msg = "design must be created with hrfdecoder_design()")

  # Default to blocked CV when block_var present
  if (is.null(crossval) && !is.null(design$block_var)) {
    # Warn if only a single block/run is present since blocked CV may fail
    if (length(unique(design$block_var)) < 2) {
      warning("hrfdecoder_model: only one run/block detected; consider providing a custom crossval or using multiple runs.")
    }
    crossval <- blocked_cross_validation(design$block_var)
  }

  # Validate aggregation window
  if (!is.numeric(window) || length(window) != 2L || any(!is.finite(window))) {
    warning("hrfdecoder_model: invalid 'window'; using default c(4, 8).")
    window <- c(4, 8)
  } else {
    window <- sort(as.numeric(window))
  }

  # Validate lambda parameters
  if (!is.finite(lambda_W) || lambda_W < 0) {
    warning("hrfdecoder_model: invalid 'lambda_W'; using default 10.")
    lambda_W <- 10
  }
  if (!is.finite(lambda_HRF) || lambda_HRF < 0) {
    warning("hrfdecoder_model: invalid 'lambda_HRF'; using default 1.")
    lambda_HRF <- 1
  }
  if (!is.finite(lambda_smooth) || lambda_smooth < 0) {
    warning("hrfdecoder_model: invalid 'lambda_smooth'; using default 5.")
    lambda_smooth <- 5
  }

  # Optionally validate HRF basis if provided but fmrihrf is unavailable
  if (!is.null(basis) && !requireNamespace("fmrihrf", quietly = TRUE)) {
    warning("hrfdecoder_model: 'fmrihrf' not available; ignoring provided basis and using default HRF.")
    basis <- NULL
  }

  # Choose a performance function: default to multiclass metrics on the aggregated
  # classification_result that merge_results will construct.
  perf_fun <- if (!is.null(performance) && is.function(performance)) {
    get_custom_perf(performance, design$split_groups)
  } else {
    get_multiclass_perf(design$split_groups, class_metrics = TRUE)
  }

  create_model_spec(
    "hrfdecoder_model",
    dataset = dataset,
    design  = design,
    lambda_W = lambda_W,
    lambda_HRF = lambda_HRF,
    lambda_smooth = lambda_smooth,
    basis = basis,
    window = window,
    nonneg = nonneg,
    max_iter = max_iter,
    tol = tol,
    crossval = crossval,
    performance = perf_fun,
    compute_performance = TRUE,
    return_predictions = return_predictions,
    return_fits = return_fits,
    metrics = metrics,
    primary_metric = primary_metric
  )
}


#' Provide a dummy y_train for fold construction (length = #TRs)
#'
#' internal_crossval() uses y_train only to create CV folds; the actual decoder
#' ignores it. We return a simple sequence matching the number of rows (TRs).
#'
#' @param obj hrfdecoder_model specification object
#' @return Integer vector of length equal to number of TRs (observations) in the dataset
#' @export
y_train.hrfdecoder_model <- function(obj) {
  seq_len(nobs(obj$dataset))
}


#' Train per ROI/fold using hrfdecoder
#'
#' This method trains the continuous-time decoder on TR-level data from the
#' training fold. Unlike standard MVPA, it does not use the `y` parameter
#' for training - the actual decoding targets come from the event model and
#' events stored in the design object.
#'
#' @param obj hrfdecoder_model spec
#' @param train_dat Tibble/matrix of TR x features (ROI data for train rows)
#' @param y IGNORED. This is a dummy TR sequence from CV fold construction.
#'   Actual targets come from obj$design$event_model and obj$design$events.
#' @param sl_info Searchlight info list with center ids (optional)
#' @param cv_spec Cross-validation spec (optional)
#' @param indices Global ROI indices
#' @param ... Additional arguments (unused)
#' @return A list with class "hrfdecoder_fit_wrap" containing the fitted model,
#'   searchlight info, and ROI indices
#' @export
train_model.hrfdecoder_model <- function(obj, train_dat, y, sl_info, cv_spec, indices, ...) {
  if (!requireNamespace("hrfdecoder", quietly = TRUE)) {
    stop("hrfdecoder_model: The 'hrfdecoder' package is required at runtime. Please install it.")
  }

  X <- as.matrix(train_dat)

  # Pull extras from the design
  evm <- obj$design$event_model
  events <- obj$design$events
  if (is.null(evm) || is.null(events)) {
    stop("hrfdecoder_model: design must contain $event_model and $events for aggregation.")
  }

  # Basis: allow NULL to defer to hrfdecoder defaults; if provided, pass through
  basis <- obj$basis

  # Optional baseline component if present
  baseline <- NULL
  if (!is.null(obj$design$baseline_model)) baseline <- obj$design$baseline_model
  if (is.null(baseline) && !is.null(evm$baseline)) baseline <- evm$baseline

  fit <- tryCatch({
    hrfdecoder::fit_hrfdecoder(
      Y = X,
      ev_model = evm,
      base_model = baseline,
      hrf = basis,
      lambda_W = obj$lambda_W,
      lambda_HRF = obj$lambda_HRF,
      lambda_smooth = obj$lambda_smooth,
      nonneg = obj$nonneg,
      max_iter = obj$max_iter,
      tol = obj$tol
    )
  }, error = function(e) {
    stop(paste0("fit_hrfdecoder failed: ", e$message))
  })

  structure(list(fit = fit, sl_info = sl_info, indices = indices), class = "hrfdecoder_fit_wrap")
}


#' Format per-fold results: predict TR-level soft labels, aggregate to events
#'
#' This method performs TR-level prediction on the test fold, then aggregates
#' the continuous predictions to event-level using the time window specified
#' in the model. Ground truth labels are extracted from the events table during
#' aggregation.
#'
#' @param obj hrfdecoder_model specification object
#' @param result Output from train_model.hrfdecoder_model containing the fitted model
#' @param error_message Optional error message if training failed
#' @param context List containing test data and other fold information from crossval
#' @param ... Additional arguments (unused)
#' @return A tibble with one row containing class predictions, probabilities,
#'   true labels, test indices, optional fit object, error status, and error message
#' @export
format_result.hrfdecoder_model <- function(obj, result, error_message = NULL, context, ...) {
  if (!is.null(error_message)) {
    return(tibble::tibble(
      class = list(NULL),
      probs = list(NULL),
      y_true = list(context$ytest),
      fit = list(NULL),
      error = TRUE,
      error_message = error_message
    ))
  }

  if (!requireNamespace("hrfdecoder", quietly = TRUE)) {
    return(tibble::tibble(
      class = list(NULL), probs = list(NULL), y_true = list(context$ytest),
      fit = list(NULL), error = TRUE,
      error_message = "Missing required package 'hrfdecoder'"
    ))
  }

  # Predict on this fold's TEST rows
  Xtest <- as.matrix(tibble::as_tibble(context$test, .name_repair = .name_repair))
  Ptest <- tryCatch({
    hrfdecoder::predict_hrfdecoder(
      object = result$fit,
      Y_test = Xtest,
      mode = "tr"
    )
  }, error = function(e) {
    attr(e, "message") <- paste0("predict_hrfdecoder failed: ", e$message)
    e
  })

  if (inherits(Ptest, "error")) {
    return(tibble::tibble(
      class = list(NULL), probs = list(NULL), y_true = list(context$ytest),
      fit = list(NULL), error = TRUE, error_message = attr(Ptest, "message")
    ))
  }

  # Aggregate TR-level soft labels to event-level probabilities
  events <- obj$design$events
  # Subset to event conditions if background column exists
  conditions <- tryCatch(result$fit$conditions, error = function(...) NULL)
  if (!is.null(conditions) && !is.null(colnames(Ptest))) {
    common <- intersect(colnames(Ptest), conditions)
    if (length(common) == length(conditions)) {
      Ptest <- Ptest[, conditions, drop = FALSE]
    }
  }
  agg <- tryCatch({
    hrfdecoder::aggregate_events(
      P = Ptest,
      events = events,
      TR = result$fit$settings$TR,
      conditions = if (!is.null(conditions)) conditions else colnames(Ptest),
      window = obj$window,
      hrf = tryCatch(result$fit$hrf, error = function(...) NULL),
      normalize = TRUE
    )
  }, error = function(e) {
    attr(e, "message") <- paste0("aggregate_events failed: ", e$message)
    e
  })

  if (inherits(agg, "error")) {
    return(tibble::tibble(
      class = list(NULL), probs = list(NULL), y_true = list(context$ytest),
      fit = list(NULL), error = TRUE, error_message = attr(agg, "message")
    ))
  }

  probs <- as.matrix(agg$probs)
  observed <- agg$y_true
  if (is.null(colnames(probs))) {
    stop("hrfdecoder_model: aggregated probabilities must have column names for classes.")
  }

  # Check for events with zero probabilities (empty aggregation windows)
  zero_prob_events <- rowSums(probs) == 0
  if (any(zero_prob_events)) {
    warning(sprintf(
      "hrfdecoder_model: %d event(s) have zero probabilities (aggregation window may fall outside available TRs)",
      sum(zero_prob_events)
    ))
  }

  # Ensure observed is a factor with the same level set ordering
  observed <- factor(as.character(observed), levels = colnames(probs))
  pred_class <- factor(colnames(probs)[max.col(probs, ties.method = "first")], levels = colnames(probs))

  plist <- list(
    class    = list(pred_class),
    probs    = list(probs),
    y_true   = list(observed),
    test_ind = list(as.integer(context$test)),
    fit      = list(if (isTRUE(obj$return_fits)) result$fit else NULL),
    error    = FALSE,
    error_message = "~"
  )
  tibble::as_tibble(plist, .name_repair = .name_repair)
}


#' Merge fold results within an ROI to a classification_result and compute performance
#'
#' This method concatenates event-level predictions across all CV folds, filters
#' out zero-probability events (which occur when events fall outside the test fold),
#' and builds a standard classification_result for metric computation.
#'
#' @param obj hrfdecoder_model specification object
#' @param result_set Tibble containing format_result outputs from all folds
#' @param indices ROI/searchlight voxel indices
#' @param id Unique identifier for this ROI/searchlight sphere
#' @param ... Additional arguments (unused)
#' @return A tibble with one row containing the classification_result, performance
#'   metrics, optional fast metrics, primary metric name/value, and error status
#' @export
merge_results.hrfdecoder_model <- function(obj, result_set, indices, id, ...) {
  if (any(result_set$error)) {
    emsg <- result_set$error_message[which(result_set$error)[1]]
    return(tibble::tibble(
      result = list(NULL), indices = list(indices),
      performance = list(NULL), id = id,
      error = TRUE, error_message = emsg
    ))
  }

  # Concatenate per-fold event-level outputs
  probs_list <- lapply(result_set$probs, as.matrix)
  probs <- do.call(rbind, probs_list)
  y_true <- do.call(c, result_set$y_true)
  predicted <- do.call(c, result_set$class)

  # Filter out rows with zero total probability (events not belonging to this fold)
  valid_mask <- rowSums(probs) > 1e-10
  if (!any(valid_mask)) {
    return(tibble::tibble(
      result = list(NULL), indices = list(indices),
      performance = list(NULL), id = id,
      error = TRUE,
      error_message = "hrfdecoder_model: All events have zero probabilities across all folds"
    ))
  }
  probs <- probs[valid_mask, , drop = FALSE]
  y_true <- y_true[valid_mask]
  predicted <- predicted[valid_mask]

  # Build a standard classification_result for downstream metrics
  cres <- classification_result(y_true, predicted, probs, testind = NULL, test_design = NULL, predictor = NULL)

  # Legacy performance (named vector) used by rMVPA searchlight combiners
  perf <- if (isTRUE(obj$compute_performance)) compute_performance(obj, cres) else NULL

  # Optional fast metrics via hrfdecoder (kept separate to avoid breaking callers)
  metrics_tbl <- NULL
  if (requireNamespace("hrfdecoder", quietly = TRUE)) {
    # Select metrics requested or use defaults
    req_metrics <- obj$metrics %||% c("auc_ovo", "accuracy")
    metrics_tbl <- tryCatch(
      hrfdecoder::hrf_metrics(probs, y_true, metrics = req_metrics),
      error = function(...) NULL
    )
  }

  # Primary metric/value surfaced at top-level for easy display
  pm <- obj$primary_metric %||% if (!is.null(metrics_tbl)) metrics_tbl$metric[1] else NA_character_
  pv <- if (!is.null(metrics_tbl) && !is.na(pm)) {
    metrics_tbl$value[match(pm, metrics_tbl$metric)]
  } else {
    # Fallback to accuracy if possible, otherwise NA
    mean(predicted == y_true)
  }

  tibble::tibble(
    result = list(cres),
    indices = list(indices),
    performance = list(perf),
    metrics = list(metrics_tbl),
    primary_metric = pm,
    primary_value = pv,
    id = id,
    error = FALSE,
    error_message = "~"
  )
}


#' Compute performance for hrfdecoder_model
#'
#' This method calls the performance function stored in the model specification
#' to compute metrics on the classification_result produced by merge_results.
#'
#' @param obj hrfdecoder_model specification object
#' @param result A classification_result object from merge_results
#' @return Named numeric vector of performance metrics (e.g., Accuracy, AUC)
#' @export
compute_performance.hrfdecoder_model <- function(obj, result) {
  obj$performance(result)
}
