#' Get Unique Region IDs
#'
#' Extract unique region IDs from a region mask, handling both volume and surface data.
#'
#' @param region_mask A region mask object (NeuroVol or NeuroSurface)
#' @param ... Additional arguments passed to methods
#'
#' @return A sorted vector of unique region IDs greater than 0
#' @keywords internal
get_unique_regions <- function(region_mask, ...) {
  UseMethod("get_unique_regions")
}

#' Strip Dataset from Model Specification
#'
#' Removes the potentially large dataset component from a model specification object
#' to avoid copying it during parallel processing.
#'
#' @param obj The model specification object.
#' @param ... Additional arguments.
#' @return The model specification object with the `dataset` element removed or set to NULL.
#' @note For internal parallel dispatch, \code{as_worker_spec()} is now
#'   preferred. \code{strip_dataset} remains exported for backward
#'   compatibility and external use.
#' @rdname strip_dataset-methods
#' @examples
#' \donttest{
#'   ds <- gen_sample_dataset(D = c(4, 4, 4), nobs = 20)
#'   mdl <- load_model("sda_notune")
#'   mspec <- mvpa_model(mdl, ds$dataset, ds$design, "classification")
#'   stripped <- strip_dataset(mspec)
#'   is.null(stripped$dataset)
#' }
#' @export
strip_dataset <- function(obj, ...) {
  UseMethod("strip_dataset")
}

#' Select Features
#'
#' Given a \code{feature_selection} specification object and a dataset, returns the set of selected features as a binary vector.
#'
#' @param obj The \code{feature_selection} object specifying the feature selection method and its parameters.
#' @param X The dataset containing the training features. This can be a matrix or a \code{ROIVolume} or \code{ROISurface} object.
#' @param Y The dependent variable as a factor or numeric variable.
#' @param ... Additional arguments to be passed to the method-specific function.
#' 
#' @return A logical vector indicating the columns of \code{X} matrix that were selected.
#' 
#' @examples 
#' fsel <- feature_selector("FTest", "top_k", 2)
#' coords <- rbind(c(1,1,1), c(2,2,2), c(3,3,3))
#' space <- neuroim2::NeuroSpace(c(10,10,10))
#' roi_data <- matrix(rnorm(100*3), 100, 3)
#' ROI <- neuroim2::ROIVec(space, coords=coords, roi_data)
#' Y <- factor(rep(c("a", "b"), each=50))
#' featureMask <- select_features(fsel, neuroim2::values(ROI), Y)
#' sum(featureMask) == 2
#'
#' fsel2 <- feature_selector("FTest", "top_p", .1)
#' featureMask <- select_features(fsel2, neuroim2::values(ROI), Y)
#'
#' @rdname select_features-methods
#' @export
select_features <- function(obj, X, Y, ...) {
  UseMethod("select_features")
}

#' Format Result Object
#'
#' @param obj The result object to be formatted.
#' @param result The result object to be formatted.
#' @param error_message An optional error message.
#' @param context Optional contextual metadata (e.g., ROI details) supplied to method implementations.
#' @param ... Additional arguments to be passed to the method-specific function.
#' @return A formatted result object, typically a tibble row.
#'
#' @examples
#' \donttest{
#' dataset <- gen_sample_dataset(D = c(6, 6, 6), nobs = 20,
#'                               response_type = "categorical",
#'                               data_mode = "image")
#' cval <- blocked_cross_validation(dataset$design$block_var)
#' model <- load_model("sda_notune")
#' mspec <- mvpa_model(model, dataset$dataset, dataset$design,
#'                     "classification", crossval = cval)
#' # Typically called internally during processing
#' format_result(mspec, result = NULL, error_message = "example",
#'               context = list(test = 1, ytest = factor("a")))
#' }
#' @export
format_result <- function(obj, result, error_message, ...) {
  UseMethod("format_result")
}


#' Merge Multiple Results
#'
#' @param obj The base object containing merge specifications
#' @param result_set List of results to be merged
#' @param indices List of indices corresponding to each result
#' @param id Identifier for the merged result
#' @param ... Additional arguments passed to specific merge methods
#'
#' @return A merged result object containing:
#' \itemize{
#'   \item Combined results from all input objects
#'   \item Associated indices
#'   \item Merged metadata
#' }
#'
#' @examples
#' \donttest{
#' ds <- gen_sample_dataset(D = c(5, 5, 5), nobs = 20, nlevels = 2)
#' model <- load_model("sda_notune")
#' mspec <- mvpa_model(
#'   model   = model,
#'   dataset = ds$dataset,
#'   design  = ds$design,
#'   model_type = "classification"
#' )
#'
#' # Construct a minimal result_set resembling a single CV fold:
#' obs  <- ds$design$y_train
#' levs <- levels(obs)
#' prob_mat <- matrix(
#'   1 / length(levs),
#'   nrow = length(obs),
#'   ncol = length(levs),
#'   dimnames = list(NULL, levs)
#' )
#'
#' result_set <- tibble::tibble(
#'   test_ind      = list(seq_along(obs)),
#'   probs         = list(prob_mat),
#'   error         = FALSE,
#'   error_message = "~"
#' )
#'
#' merge_results(mspec, result_set, indices = list(1:10), id = 1)
#' }
#'
#' @export
merge_results <- function(obj, result_set, indices, id, ...) {
  UseMethod("merge_results")
}

#' Run Future
#'
#' Run a future-based computation defined by the object and frame.
#'
#' @param obj An object specifying the computation.
#' @param frame A data frame or environment containing data for the computation.
#' @param processor A function or object specifying how to process the frame.
#' @param ... Additional arguments passed to the method-specific function.
#'
#' @return A tibble containing the processed results.
#'
#' @examples
#' frame <- tibble::tibble(
#'   .id = 1:2,
#'   rnum = c("roi1", "roi2"),
#'   roi = list(1:3, 4:5),
#'   size = c(3, 2)
#' )
#' mod_spec <- list(process_roi = function(mod_spec, roi, rnum, ...) {
#'   tibble::tibble(
#'     result = list(mean(roi)),
#'     indices = list(seq_along(roi)),
#'     performance = list(NULL),
#'     id = rnum
#'   )
#' })
#' run_future(mod_spec, frame, NULL)
#'
#' @rdname run_future-methods
#' @keywords internal
#' @export
run_future <- function(obj, frame, processor, ...) {
  UseMethod("run_future")
}

#' Fit a Model on a Single ROI
#'
#' The primary dispatch point for per-ROI analysis in the new architecture.
#' Each model type implements this method to encapsulate its complete analysis
#' pipeline (including any internal cross-validation).
#'
#' Models that implement \code{fit_roi} are automatically preferred by
#' \code{\link{process_roi.default}} over the legacy dispatch chain
#' (\code{internal_crossval} / \code{external_crossval} / \code{train_model}).
#'
#' @param model The model specification object.
#' @param roi_data A list with components:
#'   \describe{
#'     \item{train_data}{Numeric matrix (observations x features)}
#'     \item{test_data}{Numeric matrix or NULL}
#'     \item{indices}{Integer vector of voxel/vertex indices}
#'     \item{train_roi}{The raw ROI object (for models needing spatial info)}
#'     \item{test_roi}{The raw test ROI object or NULL}
#'   }
#' @param context A list with components:
#'   \describe{
#'     \item{design}{The design object}
#'     \item{cv_spec}{Cross-validation specification or NULL}
#'     \item{id}{ROI identifier (center voxel ID or region number)}
#'     \item{center_global_id}{Global index of center voxel or NA}
#'   }
#' @param ... Additional model-specific arguments.
#' @return A \code{\link{roi_result}} object.
#'
#' @seealso \code{\link{roi_result}}, \code{\link{output_schema}},
#'   \code{\link{process_roi}}
#'
#' @section Architecture TODO:
#' \code{fit_roi} is currently scoped to ROI/searchlight iteration.
#' Whole-brain global analysis uses \code{\link{run_global}} via a separate
#' \code{cv_run_global -> train_model} path. A future refactor may unify these
#' flows under a single fit contract.
#'
#' @examples
#' \donttest{
#'   # fit_roi is typically called internally by process_roi.default.
#'   # To implement for a new model class:
#'   # fit_roi.my_model <- function(model, roi_data, context, ...) {
#'   #   metrics <- c(accuracy = 0.85)
#'   #   roi_result(metrics, roi_data$indices, context$id)
#'   # }
#' }
#' @export
fit_roi <- function(model, roi_data, context, ...) {
  UseMethod("fit_roi")
}

#' Declare the Output Metric Schema for a Model
#'
#' Returns a named list describing the metrics that \code{\link{fit_roi}}
#' produces for this model type. Used by the generic combiner (Phase 2) to
#' build the correct number of output maps without model-specific
#' \code{combine_*} functions.
#'
#' @param model A model specification object.
#' @return A named list where names are metric names and values describe
#'   the type: \code{"scalar"} (one value per ROI) or \code{"vector[N]"}
#'   (N values per ROI). Returns \code{NULL} for models using the legacy
#'   combiner path.
#'
#' @examples
#' \donttest{
#'   # Default returns NULL (legacy path)
#'   ds <- gen_sample_dataset(c(4,4,4), 20)
#'   mdl <- load_model("sda_notune")
#'   mspec <- mvpa_model(mdl, ds$dataset, ds$design, "classification")
#'   output_schema(mspec)  # NULL
#' }
#' @export
output_schema <- function(model) {
  UseMethod("output_schema")
}

#' @export
output_schema.default <- function(model) {
  NULL
}

#' Check if a specific fit_roi method exists for a model object
#'
#' Walks the S3 class chain looking for a registered \code{fit_roi} method.
#' Returns \code{FALSE} if only the implicit default dispatch would apply.
#'
#' @param obj A model specification object.
#' @return Logical.
#' @keywords internal
#' @noRd
.has_fit_roi <- function(obj) {
  for (cls in class(obj)) {
    if (!is.null(utils::getS3method("fit_roi", cls, optional = TRUE))) {
      return(TRUE)
    }
  }
  FALSE
}

#' Process ROI
#'
#' Process a region of interest (ROI) and return the formatted results.
#'
#' @param mod_spec The model specification object.
#' @param roi The region of interest data.
#' @param rnum A numeric or string identifier for the ROI.
#' @param ... Additional arguments passed to the method-specific function.
#'
#' @return A tibble row containing the performance metrics for the ROI.
#'
#' @examples
#' \donttest{
#'   ds <- gen_sample_dataset(c(4, 4, 4), 20, blocks = 2)
#'   cv <- blocked_cross_validation(ds$design$block_var)
#'   mdl <- load_model("sda_notune")
#'   spec <- mvpa_model(
#'     model = mdl,
#'     dataset = ds$dataset,
#'     design = ds$design,
#'     model_type = "classification",
#'     crossval = cv
#'   )
#'   vox <- sample(which(ds$dataset$mask > 0), 30)
#'   samp <- data_sample(ds$dataset, vox)
#'   roi_obj <- as_roi(samp, ds$dataset)
#'   process_roi(spec, roi_obj, 1)
#' }
#'
#' @rdname process_roi-methods
#' @export
process_roi <- function(mod_spec, roi, rnum, ...) {
  UseMethod("process_roi")
}

#' Default Process ROI Method
#'
#' @param center_global_id Optional global ID of the center voxel. Defaults to NA.
#'
#' @rdname process_roi-methods
#' @export
process_roi.default <- function(mod_spec, roi, rnum, center_global_id = NA, ...) {
  # --- fit_roi shim (Phase 1) ---
  # If the model class provides a fit_roi method, prefer it over the

  # legacy 3-way dispatch. The roi_result is converted to the tibble
  # format that existing combiners expect.
  if (.has_fit_roi(mod_spec)) {
    roi_data <- list(
      train_data = as.matrix(neuroim2::values(roi$train_roi)),
      test_data = if (!is.null(roi$test_roi)) as.matrix(neuroim2::values(roi$test_roi)) else NULL,
      indices = neuroim2::indices(roi$train_roi),
      train_roi = roi$train_roi,
      test_roi = roi$test_roi
    )
    context <- list(
      design = mod_spec$design,
      cv_spec = if (!is.null(mod_spec$crossval)) mod_spec$crossval else mod_spec$cv_spec,
      id = rnum,
      center_global_id = center_global_id
    )
    res <- tryCatch(
      fit_roi(mod_spec, roi_data, context, ...),
      error = function(e) {
        futile.logger::flog.warn(
          "process_roi: fit_roi failed for ROI %s: %s", rnum, conditionMessage(e))
        roi_result(
          metrics = NULL,
          indices = neuroim2::indices(roi$train_roi),
          id = rnum,
          error = TRUE,
          error_message = conditionMessage(e)
        )
      }
    )
    return(roi_result_to_tibble(res))
  }

  # --- Legacy dispatch ---
  warning("process_roi legacy dispatch is deprecated. Implement a fit_roi.* method for your model class.",
          call. = FALSE)
  # Capture additional arguments to pass down
  dots <- list(...)
  if (!is.null(mod_spec$process_roi)) {
    # Pass center_global_id and dots to user's custom processor
    do.call(mod_spec$process_roi, c(list(mod_spec, roi, rnum, center_global_id = center_global_id), dots))
  } else if (has_test_set(mod_spec)) {
    # Pass center_global_id and dots to external_crossval
    do.call(external_crossval, c(list(mod_spec, roi, rnum, center_global_id = center_global_id), dots))
  } else if (has_crossval(mod_spec)) {
    # Pass center_global_id and dots to internal_crossval
    do.call(internal_crossval, c(list(mod_spec, roi, rnum, center_global_id = center_global_id), dots))
  } else {
    # Pass center_global_id and dots to process_roi_default
    do.call(process_roi_default, c(list(mod_spec, roi, rnum, center_global_id = center_global_id), dots))
  }
}

#' Default ROI Processing Helper
#'
#' @param mod_spec The model specification object.
#' @param roi The ROI containing training data.
#' @param rnum The region number or identifier.
#' @param center_global_id Optional global ID of the center voxel. Defaults to NA.
#' @param ... Additional arguments passed to specific methods.
#' @keywords internal
#' @noRd
#' @importFrom neuroim2 indices values
#' @importFrom tibble as_tibble tibble
#' @importFrom futile.logger flog.warn
process_roi_default <- function(mod_spec, roi, rnum, center_global_id = NA, ...) {
  # This helper is called by process_roi.default for models 
  # that don't use internal cross-validation.
  # It runs train_model and then passes the result to merge_results
  # for final performance computation and formatting.

  xtrain_mat <- as.matrix(neuroim2::values(roi$train_roi))
  xtrain <- if (.matrix_first_roi_enabled()) {
    xtrain_mat
  } else {
    tibble::as_tibble(xtrain_mat, .name_repair=.name_repair)
  }
  ind <- indices(roi$train_roi)
  
  # Determine center_local_id based on center_global_id
  center_local_id <- NA
  if (!is.na(center_global_id)) {
      center_local_id <- match(center_global_id, ind)
      if (is.na(center_local_id)) {
          stop(paste0("process_roi_default: Provided center_global_id ", center_global_id, 
                      " not found within the voxel indices for this ROI/searchlight (id: ", rnum, ")."))
      }
  }
  
  # Prepare sl_info
  sl_info <- list(center_local_id = center_local_id, center_global_id = center_global_id)
  
  # Run train_model, passing sl_info
  # Assuming train_model methods will accept sl_info if needed
  # Also pass other ... args
  dots <- list(...)
  train_args <- c(list(mod_spec, xtrain, y = NULL, indices = ind, sl_info = sl_info), dots)
  if (is.null(train_args$cv_spec)) {
    if (!is.null(mod_spec$cv_spec)) {
      train_args$cv_spec <- mod_spec$cv_spec
    } else if (!is.null(mod_spec$crossval)) {
      train_args$cv_spec <- mod_spec$crossval
    }
  }
  train_args$quiet_error <- TRUE
  train_result_obj <- try(do.call(train_model, train_args), silent = TRUE)
  
  # Prepare a result set structure for merge_results
  if (inherits(train_result_obj, "try-error")) {
    # If training failed, create an error result set for merge_results
    error_msg <- attr(train_result_obj, "condition")$message
    result_set <- tibble::tibble(
      result = list(NULL), # No result from train_model
      error = TRUE,
      error_message = ifelse(is.null(error_msg), "Unknown training error", error_msg)
      # We don't need to mimic all columns internal_crossval might produce,
      # only what merge_results requires for error handling.
    )
     futile.logger::flog.warn("ROI %s: train_model failed: %s", rnum, error_msg)
     
  } else {
    # If training succeeded, create a success result set for merge_results
    # Store the *output* of train_model in the 'result' column. 
    # merge_results.vector_rsa_model expects the scores vector here.
     result_set <- tibble::tibble(
       result = list(train_result_obj), # Store train_model output here
       error = FALSE,
       error_message = "~"
       # merge_results will compute the 'performance' column.
     )
  }
  
  # Call merge_results to compute final performance and format the output tibble
  # merge_results handles both success and error cases based on result_set$error
  final_result <- merge_results(mod_spec, result_set, indices=ind, id=rnum)
  return(final_result)
}

#'
#' Train a classification, regression, or representational model.
#'
#' This is a generic function that trains a model based on the provided
#' model specification object. Different model types will have different
#' methods implemented with specific parameters.
#'
#' @param obj The model specification object.
#' @param ... Additional arguments to be passed to the method-specific function.
#'
#' @return A trained model object. The exact return value depends on the specific
#'   method implementation.
#'
#' @examples
#' \donttest{
#'   # Generate a small sample dataset for classification
#'   dset_info <- gen_sample_dataset(
#'     D = c(8, 8, 8),
#'     nobs = 20,
#'     response_type = "categorical",
#'     data_mode = "image",
#'     nlevels = 2
#'   )
#'
#'   # Create a cross-validation specification
#'   cval <- blocked_cross_validation(dset_info$design$block_var)
#'
#'   # Load a simple classifier
#'   sda_model <- load_model("sda_notune")
#'
#'   # Create an MVPA model specification
#'   mspec <- mvpa_model(
#'     model = sda_model,
#'     dataset = dset_info$dataset,
#'     design = dset_info$design,
#'     model_type = "classification",
#'     crossval = cval
#'   )
#'
#'   # Extract voxel-wise data as a numeric matrix for training
#'   vox <- which(dset_info$dataset$mask > 0)
#'   X   <- neuroim2::series(dset_info$dataset$train_data, vox)
#'
#'   # Train the model on the voxel matrix
#'   fit <- train_model(
#'     mspec,
#'     X,
#'     dset_info$design$y_train,
#'     indices = vox
#'   )
#' }
#' @export
train_model <- function(obj,...) {
  UseMethod("train_model")
}

#' Training Labels/Response Extraction
#'
#' Extract the training labels or response variable from an object.
#'
#' @param obj The object from which to extract the training response variable.
#' @return The training response variable (factor for classification, numeric for regression).
#'
#' @rdname y_train-methods
#' @examples
#' ds <- gen_sample_dataset(D = c(4, 4, 4), nobs = 10)
#' y_train(ds$design)
#' @export
y_train <- function(obj) {
  UseMethod("y_train")
}

#' Test Labels/Response Extraction
#'
#' Extract the test labels or response variable from an object.
#'
#' @param obj The object from which to extract the test response variable.
#' @return The test response variable.
#'
#' @rdname y_test-methods
#' @examples
#' ds <- gen_sample_dataset(D = c(4, 4, 4), nobs = 10, external_test = TRUE)
#' y_test(ds$design)
#' @export
y_test <- function(obj) {
  UseMethod("y_test")
}

#' Test Design Extraction
#'
#' Return the design table associated with the test set from an object.
#'
#' @param obj The object from which to extract the test design table.
#' @return A data frame containing the test set design variables.
#'
#' @rdname test_design-methods
#' @examples
#' ds <- gen_sample_dataset(D = c(4, 4, 4), nobs = 10, external_test = TRUE)
#' test_design(ds$design)
#' @export
test_design <- function(obj) {
  UseMethod("test_design")
}

#' Fit Model
#'
#' Fit a classification or regression model.
#'
#' @param obj A model fitting object.
#' @param roi_x An ROI containing the training data.
#' @param y The response vector.
#' @param wts A set of case weights.
#' @param param Tuning parameters.
#' @param lev Factor levels (for classification).
#' @param last Logical indicating if this is the last iteration.
#' @param classProbs Logical indicating if class probabilities should be returned.
#' @param ... Additional arguments to be passed to the method-specific function.
#'
#' @examples
#' \donttest{
#'   ds <- gen_sample_dataset(
#'     D = c(6, 6, 6), nobs = 20,
#'     response_type = "categorical",
#'     data_mode = "image", nlevels = 2
#'   )
#'   mdl <- load_model("sda_notune")
#'   mspec <- mvpa_model(
#'     model = mdl,
#'     dataset = ds$dataset,
#'     design  = ds$design,
#'     model_type = "classification"
#'   )
#'   grid <- tune_grid(mspec, ds$dataset$train_data, ds$design$y_train, len = 1)
#'   fit  <- fit_model(mspec, ds$dataset$train_data,
#'                    ds$design$y_train, wts = NULL, param = grid)
#' }
#'
#' @rdname fit_model-methods
#' 
#' @export
fit_model <- function(obj, roi_x, y, wts, param, lev=NULL, last=FALSE, classProbs=FALSE, ...) {
  UseMethod("fit_model")
}

#' Extract Tuning Grid
#'
#' Returns the parameter grid used to tune a model.
#'
#' @param obj A model or model specification.
#' @param x Training data.
#' @param y Response variable.
#' @param len Number of parameter sets to generate.
#'
#' @return A data frame of tuning parameter combinations.
#' @rdname tune_grid-methods
#' @export
#'
#' @examples
#' \donttest{
#' ds  <- gen_sample_dataset(D = c(5, 5, 5), nobs = 10)
#' mdl <- load_model("sda_notune")
#' mspec <- mvpa_model(mdl, ds$dataset, ds$design, model_type = "classification")
#' tune_grid(mspec, ds$dataset$train_data, ds$design$y_train, len = 1)
#' }
tune_grid <- function(obj, x, y, len) {
  UseMethod("tune_grid")
}

#' Test Set Availability
#'
#' Determine whether the object contains a separate test set.
#'
#' @param obj Object to query.
#'
#' @return Logical indicating if a test set exists.
#' @rdname has_test_set-methods
#' @export
#'
#' @examples
#' ds <- gen_sample_dataset(D = c(4, 4, 4), nobs = 10, external_test = TRUE)
#' has_test_set(ds$design)
has_test_set <- function(obj) {
  UseMethod("has_test_set")
}

#' Cross-Validation Availability
#'
#' Determine whether cross-validation is specified for the object.
#'
#' @param obj Model specification object.
#'
#' @return Logical indicating if cross-validation will be performed.
#' @rdname has_crossval-methods
#' @export
#'
#' @examples
#' ds <- gen_sample_dataset(D = c(4, 4, 4), nobs = 10)
#' cval <- blocked_cross_validation(ds$design$block_var)
#' mdl <- load_model("sda_notune")
#' mspec <- mvpa_model(mdl, ds$dataset, ds$design,
#'                     "classification", crossval = cval)
#' has_crossval(mspec)
has_crossval <- function(obj) {
  UseMethod("has_crossval")
}

#' @export
has_crossval.default <- function(obj) {
  FALSE
}

#' Compute Performance Metrics
#'
#' Generic function to compute performance metrics from result objects.
#'
#' @param x Result object from a classification or regression analysis.
#' @param ... Additional arguments passed to methods.
#'
#' @return Named numeric vector of performance metrics.
#' @rdname performance-methods
#' @export
#'
#' @examples
#' cres <- binary_classification_result(
#'   observed  = factor(c("a", "b")),
#'   predicted = factor(c("a", "b")),
#'   probs     = matrix(c(0.8, 0.2, 0.3, 0.7), ncol = 2,
#'                      dimnames = list(NULL, c("a", "b")))
#' )
#' performance(cres)
performance <- function(x, ...) {
  UseMethod("performance")
}

#' Compute Performance for an Object
#'
#' Delegates calculation of performance metrics to the appropriate method.
#'
#' @param obj Model specification or object capable of computing performance.
#' @param result The classification/regression result to evaluate.
#'
#' @return Named numeric vector of performance metrics.
#' @rdname compute_performance-methods
#' @export
#'
#' @examples
#' cres <- binary_classification_result(
#'   observed  = factor(c("a", "b")),
#'   predicted = factor(c("a", "b")),
#'   probs     = matrix(c(0.8, 0.2, 0.3, 0.7), ncol = 2,
#'                      dimnames = list(NULL, c("a", "b")))
#' )
#' dummy <- list(performance = performance)
#' class(dummy) <- "mvpa_model"
#' compute_performance(dummy, cres)
compute_performance <- function(obj, result) {
  UseMethod("compute_performance")
}

#' Merge Multiple Classification/Regression Results
#'
#' This function merges two or more classification/regression result objects.
#'
#' @param x The first classification/regression result object.
#' @param ... Additional classification/regression result objects.
#'
#' @return A single merged classification/regression result object.
#' @examples
#' cres1 <- binary_classification_result(
#'   observed = factor(c("a","b")),
#'   predicted = factor(c("a","b")),
#'   probs = matrix(c(.8,.2,.3,.7), ncol=2, dimnames=list(NULL,c("a","b")))
#' )
#' cres2 <- binary_classification_result(
#'   observed = factor(c("b","a")),
#'   predicted = factor(c("b","a")),
#'   probs = matrix(c(.2,.8,.7,.3), ncol=2, dimnames=list(NULL,c("a","b")))
#' )
#' merged <- merge_classif_results(cres1, cres2)
#' @export
merge_classif_results <- function(x, ...) {
  UseMethod("merge_classif_results")
}

#' Get Multiple Data Samples
#'
#' Extract multiple data samples based on a list of voxel/index sets from a dataset object.
#'
#' @param obj The input dataset object.
#' @param vox_list A list of vectors containing voxel indices to extract.
#'
#' @return A list of data samples.
#' @examples
#' \donttest{
#'   ds <- gen_sample_dataset(c(5,5,5), 20, blocks=2)
#'   vox_list <- list(1:10, 11:20)
#'   samples <- get_samples(ds$dataset, vox_list)
#' }
#' @export
get_samples <- function(obj, vox_list) {
  UseMethod("get_samples")
}

#' Extract Sample from Dataset
#'
#' Extract a sample from a given dataset object.
#'
#' @param obj The input dataset object.
#' @param vox The voxel indices/coordinates.
#' @param ... Additional arguments to methods.
#'
#' @return A sample extracted from the dataset.
#' @export
data_sample <- function(obj, vox, ...) {
  UseMethod("data_sample")
}

#' Convert object to ROI
#'
#' Convert the provided object into an ROIVolume or ROISurface object.
#'
#' @param obj The object to be converted.
#' @param data The associated data object.
#' @param ... Additional arguments passed to methods.
#' @return An ROIVolume or ROISurface object.
#' @examples
#' \donttest{
#'   ds <- gen_sample_dataset(c(5,5,5), 20, blocks=2)
#'   vox <- sample(which(ds$dataset$mask > 0), 10)
#'   samp <- data_sample(ds$dataset, vox)
#'   roi <- as_roi(samp, ds$dataset)
#' }
#' @export
as_roi <- function(obj, data, ...) {
  UseMethod("as_roi")
}

#' Generate Searchlight Iterator
#'
#' Generate a searchlight iterator suitable for given data.
#'
#' @param obj The input dataset object.
#' @param ... Additional arguments to methods.
#'
#' @return A searchlight iterator object.
#' @examples
#' \donttest{
#'   ds <- gen_sample_dataset(c(5,5,5), 20)
#'   sl <- get_searchlight(ds$dataset, radius=4)
#' }
#' @export
get_searchlight <- function(obj, ...) {
  UseMethod("get_searchlight")
}

#' Wrap Output
#'
#' Wrap output values into a desired format.
#'
#' @param obj The object used to determine the wrapping method.
#' @param vals The values to be wrapped.
#' @param ... Additional arguments passed to methods.
#' @return A wrapped output object.
#' @keywords internal
wrap_output <- function(obj, vals, ...) {
  UseMethod("wrap_output")
}

#' Merge Predictions
#'
#' Combine predictions from multiple models on the same test set.
#'
#' @param obj1 The first object containing predictions.
#' @param rest Other objects containing predictions.
#' @param ... Additional arguments. Methods for this generic may implement specific arguments
#'   such as `weights` to control how predictions are combined.
#' @return A combined object with merged predictions.
#' @examples
#' \donttest{
#' p1 <- structure(
#'   list(class = c("a","b"),
#'        probs = matrix(c(.8,.2,.3,.7), ncol=2, dimnames=list(NULL,c("a","b")))),
#'   class = c("classification_prediction", "prediction", "list")
#' )
#' p2 <- structure(
#'   list(class = c("b","a"),
#'        probs = matrix(c(.4,.6,.6,.4), ncol=2, dimnames=list(NULL,c("a","b")))),
#'   class = c("classification_prediction", "prediction", "list")
#' )
#' merge_predictions(p1, list(p2))
#' }
#' @export
merge_predictions <- function(obj1, rest, ...) {
  UseMethod("merge_predictions")
}

#' Extract Row-wise Subset of a Result
#'
#' Extract a subset of rows from a classification/regression result object.
#'
#' @param x The input result object.
#' @param indices Row indices to extract.
#'
#' @return A new result object with the specified rows.
#' @examples
#' cres <- binary_classification_result(
#'   observed = factor(c("a","b","a")),
#'   predicted = factor(c("a","b","b")),
#'   probs = matrix(c(.8,.2,.4,.3,.7,.6), ncol=2, dimnames=list(NULL,c("a","b")))
#' )
#' sub_result(cres, 1:2)
#' @export
sub_result <- function(x, indices) {
  UseMethod("sub_result")
}

#' Get Number of Observations
#'
#' Retrieve the number of observations in an object.
#'
#' @param x The input object.
#' @return The number of observations.
#' @examples
#' ds <- gen_sample_dataset(c(5,5,5), 20)
#' nobs(ds$dataset)
#' @export
nobs <- function(x) {
  UseMethod("nobs")
}

#' Probability of Observed Class
#'
#' Extract the predicted probability for the observed class.
#'
#' @param x The object from which to extract the probability.
#' @return A vector of predicted probabilities.
#' @examples
#' cres <- binary_classification_result(
#'   observed = factor(c("a","b")),
#'   predicted = factor(c("a","b")),
#'   probs = matrix(c(.8,.2,.3,.7), ncol=2, dimnames=list(NULL,c("a","b")))
#' )
#' prob_observed(cres)
#' @export
prob_observed <- function(x) {
  UseMethod("prob_observed")
}

#' Number of Response Categories
#'
#' Get the number of response categories or levels.
#'
#' @param x The object from which to extract the number of categories.
#' @return The number of response categories.
#' @examples
#' des <- mvpa_design(
#'   train_design = data.frame(y = factor(c("a","b","a","b"))),
#'   y_train = factor(c("a","b","a","b")),
#'   block_var = factor(c("1","1","2","2"))
#' )
#' nresponses(des)
#' @export
nresponses <- function(x) {
  UseMethod("nresponses")
}

#' Predict Model Output
#'
#' Generic function to predict outcomes from a fitted model object using new data.
#'
#' @param object A fitted model object for which a prediction method is defined.
#' @param fit The fitted model object, often returned by `train_model`.
#'              (Note: For some models, `object` itself might be the fit).
#' @param newdata New data (e.g., a matrix or data frame) for which to make predictions.
#'                The structure should be compatible with what the model was trained on.
#' @param ... Additional arguments passed to specific prediction methods.
#'
#' @return Predictions whose structure depends on the specific method (e.g., a vector,
#'   matrix, or data frame).
#' @examples
#' \dontrun{
#'   preds <- predict_model(model_spec, fitted_model, new_data)
#' }
#' @export
predict_model <- function(object, fit, newdata, ...) {
  UseMethod("predict_model")
}

#' @keywords internal
#' @noRd
.name_repair = ~ vctrs::vec_as_names(..., repair = "unique", quiet = TRUE)

#' @keywords internal
#' @noRd
.matrix_first_roi_enabled <- function() {
  TRUE
}

#' Run Searchlight Analysis
#'
#' Execute a searchlight analysis using multivariate pattern analysis.
#'
#' @param model_spec A \code{mvpa_model} instance containing the model specifications
#' @param radius The searchlight radius in millimeters
#' @param method The type of searchlight, either 'randomized' or 'standard'
#' @param niter The number of searchlight iterations (used only for 'randomized' method)
#' @param backend Execution backend: \code{"default"} (standard pipeline),
#'   \code{"shard"} (shared-memory backend), or \code{"auto"} (try shard and
#'   fall back to default).
#' @param ... Extra arguments passed to specific searchlight methods. Currently supported:
#'   \itemize{
#'     \item \code{batch_size}: Integer specifying how many searchlights (ROIs) are grouped
#'           into a single batch for processing by \code{\link{mvpa_iterate}}. Each batch is
#'           launched sequentially, while ROIs within a batch are processed in parallel (using
#'           the active \pkg{future}/\pkg{furrr} plan). The default is 10\% of the total number
#'           of searchlights. Smaller values lower peak memory usage and can improve load
#'           balancing, at the cost of more scheduling overhead; larger values reduce overhead
#'           but require more memory per worker. As a rule of thumb, start with the default,
#'           decrease \code{batch_size} if you hit memory limits, and increase it if you have
#'           many ROIs, ample RAM, and see low CPU utilization.
#'   }
#'
#' @section Progress reporting:
#' When \code{verbose = TRUE} is passed via \code{...} and the
#' \pkg{progressr} package is installed, real-time per-ROI progress updates
#' are emitted from parallel workers. To enable a progress bar:
#' \preformatted{
#'   library(progressr)
#'   handlers(global = TRUE)               # once per session
#'   result <- run_searchlight(mspec, radius = 8, verbose = TRUE)
#' }
#' Without \pkg{progressr}, only coarse batch-level log messages are shown.
#'
#' @return A named list of \code{NeuroVol} objects containing performance metrics (e.g., AUC) at each voxel location
#'
#' @examples
#' \donttest{
#'   # Generate sample dataset with categorical response
#'   dataset <- gen_sample_dataset(
#'     D = c(8,8,8),           # 8x8x8 volume
#'     nobs = 100,             # 100 observations
#'     response_type = "categorical",
#'     data_mode = "image",
#'     blocks = 3,             # 3 blocks for cross-validation
#'     nlevels = 2             # binary classification
#'   )
#'   
#'   # Create cross-validation specification using blocks
#'   cval <- blocked_cross_validation(dataset$design$block_var)
#'   
#'   # Load the SDA classifier (Shrinkage Discriminant Analysis)
#'   model <- load_model("sda_notune")
#'   
#'   # Create MVPA model
#'   mspec <- mvpa_model(
#'     model = model,
#'     dataset = dataset$dataset,
#'     design = dataset$design,
#'     model_type = "classification",
#'     crossval = cval
#'   )
#'   
#'   # Run searchlight analysis
#'   results <- run_searchlight(
#'     mspec,
#'     radius = 8,            # 8mm radius
#'     method = "standard"    # Use standard searchlight
#'   )
#'   
#'   # Run with custom batch size for memory management
#'   # results <- run_searchlight(
#'   #   mspec,
#'   #   radius = 8,
#'   #   method = "standard",
#'   #   batch_size = 500      # Process 500 searchlights per batch
#'   # )
#'   
#'   # Results contain performance metrics
#'   # Access them with results$performance
#' }
#'
#' @export
run_searchlight <- function(model_spec, radius, method = c("standard", "randomized", "resampled"),
                            niter = 4, backend = c("default", "shard", "auto"), ...) {
  UseMethod("run_searchlight")
}

#' Region of Interest Based MVPA Analysis
#'
#' Run a separate MVPA analysis for multiple disjoint regions of interest.
#'
#' @param model_spec A \code{mvpa_model} instance containing the model specifications
#' @param region_mask A \code{NeuroVol} or \code{NeuroSurface} object where each region is identified by a unique integer
#' @param coalesce_design_vars If \code{TRUE}, merges design variables into the prediction table (if present and generated). Default is \code{FALSE}.
#' @param processor An optional custom processor function for each region (ROI). If NULL (default), behavior depends on the \code{model_spec} class.
#' @param verbose If \code{TRUE}, print progress messages during iteration (default is \code{FALSE}).
#' @param backend Execution backend: \code{"default"} (standard pipeline),
#'   \code{"shard"} (shared-memory backend), or \code{"auto"} (try shard and
#'   fall back to default).
#' @param ... Extra arguments passed to specific regional analysis methods (e.g., `return_fits`, `compute_performance`).
#'
#' @section Progress reporting:
#' When \code{verbose = TRUE} is passed via \code{...} and the
#' \pkg{progressr} package is installed, real-time per-ROI progress updates
#' are emitted from parallel workers. To enable a progress bar:
#' \preformatted{
#'   library(progressr)
#'   handlers(global = TRUE)               # once per session
#'   result <- run_regional(mspec, region_mask, verbose = TRUE)
#' }
#' Without \pkg{progressr}, only coarse batch-level log messages are shown.
#'
#' @return A \code{regional_mvpa_result} object (list) containing:
#'   \item{performance_table}{A tibble of performance metrics for each region (if computed).}
#'   \item{prediction_table}{A tibble with detailed predictions for each observation/region (if generated).}
#'   \item{pooled_prediction_table}{Optional pooled trial-level prediction table when pooling is requested.}
#'   \item{pooled_performance}{Optional pooled performance metrics when pooling is requested.}
#'   \item{vol_results}{A list of volumetric maps representing performance metrics across space (if computed).}
#'   \item{fits}{A list of fitted model objects for each region (if requested via `return_fits=TRUE`).}
#'   \item{model_spec}{The original model specification object provided.} # Note: Original documentation said 'performance', clarified here.
#'
#' @examples
#' \donttest{
#'   # Generate sample dataset (3D volume with categorical response)
#'   dataset <- gen_sample_dataset(
#'     D = c(10,10,10),       # Small 10x10x10 volume
#'     nobs = 100,            # 100 observations
#'     nlevels = 3,           # 3 classes
#'     response_type = "categorical",
#'     data_mode = "image",
#'     blocks = 3             # 3 blocks for cross-validation
#'   )
#'   
#'   # Create region mask with 5 ROIs
#'   region_mask <- neuroim2::NeuroVol(
#'     sample(1:5, size=length(dataset$dataset$mask), replace=TRUE),
#'     neuroim2::space(dataset$dataset$mask)
#'   )
#'   
#'   # Create cross-validation specification
#'   cval <- blocked_cross_validation(dataset$design$block_var)
#'   
#'   # Load SDA classifier (Shrinkage Discriminant Analysis)
#'   model <- load_model("sda_notune")
#'   
#'   # Create MVPA model
#'   mspec <- mvpa_model(
#'     model = model,
#'     dataset = dataset$dataset,
#'     design = dataset$design,
#'     model_type = "classification",
#'     crossval = cval,
#'     return_fits = TRUE    # Return fitted models
#'   )
#'   
#'   # Run regional analysis
#'   results <- run_regional(mspec, region_mask)
#'   
#'   # Access results
#'   head(results$performance_table)     # Performance metrics
#'   head(results$prediction_table)      # Predictions
#'   first_roi_fit <- results$fits[[1]]  # First ROI's fitted model
#' }
#'
#' @rdname run_regional-methods
#' @export
run_regional <- function(model_spec, region_mask, backend = c("default", "shard", "auto"), ...) {
  UseMethod("run_regional")
}

#' Cross-validation samples
#'
#' Apply a cross-validation scheme to split the data into training and testing sets.
#'
#' @param obj A cross-validation control object.
#' @param data A data frame containing the predictors.
#' @param y A vector containing the response variable.
#' @param ... Extra arguments passed to the specific cross-validation methods (e.g., `id` for custom cross-validation).
#'
#' @return A tibble containing training and testing sets for each fold.
#'
#' @examples
#' cval <- kfold_cross_validation(len = 20, nfolds = 4)
#' dat  <- as.data.frame(matrix(rnorm(20 * 2), 20, 2))
#' y    <- factor(rep(letters[1:4], 5))
#' crossval_samples(cval, dat, y)
#' @export
crossval_samples <- function(obj, data, y,...) {
  UseMethod("crossval_samples")
}

#' Generic Pairwise Distance Computation
#'
#' Compute pairwise distances between rows of X using a specified distance object.
#'
#' @param obj A distance object specifying the distance measure.
#' @param X A numeric matrix of data points (rows = samples).
#' @param ... Additional arguments passed to methods.
#'
#' @return A matrix or dist object of pairwise distances.
#' @keywords internal
#' @noRd
pairwise_dist <- function(obj, X,...) {
  UseMethod("pairwise_dist")
}

#' Filter Region of Interest (ROI)
#'
#' Filter an ROI by removing columns with missing values or zero std dev.
#'
#' @param roi A list containing the train and test ROI data.
#' @param ... Additional arguments passed to methods.
#'
#' @return A list with filtered train and test ROI data.
#' @keywords internal
#' @export
filter_roi <- function(roi, ...) {
  UseMethod("filter_roi", roi$train_roi)
}

#' Get the Number of Folds
#'
#' An S3 generic method to retrieve the number of folds from a cross-validation specification object.
#'
#' @param obj A cross-validation specification object (e.g., inheriting from `cross_validation`).
#' @param ... Additional arguments passed to methods.
#' @return An integer representing the number of folds.
#' @examples
#' cval <- kfold_cross_validation(len = 20, nfolds = 4)
#' get_nfolds(cval)
#' @export
get_nfolds <- function(obj, ...) {
  UseMethod("get_nfolds")
}

#' Get Training Indices for a Fold
#'
#' An S3 generic method to retrieve the training sample indices for a specific fold 
#' from a cross-validation specification object.
#'
#' @param obj A cross-validation specification object (e.g., inheriting from `cross_validation`).
#' @param fold_num An integer specifying the fold number for which to retrieve training indices.
#' @param ... Additional arguments passed to methods.
#' @return An integer vector of training indices.
#' @examples
#' cval <- kfold_cross_validation(len = 20, nfolds = 4)
#' train_indices(cval, 1)
#' @export
train_indices <- function(obj, fold_num, ...) {
  UseMethod("train_indices")
}


#' Data sample
#'
#' Construct a light-weight data sample descriptor for a set of features (e.g., voxels)
#' associated with a dataset. Methods define how to interpret and convert this descriptor.
#'
#' @param obj An object (typically a dataset) from which to draw a sample.
#' @param vox The feature indices or coordinates defining the sample.
#' @param ... Additional arguments passed to methods.
#'
#' @return An object of class `data_sample` describing the sample; methods provide
#'   conversions to ROI structures or data frames as needed.
#' @examples
#' \dontrun{
#'   ds <- gen_sample_dataset(c(5,5,5), 20)
#'   vox <- sample(which(ds$dataset$mask > 0), 10)
#'   samp <- data_sample(ds$dataset, vox)
#' }
#' @name data_sample
NULL

#' Get Center IDs for Searchlight Iteration
#'
#' Returns the set of center IDs over which a searchlight analysis iterates.
#' For volumetric datasets this is the set of nonzero mask voxels; for clustered
#' datasets it is the sequence of cluster indices.
#'
#' @param dataset The dataset object.
#' @param ... Additional arguments.
#' @return Integer vector of center IDs.
#' @examples
#' ds <- gen_sample_dataset(c(5,5,5), 20)
#' ids <- get_center_ids(ds$dataset)
#' length(ids)
#' @export
get_center_ids <- function(dataset, ...) UseMethod("get_center_ids")

#' Build Spatial Output Map for a Single Metric
#'
#' Constructs a spatial object (e.g., \code{NeuroVol}, \code{NeuroSurface}) from
#' a numeric vector of per-center metric values and the corresponding center IDs.
#'
#' @param dataset The dataset object.
#' @param metric_vector Numeric vector of metric values (one per center).
#' @param ids Integer vector of center IDs corresponding to \code{metric_vector}.
#' @param ... Additional arguments.
#' @return A spatial object (\code{NeuroVol}, \code{NeuroSurface}, etc.)
#' @keywords internal
build_output_map <- function(dataset, metric_vector, ids, ...) UseMethod("build_output_map")

#' Get Searchlight Scope
#'
#' Returns the analysis scope for a dataset, controlling how \code{mvpa_iterate}
#' handles center IDs and whether \code{pobserved} maps are constructed.
#'
#' @param dataset The dataset object.
#' @param ... Additional arguments.
#' @return Character string: \code{"searchlight"} or \code{"regional"}.
#' @keywords internal
searchlight_scope <- function(dataset, ...) UseMethod("searchlight_scope")

#' Extract Raw Model Weights
#'
#' Extract the raw weight matrix from a fitted model object. Returns a
#' P x D numeric matrix (features x discriminant directions).
#'
#' @param object A fitted model object (e.g., from \code{sda}, \code{glmnet}).
#' @param ... Additional arguments passed to methods.
#' @return A numeric matrix of dimension P x D.
#' @examples
#' \donttest{
#'   ds <- gen_sample_dataset(c(5,5,5), 20, nlevels=2)
#'   mdl <- load_model("sda_notune")
#'   mspec <- mvpa_model(mdl, ds$dataset, ds$design, "classification")
#'   vox <- which(ds$dataset$mask > 0)
#'   X <- neuroim2::series(ds$dataset$train_data, vox)
#'   fit <- train_model(mspec, X, ds$design$y_train, indices=vox)
#'   w <- extract_weights(fit)
#' }
#' @export
extract_weights <- function(object, ...) UseMethod("extract_weights")

#' Extract Full Feature Matrix from a Dataset
#'
#' Returns the full T x P feature matrix (observations x features) for
#' whole-brain analyses. For clustered datasets this is the parcel time-series;
#' for image datasets this is the voxel series under the mask.
#'
#' @param dataset An \code{mvpa_dataset} or a plain matrix.
#' @param ... Additional arguments passed to methods.
#' @return A numeric matrix of dimension T x P.
#' @examples
#' ds <- gen_sample_dataset(c(5,5,5), 20)
#' X <- get_feature_matrix(ds$dataset)
#' dim(X)
#' @export
get_feature_matrix <- function(dataset, ...) UseMethod("get_feature_matrix")

#' Region Importance via Random Subset Comparison
#'
#' Assess each feature's (region/parcel/voxel) contribution to model performance
#' by comparing cross-validated accuracy when the feature is included vs. excluded
#' across many random feature subsets.
#'
#' @details
#' \strong{Interpretation.}
#' \code{region_importance} measures each feature's \emph{marginal predictive
#' contribution}: how much cross-validated performance improves, on average,
#' when the feature is included versus excluded across random feature subsets.
#' This is an approximation of Shapley values and is purely a \strong{decoding}
#' (backward) measure -- it never examines model weights.
#'
#' \strong{Haufe et al. (2014) considerations.}
#' Because this method works through out-of-sample performance deltas rather
#' than weight interpretation, it sidesteps the core failure mode described by
#' Haufe et al. (2014), where backward-model weights are misinterpreted as
#' activation patterns.
#'
#' However, \emph{suppressor variables} -- features that improve classification
#' by cancelling correlated noise rather than carrying signal -- will receive
#' positive importance, since they genuinely improve generalization.
#' This is correct from a decoding perspective ("which features help classify?")
#' but potentially misleading from a neuroscience perspective ("where does the
#' signal originate?").
#'
#' \strong{Complementary methods.}
#' For forward-model (activation-pattern) interpretation of linear models,
#' use \code{\link{haufe_importance}} directly or rely on the
#' \code{importance_vector} returned by \code{\link{run_global}}, which uses
#' \code{\link{model_importance}} internally.  For non-linear models such as
#' random forests, \code{region_importance} is the recommended approach.
#'
#' @param model_spec An \code{mvpa_model} specification.
#' @param ... Additional arguments passed to methods.
#' @return A \code{region_importance_result} object.
#' @examples
#' \donttest{
#'   ds <- gen_sample_dataset(c(5,5,5), 40, nlevels=2, blocks=3)
#'   cval <- blocked_cross_validation(ds$design$block_var)
#'   mdl <- load_model("sda_notune")
#'   mspec <- mvpa_model(mdl, ds$dataset, ds$design,
#'     "classification", crossval=cval)
#'   imp <- region_importance(mspec, n_regions=5, n_subsets=10)
#' }
#' @export
region_importance <- function(model_spec, ...) UseMethod("region_importance")

#' Per-Feature Model Importance
#'
#' Generic function that extracts a per-feature importance vector from a fitted
#' model object.
#'
#' @details
#' The importance measure returned depends on the model class:
#'
#' \describe{
#'   \item{\strong{SDA / glmnet (linear models)}}{Computes Haufe et al. (2014)
#'     \emph{forward-model} activation patterns: A = Sigma_x * W * inv(W' Sigma_x W).
#'     The returned vector is the L2 norm (or custom \code{summary_fun}) of the
#'     rows of A.
#'     This is the recommended importance measure for neuroscience interpretation
#'     because it reflects where the neural signal originates, not merely which
#'     features carry discriminative weight.}
#'   \item{\strong{randomForest}}{Returns MeanDecreaseGini (or MeanDecreaseAccuracy
#'     when available).
#'     This is a \strong{backward} (decoding) measure and is \strong{not}
#'     Haufe-safe: suppressor variables that reduce node impurity without
#'     carrying neural signal will receive high importance.
#'     Use \code{\link{region_importance}} for a slower but more robust
#'     backward measure that is bounded by out-of-sample performance, or
#'     restrict neuroscience interpretation to linear models with Haufe-based
#'     importance.}
#'   \item{\strong{default}}{Returns \code{NULL}, signaling that no importance
#'     is available for the model class.}
#' }
#'
#' @param object A fitted model object (e.g., from \code{sda}, \code{glmnet},
#'   \code{randomForest}).
#' @param X_train The training data matrix (T x P) used to fit \code{object}.
#'   Required for Haufe-based methods; ignored by tree-based methods.
#' @param ... Additional arguments passed to methods.
#' @return A numeric vector of length P with per-feature importance scores,
#'   or \code{NULL} if importance is not available for this model class.
#' @examples
#' \donttest{
#'   ds <- gen_sample_dataset(c(5,5,5), 40, nlevels=2)
#'   mdl <- load_model("sda_notune")
#'   mspec <- mvpa_model(mdl, ds$dataset, ds$design, "classification")
#'   vox <- which(ds$dataset$mask > 0)
#'   X <- neuroim2::series(ds$dataset$train_data, vox)
#'   fit <- train_model(mspec, X, ds$design$y_train, indices=vox)
#'   imp <- model_importance(fit, X)
#' }
#' @export
model_importance <- function(object, X_train, ...) UseMethod("model_importance")

#' Run Global (Whole-Brain) MVPA Analysis
#'
#' Train a single classifier on ALL features (parcels or voxels) with
#' cross-validation, and compute per-feature importance via Haufe et al.
#' (2014) activation patterns.
#'
#' @section Architecture TODO:
#' \code{run_global} currently uses a dedicated global CV/training pipeline
#' rather than dispatching through \code{\link{fit_roi}}. This is intentional
#' for now; future cleanup may unify global and ROI fitting interfaces.
#'
#' @param model_spec An \code{mvpa_model} specification.
#' @param ... Additional arguments passed to methods.
#' @return A \code{global_mvpa_result} object.
#' @examples
#' \donttest{
#'   ds <- gen_sample_dataset(c(5,5,5), 40, nlevels=2, blocks=3)
#'   cval <- blocked_cross_validation(ds$design$block_var)
#'   mdl <- load_model("sda_notune")
#'   mspec <- mvpa_model(mdl, ds$dataset, ds$design,
#'     "classification", crossval=cval)
#'   result <- run_global(mspec)
#' }
#' @export
run_global <- function(model_spec, ...) UseMethod("run_global")
