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
#' @param ... Additional arguments to be passed to the method-specific function.
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
#' @noRd
run_future <- function(obj, frame, processor, ...) {
  UseMethod("run_future")
}

#' Process ROI
#'
#' Process a region of interest (ROI) with possible cross-validation.
#'
#' @param mod_spec The model specification object.
#' @param roi The region of interest data.
#' @param rnum A numeric or string identifier for the ROI.
#' @param ... Additional arguments passed to the method-specific function.
#' @keywords internal
#' @noRd
process_roi <- function(mod_spec, roi, rnum, ...) {
  UseMethod("process_roi")
}

#' Default Process ROI Method
#'
#' Default implementation for processing ROIs when no custom ROI processor is provided.
#'
#' @inheritParams process_roi
#' @keywords internal
#' @noRd
process_roi.default <- function(mod_spec, roi, rnum, center_global_id = NA, ...) {
  # Capture additional arguments to pass down
  dots <- list(...)
  if (!is.null(mod_spec$process_roi)) {
    # Pass center_global_id and dots to user's custom processor
    do.call(mod_spec$process_roi, c(list(mod_spec, roi, rnum, center_global_id = center_global_id), dots))
  } else if (has_test_set(mod_spec)) { # Changed from mod_spec to mod_spec$dataset
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

  xtrain <- tibble::as_tibble(neuroim2::values(roi$train_roi), .name_repair=.name_repair)
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
  train_result_obj <- try(do.call(train_model, c(list(mod_spec, xtrain, y = NULL, indices=ind, sl_info = sl_info), dots)))
  
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
#' @export
train_model <- function(obj,...) {
  UseMethod("train_model")
}

#' Training Labels/Response Extraction
#'
#' Extract the training labels or response variable from an object.
#'
#' @param obj The object from which to extract the training response variable.
#' @export
y_train <- function(obj) {
  UseMethod("y_train")
}

#' Test Labels/Response Extraction
#'
#' Extract the test labels or response variable from an object.
#'
#' @param obj The object from which to extract the test response variable.
#' @export
y_test <- function(obj) {
  UseMethod("y_test")
}

#' Test Design Extraction
#'
#' Return the design table associated with the test set from an object.
#'
#' @param obj The object from which to extract the test design table.
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
#' @export
fit_model <- function(obj, roi_x, y, wts, param, lev=NULL, last=FALSE, classProbs=FALSE, ...) {
  UseMethod("fit_model")
}

#' Tune Grid Extraction
#'
#' Extract the parameter grid to optimize for a model.
#'
#' @param obj The model object.
#' @param x The training data.
#' @param y The response vector.
#' @param len The number of elements in the tuning grid.
tune_grid <- function(obj, x,y,len) {
  UseMethod("tune_grid")
}

#' Test Set Availability
#'
#' Check if an object has a test set available.
#'
#' @param obj The object to check for a test set.
#' @export
has_test_set <- function(obj) {
  UseMethod("has_test_set")
}

#' Requires cross-validation to be performed
#' @param obj The model object.
has_crossval <- function(obj) {
  UseMethod("has_crossval")
}

#' @export
has_crossval.default <- function(obj) {
  FALSE
}

#' Compute Performance Metrics
#'
#' Compute performance metrics (accuracy, AUC, RMSE, etc.) for classification/regression results.
#'
#' @param x The classification/regression result object to evaluate.
#' @param ... Additional arguments passed to method-specific performance functions.
#'
#' @return A list of performance metrics.
#' @export
performance <- function(x,...) {
  UseMethod("performance")
}

#' Compute Performance for an Object
#'
#' Delegates calculation of performance metrics to the appropriate method.
#'
#' @param obj The input object.
#' @param result The classification/regression result object to evaluate.
#'
#' @return A list of performance metrics.
#' @export
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
#'
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
#' @keywords internal
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
#' @export
predict_model <- function(object, fit, newdata, ...) {
  UseMethod("predict_model")
}

#' @keywords internal
#' @noRd
.name_repair = ~ vctrs::vec_as_names(..., repair = "unique", quiet = TRUE)

#' Run Searchlight Analysis
#'
#' Execute a searchlight analysis using multivariate pattern analysis.
#'
#' @param model_spec A \code{mvpa_model} instance containing the model specifications
#' @param radius The searchlight radius in millimeters
#' @param method The type of searchlight, either 'randomized' or 'standard'
#' @param niter The number of searchlight iterations (used only for 'randomized' method)
#' @param ... Extra arguments passed to specific searchlight methods
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
#'   # Results contain performance metrics
#'   # Access them with results$performance
#' }
#'
#' @export
run_searchlight <- function(model_spec, radius, method = c("standard", "randomized"), niter = NULL, ...) {
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
#' @param ... Extra arguments passed to specific regional analysis methods (e.g., `return_fits`, `compute_performance`).
#'
#' @return A \code{regional_mvpa_result} object (list) containing:
#'   \item{performance_table}{A tibble of performance metrics for each region (if computed).}
#'   \item{prediction_table}{A tibble with detailed predictions for each observation/region (if generated).}
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
#'   region_mask <- NeuroVol(
#'     sample(1:5, size=length(dataset$dataset$mask), replace=TRUE),
#'     space(dataset$dataset$mask)
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
#'   head(results$performance)           # Performance metrics
#'   head(results$prediction_table)      # Predictions
#'   first_roi_fit <- results$fits[[1]]  # First ROI's fitted model
#' }
#'
#' @rdname run_regional-methods
#' @export
run_regional <- function(model_spec, region_mask, ...) {
  UseMethod("run_regional")
}

#' crossval_samples
#'
#' Apply a cross-validation scheme to split the data into training and testing sets.
#'
#' @param obj A cross-validation control object.
#' @param data A data frame containing the predictors.
#' @param y A vector containing the response variable.
#' @param ... Extra arguments passed to the specific cross-validation methods.
#'
#' @return A tibble containing training and testing sets for each fold.
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
#' @export
train_indices <- function(obj, fold_num, ...) {
  UseMethod("train_indices")
}



