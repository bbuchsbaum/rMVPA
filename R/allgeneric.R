


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
#' ROI <- neuroim2::ROIVec(neuroim2::NeuroSpace(c(10,10,10)), coords=coords, matrix(rnorm(100*3), 100, 3))
#' Y <- factor(rep(c("a", "b"), each=50))
#' featureMask <- select_features(fsel, neuroim2::values(ROI), Y)
#' sum(featureMask) == 2
#'
#' fsel2 <- feature_selector("FTest", "top_p", .1)
#' featureMask <- select_features(fsel2, neuroim2::values(ROI), Y)
#' 
#' @export
select_features <- function(obj, X, Y, ...) {
  UseMethod("select_features")
}



#' Train Model
#'
#' Train a classification or regression model.
#'
#' @param obj The model specification object.
#' @param ... Additional arguments to be passed to the method-specific function.
#'
#' @export
#' 
train_model <- function(obj,...) {
  UseMethod("train_model")
}



#' Training Labels/Response Extraction
#'
#' Extract the training labels or response variable from an object.
#'
#' @param obj The object from which to extract the training response variable.
#' @export
#' 
y_train <- function(obj) {
  UseMethod("y_train")
}


#' Test Labels/Response Extraction
#'
#' Extract the test labels or response variable from an object.
#'
#' @param obj The object from which to extract the test response variable.
#' @export
#' 
y_test <- function(obj) {
  UseMethod("y_test")
}

#' Test Design Extraction
#'
#' Return the design table associated with the test set from an object.
#'
#' @param obj The object from which to extract the test design table.
#' @export
#' 
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
#' @param lev Unused.
#' @param last Unused.
#' @param classProbs Unused.
#' @param ... Additional arguments to be passed to the method-specific function.
#' 
#' @export
fit_model <- function(obj, roi_x, y, wts, param, lev, last, classProbs, ...) {
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
#' @export
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


#' performance
#' 
#' Compute appropriate performance metrics such as accuracy/AUC/RMSE from a classification/regression result
#' 
#' @param x the result to evaluate performance of
#' @param ... extra args
#' @export
#' @rdname performance-methods
performance <- function(x,...) {
  UseMethod("performance")
}


#' compute_performance
#' 
#' Delegate calculation of performance metrics
#' 
#' @param obj the object
#' @param result the 'result' to evaluate
compute_performance <- function(obj, result) {
  UseMethod("compute_performance")
}

#' merge_results
#' 
#' merge two classification/regression results
#' 
#' @param x the first result
#' @param ... the rest
#' @export
merge_results <- function(x, ...) {
  UseMethod("merge_results")
}


#' get_samples
#' 
#' @param obj the object
#' @param vox_list the list of voxel/index sets
#' @export
get_samples <- function(obj, vox_list) {
  UseMethod("get_samples")
}

#' data_sample
#' 
#' @param obj the object
#' @param vox the voxel indices
#' @export
data_sample <- function(obj, vox) {
  UseMethod("data_sample")
}


#' extract_sample
#' 
#' @param obj the object
#' @param ... extra args
extract_sample <- function(obj,...) {
  UseMethod("extract_sample")
}


#' Convert an object to an ROIVolume or ROISurface object
#'
#' The as_roi function is a generic function that converts the provided object
#' into an ROIVolume or ROISurface object, given the associated data object and
#' any additional arguments.
#'
#' @param obj The object to be converted.
#' @param data The associated data object used in the conversion process.
#' @param ... Additional arguments to be passed to the function.
#' @return An ROIVolume or ROISurface object.
#' @keywords internal
as_roi <- function(obj, data, ...) {
  UseMethod("as_roi")
}


#' get_searchlight
#' 
#' generate a searchlight iterator appropriate for a given input dataset (e.g. volumetric or surface).
#' 
#' @param obj the object
#' @param ... extra args
#' @export
get_searchlight <- function(obj, ...) {
  UseMethod("get_searchlight")
}




#' Wrap output values into a desired format
#'
#' The wrap_output function is a generic function that takes an object,
#' values to be wrapped, and additional arguments to create a wrapped
#' output based on the provided object and values.
#'
#' @param obj The object used to determine the appropriate method for wrapping the output.
#' @param vals The values to be wrapped.
#' @param ... Additional arguments to be passed to the function.
#' @return A wrapped output based on the provided object and values.
#' @keywords internal
wrap_output <- function(obj, vals, ...) {
  UseMethod("wrap_output")
}

#' Merge predictions from multiple models
#'
#' Combine predictions from several models applied to the same test set.
#'
#' @param obj1 The first object containing predictions.
#' @param rest The rest of the objects containing predictions.
#' @param ... Additional arguments to be passed to the function.
#' @return A combined object containing merged predictions from multiple models.
#' @export
merge_predictions <- function(obj1, rest, ...) {
  UseMethod("merge_predictions")
}


#' sub_result
#' 
#' exctract a row-wise subset of a classification/regression result object.
#' 
#' @param x the result object to subset
#' @param indices the row indices
sub_result <- function(x, indices) {
  UseMethod("sub_result")
}

#' nobs
#' 
#' get number of observations 
#' 
#' @param x the object
#' @export
nobs <- function(x) {
  UseMethod("nobs")
}

#' Probability of observed class
#'
#' Extract the predicted probability for the observed class.
#'
#' @param x The object from which to extract the predicted probability for the observed class.
#' @return A vector of predicted probabilities for the observed class.
#' @seealso \code{\link{prob_observed.binary_classification_result}}, \code{\link{prob_observed.multiway_classification_result}}
#' @export
prob_observed <- function(x) {
  UseMethod("prob_observed")
}



#' Number of response categories
#'
#' Get the number of response categories or category levels for the dependent variable.
#'
#' @param x The object from which to extract the number of response categories.
#' @return The number of response categories or category levels for the dependent variable.
#' @export
nresponses <- function(x) {
  UseMethod("nresponses")
}


.name_repair = ~ vctrs::vec_as_names(..., repair = "unique", quiet = TRUE)

#' run_searchlight
#'
#' Execute a searchlight analysis.
#'
#' @param model_spec A \code{mvpa_model} instance containing the model specifications.
#' @param radius The searchlight radius in millimeters.
#' @param method The type of searchlight, either 'randomized' or 'standard'.
#' @param niter The number of searchlight iterations (used only for the 'randomized' method).
#' @param ... Extra arguments passed to the specific searchlight methods.
#' @return A named list of \code{NeuroVol} objects, where each element contains a performance metric (e.g. AUC) at every voxel location.
#' @seealso \code{\link{run_searchlight.randomized}}, \code{\link{run_searchlight.standard}}
#' @examples
#' # TODO: Add an example
#' @export
run_searchlight <- function(model_spec, radius, method, niter,...) {
  UseMethod("run_searchlight")
}

#' Region of interest based MVPA analysis
#'
#' Run a separate MVPA analysis for multiple disjoint regions of interest.
#'
#' @param model_spec A \code{mvpa_model} instance containing the model specifications.
#' @param region_mask A \code{NeuroVol} or \code{NeuroSurface} object where each region is identified by a unique integer.
#'        Every non-zero set of positive integers will be used to define a set of voxels for classification analysis.
#' @param ... Extra arguments passed to the specific regional analysis methods.
#' @return A named list of results for each region of interest.
#' @examples
#' # TODO: Add an example
#' @export      
run_regional <- function(model_spec, region_mask, ...) {
  UseMethod("run_regional")
}


#' crossval_samples
#'
#' A generic function that applies a cross-validation scheme to split the data into training and testing sets.
#' It is used along with cross-validation control objects and S3 implementation functions to perform the cross-validation process.
#'
#' @param obj A cross-validation control object.
#' @param data A data frame containing the predictor variables.
#' @param y A vector containing the response variable.
#' @param ... Extra arguments passed to the specific cross-validation methods.
#' @return A tibble containing the training and testing sets for each fold, as well as the response variables for both sets.
#' @seealso \code{\link{crossval_samples.sequential_blocked_cross_validation}},
#'   \code{\link{crossval_samples.kfold_cross_validation}},
#'   \code{\link{crossval_samples.blocked_cross_validation}},
#'   \code{\link{crossval_samples.bootstrap_blocked_cross_validation}},
#'   \code{\link{crossval_samples.custom_cross_validation}},
#'   \code{\link{crossval_samples.twofold_blocked_cross_validation}}
#' @examples
#' # Example with k-fold cross-validation
#' cval <- kfold_cross_validation(len=100, nfolds=10)
#' samples <- crossval_samples(cval, data=as.data.frame(matrix(rnorm(100*10), 100, 10)), y=rep(letters[1:5],20))
#' stopifnot(nrow(samples) == 10)
#' @export
crossval_samples <- function(obj, data, y,...) { 
  UseMethod("crossval_samples") 
}




