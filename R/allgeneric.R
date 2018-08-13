


#' select_features
#' 
#' Given a \code{feature_selection} specification object and a dataset return the set of selected features as a binary vector.
#' 
#' @param obj the \code{feature_selection} object
#' @param X region of interest containing training features, a class of type \code{ROIVolume} or \code{ROISurface} or \code{matrix}
#' @param Y the dependent variable as a \code{factor} or \code{numeric} variable.
#' @param ... additional arguments
#' @return a \code{logical} vector indicating the columns of \code{X} matrix that were selected.
#' @examples 
#' 
#' fsel <- feature_selector("FTest", "top_k", 2)
#' coords <- rbind(c(1,1,1), c(2,2,2), c(3,3,3))
#' 
#' ROI <- neuroim::ROIVector(neuroim::BrainSpace(c(10,10,10)), coords=coords, matrix(rnorm(100*3), 100, 3))
#' Y <- factor(rep(c("a", "b"), each=50))
#' featureMask <- select_features(fsel, ROI@data, Y)
#' sum(featureMask) == 2
#' 
#' fsel2 <- feature_selector("FTest", "top_p", .1)
#' featureMask <- select_features(fsel2, ROI@data, Y)
#' 
#' @export
select_features <- function(obj, X, Y, ...) {
  UseMethod("select_features")
}



#' train_model
#' 
#' train a classification or regression model
#' @param obj the model specification
#' @param ... extra args
#' @export
train_model <- function(obj,...) {
  UseMethod("train_model")
}



#' y_train
#' 
#' extract the training labels/response
#' 
#' @param obj the object to extract training response variable from.
#' @export
y_train <- function(obj) {
  UseMethod("y_train")
}


#' y_test
#' 
#' extract the test labels/response.
#' @param obj the object to extract test response variable from.
#' @export
y_test <- function(obj) {
  UseMethod("y_test")
}

#' fit_model
#' 
#' fit a classification or regression model
#' 
#' @param obj a model fitting object
#' @param roi_x an ROI containing the training data
#' @param y the response vector
#' @param wts a set of case weights
#' @param param tuning parameters
#' @param lev unused
#' @param last unused
#' @param classProbs unused 
#' @param ... extra args
fit_model <- function(obj, roi_x, y, wts, param, lev, last, classProbs, ...) {
  UseMethod("fit_model")
}


#' tune_grid
#' 
#' extract the parameter grid to optimize.
#' 
#' @param obj the model object
#' @param x the training data
#' @param y the response vector
#' @param len the number of elements in the tuning grid
#' @export
tune_grid <- function(obj, x,y,len) {
  UseMethod("tune_grid")
}

#' has_test_set
#' 
#' @param obj the object
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


#' as_roi
#' 
#' convert to an \code{ROIVolume} or \code{ROISurface} object
#' 
#' @param obj the object to convert
#' @param ... extra args
as_roi <- function(obj,...) {
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




#' wrap_output
#' 
#' @param obj the object
#' @param vals the values to wrap
#' @param ... extra args
wrap_output <- function(obj, vals, ...) {
  UseMethod("wrap_output")
}

#' merge_predictions
#' 
#' combine predictions from several models applied to the same test set.
#' 
#' @param obj1 the first object
#' @param rest the rest of the objects
#' @param ... extra args
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

#' nresponses
#' 
#' get number of responses or category levels for the dependent variable
#' 
#' @param x the object
#' @export
nresponses <- function(x) {
  UseMethod("nresponses")
}


#' run_searchlight
#' 
#' execute a searchlight analysis
#' 
#' @param model_spec a \code{mvpa_model} instance.
#' @param radius the searchlight radus in millimeters.
#' @param method the type of searchlight (randomized or standard)
#' @param niter the number of searchlight iterations (used only for 'randomized' method)
#' @param ... extra args
#' @return a named list of \code{BrainVolume} objects, where each element contains a performance metric (e.g. AUC) at every voxel location.
#' @export
run_searchlight <- function(model_spec, radius, method, niter,...) {
  UseMethod("run_searchlight")
}


#' crossval_samples
#' 
#' @param obj the cross-validation control object
#' @param data the data to split up
#' @param y the response variable
#' @param ... extra args
#' @export
crossval_samples <- function(obj, data, y,...) { 
  UseMethod("crossval_samples") 
}




