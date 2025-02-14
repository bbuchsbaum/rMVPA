#' @keywords internal
"_PACKAGE"

#' Extending the MVPA Framework: Creating a New Analysis Type
#'
#' This documentation describes how to implement a **new MVPA analysis type** within the package. 
#' The MVPA framework here is designed to be extensible via S3 generics. By providing the necessary 
#' methods and classes, you can integrate custom analyses (beyond standard classification, regression, or RSA).
#'
#' @section Overview:
#' 1. **Define a Model Specification**  
#'    - Create an S3 class inheriting from \code{model_spec}, e.g. \code{"myanalysis_model"}.
#'    - This class should hold all the **parameters** or **configuration** needed for your analysis 
#'      (e.g., hyperparameters, references, design objects).
#'    - Construct it using a function like \code{create_model_spec("myanalysis_model", ...)}.
#'
#' 2. **Implement Required Generics**  
#'    - \strong{\code{train_model.myanalysis_model(obj, train_dat, ...)}:}  
#'      Defines how to \emph{"train"} or \emph{execute} your analysis on a subset of data 
#'      (e.g., ROI voxel intensities).  
#'      - This might involve fitting a specialized model, computing a metric, or doing an 
#'        unsupervised transform.  
#'      - Return a structured object (often a named list or model fit) that can be further 
#'        handled by the framework.
#'    - \strong{\code{process_roi.myanalysis_model(mod_spec, roi, rnum, ...)} (Optional):}  
#'      If you need custom \emph{per-ROI} logic beyond the default pipeline, implement a 
#'      \code{process_roi} method for your class. This method is called automatically by
#'      \code{run_regional(...)} for each region.  
#'      - By default, the framework calls \code{train_model} on the ROI's data, but if you 
#'        want to skip or alter that flow, provide a custom \code{process_roi} method.
#'    - (Optional) \strong{\code{predict_model.myanalysis_model(...)}} or \strong{\code{evaluate_model.myanalysis_model(...)}}:  
#'      If your analysis involves a separate test step or specialized evaluation, you can 
#'      define these S3 methods.  
#'    - (Optional) \strong{\code{run_future.myanalysis_model(...)}}, \strong{\code{merge_results.myanalysis_model(...)}}, etc.:  
#'      If you need custom behavior for parallelization or merging partial results, override 
#'      these generics as well.
#'
#' 3. **Use `run_regional()` or `run_searchlight()`**  
#'    - Once your new model class is defined, you can pass an object of that class to 
#'      \code{run_regional(...)} or \code{run_searchlight(...)} just like a built-in analysis.  
#'    - The framework automatically dispatches calls to your \code{train_model.myanalysis_model} 
#'      (and optionally your \code{process_roi.myanalysis_model}), returning a result object 
#'      consistent with other MVPA analyses.
#'
#' @section Why These Pieces Are Needed:
#' - **Model Specification**: Provides the blueprint for how the pipeline should handle your analysis 
#'   (what data you need, how you store intermediate results, etc.).
#' - **Train Model Generic**: Tells the pipeline what to \emph{do} with each chunk of data 
#'   (ROI subset, searchlight voxel set, cross-validation fold). 
#'   This is the core of your analysis.
#' - **\code{process_roi}**: If your analysis requires more complex per-region logic than the default 
#'   pipeline (e.g., skipping certain steps, altering the data structure), implementing 
#'   \code{process_roi.myanalysis_model} gives you full control.
#'
#' @section Example Code Skeleton:
#' \preformatted{
#' # Minimal template for a new analysis model
#'
#' myanalysis_model <- function(dataset, design, param1=NULL, param2=5, ...) {
#'   # Create an object of class "myanalysis_model"
#'   create_model_spec(
#'     "myanalysis_model",
#'     dataset  = dataset,
#'     design   = design,
#'     param1   = param1,
#'     param2   = param2,
#'     ...
#'   )
#' }
#'
#' # S3 method to "train" your analysis
#' train_model.myanalysis_model <- function(obj, train_dat, y, indices, ...) {
#'   # train_dat = subset of voxel intensities or features
#'   # y = optional labels/response if relevant
#'   # indices = location info
#'
#'   # 1) Do your computations (model fitting, metrics, transforms, etc.)
#'   fit_object <- your_method(train_dat, some_params=obj$param1)
#'
#'   # 2) Return a structure. The pipeline may store or evaluate it further.
#'   list(
#'     fit = fit_object,
#'     # any other stuff you want to keep
#'   )
#' }
#'
#' # Optional custom ROI processor
#' process_roi.myanalysis_model <- function(mod_spec, roi, rnum, ...) {
#'   # If you want to skip or augment the default cross-validation / train_model approach
#'   # in run_regional, define your logic here. 
#'   # Otherwise, the default calls train_model(...) on the region's data.
#' }
#'
#' # You can now call:
#' # new_model <- myanalysis_model(dataset, design, param1="xyz")
#' # results_regional <- run_regional(new_model, region_mask)
#' # results_searchlight <- run_searchlight(new_model, radius=8, method="standard")
#' }
#'
#' @name new-analysis-overview
NULL