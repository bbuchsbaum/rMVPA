#' @keywords internal
"_PACKAGE"

#' Extending rMVPA: Plugin Analysis Models
#'
#' rMVPA can be extended from external packages by defining S3 model classes
#' that inherit from \code{model_spec}. The preferred extension path is the
#' modern \code{fit_roi}/\code{output_schema} contract.
#'
#' @section Recommended Contract:
#' 1. **Constructor**  
#'    Build a model object with \code{\link{create_model_spec}}, assigning
#'    your model class name (for example \code{"my_model"}).
#'
#' 2. **\code{fit_roi.my_model} (required)**  
#'    Implement per-ROI work and return \code{\link{roi_result}}.
#'
#' 3. **\code{output_schema.my_model} (recommended)**  
#'    Return a named list with values \code{"scalar"} or \code{"vector[N]"}.
#'    This enables schema-driven map construction in searchlight analyses.
#'
#' 4. **Optional overrides**  
#'    - \code{y_train.my_model} / \code{y_test.my_model} if labels are not
#'      available from the attached design object.  
#'    - \code{strip_dataset.my_model} for advanced parallel memory control.
#'
#' @section Plugin Test Helpers:
#' - \code{\link{mock_roi_data}} builds lightweight ROI payloads for
#'   unit tests.
#' - \code{\link{mock_context}} builds the standard \code{fit_roi}
#'   context list.
#' - \code{\link{validate_plugin_model}} runs a contract check on one ROI
#'   and validates schema/metric agreement.
#' - \code{\link{validate_model_spec}} performs a plugin-readiness lint with
#'   optional dry-run execution.
#'
#' @section Two Registry Types:
#' rMVPA exposes two different extension registries:
#' \describe{
#'   \item{Classifier registry}{Use \code{\link{register_mvpa_model}} to add
#'   a new low-level estimator (fit/predict/prob) that plugs into the
#'   built-in \code{mvpa_model} analysis type.}
#'   \item{Analysis-type registry}{Define a new S3 model class (via
#'   \code{\link{create_model_spec}} + \code{fit_roi.<class>}) when you need
#'   a new per-ROI analysis workflow, not just a new estimator backend.}
#' }
#'
#' @section Notes:
#' - \code{\link{process_roi.default}} automatically prefers \code{fit_roi}
#'   when a method exists for your class.
#' - Legacy \code{train_model}/\code{process_roi} extension paths are still
#'   available for backward compatibility, but new plugins should use
#'   \code{fit_roi}.
#'
#' @section Minimal Example:
#' \preformatted{
#' my_model <- function(dataset, design, alpha = 0.1, ...) {
#'   create_model_spec(
#'     "my_model",
#'     dataset = dataset,
#'     design = design,
#'     alpha = alpha,
#'     ...
#'   )
#' }
#'
#' fit_roi.my_model <- function(model, roi_data, context, ...) {
#'   score <- mean(roi_data$train_data)
#'   roi_result(
#'     metrics = c(score = score),
#'     indices = roi_data$indices,
#'     id = context$id
#'   )
#' }
#'
#' output_schema.my_model <- function(model) {
#'   list(score = "scalar")
#' }
#'
#' # optional only when you do not rely on design defaults:
#' # y_train.my_model <- function(obj) y_train(obj$design)
#' }
#'
#' @name new-analysis-overview
NULL
