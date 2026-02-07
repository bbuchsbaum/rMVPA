#' Register a Custom MVPA Model
#'
#' Adds a user-defined model specification to the rMVPA model registry (MVPAModels).
#'
#' @param name A character string, the unique name for the model.
#' @param model_spec A list containing the model specification. It must include
#'   elements: `type` ("Classification" or "Regression"), `library` (character vector
#'   of required packages for the *model itself*, not for rMVPA's wrappers), 
#'   `label` (character, usually same as name), `parameters`
#'   (data.frame of tunable parameters: parameter, class, label), `grid` (function to
#'   generate tuning grid, takes x, y, len args), `fit` (function), `predict` (function), 
#'   and `prob` (function for classification, takes modelFit, newdata; should return matrix/df with colnames as levels).
#' @export
#' @examples
#' \dontrun{
#' # Example of how a user might define an e1071 SVM spec
#' my_svm_spec <- list(
#'   type = "Classification", library = "e1071", label = "my_svm",
#'   parameters = data.frame(parameter = "cost", class = "numeric", label = "Cost (C)"),
#'   # grid should return a data.frame with columns matching 'parameter' names in 'parameters'
#'   grid = function(x, y, len = NULL) { 
#'      data.frame(cost = if (is.null(len) || len == 1) 1 else 10^seq(-2, 2, length.out = len))
#'   },
#'   # fit function receives: x, y, wts (weights), param (current params from grid), 
#'   # lev (levels of y), last (unused), weights (unused), classProbs (unused by e1071::svm)
#'   fit = function(x, y, wts, param, lev, last, weights, classProbs, ...) {
#'      e1071::svm(x, y, cost = param$cost, probability = TRUE, ...) # Ensure probability=TRUE for prob
#'   },
#'   # predict function receives: modelFit (output of $fit), newdata
#'   predict = function(modelFit, newdata, ...) {
#'      predict(modelFit, newdata, ...)
#'   },
#'   # prob function receives: modelFit, newdata
#'   # Should return a matrix/df with columns named as in levels(y)
#'   prob = function(modelFit, newdata, ...) {
#'     pred_obj <- predict(modelFit, newdata, probability = TRUE)
#'     attr(pred_obj, "probabilities") 
#'   }
#' )
#' register_mvpa_model("my_svm", my_svm_spec)
#' # Now load_model("my_svm") would work.
#' }
register_mvpa_model <- function(name, model_spec) {
  required_elements <- c("type", "library", "label", "parameters", "grid", "fit", "predict", "prob")
  if (!all(required_elements %in% names(model_spec))) {
    stop("model_spec is missing one or more required elements: ", 
         paste(setdiff(required_elements, names(model_spec)), collapse=", "))
  }
  if (!is.data.frame(model_spec$parameters) || 
      !all(c("parameter", "class", "label") %in% names(model_spec$parameters))) {
    stop("'model_spec$parameters' must be a data.frame with columns: parameter, class, label.")
  }
  MVPAModels[[name]] <- model_spec
  invisible(NULL)
}