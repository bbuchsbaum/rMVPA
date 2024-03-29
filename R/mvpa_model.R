

#' @keywords internal
get_multiclass_perf <- function(split_list=NULL, class_metrics=TRUE) {
  function(result) {
    performance(result, split_list, class_metrics)
  }
}

#' @keywords internal
get_binary_perf <- function(split_list=NULL) {
  function(result) {
    performance(result, split_list)
  }
}

#' @keywords internal
get_regression_perf <- function(split_list=NULL) {

  function(result) {
    performance(result, split_list)
  }
}


#' @keywords internal
get_custom_perf <- function(fun, split_list) {
  function(result) {
    custom_performance(result, fun, split_list)
  }
}


#' @keywords internal
compute_performance.mvpa_model <- function(obj, result) {
  obj$performance(result)
}

#' @export
tune_grid.mvpa_model <- function(obj, x,y,len=1) {
  if (is.null(obj$tune_grid)) {
    obj$model$grid(x,y,len)
  } else {
    obj$tune_grid
  }
}


#' @export
select_features.mvpa_model <- function(obj, X, Y,...) {
  if (!is.null(obj$feature_selector)) {
    select_features(obj$feature_selector, X, Y,...)
  } else {
    ## todo check 'length' function in ROIVector
    rep(TRUE, length(X))
  }
}


#' @export
crossval_samples.mvpa_model <- function(obj,data,y,...) { 
  crossval_samples(obj$crossval,data,y) 
}

#' @export
has_test_set.mvpa_model <- function(obj) {
  !is.null(obj$design$y_test) 
}

#' @export
y_train.mvpa_model <- function(obj) y_train(obj$design)


#' @export
y_test.mvpa_model <- function(obj) y_test(obj$design)


#' Create an MVPA Model
#'
#' Create an MVPA model based on a caret-based classification or regression model.
#'
#' @param model A caret-based classification or regression model.
#' @param dataset An `mvpa_dataset` instance.
#' @param design An `mvpa_design` instance.
#' @param model_type A character string indicating the problem type: "classification" or "regression".
#' @param crossval An optional `cross_validation` instance.
#' @param feature_selector An optional `feature_selector` instance.
#' @param tune_grid An optional parameter tuning grid as a `data.frame`.
#' @param tune_reps The number of replications used during parameter tuning. Only relevant if `tune_grid` is supplied.
#' @param performance An optional custom function for computing performance metrics.
#' @param class_metrics A logical flag indicating whether to compute performance metrics for each class.
#'
#' @export
#'
#' @details
#'
#' If `performance` is supplied, it must be a function that takes one argument and returns a named list of scalar values. 
#' The argument the function takes is a class deriving from `classification_result` appropriate for the problem at hand.
#' See example below.
#'
#' @examples
#'
#' mod <- load_model("sda")
#' traindat <- neuroim2::NeuroVec(array(rnorm(6*6*6*100), c(6,6,6,100)), neuroim2::NeuroSpace(c(6,6,6,100)))
#' mask <- neuroim2::LogicalNeuroVol(array(rnorm(6*6*6)>-.2, c(6,6,6)), neuroim2::NeuroSpace(c(6,6,6)))
#'
#' mvdset <- mvpa_dataset(traindat,mask=mask)
#' design <- data.frame(fac=rep(letters[1:4], 25), block=rep(1:10, each=10))
#' cval <- blocked_cross_validation(design$block)
#' mvdes <- mvpa_design(design, ~ fac, block_var=~block)
#'
#' custom_perf <- function(result) {
#'   c(accuracy=sum(result$observed == result$predicted)/length(result$observed))
#' }
#' mvpmod <- mvpa_model(mod, dataset=mvdset, design=mvdes, crossval=cval, performance=custom_perf)
#' ret <- run_searchlight(mvpmod)
#' stopifnot("accuracy" %in% names(ret))
mvpa_model <- function(model, 
                       dataset,
                       design,
                       model_type=c("classification", "regression"), 
                       crossval=NULL, 
                       feature_selector=NULL, 
                       tune_grid=NULL, 
                       tune_reps=15,
                       performance=NULL,
                       class_metrics=TRUE) {
  
  assert_that(!is.null(model$fit))
  
  assert_that(inherits(design, "mvpa_design"))
  assert_that(inherits(dataset, "mvpa_dataset"))
  assert_that(is.logical(class_metrics))
  
  if (is.null(dataset$test_data) && !is.null(design$y_test)) {
    stop("mvpa_model: design has `y_test` dataset must have `test_data`")
  }
  
  if (!is.null(dataset$test_data) && is.null(design$y_test)) {
    stop("mvpa_model: if dataset has `test_data` design must have `y_test`")
  }
  
  perf <- if (!is.null(performance) && is.function(performance)) {
    #assert_that(is.function(performance)) 
    get_custom_perf(performance, design$split_groups)
  } else if (is.numeric(design$y_train)) {
    get_regression_perf(design$split_groups)
  } else if (length(levels(design$y_train)) > 2) {
    get_multiclass_perf(design$split_groups, class_metrics)
  } else if (length(levels(design$y_train)) == 2) {
    get_binary_perf(design$split_groups)
  } else {
    flog.info("class of design$ytrain: %s", class(design$ytrain))
    stop("performance method not found")
  }
  

  model_type <- match.arg(model_type)
  
  if (is.null(crossval) && !is.null(design$block_var)) {
    crossval <- blocked_cross_validation(design$block_var)
  }
  
  assertthat::assert_that(!is.null(crossval))
  
  ## TODO check that crossval is compatible with design
  
  ## TODO check that mvpa_design is compatible with mvpa_dataset (n training obs == length(y_train))
  #assert_that(length(y_train(design)) == )
  ret <- list(model=model,
              dataset=dataset,
              design=design,
              model_type=model_type,
              model_name=model$label,
              tune_grid=tune_grid,
              tune_reps=tune_reps,
              feature_selector=feature_selector,
              crossval=crossval,
              performance=perf)
  
  class(ret) <- "mvpa_model"
  ret
  
}



#' @export
#' @method print mvpa_model
print.mvpa_model <- function(x,...) {
  cat("mvpa_model object. \n")
  cat("model: ", x$model_name, "\n")
  cat("model type: ", x$model_type, "\n")
  if (!is.null(x$tune_grid)) {
    cat("tune_reps: ", x$tune_reps, "\n")
    cat("tune_grid: ", "\n")
    print(x$tune_grid)
  }
  
  print(x$crossval)
  print(x$dataset)
  print(x$design)
}
  
  



