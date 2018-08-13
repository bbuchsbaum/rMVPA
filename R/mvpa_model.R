
get_multiclass_perf <- function(split_list=NULL, class_metrics=TRUE) {
  function(result) {
    performance(result, split_list, class_metrics)
  }
}

get_binary_perf <- function(split_list=NULL) {
  function(result) {
    performance(result, split_list)
  }
}

get_regression_perf <- function(split_list=NULL) {

  function(result) {
    performance(result, split_list)
  }
}

get_custom_perf <- function(fun, split_list) {
  function(result) {
    custom_performance(result, fun, split_list)
  }
}


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


#' mvpa_model
#' 
#' @param model a caret-based classification or regression model
#' @param dataset a \code{mvpa_dataset} instance
#' @param design a \code{mvpa_design} instance
#' @param model_type a \code{character} string indicating problem type: 'classification' or 
#' regression'
#' @param crossval a \code{cross_validation} instance
#' @param feature_selector an optional \code{feature_selector} instance
#' @param tune_grid an optional parameter tuning grid as a \code{data.frame}
#' @param tune_reps the number of replications used during paramter tuning. Only relevant if \code{tune_grid} supplied.
#' @param performance an optional custom function for computing performance metrics.
#' @param class_metrics \code{logical} flag indicating whether to compute performance metrics for each class.
#' @export
mvpa_model <- function(model, 
                       dataset,
                       design,
                       model_type=c("classification", "regression"), 
                       crossval, 
                       feature_selector=NULL, 
                       tune_grid=NULL, 
                       tune_reps=5,
                       performance=NULL,
                       class_metrics=TRUE) {
  
  assert_that(!is.null(model$fit))
  
  assert_that(inherits(design, "mvpa_design"))
  assert_that(inherits(dataset, "mvpa_dataset"))
  assert_that(is.logical(class_metrics))
  
  perf <- if (!is.null(performance)) {
    assert_that(is.function(performance)) 
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
print.mvpa_model <- function(x) {
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
  
  



