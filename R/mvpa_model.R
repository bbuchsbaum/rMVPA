#' @keywords internal
#' @noRd
wrap_result <- function(result_table, design, fit=NULL) {
  
  observed <- y_test(design)
  
  ## It could happen that not all design rows are actually tested, which is why we find the unique set of test indices
  testind <- unique(sort(unlist(result_table$test_ind)))
  
  if (is.factor(observed)) {
    prob <- matrix(0, length(testind), length(levels(observed)))
    colnames(prob) <- levels(observed)
    
    for (i in seq_along(result_table$probs)) {
      p <- as.matrix(result_table$probs[[i]])
      tind <- match(result_table$test_ind[[i]], testind)
      prob[tind,] <- prob[tind,] + p
    }
    
    ## probs must sum to one, can divide by sum.
    prob <- t(apply(prob, 1, function(vals) vals / sum(vals)))
    maxid <- max.col(prob)
    pclass <- levels(observed)[maxid]
    
    ## storing observed, testind, test_design 
    classification_result(observed[testind], pclass, prob, testind=testind, design$test_design, fit)
  } else {
    
    testind <- unique(sort(unlist(result_table$test_ind)))
    preds <- numeric(length(testind))
    
    for (i in seq_along(result_table$preds)) {
      #tind <- result_table$test_ind[[i]]
      tind <- match(result_table$test_ind[[i]], testind)
      preds[tind] <- preds[tind] + result_table$preds[[i]]
    }
    
    ## TODO check me
    counts <- table(sort(unlist(result_table$test_ind)))
    # Ensure counts aligns with testind order and convert to numeric vector
    preds <- preds/as.numeric(counts[as.character(testind)])
    regression_result(observed, preds, testind=testind, test_design=design$test_design, fit)
  }
}



#' @rdname merge_results-methods
#' @export
merge_results.mvpa_model <- function(obj, result_set, indices, id, ...) {
  # Check if any errors occurred during the process
  if (any(result_set$error)) {
    emessage <- result_set$error_message[which(result_set$error)[1]]
    tibble::tibble(result=list(NULL), indices=list(indices), performance=list(NULL), error=TRUE, error_message=emessage)
  } else {
    # If no errors, wrap the result and compute performance if required
    cres <- if (obj$return_fit) {
      predictor <- weighted_model(result_set$fit)
      wrap_result(result_set, obj$design, predictor)
    } else {
      wrap_result(result_set, obj$design)
    }
    
    if (obj$compute_performance) {
      tibble::tibble(result=list(cres), indices=list(indices),
                     performance=list(compute_performance(obj, cres)),
                     id=id, error=FALSE, error_message="~")
    } else {
      tibble::tibble(result=list(cres), indices=list(indices),
                     performance=list(NULL),
                     id=id, error=FALSE, error_message="~")
    }
  }
 
}


#' @rdname format_result
#' @method format_result mvpa_model
#' @export
#' @importFrom stats predict
format_result.mvpa_model <- function(obj, result, error_message=NULL, context, ...) {
  if (!is.null(error_message)) {
    return(tibble::tibble(class=list(NULL), probs=list(NULL), y_true=list(context$ytest),
                          fit=list(NULL), error=TRUE, error_message=error_message))
  } else {
    pred <- predict(result, tibble::as_tibble(context$test, .name_repair=.name_repair), NULL)
    plist <- lapply(pred, list)
    plist$y_true <- list(context$ytest)
    plist$test_ind <- list(as.integer(context$test))
    plist$fit <- if (obj$return_fit) list(result) else list(NULL)
    plist$error <- FALSE
    plist$error_message <- "~"
    tibble::as_tibble(plist, .name_repair=.name_repair)
  }
}



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
#' @rdname compute_performance-methods
#' @export
compute_performance.mvpa_model <- function(obj, result) {
  obj$performance(result)
}

#' @export
tune_grid.mvpa_model <- function(obj, x,y,len=1) {
  if (is.null(obj$tune_grid)) {
    obj$tune_grid <- obj$model$grid(x,y,len)
    obj$tune_grid
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
#' @rdname crossval_samples
crossval_samples.mvpa_model <- function(obj,data,y,...) {
  crossval_samples(obj$crossval,data,y)
}

#' @export
has_test_set.mvpa_model <- function(obj) {
  # Use the stored flag rather than checking for existence of design$y_test
  isTRUE(obj$has_test_set)
}

#' @export
has_test_set.model_spec <- function(obj) {
  # Use the stored flag rather than checking for existence of design$y_test
  isTRUE(obj$has_test_set)
}



#' @rdname y_train-methods
#' @export
y_train.mvpa_model <- function(obj) y_train(obj$design)


#' @rdname y_test-methods
#' @export
y_test.mvpa_model <- function(obj) y_test(obj$design)


#' @param name A character string indicating the name of the model.
#' @param dataset An `mvpa_dataset` instance.
#' @param design An `mvpa_design` instance.
#' @param return_predictions A \code{logical} indicating whether to return row-wise predictions for each voxel set (defaults to TRUE).
#' @param tune_reps The number of replications used during parameter tuning. Only relevant if `tune_grid` is supplied.
#' @param ... Additional arguments to be passed to the model, including `has_test_set` flag.
#' @noRd
create_model_spec <- function(name, dataset, design, return_predictions=FALSE, 
                              compute_performance=FALSE, tune_reps=FALSE, ...) {
  # Automatically detect test set presence from dataset and design if not passed as parameter
  args <- list(...)
  if (is.null(args$has_test_set)) {
    args$has_test_set <- !is.null(dataset$test_data) && !is.null(design$y_test)
  }
  
  ret <- list(name=name, dataset=dataset, design=design, 
              return_predictions=return_predictions, compute_performance=compute_performance, 
              tune_reps=tune_reps)
  
  # Add the remaining arguments
  ret <- c(ret, args)
  
  class(ret) <- c(name, "model_spec", "list")
  ret
}

#' Create an MVPA Model
#'
#' Create an MVPA model based on a classification or regression model from the MVPAModels registry.
#'
#' @param model A character string naming a model from the MVPAModels registry, or a custom model specification list.
#' @param dataset An `mvpa_dataset` instance.
#' @param design An `mvpa_design` instance.
#' @param model_type A character string indicating the problem type: "classification" or "regression".
#' @param crossval An optional `cross_validation` instance.
#' @param feature_selector An optional `feature_selector` instance.
#' @param tune_grid An optional parameter tuning grid as a `data.frame`.
#' @param tune_reps The number of replications used during parameter tuning. Only relevant if `tune_grid` is supplied.
#' @param performance An optional custom function for computing performance metrics.
#' @param class_metrics A logical flag indicating whether to compute performance metrics for each class.
#' @param compute_performance A \code{logical} indicating whether to compute and store performance measures for each voxel set (defaults to TRUE).
#' @param return_predictions A \code{logical} indicating whether to return row-wise predictions for each voxel set (defaults to TRUE).
#' @param return_fits A \code{logical} indicating whether to return the model fit for each voxel set (defaults to FALSE).
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
#' arr_data <- array(rnorm(6*6*6*100), c(6,6,6,100))
#' sp <- neuroim2::NeuroSpace(c(6,6,6,100))
#' traindat <- neuroim2::NeuroVec(arr_data, sp)
#' mask <- neuroim2::LogicalNeuroVol(array(rnorm(6*6*6)>-.2, c(6,6,6)), neuroim2::NeuroSpace(c(6,6,6)))
#'
#' mvdset <- mvpa_dataset(traindat,mask=mask)
#' design <- data.frame(fac=rep(letters[1:4], 25), block=rep(1:10, each=10))
#' cval <- blocked_cross_validation(design$block)
#' mvdes <- mvpa_design(design, y_train = ~ fac, block_var = ~ block)
#'
#' custom_perf <- function(result) {
#'   c(accuracy=sum(result$observed == result$predicted)/length(result$observed))
#' }
#' mvpmod <- mvpa_model(mod, dataset=mvdset, design=mvdes, crossval=cval, performance=custom_perf)
#' ret <- run_searchlight(mvpmod)
#' stopifnot("accuracy" %in% ret$metrics)
mvpa_model <- function(model, 
                       dataset,
                       design,
                       model_type=c("classification", "regression"), 
                       crossval=NULL, 
                       feature_selector=NULL, 
                       tune_grid=NULL, 
                       tune_reps=15,
                       performance=NULL,
                       class_metrics=TRUE,
                       compute_performance=TRUE,
                       return_predictions=TRUE,
                       return_fits=FALSE) {
                       
  
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
  
  # Determine if we have a test set based on both dataset and design
  has_test <- !is.null(design$y_test) && !is.null(dataset$test_data)
  
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
  ret <- create_model_spec("mvpa_model", 
              dataset=dataset,
              design=design,
              model=model,
              model_type=model_type,
              model_name=model$label,
              tune_grid=tune_grid,
              tune_reps=tune_reps,
              feature_selector=feature_selector,
              crossval=crossval,
              performance=perf,
              compute_performance=compute_performance,
              return_predictions=return_predictions,
              return_fits=return_fits,
              has_test_set=has_test)  # Add flag for test set presence
  
  ret
  
}


#' @export
has_crossval.mvpa_model <- function(obj) {
  !is.null(obj$crossval)
}

#' @export
has_crossval.model_spec <- function(obj) {
  !is.null(obj$crossval)
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
  

#' @rdname strip_dataset-methods
#' @export
strip_dataset.default <- function(obj, ...) {
  if (!is.null(obj$dataset)) {
    futile.logger::flog.debug("Stripping dataset from model specification.")
    obj$dataset <- NULL
  } else {
    futile.logger::flog.debug("Dataset already NULL or missing in model specification.")
  }
  obj
}

#' @rdname strip_dataset-methods
#' @export
strip_dataset.mvpa_model <- function(obj, ...) {
  obj$dataset <- NULL
  obj
}
  
  
