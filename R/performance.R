


#' Calculate the Predicted Class from Probability Matrix
#'
#' This function calculates the predicted class from a matrix of predicted probabilities. The class with the highest probability is selected as the predicted class.
#'
#' @param prob A matrix of predicted probabilities with column names indicating the classes.
#' @return A vector of predicted classes corresponding to the highest probability for each row in the input matrix.
#' @examples
#' prob <- matrix(c(0.2, 0.8,
#'                  0.6, 0.4),
#'                nrow = 2, byrow = TRUE,
#'                dimnames = list(NULL, c("A", "B")))
#' predicted_class(prob)
#' @export
predicted_class <- function(prob) {
  maxid <- max.col(prob, ties.method="random")
  pclass <- colnames(prob)[maxid]
}

#' Calculate Performance Metrics for Regression Result
#'
#' This function calculates performance metrics for a regression result object, including R-squared, Root Mean Squared Error (RMSE), and Spearman correlation.
#'
#' @param x A \code{regression_result} object.
#' @param split_list Split results by indexed sub-groups (not supported for regression analyses yet).
#' @param ... extra args (not used).
#' @return A named vector with the calculated performance metrics: R-squared, RMSE, and Spearman correlation.
#' @details
#' The function calculates the following performance metrics for the given regression result object:
#' - R-squared: proportion of variance in the observed data that is predictable from the fitted model.
#' - RMSE: root mean squared error, a measure of the differences between predicted and observed values.
#' - Spearman correlation: a measure of the monotonic relationship between predicted and observed values.
#' @seealso \code{\link{regression_result}}
#' @export
#' @importFrom yardstick rsq_vec rmse_vec
#' @importFrom stats cor
performance.regression_result <- function(x, split_list=NULL,...) {
  if (!is.null(split_list)) {
    ## TODO: add support
    stop("split_by not supported for regression analyses yet.")
  }
  
  obs <- x$observed
  pred <- x$predicted
  
  # Ensure pred is a numeric vector without table attributes
  pred <- as.numeric(pred)
  
  res_rsq <- yardstick::rsq_vec(truth = obs, estimate = pred)
  res_rmse <- yardstick::rmse_vec(truth = obs, estimate = pred)
  res_spearcor <- tryCatch(stats::cor(obs, pred, method="spearman", use="pairwise.complete.obs"), error = function(e) NA_real_)
  
  c(R2=res_rsq, RMSE=res_rmse, spearcor=res_spearcor)
}


#' Apply Custom Performance Metric to Prediction Result
#'
#' This function applies a user-supplied performance metric to a prediction result object.
#'
#' @param x The prediction result object.
#' @param custom_fun The function used to compute performance metrics, i.e., \code{custom_fun(x)}.
#' @param split_list An optional named list of splitting groups. If provided, the performance metric will be computed for each group and returned as a named vector.
#' @return A named vector with the calculated custom performance metric(s).
#' @details
#' The function allows users to apply a custom performance metric to a prediction result object.
#' If a split list is provided, the performance metric will be computed for each group separately, and the results will be returned as a named vector.
#' @examples
#' cres <- binary_classification_result(
#'   observed  = factor(c("A", "B")),
#'   predicted = factor(c("A", "A")),
#'   probs = matrix(c(0.9, 0.1,
#'                    0.6, 0.4),
#'                  ncol = 2, byrow = TRUE,
#'                  dimnames = list(NULL, c("A", "B")))
#' )
#' acc_fun <- function(x) mean(x$observed == x$predicted)
#' custom_performance(cres, acc_fun)
#' @export
custom_performance <- function(x, custom_fun, split_list=NULL) {
  if (is.null(split_list)) {
    custom_fun(x)
  } else {
    total <- custom_fun(x)
    subtots <- unlist(lapply(names(split_list), function(tag) {
      ind <- split_list[[tag]]
      ret <- custom_fun(sub_result(x, ind))
      names(ret) <- paste0(names(ret), "_", tag)
      ret
    }))
    
    c(total, subtots)
  }
  
}

#' @rdname merge_results-methods
#' @method merge_results binary_classification_result
#' @export
merge_results.binary_classification_result <- function(x,...) {
  rlist <- list(x,...)
  probs <- Reduce("+", lapply(rlist, function(x) x$probs))/length(rlist)
  
  mc <- max.col(probs)
  predicted <- levels(x$observed)[mc]
  binary_classification_result(observed=x$observed, predicted=predicted, probs=probs, testind=x$testind, 
                               test_design=x$test_design, predictor=x$predictor)
}

#' @rdname merge_results-methods
#' @method merge_results regression_result
#' @export
merge_results.regression_result <- function(x,...) {
  rlist <- list(x,...)
  pred <- Reduce("+", lapply(rlist, function(x) x$predicted))/length(rlist)
  regression_result(observed=x$observed, predicted=pred, testind=x$testind, 
                               test_design=x$test_design, predictor=x$predictor)
}



#' @export
prob_observed.binary_classification_result <- function(x) {
  x$probs[cbind(seq(1,nrow(x$probs)),as.integer(x$observed))]
}

#' @export
prob_observed.multiway_classification_result <- function(x) {
  x$probs[cbind(seq(1,nrow(x$probs)),as.integer(x$observed))]
}

#' @export
prob_observed.regression_result <- function(x) {
  # Regression results don't have probabilities, return NULL
  # This allows the searchlight combiner to skip prob_observed for regression
  NULL
}

#' @rdname merge_results-methods
#' @method merge_results multiway_classification_result
#' @export
merge_results.multiway_classification_result <- function(x,...) {
  
  rlist <- list(x,...)
  #ds <- sapply(rlist, function(x) nrow(x$probs))
  
  probs <- Reduce("+", lapply(rlist, function(x) x$probs))/length(rlist)
  mc <- max.col(probs)
  predicted <- levels(x$observed)[mc]
  
  multiway_classification_result(observed=x$observed, predicted=predicted, probs=probs, 
                                 testind=x$testind,  test_design=x$test_design, predictor=x$predictor)
}

#' @export
performance.binary_classification_result <- function(x, split_list=NULL,...) {
  stopifnot(length(x$observed) == length(x$predicted))
  
  if (is.null(split_list)) {
    ret <- binary_perf(x$observed, x$predicted, x$probs)
  } else {
    total <- binary_perf(x$observed, x$predicted, x$probs)
    
    subtots <- unlist(lapply(names(split_list), function(tag) {
      ind <- split_list[[tag]]
      if (!is.null(x$testind)) {
        ind <- which(x$testind %in% ind)
      }
      ret <- binary_perf(x$observed[ind], x$predicted[ind], x$probs[ind,])
      names(ret) <- paste0(names(ret), "_", tag)
      ret
    }))
    
    ret <- c(total, subtots)
  }
}


#' @export
performance.multiway_classification_result <- function(x, split_list=NULL, class_metrics=FALSE,...) {
  stopifnot(length(x$observed) == length(x$predicted))

  if (is.null(split_list)) {
    multiclass_perf(x$observed, x$predicted, x$probs, class_metrics)
  } else {
    total <- multiclass_perf(x$observed, x$predicted, x$probs, class_metrics)
    subtots <- unlist(lapply(names(split_list), function(tag) {
      ind <- split_list[[tag]]
      
      if (!is.null(x$testind)) {
        ind <- which(x$testind %in% ind)
      }
      
      ret <- multiclass_perf(x$observed[ind], x$predicted[ind], x$probs[ind,], class_metrics)
      names(ret) <- paste0(names(ret), "_", tag)
      ret
    }))
    
    c(total, subtots)
    
  }
  
}

#' @keywords internal
#' @noRd
combinedACC <- function(Pred, Obs) {
  levs <- levels(as.factor(Obs))
  maxind <- apply(Pred, 1, which.max)
  pclass <- levs[maxind]
  sum(pclass == Obs)/length(pclass)
  
}


#' @keywords internal
#' @importFrom yardstick accuracy_vec roc_auc_vec
binary_perf <- function(observed, predicted, probs) {
  # Ensure observed is a factor with levels in a consistent order
  # and probs columns match this order.
  lvls <- levels(observed)
  if (length(lvls) != 2) stop("binary_perf expects 2 levels in observed.")
  
  # Ensure predicted is a factor with the same levels as observed
  predicted_factor <- factor(predicted, levels = lvls)
  
  # Assuming probs has columns named after levels(observed) or in the same order.
  # And positive class is the second level.
  prob_positive_class <- if (ncol(probs) == 2) probs[, lvls[2]] else probs[,1] # Adapt if probs is single col

  res_acc <- yardstick::accuracy_vec(truth = observed, estimate = predicted_factor)
  res_auc <- tryCatch(
     yardstick::roc_auc_vec(truth = observed, estimate = prob_positive_class, event_level = "second"),
     error = function(e) NA_real_
  )
  
  # Note: The original code subtracted 0.5 from AUC. This is unusual but preserved for compatibility
  c(Accuracy = res_acc, AUC = res_auc - 0.5) # yardstick returns numeric, ensure it's named
}

#' @keywords internal
#' @importFrom yardstick accuracy_vec roc_auc_vec
multiclass_perf <- function(observed, predicted, probs, class_metrics=FALSE) {
  lvls <- levels(observed)
  predicted_factor <- factor(predicted, levels = lvls)

  acc <- yardstick::accuracy_vec(truth = observed, estimate = predicted_factor)
  
  # Calculate per-class AUC using one-vs-rest approach (matching original logic)
  aucres <- sapply(seq_along(lvls), function(i) {
    lev <- lvls[i]
    pos <- observed == lev
    pclass <- probs[,i]
    pother <- rowMeans(probs[,-i, drop=FALSE])
    # Original uses pclass - pother as the score
    score <- pclass - pother
    binary_truth <- factor(ifelse(pos, "positive", "negative"), levels = c("negative", "positive"))
    tryCatch(
      yardstick::roc_auc_vec(truth = binary_truth, estimate = score, event_level = "second") - 0.5,
      error = function(e) NA_real_
    )
  })
  
  names(aucres) <- paste0("AUC_", colnames(probs))
  
  metrics <- c(Accuracy = acc, AUC = mean(aucres, na.rm=TRUE))
  
  if (class_metrics) {
    c(metrics, aucres)
  } else {
    metrics
  }
}
  




