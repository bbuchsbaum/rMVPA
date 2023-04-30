


#' Calculate the Predicted Class from Probability Matrix
#'
#' This function calculates the predicted class from a matrix of predicted probabilities. The class with the highest probability is selected as the predicted class.
#'
#' @param prob A matrix of predicted probabilities with column names indicating the classes.
#' @return A vector of predicted classes corresponding to the highest probability for each row in the input matrix.
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
#' @return A named vector with the calculated performance metrics: R-squared, RMSE, and Spearman correlation.
#' @details
#' The function calculates the following performance metrics for the given regression result object:
#' - R-squared: proportion of variance in the observed data that is predictable from the fitted model.
#' - RMSE: root mean squared error, a measure of the differences between predicted and observed values.
#' - Spearman correlation: a measure of the monotonic relationship between predicted and observed values.
#' @seealso \code{\link{regression_result}}
#' @export
performance.regression_result <- function(x, split_list,...) {
  if (!is.null(split_list)) {
    ## TODO: add support
    stop("split_by not supported for regression analyses yet.")
  }
  
  #browser()
  R2 <- 1 - sum((x$observed - x$predicted)^2)/sum((x$observed-mean(x$observed))^2)
  rmse <- sqrt(mean((x$observed-x$predicted)^2))
  rcor <- cor(x$observed, x$predicted, method="spearman")
  c(R2=R2, RMSE=rmse, spearcor=rcor)
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

#' @export
merge_results.binary_classification_result <- function(x,...) {
  rlist <- list(x,...)
  probs <- Reduce("+", lapply(rlist, function(x) x$probs))/length(rlist)
  
  mc <- max.col(probs)
  predicted <- levels(x$observed)[mc]
  binary_classification_result(observed=x$observed, predicted=predicted, probs=probs, testind=x$testind, 
                               test_design=x$test_design, predictor=x$predictor)
}

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
combinedAUC <- function(Pred, Obs) {
  Obs <- as.factor(Obs)
  mean(sapply(1:ncol(Pred), function(i) {
    lev <- levels(Obs)[i]
    pos <- Obs == lev
    pclass <- Pred[,i]
    pother <- rowMeans(Pred[,-i,drop=FALSE])
    Metrics::auc(as.numeric(pos), pclass - pother)-.5
  }))
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
binary_perf <- function(observed, predicted, probs) {
  obs <- as.character(observed)
  ncorrect <- sum(obs == predicted)
  ntotal <- length(obs)
  #maxClass <- max(table(obs))
  
  #out <- binom.test(ncorrect,
  #                  ntotal,
  #                  p = maxClass/ntotal,
  #                  alternative = "greater")
  
  
  #c(ZAccuracy=-qnorm(out$p.value), Accuracy=ncorrect/ntotal, AUC=Metrics::auc(observed == levels(observed)[2], probs[,2])-.5)
  c(Accuracy=ncorrect/ntotal, AUC=Metrics::auc(observed == levels(observed)[2], probs[,2])-.5)
  
}

#' @keywords internal
multiclass_perf <- function(observed, predicted, probs, class_metrics=FALSE) {
  
  obs <- as.character(observed)
  ntotal <- length(obs)
 
  aucres <- sapply(1:ncol(probs), function(i) {
    lev <- try(levels(observed)[i])
    pos <- obs == lev
    pclass <- probs[,i]
    pother <- rowMeans(probs[,-i, drop=FALSE])
    Metrics::auc(as.numeric(pos), pclass - pother)-.5
  })
  
  names(aucres) <- paste0("AUC_", colnames(probs))
  
  
  if (class_metrics) {
    c(Accuracy=sum(obs == as.character(predicted))/length(obs), AUC=mean(aucres, na.rm=TRUE), aucres)
  } else {
    c(Accuracy=sum(obs == as.character(predicted))/length(obs), AUC=mean(aucres, na.rm=TRUE))
  }
}
  




