


#' predicted_class
#' @param prob a matrix of predicted probabilities with column names indicating the classes
predicted_class <- function(prob) {
  maxid <- max.col(prob, ties.method="random")
  pclass <- colnames(prob)[maxid]
}


#' @param x the result object
#' @param split_list split results by indexed sub-groups
#' @export
#' @rdname performance-methods
performance.regression_result <- function(x, split_list) {
  if (!is.null(split_list)) {
    ## TODO: add support
    stop("split_by not supported for regression analyses yet.")
  }
  R2 <- 1 - sum((x$observed - x$predicted)^2)/sum((x$observed-mean(x$observed))^2)
  rmse <- sqrt(mean((x$observed-x$predicted)^2))
  rcor <- cor(x$observed, x$predicted, method="spearman")
  c(R2=R2, RMSE=rmse, spearcor=rcor)
}


#' custom_performance
#' 
#' applies a user-supplied performance metric to a prediction result
#' 
#' @param x the result object
#' @param custom_fun the function used to compute performance metrics, i.e. \code{custom_fun(x)}
#' @param split_list the splitting groups
#' @export
custom_performance <- function(x, custom_fun, split_list=NULL) {
  if (is.null(split_list)) {
    custom_fun(x)
  } else {
    total <- custom_fun(x)
    subtots <- unlist(lapply(names(split_list), function(tag) {
      ind <- split_list[[tag]]
      ret <- custom_fun(subResult(x, ind))
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
merge_results.multiway_classification_result <- function(x,...) {
  
  rlist <- list(x,...)
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
  
  #ncorrect <- sum(obs == predicted)
  ntotal <- length(obs)
  #maxClass <- max(table(obs))
  
  #out <- binom.test(ncorrect,
  #                  ntotal,
  #                  p = maxClass/ntotal,
  #                  alternative = "greater")
  
 
  aucres <- sapply(1:ncol(probs), function(i) {
    lev <- try(levels(observed)[i])
    pos <- obs == lev
    pclass <- probs[,i]
    pother <- rowMeans(probs[,-i, drop=FALSE])
    Metrics::auc(as.numeric(pos), pclass - pother)-.5
  })
  
  names(aucres) <- paste0("AUC_", colnames(probs))
  
  
  if (class_metrics) {
    #c(ZAccuracy=-qnorm(out$p.value), Accuracy=sum(obs == as.character(predicted))/length(obs), AUC=mean(aucres), aucres)
    c(Accuracy=sum(obs == as.character(predicted))/length(obs), AUC=mean(aucres), aucres)
  } else {
    #c(ZAccuracy=-qnorm(out$p.value), Accuracy=sum(obs == as.character(predicted))/length(obs), AUC=mean(aucres))
    c(Accuracy=sum(obs == as.character(predicted))/length(obs), AUC=mean(aucres))
  }
}
  




