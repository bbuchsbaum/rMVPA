
#' predicted_class
#' @param prob a matrix of predicted probabilities with column names indicating the classes
predicted_class <- function(prob) {
  maxid <- apply(prob, 1, which.max)
  pclass <- colnames(prob)[maxid]
}


#' @export
#' @param x
#' @param splitList
#' @param classMetrics
performance.SimilarityResult <- function(x, splitList=NULL, classMetrics=FALSE) {  
  #simAll <- x$simWithinTable$sim
  #names(simAll) <- paste0("sim_", x$simWithinTable$label)
  #c(simWithin=x$sWithin, simDiff=x$sWithin - x$sBetween, simAll)
}

.twowayPerf <- function(observed, predicted, probs) {
  ncorrect <- sum(observed == predicted)
  ntotal <- length(observed)
  maxClass <- max(table(observed))
  
  out <- binom.test(ncorrect,
                    ntotal,
                    p = maxClass/ntotal,
                    alternative = "greater")
  
  c(ZAccuracy=-qnorm(out$p.value), Accuracy=ncorrect/ntotal, AUC=Metrics::auc(observed == levels(observed)[2], probs[,2])-.5)
  
}

#' @export
#' @param x
#' @param splitList
#' @param classMetrics
performance.RegressionResult <- function(x, splitList, classMetrics=FALSE) {
  R2 <- 1 - sum((x$observed - x$predicted)^2)/sum((x$observed-mean(x$observed))^2)
  rmse <- sqrt(mean((x$observed-x$predicted)^2))
  rcor <- cor(x$observed, x$predicted, method="spearman")
  c(R2=R2, RMSE=rmse, spearcor=rcor)
}

#' @export
#' 
customPerformance <- function(x, customFun, splitList=NULL) {
  if (is.null(splitList)) {
    customFun(x)
  } else {
    total <- customFun(x)
    subtots <- unlist(lapply(names(splitList), function(tag) {
      ind <- splitList[[tag]]
      ret <- customFun(subResult(x, ind))
      names(ret) <- paste0(names(ret), "_", tag)
      ret
    }))
    
    c(total, subtots)
  }
  
}

#' @export
merge_results.TwoWayClassificationResult <- function(x,y) {
  probs <- (x$probs + y$probs)/2
  TwoWayClassificationResult(observed=x$observed, predicted=NULL, probs=probs, testDesign=x$testDesign, predictor=x$predictor)
}

#' @export
merge_results.MultiWayClassificationResult <- function(x,y) {
  probs <- (x$probs + y$probs)/2
  MultiWayClassificationResult(observed=x$observed, predicted=NULL, probs=probs, testDesign=x$testDesign, predictor=x$predictor)
}

#' @export
performance.TwoWayClassificationResult <- function(x, splitList=NULL, classMetrics=FALSE, customFun=NULL) {
  if (is.null(splitList)) {
    ret <- .twowayPerf(x$observed, x$predicted, x$probs)
  } else {
    total <- .twowayPerf(x$observed, x$predicted, x$probs)
    
    subtots <- unlist(lapply(names(splitList), function(tag) {
      ind <- splitList[[tag]]
      ret <- .twowayPerf(x$observed[ind], x$predicted[ind], x$probs[ind,])
      names(ret) <- paste0(names(ret), "_", tag)
      ret
    }))
    
    ret <- c(total, subtots)
  }
 
}

#' @export
performance.MultiWayClassificationResult <- function(x, splitList=NULL, classMetrics=FALSE) {
  stopifnot(length(x$observed) == length(x$predicted))
  
  if (is.null(splitList)) {
    .multiwayPerf(x$observed, x$predicted, x$probs, classMetrics)
  } else {
    total <- .multiwayPerf(x$observed, x$predicted, x$probs, classMetrics)
    subtots <- unlist(lapply(names(splitList), function(tag) {
      ind <- splitList[[tag]]
      ret <- .multiwayPerf(x$observed[ind], x$predicted[ind], x$probs[ind,], classMetrics)
      names(ret) <- paste0(names(ret), "_", tag)
      ret
    }))
    
    c(total, subtots)
    
  }
  
}

combinedAUC <- function(Pred, Obs) {
  mean(sapply(1:ncol(Pred), function(i) {
    lev <- levels(Obs)[i]
    pos <- Obs == lev
    pclass <- Pred[,i]
    pother <- rowMeans(Pred[,-i,drop=FALSE])
    Metrics::auc(as.numeric(pos), pclass - pother)-.5
  }))
}

combinedACC <- function(Pred, Obs) {
  levs <- levels(as.factor(Obs))
  maxind <- apply(Pred, 1, which.max)
  pclass <- levs[maxind]
  sum(pclass == Obs)/length(pclass)
  
}

.multiwayPerf <- function(observed, predicted, probs, classMetrics=FALSE) {
  obs <- as.character(observed)
  
  ncorrect <- sum(obs == predicted)
  ntotal <- length(obs)
  maxClass <- max(table(obs))
  
  out <- binom.test(ncorrect,
                    ntotal,
                    p = maxClass/ntotal,
                    alternative = "greater")
  
 
  aucres <- sapply(1:ncol(probs), function(i) {
    lev <- levels(observed)[i]
    pos <- obs == lev
    pclass <- probs[,i]
    pother <- rowMeans(probs[,-i, drop=FALSE])
    Metrics::auc(as.numeric(pos), pclass - pother)-.5
  })
  
  names(aucres) <- paste0("AUC_", colnames(probs))
  
  
  if (classMetrics) {
    c(ZAccuracy=-qnorm(out$p.value), Accuracy=sum(obs == as.character(predicted))/length(obs), AUC=mean(aucres), aucres)
  } else {
    c(ZAccuracy=-qnorm(out$p.value), Accuracy=sum(obs == as.character(predicted))/length(obs), AUC=mean(aucres))
  }
}
  




