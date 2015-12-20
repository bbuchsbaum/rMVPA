#' @export
performance <- function(x,...) {
  UseMethod("performance")
}

#' @export
performance.SimilarityResult <- function(x, splitList=NULL, classMetrics=FALSE) {  
  simAll <- x$simWithinTable$sim
  names(simAll) <- paste0("sim_", x$simWithinTable$label)
  c(simWithin=x$sWithin, simDiff=x$sWithin - x$sBetween, simAll)
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
performance.TwoWayClassificationResult <- function(x,splitList=NULL, classMetrics=FALSE) {
  if (is.null(splitList)) {
    .twowayPerf(x$observed, x$predicted, x$probs)
  } else {
    total <- .twowayPerf(x$observed, x$predicted, x$probs)
    subtots <- unlist(lapply(names(splitList), function(tag) {
      ind <- splitList[[tag]]
      ret <- .twowayPerf(x$observed[ind], x$predicted[ind], x$probs[ind,])
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
    c(ZAccuracy=-qnorm(out$p.value), Accuracy=sum(obs == as.character(predicted))/length(obs), Combined_AUC=mean(aucres), aucres)
  } else {
    c(ZAccuracy=-qnorm(out$p.value), Accuracy=sum(obs == as.character(predicted))/length(obs), Combined_AUC=mean(aucres))
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


