#' @export
performance <- function(x,...) {
  UseMethod("performance")
}



#' @export
performance.TwoWayClassificationResult <- function(x,...) {
  
  ncorrect <- sum(x$observed == x$predicted)
  ntotal <- length(x$observed)
  maxClass <- max(table(x$observed))
  
  out <- binom.test(ncorrect,
                    ntotal,
                    p = maxClass/ntotal,
                    alternative = "greater")
  
  c(ZAccuracy=-qnorm(out$p.value), Accuracy=ncorrect/ntotal, AUC=Metrics::auc(x$observed == levels(x$observed)[2], x$probs[,2])-.5)
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

#' @export
performance.MultiWayClassificationResult <- function(x,...) {
  obs <- as.character(x$observed)
  
  ncorrect <- sum(obs == x$predicted)
  ntotal <- length(obs)
  maxClass <- max(table(obs))
  
  out <- binom.test(ncorrect,
                    ntotal,
                    p = maxClass/ntotal,
                    alternative = "greater")
  
  aucres <- sapply(1:ncol(x$prob), function(i) {
    lev <- levels(x$observed)[i]
    pos <- obs == lev
    pclass <- x$prob[,i]
    pother <- rowMeans(x$prob[,-i])
    Metrics::auc(as.numeric(pos), pclass - pother)-.5
  })
  
  names(aucres) <- paste0("AUC_", colnames(x$prob))
  
  c(ZAccuracy=-qnorm(out$p.value), Accuracy=sum(obs == as.character(x$predicted))/length(obs), Combined_AUC=mean(aucres), aucres)
}


