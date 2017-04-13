#' Greedy optimization of the area under the curve
#' @description This algorithm optimizes the area under the curve for classification models
#' @param X the matrix of predictors
#' @param Y the dependent variable
#' @param iter an integer for the number of iterations
#' @param metricFun
#' @return A numeric of the weights for each model.
#' @details If the optimization fails to produce an error term better than the best
#' component model, a message is returned and the best optimization after N iterations
#' is returned.
#' @export
#' @examples
#' x <- matrix(runif(10), ncol=2)
#' y <- sample(c('Y', 'N'), 5, replace=TRUE)
#' greedyOptAUC(x, y)
greedyTwoClassOptim <- function(X, Y, iter = 100L, metricFUN="AUC"){ 
  
  #metricFun <- Metrics::auc
  metricFun <- switch(metric,
                      AUC=Metrics::auc)
  #ACC=combinedACC)
  
  N           <- ncol(X)
  weights     <- rep(0L, N)
  pred        <-  0 * X
  sum.weights <- 0L
  stopper     <- max(apply(X, 2, function(p) metricFun(Y, p)))
  
  while(sum.weights < iter) { 
    sum.weights   <- sum.weights + 1L
    pred          <- (pred + X) * (1L / sum.weights)
    errors        <- apply(pred, 2, function(p) metricFun(Y, p))
    best          <- which.max(errors)
    weights[best] <- weights[best] + 1L
    pred          <- pred[,best] * sum.weights
    maxtest       <- max(errors)
    
    
  }
  
  print(paste("iter:", sum.weights, "best:", maxtest))
  
  if(stopper > maxtest){
    testresult <- round(maxtest/stopper, 5) * 100
    wstr <- paste0("Optimized weights not better than best model. Ensembled result is ",
                   testresult, "%", " of best model AUC. Try more iterations.")
    message(wstr)
  }
  
  return(weights/sum(weights))
}


#' Greedy optimization of the area under the curve in the multiclass setting
#' @description This algorithm optimizes the area under the curve for multiclass classification models
#' @param PredList a list of prediction matrices
#' @param Y the dependent variable
#' @param iter an integer for the number of iterations
#' @param metric
#' @return A numeric of the weights for each model.
#' @details If the optimization fails to produce an error term better than the best
#' component model, a message is returned and the best optimization after N iterations
#' is returned.
#' @export
#' @examples
#' x <- matrix(runif(10), ncol=2)
#' y <- sample(c('Y', 'N'), 5, replace=TRUE)
#' greedOptAUC(x, y)
greedyMultiClassOptim <- function(PredList, Y, iter = 100L, metric="AUC", oneByOne=FALSE) { 
  if(is.character(Y)){
    Y <- factor(Y)
  }
  
  stopifnot(is.factor(Y))
  
  metricFun <- switch(metric,
                      AUC=combinedAUC,
                      ACC=combinedACC)
  
  
  if (oneByOne) {
    levs <- levels(Y)
    wts <- lapply(seq_along(levs), function(i) {
      lev <- levels(Y)[i]
      Y0 <- ifelse(Y == lev, 1, 0)
      X0 <- do.call(cbind, lapply(PredList, function(p) p[,i]))
      greedyTwoClassOptim(X0, Y0)
    })
    rowMeans(do.call(cbind, wts))
    
  } else {
    
    N           <- length(PredList)
    weights     <- rep(0L, N)
    pred        <- matrix(0, length(Y), ncol(PredList[[1]]))
    sum.weights <- 0L
    stopper     <- max(unlist(lapply(PredList, function(p) metricFun(p, Y))))
    
    while(sum.weights < iter) { 
      sum.weights   <- sum.weights + 1L
      pred          <- lapply(1:length(PredList), function(i) (PredList[[i]] + pred) * 1L/sum.weights)
      errors        <- unlist(lapply(pred, function(p) metricFun(p, Y)))
      best          <- which.max(errors)
      weights[best] <- weights[best] + 1L
      pred          <- pred[[best]] * sum.weights
      maxtest       <- max(errors)
    }
    
    print(paste("iter:", sum.weights, "best:", maxtest))
    
    if(stopper > maxtest){
      testresult <- round(maxtest/stopper, 5) * 100
      wstr <- paste0("Optimized weights not better than best model. Ensembled result is ",
                     testresult, "%", " of best model AUC. Try more iterations.")
      message(wstr)
    }
    
    weights/sum(weights)
  }
}




#' @export
BaggedMultiClassOptAUC <- function(PredList, Y, iter = 20, bagIter=10, bagFrac=.5) {
  wtlist <- lapply(1:bagIter, function(i) {
    sam <- sort(sample(seq_along(PredList), bagFrac * length(PredList))) 
    wts <- greedyMultiClassOptim(PredList[sam], Y, iter=iter, metric="AUC")
    ret <- numeric(length(PredList))
    ret[sam] <- wts
    ret
  })
  
  wts <- rowSums(do.call(cbind, wtlist))
  
}

#' @export
BaggedMultiClassOptACC <- function(PredList, Y, iter = 20, bagIter=10, bagFrac=.5) {
  wtlist <- lapply(1:bagIter, function(i) {
    sam <- sort(sample(seq_along(PredList), bagFrac * length(PredList))) 
    wts <- greedyMultiClassOptim(PredList[sam], Y, iter=iter, metric="ACC")
    ret <- numeric(length(PredList))
    ret[sam] <- wts
    ret
  })
  
  wts <- rowSums(do.call(cbind, wtlist))
  
}

#' @export
BaggedTwoClassOptAUC <- function(X, Y, iter = 20L, bagIter=20, bagFrac=.5) {
  wtlist <- lapply(1:bagIter, function(i) {
    sam <- sort(sample(1:ncol(X), bagFrac * ncol(X))) 
    wts <- greedyTwoClassOptim(X[,sam], Y, iter=30)
    ret <- numeric(ncol(X))
    ret[sam] <- wts
    ret
  })
  
  wts <- rowSums(do.call(cbind, wtlist))
  
}


#' @export
BaggedTwoClassOptACC <- function(X, Y, iter = 20, bagIter=20, bagFrac=.5) {
  wtlist <- lapply(1:bagIter, function(i) {
    sam <- sort(sample(1:ncol(X), bagFrac * ncol(X))) 
    wts <- greedyTwoClassOptim(X[,sam], Y, iter=iter)
    ret <- numeric(ncol(X))
    ret[sam] <- wts
    ret
  })
  
  wts <- rowSums(do.call(cbind, wtlist))
  
}
