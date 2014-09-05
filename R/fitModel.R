


corsimFit <- function(x, y) {
  list(conditionMeans = splitReduce(x, y, mean))
}

predict.corSimFit <- function(modelFit, newData) {
  res <- sapply(1:nrow(newData), function(i) {
    pattern <- newData[i,]
    cor(pattern, modelFit$conditionMeans)
  })
  
}

MVPAModels <- list()
    
MVPAModels$corsim <- list(type = "Classification", 
                    library = "rMVPA", 
                    loop = NULL, 
                    parameters=data.frame(parameters="lambda", class="numeric", labels="lambda"),
                    grid=function(x, y, len = NULL) data.frame(lambda=seq(.2, .9, length.out=len)),
                    fit=function(x, y, wts, param, lev, last, weights, classProbs, ...) corsimFit(x,y),
                    predict=function(modelFit, newdata, preProc = NULL, submodels = NULL) predict(modelFit, as.matrix(newdata))$class,
                    prob=function(modelFit, newdata, preProc = NULL, submodels = NULL) {
                      predict(modelFit, as.matrix(newdata))$posterior
                    })


MVPAModels$lda_strimmer <- list(type = "Classification", 
                 library = "sda", 
                 loop = NULL, 
                 parameters=data.frame(parameters="lambda", class="numeric", labels="lambda"),
                 grid=function(x, y, len = NULL) data.frame(lambda=0),
                 fit=function(x, y, wts, param, lev, last, weights, classProbs, ...) sda(Xtrain=as.matrix(x), L=y, verbose=FALSE, ...),
                 predict=function(modelFit, newdata, preProc = NULL, submodels = NULL) predict(modelFit, as.matrix(newdata))$class,
                 prob=function(modelFit, newdata, preProc = NULL, submodels = NULL) {
                    predict(modelFit, as.matrix(newdata))$posterior
                })


MVPAModels$lda_thomaz <- list(type = "Classification", 
                        library = "sparsediscrim", 
                        loop = NULL, 
                        parameters=data.frame(parameters="prior", class="numeric", labels="prior"),
                        grid=function(x, y, len = NULL) { data.frame(prior=.5) },
                        fit=function(x, y, wts, param, lev, last, weights, classProbs, ...) { lda_thomaz(x,y, ...) },
                        predict=function(modelFit, newdata, preProc = NULL, submodels = NULL) { predict(modelFit, as.matrix(newdata))$class },
                        prob=function(modelFit, newdata, preProc = NULL, submodels = NULL) { 
                          scores <- t(predict(modelFit, newdata)$scores)
                          mc <- scores[cbind(1:nrow(scores), max.col(scores, ties.method = "first"))]
                          probs <- exp(scores - mc)
                          zapsmall(probs/rowSums(probs))
                        })


#' create an \code{MVPA} instance
#' @param trainVec
#' @param Y
#' @param mask
#' @param blockVar
MVPADataset <- function(trainVec, Y, mask, blockVar) {
  ret <- list(
    trainVec=trainVec,
    Y=Y,
    mask=mask,
    blockVar=blockVar)
  
  class(ret) <- c("MVPADataset", "list")
  ret
}


TwoWayClassificationResult <- function(observed, predicted, probs) {
  ret <- list(observed=observed,
              predicted=predicted,
              probs=probs)
  
  class(ret) <- c("TwoWayClassificationResult", "list")
  ret
  
}

MultiWayClassificationResult <- function(observed, predicted, probs) {
  ret <- list(observed=observed,
              predicted=predicted,
              probs=probs)
  
  class(ret) <- c("MultiWayClassificationResult", "list")
  ret
  
}

performance <- function(x,...) {
  UseMethod("performance")
}



performance.TwoWayClassificationResult <- function(x) {
  c(acc=sum(x$observed == x$predicted)/length(x$observed), acc=Metrics::auc(x$observed == levels(x$observed)[2], x$probs[,2]))
}

performance.MultiWayClassificationResult <- function(x) {
  c(accuracy=sum(x$observed == x$predicted)/length(x$observed))
}


#' create a classification model
createModel <- function(method) {
  ret <- list(method=method)
  class(ret) <- c(method, "list")
  ret
}


trainModel <- function(x, ...) {
  UseMethod("trainModel")
}


crossval <- function(X, Y, foldSplit, method, ncores=2, tuneGrid=NULL, tuneLength=1) {
  
  if (is.null(tuneGrid) || tuneLength == 1 || nrow(tuneGrid) == 1) {
    ctrl <- caret::trainControl("none", verboseIter=TRUE, classProb=TRUE)
  } else {
    ctrl <- caret::trainControl("cv", verboseIter=TRUE, classProb=TRUE)
  }
  
  res <- parallel::mclapply(foldSplit, function(fidx) {
    Xtrain <- X[-fidx,]
    Ytrain <- Y[-fidx]
    Xtest <- X[fidx,]   
    if (!is.null(tuneGrid)) {
      fit <- caret::train(Xtrain, Ytrain, method=method, trControl=ctrl, tuneGrid=tuneGrid)
    }
    else {
      fit <- caret::train(Xtrain, Ytrain, method=method, trControl=ctrl, tuneLength=tuneLength)
    }
    cbind(class=predict(fit, newdata=Xtest), predict(fit, newdata=Xtest, type="prob"))
  }, mc.cores=ncores)
  

  
  ret <- do.call(rbind, res)
  if (is.factor(Y) && length(levels(Y)) == 2) {
    TwoWayClassificationResult(Y, ret[,1], ret[,2:3])
  } else if (is.factor(Y) && length(levels(Y)) > 2) {
    MultiWayClassificationResult(Y, ret[,1], ret[,2:ncol(ret)])
  }
}

fitModel <- function(model, dset, voxelGrid, ncores=2, tuneGrid=NULL, tuneLength=NULL) {
  M <- series(dset$trainVec, voxelGrid)
  result <- crossval(M, dset$Y, split(1:length(dset$blockVar), dset$blockVar), model, ncores, tuneGrid,tuneLength)
}

trainModel.default <- function(model, dset, voxelGrid, ncores=2) {
  M <- series(dset$trainVec, voxelGrid)
  fitlist <- crossval(M, dset$Y, dset$blockVar, model)
}



trainModel.lda_strimmer <- function(dset, voxelGrid, ncores=2) {
  M <- series(dset$trainVec, voxelGrid)
  ctrl <- caret::trainControl("none")
  fit <- caret::train(x$x, x$y, method=MVPA.sda, trControl=ctrl, tuneLength=2)
}

trainModel.lda_thomaz <- function(x) {
  ctrl <- caret::trainControl("none")
  fit <- caret::train(x$x, x$y, method=MVPA.lda_thomaz, trControl=ctrl, tuneLength=0)
}
