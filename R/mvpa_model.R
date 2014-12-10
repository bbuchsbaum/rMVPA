


.noneControl <- caret::trainControl("none", verboseIter=TRUE, classProbs=TRUE)
.cvControl <- caret::trainControl("cv", verboseIter=TRUE, classProbs=TRUE)  


#' @export
invertFolds <- function(foldSplit, len) {
  index <- relist(1:length(unlist(foldSplit)), foldSplit)
  allInd <- 1:len
  index <- lapply(index, function(i) allInd[-i])
}


#' Create an exhaustive searchlight iterator
#' @param X data matrix
#' @param Y the labels
#' @param blockVar variable denoting the cross-validation folds
#' @export
MatrixFoldIterator <- function(X, Y, blockVar, trainSets=NULL, testSets=NULL) {
  
  if (is.null(trainSets) & is.null(testSets)) {
    testSets <- split(1:length(blockVar), blockVar)
    trainSets <- invertFolds(testSets, length(blockVar))
  }
  
  if (is.null(trainSets) && !is.null(testSets)) {
    stop("must supply both trainSets and testSets")
  }
  
  if (is.null(testSets) && !is.null(trainSets)) {
    stop("must supply both trainSets and testSets")
  }
  
  if (nrow(X) != length(Y)) {
    stop("X matrix must have same number of rows as Y variable")
  }
  
  
  index <- 0
  
  .getTrainSets <- function() {
    trainSets
  }
  .getTestSets <- function() {
    testSets
  }
  
  .getIndex <- function() {
    index
  }
  
  nextEl <- function() {
    if (index < length(trainSets)) { 
      index <<- index + 1
      list(Xtrain=X[trainSets[[index]], ], Ytrain=Y[trainSets[[index]]], Xtest=X[testSets[[index]],], Ytest=Y[testSets[[index]]])

    } else {
      stop('StopIteration')
    }
  }
  
  obj <- list(X=X, Y=Y,nextElem=nextEl, index=.getIndex, getTrainSets=.getTrainSets, getTestSets=.getTestSets)
  class(obj) <- c("FoldIterator", 'abstractiter', 'iter')
  obj
  
}



#' create an \code{MVPA} instance
#' @param trainVec
#' @param Y
#' @param mask
#' @param blockVar
#' @param testVec
#' @param testY
#' @param modelName
#' @param tuneGrid
#' @export
MVPADataset <- function(trainVec, Y, mask, blockVar, testVec, testY, modelName="corsim", tuneGrid=NULL) {
  
  model <- loadModel(modelName)
  
  testSets <- split(1:length(blockVar), blockVar)
  trainSets <- invertFolds(testSets, length(Y))
   
  ret <- list(
    trainVec=trainVec,
    Y=Y,
    mask=mask,
    blockVar=blockVar,
    testVec=testVec,
    testY=testY,
    model=model,
    tuneGrid=tuneGrid,
    testSets=testSets,
    trainSets=trainSets)
   
  class(ret) <- c("MVPADataset", "list")
  ret
}

NullResult <- function(vox, model, observed) {
  ret <- list(vox=vox,
              model=model,
              observed=observed)
  class(ret) <- c("NullResult", "list")
  ret
  
}


#' create an \code{TwoWayClassification} instance
#' @param observed
#' @param predicted
#' @param probs
#' @export
TwoWayClassificationResult <- function(observed, predicted, probs, predictor=NULL) {
  ret <- list(
              observed=observed,
              predicted=predicted,
              probs=as.matrix(probs),
              predictor=predictor
              )
  
  class(ret) <- c("TwoWayClassificationResult", "list")
  ret
 
}

#' create an \code{TwoWayClassification} instance
#' @param observed
#' @param predicted
#' @param probs

#' @export
MultiWayClassificationResult <- function(observed, predicted, probs, predictor=NULL) {
  ret <- list(
              observed=observed,
              predicted=predicted,
              probs=as.matrix(probs),
              predictor=predictor
              )
  
  class(ret) <- c("MultiWayClassificationResult", "list")
  ret
  
}

#' @export
classificationResult <- function(observed, predicted, probs) {
  if (length(levels(as.factor(observed))) == 2) {
    TwoWayClassificationResult(observed,predicted, probs)
  } else if (length(levels(as.factor(observed))) > 2) {
    MultiWayClassificationResult(observed,predicted, probs)
  } else {
    stop("observed data must be a factor with 2 or more levels")
  }
}


#' @export
RawModel <- function(model, Xtrain, Ytrain, Xtest, Ytest, tuneGrid) {
  fit <- model$fit(Xtrain, Ytrain, NULL, tuneGrid, classProbs=TRUE)  
  ret <- list(
              model=model,
              Xtrain=Xtrain,
              Ytrain=Ytrain,
              Xtest=Xtest,
              Ytest=Ytest,
              tuneGrid=tuneGrid,
              modelFit=fit)
  
  class(ret) <- c("RawModel", "list")
  ret
}

#' @export
RawPredictor <- function(fit, model) {
  ret <- list(fit=fit,
              model=model)
  
  class(ret) <- c("RawPredictor", "list")
  ret
}

#' @export
CaretPredictor <- function(fit) {
  ret <- list(fit=fit)
  
  class(ret) <- c("CaretPredictor", "list")
  ret
}


#' @export
MVPAVoxelPredictor <- function(predictor, voxelGrid) {
  ret <- list(predictor=predictor, voxelGrid=voxelGrid)
  
  class(ret) <- c("MVPAVoxelPredictor", "list")
  ret
}


#' @export
ListModel <- function(fits) {
  stopifnot(is.list(fits))
  ret <- fits
  class(ret) <- c("ListModel", "list")
  ret
}


#' @export
CaretModel <- function(model, Xtrain, Ytrain, Xtest, Ytest, tuneGrid, tuneControl,...) {
  fit <- caret::train(Xtrain, Ytrain, method=model, trControl=tuneControl, tuneGrid=tuneGrid, ...)
  ret <- list(
              model=model,
              Xtrain=Xtrain,
              Ytrain=Ytrain,
              Xtest=Xtest,
              Ytest=Ytest,
              tuneGrid=tuneGrid,
              tuneControl=tuneControl,
              modelFit=fit)
  
  class(ret) <- c("CaretModel", "list")
  ret
}



#' @export
evaluateModel <- function(x, newdata) {
  UseMethod("evaluateModel")
}

#' @export
asPredictor <- function(x, voxelGrid) {
  UseMethod("asPredictor")
}

#' @export
asPredictor.RawModel <- function(x) {
  RawPredictor(x$modelFit, x$model)
}

#' @export
asPredictor.CaretModel <- function(x, voxelGrid, brainSpace) {
  CaretPredictor(x$modelFit)
}


#' @export
evaluateModel.RawModel <- function(x, newdata=NULL) {
  if (is.null(newdata)) {
    newdata=x$Xtest
  }
  
  probs <- x$model$prob(x$modelFit, newdata)     
  cpred <- apply(probs,1, which.max)
  cpred <- levels(x$Ytrain)[cpred]
  list(class=cpred, probs=probs)
}

#' @export
evaluateModel.MVPAVoxelPredictor <- function(x, newdata) {
  if (is.null(newdata)) {
    stop("newdata cannot be null")
  }
  
  X <- series(newdata,x$voxelGrid)
  evaluateModel(x$predictor, X)
}


#' @export
evaluateModel.RawPredictor <- function(x, newdata=NULL) {
  if (is.null(newdata)) {
    stop("newdata cannot be null")
  }
  
  probs <- x$model$prob(x$fit, newdata=newdata)     
  cpred <- apply(probs,1, which.max)
  cpred <- levels(x$Ytrain)[cpred]
  list(class=cpred, probs=probs)
}


#' @export
evaluateModel.CaretModel <- function(x, newdata=NULL) {
  if (is.null(newdata)) {
    newdata=x$Xtest
  }
  
  list(class=predict(x$modelFit, newdata), probs=predict(x$modelFit, newdata, type="prob"))
}



#' @export
evaluateModel.ListModel <- function(x, newdata=NULL) {
  if (is.null(newdata)) {
    stop("newdata cannot be null")
  }
  
  res <- lapply(x, function(fit) {
    evaluateModel(fit, newdata)$probs
  })
  
  prob <- Reduce("+", res)/length(res)
  winner <- apply(prob, 1, which.max)
  class <- colnames(prob)[winner]
  
  list(class=class, probs=prob)
}

#' @export
evaluateModel.CaretPredictor <- function(x, newdata=NULL) {
  if (is.null(newdata)) {
    stop("newdata cannot be null")
  }
  
  list(class=predict(x$fit, newdata=newdata), probs=predict(x$fit, newdata=newdata, type="prob"))
}

#' @export
evaluateModelList <- function(modelList, newdata=NULL) {
  results <- lapply(modelList, evaluateModel, newdata)
  probMat <- do.call(rbind, lapply(results, "[[", "probs"))
  predClass <- unlist(lapply(results, "[[", "class"))
  list(class=predClass, probs=probMat)
}


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
    lev <- colnames(x$prob)[i]
    pos <- obs == lev
    pclass <- x$prob[,i]
    pother <- rowMeans(x$prob[,-i])
    Metrics::auc(as.numeric(pos), pclass - pother)-.5
  })
  
  names(aucres) <- paste0("AUC_", colnames(x$prob))

  c(ZAccuracy=-qnorm(out$p.value), Accuracy=sum(obs == as.character(x$predicted))/length(obs), Combined_AUC=mean(aucres), aucres)
}









#' @export
fitFinalModel <- function(Xtrain, Ytrain,  method, Xtest=NULL, Ytest=NULL, tuneGrid=NULL, tuneControl=NULL, ...) {  
 
  if (is.null(tuneGrid)) {
    tuneGrid <- method$grid(Xtrain, Ytrain, 1)
  }
  
  if (nrow(tuneGrid) == 1) {
    RawModel(method, Xtrain, Ytrain, Xtest, Ytest, tuneGrid)   
  } else {
    if (is.null(tuneControl)) {
      tuneControl <- .cvControl
    }
    CaretModel(method, Xtrain, Ytrain, Xtest, Ytest, tuneGrid, tuneControl, ...)
  }
}

.get_lapply <- function(ncores=1) {
  if (ncores > 1) {
    function(sets, FUN, ...) {
      parallel::mclapply(sets, FUN, ..., mc.cores=ncores)
    }
  } else {
    lapply      
  } 
}

trainModel <- function(model, Xtrain, Ytrain, Xtest, Ytest, tuneGrid, fast=TRUE, tuneControl=.noneControl) {
  ret <- if (nrow(tuneGrid) == 1 && fast) {
    RawModel(model, Xtrain, Ytrain, Xtest, Ytest, tuneGrid)        
  } else {     
    CaretModel(model, Xtrain, Ytrain, Xtest, Ytest, tuneGrid, tuneControl)
  }  
}

#' @export
cvTrainAndEvalExternal <- function(foldIterator, Xtest, Ytest, model, tuneGrid, fast=TRUE, ncores=1, returnPredictor=FALSE) {
  if (!is.factor(Y)) {
    stop("regression not supported yet") 
  }
  
  results <- if (nrow(tuneGrid) == 1) {
    model <- trainModel(model, foldIterator$X, foldIterator$Y, Xtest, Ytest, tuneGrid, fast, .noneControl)
    perf <- evaluateModel(model)
    list(perf=perf, predictor=asPredictor(model))
  } else {
    index <- invertFolds(foldIterator$getTestSets(), nrow(foldIterator$X)) 
    ctrl <- caret::trainControl("cv", verboseIter=TRUE, classProbs=TRUE, index=index)
    model <- trainModel(model, foldIterator$X, foldIterator$Y, Xtest, Ytest, tuneGrid, fast=FALSE, ctrl)
    perf <- evaluateModel(model)
    list(perf=perf, predictor=asPredictor(model))              
  } 
  
  if (returnPredictor) {
    list(class=results$class, probs=results$probs, predictor=results$predictor)
  } else {
    list(class=results$class, probs=results$probs)
  }
}

#' @export
cvTrainAndEvalInternal <- function(foldIterator, model, tuneGrid, fast=TRUE, ncores=1, returnPredictor=FALSE) {
  if (!is.factor(Y)) {
    stop("regression not supported yet") 
  }
  
  results <- foreach::foreach(fold = foldIterator, .verbose=FALSE) %do% {   
    if (nrow(tuneGrid) == 1 && fast) {
      model <- trainModel(model, fold$Xtrain, fold$Ytrain, fold$Xtest, fold$Ytest, tuneGrid, fast, .noneControl)
      evaluateModel(model)        
    } else {      
      if (nrow(tuneGrid) == 1) {
        model <- trainModel(model, fold$Xtrain, fold$Y, fold$Xtest, fold$YtestY, tuneGrid, fast=FALSE, tuneControl=.noneControl)
        evaluateModel(model)
      } else {
        index <- invertFolds(foldIterator$getTestSets()[-foldIterator$index()], nrow(fold$Xtrain)) 
        ctrl <- caret::trainControl("cv", verboseIter=TRUE, classProbs=TRUE, index=index)
        model <- trainModel(model, fold$Xtrain, fold$Ytrain, fold$Xtest, fold$Ytest, tuneGrid, fast=FALSE, tuneControl=ctrl)
        evaluateModel(model)
      }
    
    }
  }
  
  probMat <- do.call(rbind, lapply(results, "[[", "probs"))
  predClass <- unlist(lapply(results, "[[", "class"))
  
  if (returnPredictor) {
    ctrl <- if (nrow(tuneGrid) == 1) {
      .noneControl
    } else {
      caret::trainControl("cv", verboseIter=TRUE, classProbs=TRUE, index=foldIterator$getTrainSets())
    }    
    
    model <- trainModel(model, foldIterator$X, foldIterator$Y, NULL, NULL, tuneGrid, fast=FALSE, tuneControl=ctrl)
    list(class=predClass, probs=probMat, predictor=asPredictor(model))
  }
}
  
  
#   results <- 
#     .lapply(seq_along(testSets), function(blockIndex) {
#       testIndices <- testSets[[blockIndex]]
#       Xtrain <- X[-testIndices,]
#       Ytrain <- Y[-testIndices]
#       Xtest <- X[testIndices,]
#       Ytest <- Y[testIndices]
#       
#       ret <- if (nrow(tuneGrid) == 1 && fast) {
#         evaluateModel(trainModel(model, Xtrain, Y[-testIndices], X[testIndices,], Y[testIndices], tuneGrid, fast, .noneControl))        
#       } else {      
#         if (nrow(tuneGrid) == 1) {
#           evaluateModel(trainModel(model, Xtrain, Y[-testIndices], X[testIndices,], Y[testIndices], tuneGrid, fast=FALSE, tuneControl=.noneControl))
#         } else {       
#           index <- invertFolds(testSets[-blockIndex], nrow(Xtrain)) 
#           ctrl <- caret::trainControl("cv", verboseIter=TRUE, classProbs=TRUE, index=index)
#           evaluateModel(trainModel(model, Xtrain, Y[-testIndices], X[testIndices,], Y[testIndices], tuneGrid, fast=FALSE, tuneControl=ctrl))
#         }
#        
#       }
#       
#     })
#   
#   probMat <- do.call(rbind, lapply(results, "[[", "probs"))
#   predClass <- unlist(lapply(results, "[[", "class"))
#   list(class=predClass, probs=probMat)
# }
#   
#   
  
#' @export
crossValidate <- function(X, Y, trainSets, testSets, model, tuneGrid, fast=TRUE, finalFit=FALSE, ncores=1, ...) {
 
  if (!is.factor(Y)) {
    stop("regression not supported yet") 
  }
   
  .lapply <- .get_lapply(ncores)
  
  blockFits <- 
    .lapply(seq_along(testSets), function(blockIndex) {
      testIndices <- testSets[[blockIndex]]
      Xtrain <- X[-testIndices,]
      
      ret <- if (nrow(tuneGrid) == 1 && fast) {
        RawModel(model, Xtrain, Y[-testIndices], X[testIndices,], Y[testIndices], tuneGrid)        
      } else { 
          
        if (nrow(tuneGrid) == 1) {
          CaretModel(model, Xtrain, Y[-testIndices], X[testIndices,], Y[testIndices], tuneGrid, .noneControl)
        } else {
  
          index <- invertFolds(testSets[-blockIndex], nrow(Xtrain)) 
          ctrl <- caret::trainControl("cv", verboseIter=TRUE, classProbs=TRUE, index=index)
          CaretModel(model, Xtrain, Y[-testIndices], X[testIndices,], Y[testIndices], tuneGrid, ctrl)
        }
        
      }
    })
  
  
  final <- if (finalFit) {
    ## fit final model to whole data set
    ctrl <- if (nrow(tuneGrid) == 1) {
      .noneControl
    } else {
      caret::trainControl("cv", verboseIter=TRUE, classProbs=TRUE, index=trainSets)
    }    
    
    CaretModel(model, X, Y, NULL, NULL, tuneGrid, ctrl)
  } 
  
  result <- evaluateModelList(blockFits)
  list(class=result$class, prob=result$prob, observed=Y, blockFits=blockFits, finalFit=ListPredictor(blockFits))
   
}

#' @export
zeroVarianceColumns <- function(M) {
  which(apply(M, 2, sd, na.rm=TRUE) == 0)
}

#' @export
nonzeroVarianceColumns <- function(M) {
  which(apply(M, 2, sd, na.rm=TRUE) > 0)
}

#' @export
removeZeroVarianceColumns <- function(M) {
  noVariance <- which(apply(M, 2, sd, na.rm=TRUE) == 0)
  if (length(noVariance) > 0) {
    M[, hasVariance, drop=FALSE]
  } else {
    M
  }
}

#' @import neuroim
#' @export
fitMVPAModel <- function(dataset, voxelGrid, tuneLength=1, fast=TRUE, finalFit=FALSE, ncores=1) {
  
  M <- series(dataset$trainVec, voxelGrid) 
  
  if (ncol(M) < 2) {
    stop("feature matrix must have at least two columns: returning NullResult")
  }
  

  hasVariance <- which(apply(M, 2, sd, na.rm=TRUE) > 0)
  M <- M[, hasVariance, drop=FALSE]
  
  hasNA <- apply(M,1, function(x) any(is.na(x)))
  numNAs <- sum(hasNA)
  
  if (numNAs > 0) {
    stop("training data has NA values, aborting")
  } 
  
  tuneGrid <- if (is.null(dataset$tuneGrid)) {
    tuneGrid <- dataset$model$grid(M, dataset$Y, tuneLength)
  } else {
    dataset$tuneGrid
  }
  
  if (ncol(M) < 2) {
    stop("feature matrix must have at least two columns with nonzero variance")
  }
  
  voxelGrid <- voxelGrid[hasVariance, ]
  
  result <- if (is.null(dataset$testVec)) {
    crossValidate(M, dataset$Y, dataset$trainSets, dataset$testSets, dataset$model, tuneGrid, tuneLength, testVec=dataset$testVec, testY=dataset$testY, fast=fast,finalFit=finalFit, ncores=ncores)  
  } else {
    Xtest <- series(dataset$testVec, voxelGrid)
    modelFit <- if (nrow(tuneGrid) == 1) {
      fitFinalModel(M, dataset$Y,  dataset$model, Xtest, dataset$testY, tuneGrid)
    } else {      
      ctrl <- caret::trainControl("cv", verboseIter=TRUE, classProbs=TRUE, index=dataset$trainSets) 
      fitFinalModel(M, dataset$Y,  dataset$model, dataset$Xtest, dataset$testY, tuneGrid, tuneControl=ctrl)     
    }
    
    preds <- evaluateModel(modelFit)
    list(class=preds$class, prob=preds$prob, observed=dataset$testY, blockFits=NULL, finalFit=modelFit)
  }
  
  if (is.factor(dataset$Y) && length(levels(dataset$Y)) == 2) {
    TwoWayClassificationResult(voxelGrid, dataset$model, result$observed, result$class, result$prob, result$blockFits, result$finalFit)
  } else if (is.factor(dataset$Y) && length(levels(dataset$Y)) >= 3) {
    MultiWayClassificationResult(voxelGrid, dataset$model, result$observed, result$class, result$prob, result$blockFits, result$finalFit)
  } else {
    stop("regression not supported yet.")
  }
}


