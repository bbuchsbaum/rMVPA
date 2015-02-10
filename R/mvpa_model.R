


.noneControl <- caret::trainControl("none", verboseIter=TRUE, classProbs=TRUE, returnData=FALSE, returnResamp="none")
.cvControl <- caret::trainControl("cv", verboseIter=TRUE, classProbs=TRUE, returnData=FALSE, returnResamp="none")  


#' create an \code{ClassificationModel} instance
#' @param caretModel the underlying caret model object
#' @export
ClassificationModel <- function(caretModel) {
  class(caretModel) <- c("ClassificationModel", "list")
  caretModel
}

#' create an \code{EnsembleSearchlightModel} instance
#' @param caretModel
#' @export
EnsembleSearchlightModel <- function(baseLearners=list(pls=data.frame(ncomp=1:4))) {
  ret=list(baseLearners=baseLearners)
  class(ret) <- c("EnsembleSearchlightModel", "list")
  ret
}

#' create an \code{SimilarityModel} instance
#' @param type the similarity metric (pearson, kendall, spearman)
#' @export
SimilarityModel <- function(type=c("pearson", "spearman", "kendall")) {
  simFun <- if (is.function(type)) {
    type
  } else {
    type=match.arg(type)
    f <- function(x) {
      cor(x, method=type)
    }
  }
  ret <- list(simFun=simFun)
  class(ret) <- c("SimilarityModel", "list")
  ret
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
MVPADataset <- function(trainVec, Y, mask, blockVar, testVec, testY, modelName="corclass", tuneGrid=NULL, testSplitVar=NULL, testSplits=NULL) {
  
  model <- loadModel(modelName)
  
 
  testSets <- split(1:length(blockVar), blockVar)
  trainSets <- invertFolds(testSets, 1:length(blockVar))
   
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
    trainSets=trainSets,
    testSplitVar=testSplitVar,
    testSplits=testSplits)
   
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

#' @export
SimilarityResult <- function(sWithin, sBetween, simMat, simWithinTable) {
  ret <- list(sWithin=sWithin, sBetween=sBetween, simMat=simMat, simWithinTable)
  class(ret) <- c("SimilarityResult", "list")
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
classificationResult <- function(observed, predicted, probs, predictor=NULL) {
  if (length(levels(as.factor(observed))) == 2) {
    TwoWayClassificationResult(observed,predicted, probs, predictor)
  } else if (length(levels(as.factor(observed))) > 2) {
    MultiWayClassificationResult(observed,predicted, probs, predictor)
  } else {
    stop("observed data must be a factor with 2 or more levels")
  }
}


.purgeModel <- function(fit, label) {
  if (label == "Partial Least Squares") {
    fit$model = NULL
  }
  
  fit
}

#' @export
CaretModel <- function(model, Xtrain, Ytrain, Xtest, Ytest, tuneGrid, tuneControl,...) {
  
  if (model$library[1] == "gbm" && length(levels(Ytrain)) == 2) {
    fit <- caret::train(Xtrain, Ytrain, method=model, trControl=tuneControl, tuneGrid=tuneGrid, distribution="bernoulli", ...)
  } else if (model$library[1] == "gbm" && length(levels(Ytrain)) > 2) {
    fit <- caret::train(Xtrain, Ytrain, method=model, trControl=tuneControl, tuneGrid=tuneGrid, distribution="multinomial", ...)
  } else {
    fit <- caret::train(Xtrain, Ytrain, method=model, trControl=tuneControl, tuneGrid=tuneGrid, ...)
  }
  
  fit$finalModel <- .purgeModel(fit$finalModel, model$label)
  
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
RawModel <- function(model, Xtrain, Ytrain, Xtest, Ytest, tuneGrid) {
 
  fit <- model$fit(Xtrain, Ytrain, NULL, tuneGrid, lev=levels(Ytrain), classProbs=TRUE)
  fit <- .purgeModel(fit, model$label)
  
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
CalibratedPredictor <- function(predictor, calibrationX, calibrationY) {
  fits <- lapply(1:length(levels(calibrationY)), function(i) {
    print(i)
    df1=data.frame(y0=as.numeric(calibrationY == levels(calibrationY)[i]), x0=calibrationX[,i])
    #glm(y0 ~ x0, data=df1, family="binomial")
    ## scam
    gam(y0 ~ s(x0), data=df1, family=binomial())
  })
      
  ret <- list(calfits=fits, predictor=predictor)
  class(ret) <- c("CalibratedPredictor", "list")
  ret
}


#' @export
ListPredictor <- function(fits, names=1:length(fits)) {
  stopifnot(is.list(fits))
  ret <- fits
  names(ret) <- names
  
  class(ret) <- c("ListPredictor", "list")
  ret
}

WeightedPredictor <- function(fits, names=1:length(fits), weights=rep(1/length(fits), length(fits))) {
  ret <- fits
  names(ret) <- names
  attr(ret, "weights") <- weights
  class(ret) <- c("WeightedPredictor", "list")
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
  
  cpred <- apply(probs,1, function(x) {
    x[is.na(x)] <- 0
    which.max(x)    
  })
  
  cpred <- levels(x$Ytrain)[cpred]
  list(class=cpred, probs=probs)
}

#' @export
evaluateModel.CalibratedPredictor <- function(x, newdata) {
  if (is.null(newdata)) {
    stop("newdata cannot be null")
  }
  
  Preds <- evaluateModel(x$predictor, newdata)$probs
  Pcal <- do.call(cbind, lapply(1:ncol(Preds), function(i) {
    predict(x$calfits[[i]], data.frame(x0=Preds[,i]), type="response")
  }))
  
  sweep(Pcal, 1, rowSums(Pcal), "/")

}

#' @export
evaluateModel.MVPAVoxelPredictor <- function(x, newdata) {
  if (is.null(newdata)) {
    stop("newdata cannot be null")
  }
  
  X <- series(newdata,x$voxelGrid)
  
  if (is.vector(X)) {
    X <- matrix(X, nrow=1)
  }
  
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
evaluateModel.ListPredictor <- function(x, newdata=NULL) {
  if (is.null(newdata)) {
    stop("newdata cannot be null")
  }
  
  res <- lapply(x, function(fit) {
    evaluateModel(fit, newdata)
  })
  
}

#' @export
evaluateModel.WeightedPredictor <- function(x, newdata=NULL) {
  if (is.null(newdata)) {
    stop("newdata cannot be null")
  }
  
  wts <- attr(x, "weights")
  preds <- lapply(1:length(x), function(i) {
    evaluateModel(x[[i]], newdata)$probs * wts[i]
  })
  
  prob <- preds[!sapply(preds, function(x) is.null(x))]
  pfinal <- Reduce("+", prob)
  
  pclass <- names(pfinal)[apply(pfinal, 1, which.max)]
  list(class=pclass, prob=pfinal)
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

#' @export
trainModel <- function(model, Xtrain, Ytrain, Xtest, Ytest, tuneGrid, fast=TRUE, tuneControl=.noneControl) {
  ret <- if (nrow(tuneGrid) == 1 && fast) {
    RawModel(model, Xtrain, Ytrain, Xtest, Ytest, tuneGrid)        
  } else {     
    CaretModel(model, Xtrain, Ytrain, Xtest, Ytest, tuneGrid, tuneControl)
  }  
}

#' @export
#' @import foreach
crossval_external <- function(foldIterator, Xtest, Ytest, model, tuneGrid, fast=TRUE, ncores=1, returnPredictor=FALSE) {
 
  results <- if (nrow(tuneGrid) == 1) {
    fit <- trainModel(model, foldIterator$X, foldIterator$Y, Xtest, Ytest, tuneGrid, fast, .noneControl)
    perf <- evaluateModel(fit)
    list(perf=perf, predictor=asPredictor(fit))
  } else {
    index <- invertFolds(foldIterator$getTestSets(), 1:length(foldIterator$blockVar)) 
    ctrl <- caret::trainControl("cv", verboseIter=TRUE, classProbs=TRUE, index=index, returnData=FALSE, returnResamp="none")
    fit <- trainModel(model, foldIterator$X, foldIterator$Y, Xtest, Ytest, tuneGrid, fast=FALSE, ctrl)
    perf <- evaluateModel(fit)
    list(perf=perf, predictor=asPredictor(fit))              
  } 
  
  if (returnPredictor) {
    list(class=results$perf$class, probs=results$perf$probs, predictor=results$predictor)
  } else {
    list(class=results$perf$class, probs=results$perf$probs, predictor=NULL)
  }
}

#' @export
#' @import foreach
crossval_internal <- function(foldIterator, model, tuneGrid, fast=TRUE, ncores=1, returnPredictor=FALSE) {
 
  results <- foreach::foreach(fold = foldIterator, .verbose=FALSE, .packages=c(model$library)) %do% {   
    if (nrow(tuneGrid) == 1 && fast) {
      fit <- trainModel(model, fold$Xtrain, fold$Ytrain, fold$Xtest, fold$Ytest, tuneGrid, fast, .noneControl)
      evaluateModel(fit)        
    } else {      
      if (nrow(tuneGrid) == 1) {
        fit <- trainModel(model, fold$Xtrain, fold$Ytrain, fold$Xtest, fold$Ytest, tuneGrid, fast=FALSE, tuneControl=.noneControl)
        evaluateModel(fit)
      } else {
        index <- invertFolds(foldIterator$getTestSets()[-foldIterator$index()], 1:nrow(fold$Xtrain)) 
        ctrl <- caret::trainControl("cv", verboseIter=TRUE, classProbs=TRUE, index=index, returnData=FALSE, returnResamp="none")
        fit <- trainModel(model, fold$Xtrain, fold$Ytrain, fold$Xtest, fold$Ytest, tuneGrid, fast=FALSE, tuneControl=ctrl)
        evaluateModel(fit)
      }
    
    }
  }
  
  ## reorder predictions to match order of input features/labels
  ord <- foldIterator$getTestOrder()
  probMat <- do.call(rbind, lapply(results, "[[", "probs"))[ord,]
  predClass <- unlist(lapply(results, "[[", "class"))[ord]
  
  if (returnPredictor) {
    ctrl <- if (nrow(tuneGrid) == 1) {
      .noneControl
    } else {
      caret::trainControl("cv", verboseIter=TRUE, classProbs=TRUE, index=foldIterator$getTrainSets(), returnData=FALSE, returnResamp="none")
    }    
    
    fit <- trainModel(model, foldIterator$X, foldIterator$Y, NULL, NULL, tuneGrid, fast=FALSE, tuneControl=ctrl)
    list(class=predClass, probs=probMat, predictor=asPredictor(fit))
  } else {
    list(class=predClass, probs=probMat, predictor=NULL)
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
# crossValidate <- function(X, Y, trainSets, testSets, model, tuneGrid, fast=TRUE, finalFit=FALSE, ncores=1, ...) {
#  
#   if (!is.factor(Y)) {
#     stop("regression not supported yet") 
#   }
#    
#   .lapply <- .get_lapply(ncores)
#   
#   blockFits <- 
#     .lapply(seq_along(testSets), function(blockIndex) {
#       testIndices <- testSets[[blockIndex]]
#       Xtrain <- X[-testIndices,]
#       
#       ret <- if (nrow(tuneGrid) == 1 && fast) {
#         RawModel(model, Xtrain, Y[-testIndices], X[testIndices,], Y[testIndices], tuneGrid)        
#       } else { 
#           
#         if (nrow(tuneGrid) == 1) {
#           CaretModel(model, Xtrain, Y[-testIndices], X[testIndices,], Y[testIndices], tuneGrid, .noneControl)
#         } else {
#   
#           index <- invertFolds(testSets[-blockIndex], nrow(Xtrain)) 
#           ctrl <- caret::trainControl("cv", verboseIter=TRUE, classProbs=TRUE, index=index)
#           CaretModel(model, Xtrain, Y[-testIndices], X[testIndices,], Y[testIndices], tuneGrid, ctrl)
#         }
#         
#       }
#     })
#   
#   
#   final <- if (finalFit) {
#     ## fit final model to whole data set
#     ctrl <- if (nrow(tuneGrid) == 1) {
#       .noneControl
#     } else {
#       caret::trainControl("cv", verboseIter=TRUE, classProbs=TRUE, index=trainSets)
#     }    
#     
#     CaretModel(model, X, Y, NULL, NULL, tuneGrid, ctrl)
#   } 
#   
#   result <- evaluateModelList(blockFits)
#   list(class=result$class, prob=result$prob, observed=Y, blockFits=blockFits, finalFit=ListPredictor(blockFits))
#    
# }

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

# #' @import neuroim
# #' @export
# fitMVPAModel <- function(dataset, voxelGrid, tuneLength=1, fast=TRUE, finalFit=FALSE, ncores=1) {
#   
#   M <- series(dataset$trainVec, voxelGrid) 
#   
#   if (ncol(M) < 2) {
#     stop("feature matrix must have at least two columns: returning NullResult")
#   }
#   
# 
#   hasVariance <- which(apply(M, 2, sd, na.rm=TRUE) > 0)
#   M <- M[, hasVariance, drop=FALSE]
#   
#   hasNA <- apply(M,1, function(x) any(is.na(x)))
#   numNAs <- sum(hasNA)
#   
#   if (numNAs > 0) {
#     stop("training data has NA values, aborting")
#   } 
#   
#   tuneGrid <- if (is.null(dataset$tuneGrid)) {
#     tuneGrid <- dataset$model$grid(M, dataset$Y, tuneLength)
#   } else {
#     dataset$tuneGrid
#   }
#   
#   if (ncol(M) < 2) {
#     stop("feature matrix must have at least two columns with nonzero variance")
#   }
#   
#   voxelGrid <- voxelGrid[hasVariance, ]
#   
#   result <- if (is.null(dataset$testVec)) {
#     crossValidate(M, dataset$Y, dataset$trainSets, dataset$testSets, dataset$model, tuneGrid, tuneLength, testVec=dataset$testVec, testY=dataset$testY, fast=fast,finalFit=finalFit, ncores=ncores)  
#   } else {
#     Xtest <- series(dataset$testVec, voxelGrid)
#     modelFit <- if (nrow(tuneGrid) == 1) {
#       fitFinalModel(M, dataset$Y,  dataset$model, Xtest, dataset$testY, tuneGrid)
#     } else {      
#       ctrl <- caret::trainControl("cv", verboseIter=TRUE, classProbs=TRUE, index=dataset$trainSets) 
#       fitFinalModel(M, dataset$Y,  dataset$model, dataset$Xtest, dataset$testY, tuneGrid, tuneControl=ctrl)     
#     }
#     
#     preds <- evaluateModel(modelFit)
#     list(class=preds$class, prob=preds$prob, observed=dataset$testY, blockFits=NULL, finalFit=modelFit)
#   }
#   
#   if (is.factor(dataset$Y) && length(levels(dataset$Y)) == 2) {
#     TwoWayClassificationResult(voxelGrid, dataset$model, result$observed, result$class, result$prob, result$blockFits, result$finalFit)
#   } else if (is.factor(dataset$Y) && length(levels(dataset$Y)) >= 3) {
#     MultiWayClassificationResult(voxelGrid, dataset$model, result$observed, result$class, result$prob, result$blockFits, result$finalFit)
#   } else {
#     stop("regression not supported yet.")
#   }
# }


