
.noneControl <- caret::trainControl("none", verboseIter=TRUE, classProbs=TRUE, returnData=FALSE, returnResamp="none")
.cvControl <- caret::trainControl("cv", verboseIter=TRUE, classProbs=TRUE, returnData=FALSE, returnResamp="none")  
.adaptiveControl <- caret::trainControl("adaptive_cv", verboseIter=TRUE, classProbs=TRUE, returnData=FALSE, returnResamp="none")  


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
EnsembleSearchlightModel <- function(baseLearners=list(pls=data.frame(ncomp=1:5), lda_thomaz=1, sda=expand.grid(lambda=c(.01,.2,.5,.7), diagonal=FALSE))) {
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
#' @param tuneGrid
#' @export
MVPADataset <- function(trainVec, Y, mask, blockVar, testVec, testY, modelName="corclass", tuneGrid=NULL, tuneLength=1,testSplitVar=NULL,
                        testSplits=NULL, trainDesign=NULL, testDesign=NULL) {
  
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
    tuneLength=tuneLength,
    testSets=testSets,
    trainSets=trainSets,
    testSplitVar=testSplitVar,
    testSplits=testSplits,
    trainDesign=trainDesign,
    testDesign=testDesign)
   
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
  ret <- list(sWithin=sWithin, sBetween=sBetween, simMat=simMat, simWithinTable=simWithinTable)
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

#' create an \code{MultiWayClassification} instance
#' @param observed
#' @param predicted
#' @param probs
#' @export
MultiWayClassificationResult <- function(observed, predicted, probs, predictor=NULL) {
  ret <- list(
              observed=observed,
              predicted=predicted,
              probs=as.matrix(probs),
              predictor=predictor)
  
  class(ret) <- c("MultiWayClassificationResult", "list")
  ret
  
}

#' @export
classificationResult <- function(observed, predicted, probs,  predictor=NULL) {
  if (length(levels(as.factor(observed))) == 2) {
    TwoWayClassificationResult(observed,predicted, probs,  predictor)
  } else if (length(levels(as.factor(observed))) > 2) {
    MultiWayClassificationResult(observed,predicted, probs, predictor)
  } else {
    stop("observed data must be a factor with 2 or more levels")
  }
}


#' @export
CaretModel <- function(model, Xtrain, Ytrain, Xtest, Ytest, tuneGrid, tuneControl, featureSelector=NULL, parcels=NULL, ...) {
  if (!is.null(parcels)) {
    Xtrain <- groupMeans(Xtrain, 2, parcels)
  }
    
  featureMask <- if (!is.null(featureSelector)) {
    selectFeatures(featureSelector, Xtrain, Ytrain)
  } else {
    rep(TRUE, ncol(Xtrain))
  }
  
  Xtrain <- Xtrain[,featureMask]
  
  if (model$library[1] == "gbm" && length(levels(Ytrain)) == 2) {
    fit <- caret::train(Xtrain, Ytrain, method=model, trControl=tuneControl, tuneGrid=tuneGrid, distribution="bernoulli", ...)
  } else if (model$library[1] == "gbm" && length(levels(Ytrain)) > 2) {
    fit <- caret::train(Xtrain, Ytrain, method=model, trControl=tuneControl, tuneGrid=tuneGrid, distribution="multinomial", ...)
  } else {
    fit <- caret::train(Xtrain, Ytrain, method=model, trControl=tuneControl, tuneGrid=tuneGrid, ...)
  }
  
  #fit$finalModel <- cleanup(fit$finalModel)
  
  ret <- list(
    model=model,
    Xtrain=Xtrain,
    Ytrain=Ytrain,
    Xtest=Xtest,
    Ytest=Ytest,
    tuneGrid=tuneGrid,
    tuneControl=tuneControl,
    modelFit=fit,
    featureSelector=featureSelector,
    featureMask=featureMask,
    parcels=parcels)
  
  class(ret) <- c("CaretModel", "list")
  ret
}


#' @export
RawModel <- function(model, Xtrain, Ytrain, Xtest, Ytest, tuneGrid, featureSelector=NULL, parcels=NULL) {
  
  if (!is.null(parcels)) {
    message("computing parcellation", dim(Xtrain))
    Xtrain <- groupMeans(Xtrain, 2, parcels)
    message("computing parcellation", dim(Xtrain))
  }
  
  featureMask <- if (!is.null(featureSelector)) {
    selectFeatures(featureSelector, Xtrain, Ytrain)
  } else {
    rep(TRUE, ncol(Xtrain))
  }
 
  Xtrain <- Xtrain[,featureMask]
  
  fit <- model$fit(Xtrain, Ytrain, NULL, tuneGrid, lev=levels(Ytrain), classProbs=TRUE)
  #fit <- cleanup(fit)
  
  ret <- list(
              model=model,
              Xtrain=Xtrain,
              Ytrain=Ytrain,
              Xtest=Xtest,
              Ytest=Ytest,
              tuneGrid=tuneGrid,
              modelFit=fit,
              featureSelector=featureSelector,
              featureMask=featureMask,
              parcels=parcels)
              
  
  class(ret) <- c("RawModel", "list")
  ret
}

#' @export
RawPredictor <- function(fit, model, featureMask, parcels=NULL) {
  ret <- list(fit=fit,
              model=model,
              featureMask=featureMask,
              parcels=parcels)
  
  class(ret) <- c("RawPredictor", "list")
  ret
}

#' @export
CaretPredictor <- function(fit, featureMask, parcels=NULL) {
  ret <- list(fit=fit, featureMask=featureMask, parcels=parcels)
  
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


#' @export
WeightedPredictor <- function(fits, names=1:length(fits), weights=rep(1/length(fits), length(fits))) {
  stopifnot(length(weights) == length(fits))
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
  RawPredictor(x$modelFit, x$model, x$featureMask, x$parcels)
}

#' @export
asPredictor.CaretModel <- function(x) {
  CaretPredictor(x$modelFit, x$featureMask, x$parcels)
}


#' @export
evaluateModel.RawModel <- function(x, newdata=NULL) {
  if (is.null(newdata)) {
    newdata=x$Xtest
  }
  
  probs <- if (is.null(x$parcels)) {
    x$model$prob(x$modelFit, newdata[,x$featureMask])  
  } else {
    reducedData <- groupMeans(newdata, 2, x$parcels)
    x$model$prob(x$modelFit, reducedData[,x$featureMask])  
  }
  
  cpred <- apply(probs,1, function(x) {
    x[is.na(x)] <- 0
    which.max(x)    
  })
  
  cpred <- colnames(probs)[cpred]
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
  
  probs <- if (is.null(x$parcels)) {
    x$model$prob(x$fit, newdata=newdata[,x$featureMask])     
  } else {
    reducedData <- groupMeans(newdata, 2, x$parcels)
    x$model$prob(x$modelFit, reducedData[,x$featureMask])  
  }
  
  cpred <- apply(probs,1, which.max)
  cpred <- colnames(probs)[cpred]
  list(class=cpred, probs=probs)
}


#' @export
evaluateModel.CaretModel <- function(x, newdata=NULL) {
  if (is.null(newdata)) {
    newdata=x$Xtest
  }
  
  if (!is.null(x$parcels)) {
    newdata <- groupMeans(newdata, 2, x$parcels)
  }
  
  list(class=predict(x$modelFit, newdata[,x$featureMask]), probs=predict(x$modelFit, newdata[,x$featureMask], type="prob"))
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
    evaluateModel(x[[i]], newdata)$prob * wts[i]
  })
  
  prob <- preds[!sapply(preds, function(x) is.null(x))]
  pfinal <- Reduce("+", prob)
  
  
  ## TODO error ...
  ## Error in colnames(pfinal)[apply(pfinal, 1, which.max)] : 
  ##  invalid subscript type 'list'
  
  cnames <- colnames(pfinal)
  maxids <- apply(pfinal, 1, which.max)
  len <- sapply(maxids, length)
  if (any(len == 0)) {
    maxids[len == 0] <- NA
  }
  maxids <- unlist(maxids)
  
  pclass <- cnames[maxids]
  list(class=pclass, prob=pfinal)
}



#' @export
evaluateModel.CaretPredictor <- function(x, newdata=NULL) {
  if (is.null(newdata)) {
    stop("newdata cannot be null")
  }
  
  if (!is.null(x$parcels)) {
    newdata <- groupMeans(newdata, 2, x$parcels)
  }
  
  list(class=predict(x$fit, newdata=newdata[,x$featureMask]), probs=predict(x$fit, newdata=newdata[,x$featureMask], type="prob"))
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
trainModel <- function(model, Xtrain, Ytrain, Xtest, Ytest, tuneGrid, fast=TRUE, tuneControl=.noneControl, featureSelector=NULL, parcels=NULL) {
  ret <- if (nrow(tuneGrid) == 1 && fast) {
    RawModel(model, Xtrain, Ytrain, Xtest, Ytest, tuneGrid, featureSelector, parcels)        
  } else {     
    CaretModel(model, Xtrain, Ytrain, Xtest, Ytest, tuneGrid, tuneControl, featureSelector, parcels)
  }  
}

#' @export
#' @import foreach
crossval_external <- function(foldIterator, Xtest, Ytest, model, tuneGrid, fast=TRUE, ncores=1, returnPredictor=FALSE, featureSelector=NULL, parcels=NULL, ensemblePredictor=FALSE) {
  ## TODO bootstrap replications doesn't work because model is trained on full set
  ## This should be handled upstream so that "foldIterator" retruns a bootstrapped sampled version of X.
  results <- if (nrow(tuneGrid) == 1) {
    fit <- trainModel(model, foldIterator$X, foldIterator$Y, Xtest, Ytest, tuneGrid, fast, .noneControl, featureSelector, parcels)
    perf <- evaluateModel(fit)
    list(perf=perf, predictor=asPredictor(fit))
  } else {
    index <- invertFolds(foldIterator$getTestSets(), 1:length(foldIterator$blockVar)) 
    ctrl <- caret::trainControl("cv", verboseIter=TRUE, classProbs=TRUE, index=index, returnData=FALSE, returnResamp="none")
    fit <- trainModel(model, foldIterator$X, foldIterator$Y, Xtest, Ytest, tuneGrid, fast=FALSE, ctrl, featureSelector, parcels)
    perf <- evaluateModel(fit)
    list(perf=perf, predictor=asPredictor(fit))              
  } 
  
  if (returnPredictor) {
    list(class=results$perf$class, probs=results$perf$probs, predictor=results$predictor, featureMask=results$predictor$featureMask, parcels=parcels)
  } else {
    list(class=results$perf$class, probs=results$perf$probs, predictor=NULL, featureMask=results$predictor$featureMask, parcels=parcels)
  }
}

#' @export
#' @import foreach
crossval_internal <- function(foldIterator, model, tuneGrid, fast=TRUE, ncores=1, returnPredictor=FALSE, featureSelector=NULL, parcels=NULL, ensemblePredictor=FALSE) {
 
  resultList <- foreach::foreach(fold = foldIterator, .verbose=FALSE, .packages=c(model$library)) %do% {   
    if (nrow(tuneGrid) == 1 && fast) {
      fit <- trainModel(model, fold$Xtrain, fold$Ytrain, fold$Xtest, fold$Ytest, tuneGrid, fast, .noneControl, featureSelector, parcels)
      list(result=evaluateModel(fit), fit = if (returnPredictor) asPredictor(fit) else NULL, featureMask=fit$featureMask, parcels=parcels)        
    } else {      
      if (nrow(tuneGrid) == 1) {
        fit <- trainModel(model, fold$Xtrain, fold$Ytrain, fold$Xtest, fold$Ytest, tuneGrid, fast=FALSE, tuneControl=.noneControl, featureSelector)
        list(result=evaluateModel(fit), fit = if (returnPredictor) asPredictor(fit) else NULL, featureMask=fit$featureMask, parcels=parcels)    
      } else {
        index <- invertFolds(foldIterator$getTestSets()[-foldIterator$index()], 1:nrow(fold$Xtrain)) 
        ctrl <- caret::trainControl("cv", verboseIter=TRUE, classProbs=TRUE, index=index, returnData=FALSE, returnResamp="none")
        fit <- trainModel(model, fold$Xtrain, fold$Ytrain, fold$Xtest, fold$Ytest, tuneGrid, fast=FALSE, tuneControl=ctrl,featureSelector, parcels)
        list(result=evaluateModel(fit), fit = if (returnPredictor) asPredictor(fit) else NULL, featureMask=fit$featureMask, parcels=parcels)    
      } 
    }
  }
  
  results <- lapply(resultList, "[[", "result")
  
  ## reorder predictions to match order of input features/labels
  ord <- foldIterator$getTestOrder()
  probMat <- do.call(rbind, lapply(results, "[[", "probs"))[ord,]
  predClass <- unlist(lapply(results, "[[", "class"))[ord]
  
  if (returnPredictor) {
    ctrl <- if (nrow(tuneGrid) == 1) {
      .noneControl
    } else {
      caret::trainControl("cv", verboseIter=TRUE, classProbs=TRUE, index=foldIterator$getTrainSets(), returnData=FALSE, returnResamp="none", allowParallel=FALSE)
    }    
    
    #if (foldIterator$balance) {
      ## resample
    #}
    
    #wfit <- WeightedPredictor(lapply(resultList, "[[", "fit"))
    
    fit <- if (ensemblePredictor) {
      WeightedPredictor(lapply(resultList, "[[", "fit"))
    } else {
      asPredictor(trainModel(model, foldIterator$X, foldIterator$Y, NULL, NULL, tuneGrid, fast=fast, tuneControl=ctrl, featureSelector, parcels))
    }
    
    list(class=predClass, probs=probMat, predictor=fit, featureMask=fit$featureMask, parcels=parcels)
    #list(class=predClass, probs=probMat, predictor=wfit)
  } else {
    fmat <- do.call(cbind, lapply(resultList, "[[", "featureMask"))
    list(class=predClass, probs=probMat, predictor=NULL, featureMask=rowMeans(fmat), parcels=parcels)
  }
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


