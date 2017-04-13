.noneControl <- caret::trainControl("none", verboseIter=TRUE, classProbs=TRUE, returnData=FALSE, returnResamp="none", allowParallel=FALSE)

.cvControl <- caret::trainControl("cv", verboseIter=TRUE, classProbs=TRUE, returnData=FALSE, returnResamp="none",allowParallel=FALSE)  

.adaptiveControl <- caret::trainControl("adaptive_cv", verboseIter=TRUE, classProbs=TRUE, returnData=FALSE, returnResamp="none",allowParallel=FALSE)  

#' Create an \code{ClassificationResultSet} instance
#' @param blockVar the cross-validation blocking variable
#' @param resultList a list of type \code{ClassificationResult} 
#' @export
ClassificationResultSet <- function(blockVar, resultList) {
  ret <- list(
    blockVar=blockVar,
    resultList=resultList)
  
  class(ret) <- c("ClassificationResultSet", "list")
  ret
}



# wraps a caret model fit
CaretModel <- function(model, fit, Xtrain, Ytrain, Xtest, Ytest, tuneGrid, tuneControl, featureSelector=NULL, featureMask=NULL, parcels=NULL) {
  ret <- list( model=model, Xtrain=Xtrain,Ytrain=Ytrain,Xtest=Xtest,Ytest=Ytest,tuneGrid=tuneGrid,
               tuneControl=tuneControl,modelFit=fit,featureSelector=featureSelector,featureMask=featureMask,parcels=parcels)
  if (is.factor(Ytrain)) {
    class(ret) <- c("CaretModel", "list")
  } else {
    class(ret) <- c("CaretRegressionModel", "list")
  }
  ret
}


# wraps a model fit that was not bypassed caret 'train' faunction
RawModel <- function(model, fit, Xtrain, Ytrain, Xtest, Ytest, tuneGrid, featureSelector=NULL, featureMask=NULL, parcels=NULL) {
  
  ret <- list(model=model, Xtrain=Xtrain, Ytrain=Ytrain, Xtest=Xtest, Ytest=Ytest, tuneGrid=tuneGrid,modelFit=fit,
              featureSelector=featureSelector,featureMask=featureMask,parcels=parcels)
  
  if (is.factor(Ytrain)) {
    ret$modelFit$problemType <- "Classification"
    class(ret) <- c("RawModel", "list")
  } else {
    ret$modelFit$problemType <- "Regression"
    class(ret) <- c("RawRegressionModel", "list")
  }
  
  ret
}


trainModel <- function(model, ROI, Ytrain, testROI, Ytest, tuneGrid, 
                       tuneControl=.noneControl, featureSelector=NULL, parcels=NULL, ...) {
  if (!is.null(parcels)) {
    stop("parcels not supported yet")
    #message("computing parcel means", dim(ROI@data))
    #Xtrain <- groupMeans(ROI@data, 2, parcels)
  }
  
  featureMask <- if (!is.null(featureSelector)) {
    select_features(featureSelector, ROI, Ytrain)
  } else {
    rep(TRUE, length(ROI))
  }
  
  if (sum(featureMask) != length(ROI)) {
    ## subset ROI 
    ROI <- ROI[featureMask]
  } 
  
  
  
  #print(paste("Number of features: ", sum(featureMask)))
  #print(paste("length of ROI: ", length(ROI)))
  
  
  if (nrow(tuneGrid) == 1) {
    ## if model needs "structural" ROI information, then this is the "dead end" as ROI is stripped to a matrix.
    fit <- model$fit(values(ROI), Ytrain, NULL, tuneGrid, lev=levels(Ytrain), classProbs=TRUE)
    RawModel(model, fit, values(ROI), Ytrain, values(testROI), Ytest, tuneGrid, featureSelector, featureMask, parcels )
  } else {     
    fit <- caret::train(values(ROI), Ytrain, method=model, trControl=tuneControl, tuneGrid=tuneGrid, ...)
    CaretModel(model, fit, values(ROI), Ytrain, values(testROI), Ytest, tuneGrid, tuneControl, featureSelector, featureMask, parcels )
  }  
}


#' @export
RawPredictor <- function(fit, model, featureMask, parcels=NULL) {
  ret <- list(fit=fit, model=model, featureMask=featureMask,parcels=parcels)
  class(ret) <- c("RawPredictor", "list")
  ret
}

#' @export
RawRegressionPredictor <- function(fit, model, featureMask, parcels=NULL) {
  ret <- list(fit=fit, model=model, featureMask=featureMask,parcels=parcels)
  class(ret) <- c("RawRegressionPredictor", "list")
  ret
}

#' @export
CaretPredictor <- function(fit, featureMask, parcels=NULL) {
  ret <- list(fit=fit, featureMask=featureMask, parcels=parcels)
  class(ret) <- c("CaretPredictor", "list")
  ret
}

#' @export
CaretRegressionPredictor <- function(fit, featureMask, parcels=NULL) {
  ret <- list(fit=fit, featureMask=featureMask, parcels=parcels)
  class(ret) <- c("CaretRegressionPredictor", "list")
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
evaluateModel <- function(x, newdata,...) {
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
asPredictor.RawRegressionModel <- function(x) {
  RawRegressionPredictor(x$modelFit, x$model, x$featureMask, x$parcels)
}


#' @export
asPredictor.CaretModel <- function(x) {
  CaretPredictor(x$modelFit, x$featureMask, x$parcels)
}

#' @export
asPredictor.CaretRegressionModel <- function(x) {
  CaretRegressionPredictor(x$modelFit, x$featureMask, x$parcels)
}

#' @export
evaluateModel.RawRegressionModel <- function(x, newdata=NULL) {
  if (is.null(newdata)) {
    newdata=x$Xtest
  }
  
  preds <- if (is.null(x$parcels)) {
    x$model$predict(x$modelFit, newdata[,x$featureMask])  
  } else {
    reducedData <- groupMeans(newdata, 2, x$parcels)
    x$model$predict(x$modelFit, reducedData[,x$featureMask])  
  }
  
  list(preds=preds)
}


#' @export
evaluateModel.RawModel <- function(x, newdata=NULL) {
  if (is.null(newdata)) {
    newdata=x$Xtest
  }
  
  probs <- if (is.null(x$parcels)) {
    x$model$prob(x$modelFit, newdata[, x$featureMask])  
  } else {
    reducedData <- groupMeans(newdata, 2, x$parcels)
    x$model$prob(x$modelFit, reducedData[, x$featureMask])  
  }
  
  cpred <- apply(probs,1, function(x) {
    x[is.na(x)] <- 0
    which.max(x)    
  })
  
  cpred <- colnames(probs)[cpred]
  list(class=cpred, probs=probs)
}

#' @export
evaluateModel.CalibratedPredictor <- function(x, newdata,...) {
  if (is.null(newdata)) {
    stop("newdata cannot be null")
  }
  
  Preds <- evaluateModel(x$predictor, newdata,...)$probs
  Pcal <- do.call(cbind, lapply(1:ncol(Preds), function(i) {
    predict(x$calfits[[i]], data.frame(x0=Preds[,i]), type="response")
  }))
  
  sweep(Pcal, 1, rowSums(Pcal), "/")
  
}

#' @export
evaluateModel.MVPAVoxelPredictor <- function(x, newdata, subIndices=NULL) {
  if (is.null(newdata)) {
    stop("newdata cannot be null")
  }
  
  X <- if (!is.null(subIndices)) {
    series(newdata,x$voxelGrid)[subIndices,]
  } else {
    series(newdata,x$voxelGrid)
  }
  
  if (is.vector(X)) {
    X <- matrix(X, nrow=1)
  }
  
  evaluateModel(x$predictor, X)
}


#' @export
evaluateModel.RawPredictor <- function(x, newdata=NULL,...) {
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
evaluateModel.RawRegressionPredictor <- function(x, newdata=NULL) {
  if (is.null(newdata)) {
    stop("newdata cannot be null")
  }
  
  preds <- if (is.null(x$parcels)) {
    x$model$predict(x$modelFit, newdata[,x$featureMask])  
  } else {
    reducedData <- groupMeans(newdata, 2, x$parcels)
    x$model$predict(x$modelFit, reducedData[,x$featureMask])  
  }
  
  list(preds=preds)
}



#' @export
evaluateModel.CaretRegressionModel <- function(x, newdata=NULL) {
  if (is.null(newdata)) {
    newdata=x$Xtest
  }
  
  if (!is.null(x$parcels)) {
    newdata <- groupMeans(newdata, 2, x$parcels)
  }
  
  list(preds=predict(x$modelFit, newdata[,x$featureMask]))
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
evaluateModel.ListPredictor <- function(x, newdata=NULL,...) {
  if (is.null(newdata)) {
    stop("newdata cannot be null")
  }
  
  res <- lapply(x, function(fit) {
    evaluateModel(fit, newdata,...)
  })
  
}

#' @export
predict.WeightedPredictor <- function(x, newdata=NULL, ...) {
  if (is.null(newdata)) {
    stop("newdata cannot be null")
  }
  
  wts <- attr(x, "weights")
  
  preds <- lapply(1:length(x), function(i) {
    evaluateModel(x[[i]], newdata, ...)$prob * wts[i]
  })
  
  prob <- preds[!sapply(preds, function(x) is.null(x))]
  pfinal <- Reduce("+", prob)
  
  
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
evaluateModel.CaretPredictor <- function(x, newdata=NULL,...) {
  if (is.null(newdata)) {
    stop("newdata cannot be null")
  }
  
  if (!is.null(x$parcels)) {
    newdata <- groupMeans(newdata, 2, x$parcels)
  }
  
  list(class=predict(x$fit, newdata=newdata[,x$featureMask]), probs=predict(x$fit, newdata=newdata[,x$featureMask], type="prob"))
}

#' @export
evaluateModel.CaretRegressionPredictor <- function(x, newdata=NULL) {
  if (is.null(newdata)) {
    stop("newdata cannot be null")
  }
  
  if (!is.null(x$parcels)) {
    newdata <- groupMeans(newdata, 2, x$parcels)
  }
  
  list(preds=predict(x$fit, newdata=newdata[,x$featureMask]))
}

#' @export
evaluateModelList <- function(modelList, newdata=NULL) {
  results <- lapply(modelList, evaluateModel, newdata)
  probMat <- do.call(rbind, lapply(results, "[[", "probs"))
  predClass <- unlist(lapply(results, "[[", "class"))
  list(class=predClass, probs=probMat)
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


#' Carry out external cross-validation defined by a \code{FoldIterator} instance and a separate test data set.
#' @param crossVal the \code{CrossValidation} object that defines the splits of training and test samples.
#' @param Y the labels or values of dependent variable for the training data.
#' @param ROI a region of interest of type \code{ROIVolume} or \code{ROISurface} that encapsulates the training data.
#' @param testY the labels or values of dependent variable for the test data.
#' @param testROI a region of interest of type \code{ROIVolume} or \code{ROISurface} that encapsulates the test data.
#' @param subIndices an integer vector containg the subset of rows to train model on. if \code{NULL}, train on all rows.
#' @param Xtest the test data \code{matrix}
#' @param Ytest the test data labels (continuous or factor variables permitted, depending on supplied \code{model}.
#' @param model a \code{caret} model object
#' @param tuneGrid model tuning parameters stored in a \code{data.frame}
#' @param featureSelector an optional \code{FeatureSelector} instance
#' @param parcels an optional parcellation object
#' @export
#' @details this function should not typically be called by user code.
#' @import foreach
#' @export
crossval_external <- function(crossVal, Y, ROI, Ytest, testROI, subIndices=NULL, model, tuneGrid, featureSelector=NULL, parcels=NULL) {
  ## TODO bootstrap replications doesn't work because model is trained on full set
  ## This should be handled upstream so that "foldIterator" retruns a bootstrapped sampled version of X.
  
  if (!is.null(subIndices)) {
    ROI <- ROI[,subIndices]
  }
  
  results <- if (nrow(tuneGrid) == 1) {
    ## no parameter tuning required. train model of full data set.
    fit <- trainModel(model, ROI, Y, testROI, Ytest, tuneGrid, .noneControl, featureSelector, parcels)
    perf <- evaluateModel(fit)
    list(perf=perf, predictor=asPredictor(fit))
  } else {
    foldIterator <- foldIter(crossVal, Y, subIndices)
    ## parameter tuning required. Will be performed on training blocks defined in foldIterator.
    ctrl <- if (foldIterator$nfolds <= 2) {
      ## with two or fewer blocks, will perform parameter tuning using default caret cross-validation.
      caret::trainControl("cv", verboseIter=TRUE, classProbs=TRUE, returnData=FALSE, returnResamp="none")
    } else {
      ## training folds stored as list of indices and provided to caret cross-validation engine.
      index <- invertFolds(foldIterator$getTestSets(), 1:length(foldIterator$blockVar)) 
      caret::trainControl("cv", verboseIter=TRUE, classProbs=TRUE, index=index, returnData=FALSE, returnResamp="none")
    }
    
    fit <- trainModel(model, ROI, foldIterator$Y, testROI, Ytest, tuneGrid, ctrl, featureSelector, parcels)
    perf <- evaluateModel(fit)
    list(perf=perf, predictor=asPredictor(fit))              
  } 
  
  result <- list(
    prediction = list(results$perf),
    predictor  =  list(results$predictor),
    testIndices = list(1:length(Ytest)),
    featureMask = list(results$predictor$featureMask),
    subIndices = subIndices
  )
  
  result
  
}

#' Carry out internal cross-validation defined by a \code{CrossValidation} instance.
#' 
#' @param crossVal the \code{CrossValidation} object that defines the splits of training and test samples.
#' @param Y the labels or values of dependent variable.
#' @param ROI a region of interest of type \code{ROIVolume} or \code{ROISurface} that encapsulates the training data.
#' @param subIndices an integer vector containg the subset of rows to train model on.
#' @param model a \code{caret} model object
#' @param tuneGrid model tuning parameters stored in a \code{data.frame}
#' @param featureSelector an optional \code{FeatureSelector} instance
#' @param parcels an optional parcellation object
#' @export
#' @details this function should not typically be called by user code.
#' @import foreach
crossval_internal <- function(crossVal, Y, ROI, subIndices, model, tuneGrid, featureSelector=NULL, parcels=NULL) {
  foldIterator <- foldIter(crossVal, Y, subIndices)
  
  ### loop over folds
  resultList <- foreach::foreach(fold = foldIterator, .verbose=FALSE, .packages=c(model$library)) %do% {   
    if (nrow(tuneGrid) == 1) {
      ## subset ROI with training indices
      fit <- try(trainModel(model, ROI[,fold$trainIndex], fold$Ytrain, ROI[, fold$testIndex], fold$Ytest, 
                            tuneGrid,  .noneControl, featureSelector, parcels))
      list(result=evaluateModel(fit), fit = asPredictor(fit), featureMask=fit$featureMask, parcels=parcels, testIndices=fold$testIndex)        
    } else {    
      ctrl <- if (foldIterator$nfolds == 2) {
        caret::trainControl("cv", verboseIter=TRUE, classProbs=TRUE, returnData=FALSE, returnResamp="none")
      } else {
        index <- invertFolds(foldIterator$getTestSets()[-foldIterator$index()], seq_along(fold$Ytrain)) 
        caret::trainControl("cv", verboseIter=TRUE, classProbs=TRUE, index=index, returnData=FALSE, returnResamp="none")
      }
      
      fit <- trainModel(model, ROI[, fold$trainIndex], fold$Ytrain, ROI[, fold$testIndex], fold$Ytest, tuneGrid, tuneControl=ctrl, featureSelector, parcels)
      list(result=evaluateModel(fit), fit = asPredictor(fit), featureMask=fit$featureMask, parcels=parcels,testIndices=fold$testIndex)    
    } 
  }
  
  result <- list(
    prediction =  lapply(resultList, "[[", "result"),
    predictor  =  lapply(resultList, "[[", "fit"),
    testIndices = lapply(resultList, "[[", "testIndices"),
    ## redundant
    featureMask = lapply(resultList, "[[", "featureMask")
  )
  
  result
  
}



zeroVarianceColumns <- function(M) {
  which(apply(M, 2, sd, na.rm=TRUE) == 0)
}


zeroVarianceColumns2 <- function(M) {
  apply(M, 2, sd, na.rm=TRUE) == 0
}

nonzeroVarianceColumns <- function(M) {
  which(apply(M, 2, sd, na.rm=TRUE) > 0)
}

nonzeroVarianceColumns2 <- function(M) {
  apply(M, 2, sd, na.rm=TRUE) > 0
}

removeZeroVarianceColumns <- function(M) {
  noVariance <- which(apply(M, 2, sd, na.rm=TRUE) == 0)
  if (length(noVariance) > 0) {
    M[, hasVariance, drop=FALSE]
  } else {
    M
  }
}


