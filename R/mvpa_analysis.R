


colACC <- function(X, Y) {
  apply(X, 2, function(p) {
    sum(p == Y)/length(Y)
  })
}




#' @export  
matrixToVolumeList <- function(vox, mat, mask, default=NA) {
  lapply(1:ncol(mat), function(i) {
    vol <- array(default, dim(mask))   
    vol[vox] <- mat[,i]
    BrainVolume(vol, space(mask))
  })
} 

#' @export 
computePerformance <- function(result, vox) {
  perf <- t(performance(result))
  out <- cbind(vox, perf[rep(1, nrow(vox)),])   
}

###
#crossValidateROI <- 

#' @export 
mvpa_crossval <- function(dataset, vox, returnPredictor=FALSE) {
  X <- series(dataset$trainVec, vox)
  valid.idx <- nonzeroVarianceColumns(X)
  
  X <- X[,valid.idx]
  vox <- vox[valid.idx,]
  
  foldIterator <- MatrixFoldIterator(X, dataset$Y,dataset$blockVar)
  
  if (ncol(X) == 0) {
    stop("no valid columns")
  } else {
    tuneGrid <- if (!is.null(dataset$tuneGrid)) dataset$tuneGrid else dataset$model$grid(X, dataset$Y, 1)
    
    result <- if (is.null(dataset$testVec)) {
      cvres <- crossval_internal(foldIterator, dataset$model, tuneGrid, fast=TRUE, ncores=1, returnPredictor=returnPredictor)
      classificationResult(dataset$Y, as.factor(cvres$class), cvres$probs,cvres$predictor)
    } else {
      Xtest <- series(dataset$testVec, vox) 
      cvres <- crossval_external(foldIterator, Xtest, dataset$testY, dataset$model, tuneGrid, fast=TRUE, ncores=1, returnPredictor=returnPredictor)
      classificationResult(dataset$testY, as.factor(cvres$class), cvres$probs, cvres$predictor)
    }
     
    result
  }
}

.convertResultsToVolumeList <- function(res, mask) {
  invalid <- sapply(res, function(x) inherits(x, "simpleError") || is.null(x))
  validRes <- do.call(rbind, res[!invalid])
  
  if (length(validRes) == 0 || nrow(validRes) == 0) {
    print(res)
    stop("no valid results, all tasks failed")
  }
  
  vols <- matrixToVolumeList(validRes[,1:3], validRes[,4:ncol(validRes)], mask)
  names(vols) <- colnames(validRes)[4:ncol(validRes)]
  vols
}


.doStandard <- function(dataset, radius, ncores) {
  searchIter <- itertools::ihasNext(Searchlight(dataset$mask, radius)) 
  
  res <- foreach::foreach(vox = searchIter, .verbose=FALSE) %dopar% {   
    if (nrow(vox) > 1) {
      .computePerformance(mvpa_crossval(dataset, vox), vox)
    }
  }
  
  .convertResultsToVolumeList(res, dataset$mask)
}
  

.doRandomized <- function(dataset, radius, returnPredictor=FALSE) {
  searchIter <- itertools::ihasNext(RandomSearchlight(dataset$mask, radius))
  
  res <- foreach::foreach(vox = searchIter, .verbose=FALSE, .errorhandling="pass", .packages=c("rMVPA", dataset$model$library)) %do% {   
    if (nrow(vox) > 1) {  
       print(nrow(vox))
      computePerformance(mvpa_crossval(dataset,vox, returnPredictor), vox)
    }
  }
  
  .convertResultsToVolumeList(res, dataset$mask)
  
}



#runSearchLight <- function(searchIter, blockNum, dataset, model, tuneGrid) {
#    
#}


.searchEnsembleIteration <- function(dataset, trainInd, testInd, radius, model, tuneGrid) {
  searchIter <- itertools::ihasNext(RandomSearchlight(dataset$mask, radius))
  
  Ytrain <- dataset$Y[-testInd]
  blockVar <- dataset$blockVar[-testInd]
  
  res <- lapply(searchIter, function(vox) { 
    if (nrow(vox) > 2) {   
      X <- series(dataset$trainVec, vox)
      Xtrain <- X[-testInd,]    
      foldIterator <- MatrixFoldIterator(Xtrain, Ytrain, blockVar)
      cvres <- try(crossval_internal(foldIterator, model, tuneGrid, fast=TRUE, ncores=1, returnPredictor=TRUE))
      
      if (!inherits(cvres, "try-error")) {    
        cvres <- classificationResult(Ytrain, cvres$class, cvres$probs, cvres$predictor)
        attr(cvres, "vox") <- vox
        attr(cvres, "nvox") <- nrow(vox)
        param_names <- names(tuneGrid)
        attr(cvres, "radius") <- radius
        attr(cvres, "param") <- paste(sapply(1:ncol(tuneGrid), function(i) paste(param_names[i], tuneGrid[,i], sep=":")), collapse="#")
        cvres
      }
      
      cvres
    }
  })
  
  res <- res[!sapply(res, function(x) is.null(x) || inherits(x, "try-error"))]
  res

}



learners = list(
  #pls=data.frame(ncomp=1:3),
  #sda=data.frame(lambda=c(.1, .5, .9), diagonal=c(FALSE,FALSE, FALSE)),
  #corsim=expand.grid(method=c("pearson", "kendall"), robust=c(TRUE, FALSE)),
  #corsim=expand.grid(method=c("pearson", "kendall"), robust=c(TRUE, FALSE)),
  #spls=expand.grid(K=2:4, eta=c(.1,.5,.9), kappa=.5),
  #ada=expand.grid(iter=200, maxdepth=1:2,nu=.1))
  #mlpWeightDecay=expand.grid(size = 4, decay=.001),
  xgboost=expand.grid(nrounds=1:3, max.depth=c(1,2,3), eta=c(.001, .1,.8))
)


#sample_sphere <- function(vox, mask, radius) {
#  ind <- sample(1:nrow(vox), 1) 
#}


.setupModels <- function(learnerSet) {
  unlist(lapply(names(learnerSet), function(mname) {
    model <- loadModel(mname)
    params <- learnerSet[[mname]]
    lapply(1:nrow(params), function(i) {
      list(name=mname, model=model, tuneGrid=params[i,,drop=FALSE])
    }) 
  }), recursive=FALSE)
}


.summarizeResults <- function(resultList) {
  data.frame(AUC=do.call(rbind, lapply(resultList, function(x) performance(x)))[,3],
             nvox=sapply(resultList, function(x) attr(x, "nvox")),
             model=sapply(resultList, function(x) attr(x, "modelname")),
             params=sapply(resultList, function(x) attr(x, "param")),
             radius=sapply(resultList, function(x) attr(x, "radius"))) 
}

metaCombine <- function(dataset, resultList, trainIndex, blockNum, pruneFrac=1, metaLearner="spls", tuneGrid=expand.grid(K=c(1,2,3,4,5), eta=c(.2, .7, .9), kappa=.5)) {
  
  Ytrain <- dataset$Y[fold$trainIndex]
  #Ytest <- dataset$Y[fold$testIndex] 
  
  resultFrame <- .summarizeResults(resultList)
  aucorder <- order(resultFrame$AUC, decreasing=TRUE)
  nkeep <- pruneFrac * length(aucorder)
  keep <- sort(aucorder[1:nkeep])

  predmat <- do.call(cbind, lapply(resultList[keep], "[[", "probs"))
  
  foldIterator <- MatrixFoldIterator(predmat, Ytrain, dataset$blockVar[dataset$blockVar != blockNum])
   
  tuneControl <- caret::trainControl("cv", verboseIter=TRUE, classProbs=TRUE, index=foldIterator$getTrainSets(), indexOut=foldIterator$getTestSets()) 
  #modelFit <- trainModel(loadModel("spls"), predmat, Ytrain,  NULL, NULL, tuneGrid=expand.grid(K=c(1,2,3,4,5), eta=c(.2, .7, .9), kappa=.5), fast=FALSE, tuneControl)
  modelFit <- trainModel(loadModel("gbm"), predmat, Ytrain,  NULL, NULL, tuneGrid=expand.grid(interaction.depth=c(1), n.trees=c(100), shrinkage=c(.01)), fast=FALSE, tuneControl)
  modelFit <- trainModel(loadModel("rf"), predmat, Ytrain,  NULL, NULL, tuneGrid=expand.grid(ncomp=1:5), fast=FALSE, tuneControl)
  
  ensemblePredictor <- ListPredictor(lapply(resultList[keep], function(res) MVPAVoxelPredictor(res$predictor, attr(res, "vox"))))
  
  voxlist <- lapply(resultList[keep], function(el) attr(el, "vox"))
  testlist <- lapply(voxlist, function(vox) series(dataset$trainVec,vox)[fold$testIndex,])
  
  allpreds <- do.call(cbind, lapply(seq_along(resultList[keep]), function(j) {
    evaluateModel(resultList[keep][[j]]$predictor, testlist[[j]])$probs  
  }))
  
  pfinal <- evaluateModel(modelFit, as.matrix(allpreds))
  #list(predictor=CaretPredictor(modelFit), 
  
}

greedyCombine <- function(dataset, resultList, trainIndex, testIndex, calibrateProbs=FALSE, pruneFrac=1) {
  resultFrame <- summarizeResults(resultList)    
  perforder <- order(resultFrame$AUC, decreasing=TRUE)
  
  keep <- if (pruneFrac < 1) {
    sort(perforder[1:(pruneFrac * length(perforder))])
  } else {
    1:length(resultList)
  }
  
  Ytrain <- dataset$Y[trainIndex]
  Ytest <- dataset$Y[testIndex] 
  
  prunedList <- resultList[keep]
  predlist <- lapply(prunedList, "[[", "probs")
  
  baggedWeights <- if (length(levels(Ytrain)) == 2) {
    Pred <- do.call(cbind, lapply(predlist, function(x) x[,2]))
    BaggedTwoClassOptAUC(Pred, ifelse(Ytrain == levels(Ytrain)[2], 1,0))
  } else {
    baggedWeights <- BaggedMultiClassOptAUC(predlist, Ytrain)
  }
  
  baggedWeights <- baggedWeights/sum(baggedWeights)
  posWeights <- which(baggedWeights > 0)
  baggedWeights <- baggedWeights[posWeights]
   
  positiveList <- prunedList[posWeights]
  
  predictorList <- if (calibrateProbs) {
    lapply(positiveList, function(pred) {  
      CalibratedPredictor(pred$predictor, pred$probs, Ytrain)
    })
  } else {
    lapply(positiveList, "[[", "predictor")
  }
  
  voxlist <- lapply(positiveList, function(el) attr(el, "vox"))
  testlist <- lapply(voxlist, function(vox) series(dataset$trainVec,vox)[fold$testIndex,])
  
  p <- lapply(seq_along(predictorList), function(j) {
    pred <- evaluateModel(predictorList[[j]], testlist[[j]])$probs
    if (any(is.na(pred))) {
      NULL
    } else {
      pred * baggedWeights[j]
      #pred
    }
  })
  
  p <- p[!sapply(p, function(x) is.null(x))]
  pfinal <- Reduce("+", p)
  
  
}

# mvpa_searchlight_ensemble
# @param dataset a \code{MVPADataset} instance.
# @param regionMask a \code{BrainVolume} where each region is identified by a unique integer. Every non-zero set of positive integers will be used to define a set of voxels for clasisifcation analysis.
# @param ncores the number of cores for parallel processign (default is 1)
# @return a named list of \code{BrainVolume} objects, where each name indicates the performance metric and label (e.g. accuracy, AUC)
# @import itertools 
# @import foreach
# @import doParallel
# @import parallel
#' @export
mvpa_searchlight_ensemble <- function(dataset, radius=c(14,10, 6), ncores=1, learnerSet = list(pls=data.frame(ncomp=1:4)), pruneFrac=.2, combiner=c("optAUC")) {
  
  if (length(dataset$blockVar) != length(dataset$Y)) {
    stop(paste("length of 'labels' must equal length of 'cross validation blocks'", length(Y), "!=", length(blockVar)))
  }
  
  
  ## setup models and tuning parameters from specification
  models <- .setupModels(learnerSet)
  
   
  
  
  computeVoxelwiseAUC <- function(mask, AUC, radius, voxlist) {
    auc <- array(0, dim(mask))
    count <- array(0, dim(mask))
    for (i in seq_along(AUC)) {
      vox <- voxlist[[i]]
      auc[vox] <- auc[vox] + AUC[i]/radius[i]
      count[vox] <- count[vox] + 1
    }
    
    BrainVolume(auc/count, space(mask))
  }
  
  
  
  blockIterator <- FoldIterator(dataset$blockVar)
    
  allres <- lapply(blockIterator, function(fold) { 
    resultList <- unlist(lapply(radius, function(rad) {
      unlist(lapply(models, function(model) {
        searchResults <- .searchEnsembleIteration(dataset, fold$trainIndex, fold$testIndex, rad,  model$model, model$tuneGrid)     
        searchResults <- lapply(searchResults, function(sr) { 
          print(model$name)
          attr(sr, "modelname") <- model$name 
          sr
        })
      }), recursive=FALSE)
    }), recursive=FALSE)
    
    resultFrame <- summarizeResults(resultList)    
    perforder <- order(resultFrame$AUC, decreasing=TRUE)
    #voxlist <- lapply(resultList, function(el) attr(el, "vox"))
    #AUCvol <- computeVoxelwiseAUC(dataset$mask, resultFrame$AUC, resultFrame$radius, voxlist)
    
    #bestVox <- indexToGrid(dataset$mask, order(AUCvol, decreasing=TRUE)[1:25])
     
    keep <- if (length(perforder) > 50) {
      sort(perforder[1:(pruneFrac * length(perforder))])
    } else {
      which(resultFrame$AUC > 0)
    }
    
    Ytrain <- dataset$Y[fold$trainIndex]
    Ytest <- dataset$Y[fold$testIndex] 
    prunedList <- resultList[keep]
    predlist <- lapply(prunedList, "[[", "probs")
    
    baggedWeights <- if (length(levels(Ytrain)) == 2) {
      Pred <- do.call(cbind, lapply(predlist, function(x) x[,2]))
      BaggedTwoClassOptAUC(Pred, ifelse(Ytrain == levels(Ytrain)[2], 1,0))
    } else {
      #greedMultiClassOptAUC(predlist, Ytrain)
      BaggedMultiClassOptAUC(predlist, Ytrain)
    }
    baggedWeights <- baggedWeights/sum(baggedWeights)
    
     
    posWeights <- which(baggedWeights > 0)
    baggedWeights <- baggedWeights[posWeights]
    
    positiveList <- prunedList[posWeights]
    
    predictorList <- if (calibrateProbs) {
      lapply(positiveList, function(pred) {  
        CalibratedPredictor(pred$predictor, pred$probs, Ytrain)
      })
    } else {
      lapply(positiveList, "[[", "predictor")
    }
    
    voxlist <- lapply(positiveList, function(el) attr(el, "vox"))
    
    ## extract test sets
    testlist <- lapply(voxlist, function(vox) series(dataset$trainVec,vox)[fold$testIndex,])
    
    print(baggedWeights)
    ## we should now rerun all models on the feature selected subset (top 10%?) based on baggedWeights
        
     p <- lapply(seq_along(predictorList), function(j) {
       pred <- evaluateModel(predictorList[[j]], testlist[[j]])$probs
       if (any(is.na(pred))) {
         NULL
       } else {
         pred * baggedWeights[j]
         #pred
       }
     })
    
    p <- p[!sapply(p, function(x) is.null(x))]
    pfinal <- Reduce("+", p)
    
    #weightVol <- BrainVolume(array(0, dim(dataset$mask)), space(mask))
    #for (i in 1:length(voxlist)) {
    #  vox <- voxlist[[i]]
    #  weightVol[vox] <- weightVol[vox] + (baggedWeights[i] * 1/log(nrow(vox)))
    #}
    weightVol <- (weightVol - min(weightVol))/(max(weightVol)-min(weightVol))
    pfinal
  
  })
}

#' mvpa_regional
#' @param dataset a \code{MVPADataset} instance.
#' @param regionMask a \code{BrainVolume} where each region is identified by a unique integer. Every non-zero set of positive integers will be used to define a set of voxels for clasisifcation analysis.
#' @param ncores the number of cores for parallel processign (default is 1)
#' @return a named list of \code{BrainVolume} objects, where each name indicates the performance metric and label (e.g. accuracy, AUC)
#' @import itertools 
#' @import foreach
#' @import doParallel
#' @import parallel
#' @export
mvpa_regional <- function(dataset, regionMask, ncores=1, savePredictors=FALSE) {
  if (length(dataset$blockVar) != length(dataset$Y)) {
    stop(paste("length of 'labels' must equal length of 'cross validation blocks'", length(Y), "!=", length(blockVar)))
  }
  
  cl <- makeCluster(ncores)
  registerDoParallel(cl)
  
  regionSet <- sort(as.integer(unique(regionMask[regionMask > 0])))
  
  #if (ncores > 1 && length(regionSet) == 1) {
  #  mc.cores <- ncores
  #} else {
  #  mc.cores <- 1
  #}
  
  #allowParallel <- if (length(regionSet) == 1) TRUE else FALSE
  
  res <- foreach::foreach(roinum = regionSet, .verbose=TRUE, .errorhandling="pass", .packages=c("rMVPA", "MASS", "neuroim", dataset$model$library)) %dopar% {   
    idx <- which(regionMask == roinum)
    if (length(idx) > 1) {
      vox <- indexToGrid(regionMask, idx)
      result <- mvpa_crossval(dataset, vox, returnPredictor=savePredictors)
      #fit <- fitMVPAModel(dataset, vox, fast=TRUE, finalFit=FALSE, ncores=mc.cores)     
      perf <- c(ROINUM=roinum, t(performance(result))[1,])     
      predmat <- data.frame(ROI=rep(roinum, length(result$observed)), observed=result$observed, pred=result$predicted, correct=result$observed == result$predicted, prob=result$prob)
      
      perf <- structure(perf,
                predictor=result$predictor,
                predmat=predmat)
      
      perf    
    }
  }
  
  invalid <- sapply(res, function(x) inherits(x, "simpleError") || is.null(x))
  validRes <- res[!invalid]
  
  if (length(validRes) == 0) {
    print(res)
    flog.error("Classification failed for all of %s ROIs", length(regionSet))
    stop("error in mvpa_regional: aborting")
  } 
  
  perfMat <- do.call(rbind, validRes)
  predMat <- do.call(rbind, lapply(validRes, function(x) attr(x, "predmat")))
  
  outVols <- lapply(2:ncol(perfMat), function(cnum) {
    fill(regionMask, cbind(perfMat[, 1], perfMat[,cnum]))    
  })
  
  predictorList <- lapply(validRes, function(res) {
    attr(res, "predictor")
  })
  
  names(predictorList) <- regionSet[!invalid]
  names(outVols) <- colnames(perfMat)[2:ncol(perfMat)]
  list(outVols = outVols, performance=perfMat, predictorList=predictorList, predMat=predMat)

}
  
  
  
#' mvpa_searchlight
#' @param dataset a \code{MVPADataset} instance.
## @param trainVec a \code{BrainVector} instance, a 4-dimensional image where the first three dimensons are (x,y,z) and the 4th dimension is the dependent class/variable
## @param Y the dependent variable for training data. If it is a factor, then classification analysis is performed. If it is a continuous variable then regression is performed.
##        the length of \code{Y} must be the same as the length of the 4th dimension of \code{train_vec}
## @param mask a \code{BrainVolume} instance indicating the inclusion mask for voxels entering the searchlight analysis. 
## @param blockVar an \code{integer} vector indicating the blocks to be used for cross-validation. This is usually a variable indicating the scanning "run". 
##        Must be same length as \code{Y}
#' @param radius the searchlight radus in mm
## @param modelName the name of the classifcation model to be used
#' @param ncores the number of cores for parallel processign (default is 1)
#' @param method the type of searchlight (randomized, or standard)
#' @param niter the number of searchlight iterations for 'randomized' method
## @param tuneGrid parameter search grid for optimization of classifier tuning parameters
## @param testVec a \code{BrainVector} with the same spatial dimension as shape as \code{trainVec}. If supplied, this data will be held out as a test set.
## @param testY the dependent variable for test data. If supplied, this variable to evaluate classifier model trained on \code{trainVec}. 
##        \code{testY} must be the same as the length of the 4th dimension of \code{test_vec}.
#' @return a named list of \code{BrainVolume} objects, where each name indicates the performance metric and label (e.g. accuracy, AUC)
#' @import itertools 
#' @import foreach
#' @import doParallel
#' @import parallel
#' @import futile.logger
#' @export
mvpa_searchlight <- function(dataset, radius=8, method=c("randomized", "standard"), niter=4,ncores=2) {
  if (radius < 1 || radius > 100) {
    stop(paste("radius", radius, "outside allowable range (1-100)"))
  }
  
  if (length(dataset$blockVar) != length(dataset$Y)) {
    stop(paste("length of 'labels' must equal length of 'cross validation blocks'", length(dataset$Y), "!=", length(dataset$blockVar)))
  }
  
 
  cl <- makeCluster(ncores)
  registerDoParallel(cl)
  
  method <- match.arg(method)
  
  flog.info("classification model is: %s", dataset$model$label)
  
  res <- if (method == "standard") {
    .doStandard(dataset, radius, ncores)    
  } else {
    res <- parallel::mclapply(1:niter, function(i) {
      flog.info("Running randomized searchlight iteration %s", i)   
      do.call(cbind, .doRandomized(dataset, radius) )
    }, mc.cores=ncores)
   
    Xall <- lapply(1:ncol(res[[1]]), function(i) {
      X <- do.call(cbind, lapply(res, function(M) M[,i]))
      xmean <- rowMeans(X, na.rm=TRUE)
      xmean[is.na(xmean)] <- 0
      BrainVolume(xmean, space(dataset$mask))
    })
    
    names(Xall) <- colnames(res[[1]])
    Xall
    
  }
  
}