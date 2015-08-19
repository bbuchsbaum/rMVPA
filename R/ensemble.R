
### lda ensemble
### train n patches
### best best n patches via greedy search
### train lda on pairs of intra verus intraclass similarity vectors
### combine with global classifier



createModelSet <- function(modelName, ...) {
  dots <- list(...)
  tuneGrid <- expand.grid(dots)
  ret <- cbind(rep(modelName, nrow(tuneGrid)), tuneGrid)
  names(ret) <- c("model", names(tuneGrid))
  attr(ret, "caret_model") <- loadModel(modelName)
  class(ret) <- c("data.frame", "modelSet")
  ret
}

createEnsembleSpec <- function(...) {
  dots <- list(...)
  res <- unlist(lapply(dots, function(x) inherits(x, "modelSet")))
  if (!all(res)) {
    stop("all arguments must inherit from 'modelSet")
  }
  
  class(dots) <- "EnsembleSpec"
  dots
}




.setupModels <- function(learnerSet) {
  unlist(lapply(names(learnerSet), function(mname) {
    model <- loadModel(mname)
    params <- learnerSet[[mname]]
    if (is.data.frame(params)) {
      lapply(1:nrow(params), function(i) {
        list(name=mname, model=model, tuneGrid=params[i,,drop=FALSE])
      }) 
    } else {
      mgrid <- model$grid(matrix(rnorm(100*100), 100, 100),rep(1,100), params)
      lapply(1:nrow(mgrid), function(i) {
        list(name=mname, model=model, tuneGrid=mgrid[i,,drop=FALSE])
      })
    }
  }), recursive=FALSE)
}

.computeVoxelwiseAUC <- function(mask, AUC, radius, voxlist) {
  auc <- array(0, dim(mask))
  count <- array(0, dim(mask))
  for (i in seq_along(AUC)) {
    vox <- voxlist[[i]]
    auc[vox] <- auc[vox] + AUC[i]/radius[i]
    count[vox] <- count[vox] + 1
  }
  
  BrainVolume(auc/count, space(mask))
}

.summarizeResults <- function(resultList) {
  data.frame(AUC=do.call(rbind, lapply(resultList, function(x) performance(x)))[,3],
             nvox=sapply(resultList, function(x) attr(x, "nvox")),
             model=sapply(resultList, function(x) attr(x, "modelname")),
             params=sapply(resultList, function(x) attr(x, "param")))
          
}



metaCombine <- function(dataset, resultList, fold, pruneFrac=1, metaLearner="spls", tuneGrid=expand.grid(K=c(1,2,3,4,5,6), eta=c(.2, .7, .9), kappa=.5)) {
  AUC <- do.call(rbind, lapply(resultList, function(x) performance(x)))[,3]
  Ytrain <- dataset$Y[fold$trainIndex]  
  aucorder <- order(AUC, decreasing=TRUE)
  nkeep <- pruneFrac * length(aucorder)
  keep <- sort(aucorder[1:nkeep])
  
  predmat <- do.call(cbind, lapply(resultList[keep], "[[", "probs"))
  
  foldIterator <- MatrixFoldIterator(X=predmat, Y=Ytrain, blockVar=dataset$blockVar[dataset$blockVar != fold$index])
  
  tuneControl <- caret::trainControl("cv", verboseIter=TRUE, classProbs=TRUE, index=foldIterator$getTrainSets(), indexOut=foldIterator$getTestSets()) 
  modelFit <- trainModel(loadModel(metaLearner), predmat, Ytrain,  NULL, NULL, tuneGrid=tuneGrid, fast=FALSE, tuneControl)
  
  
  ## need a MetaPredictor class
  ensemblePredictor <- ListPredictor(lapply(resultList[keep], function(res) MVPAVoxelPredictor(res$predictor, attr(res, "vox"))))
  testVec <- takeVolume(dataset$trainVec, fold$testIndex,merge=TRUE)
  allpreds <- do.call(cbind, lapply(evaluateModel(ensemblePredictor, testVec), function(p) p$probs))
  pfinal <- evaluateModel(modelFit, as.matrix(allpreds))
  ## need a MetaPredicotr class
  
  
  result <- classificationResult(fold$Ytest, pfinal$class, pfinal$probs)
  
  ## may fail if we purge data from pls models.
  #vimp <- caret::varImp(modelFit$modelFit)$importance
  
  #finalWeights <- unlist(lapply(1:length(resultList), function(i) {
  #  start <- (i - 1) * length(levels(Ytrain)) + 1
  #  end <- (i - 1) * length(levels(Ytrain)) + length(levels(Ytrain))
  #  mean(colMeans(vimp[start:end,,drop=FALSE]))
  #}))
  
  #voxelList <- lapply(resultList, function(el) attr(el, "vox"))
  #weightVol <- .computeWeightVol(voxelList, finalWeights, dataset$mask)
  #AUCVol <- .computeWeightVol(voxelList, resultFrame$AUC, dataset$mask)
  list(result=result)
  
}


.computeWeightVol <- function(voxelList, finalWeights, mask) {
  weightVol <- BrainVolume(array(0, dim(mask)), space(mask))
  countVol <- BrainVolume(array(0, dim(mask)), space(mask))
  
  for (i in 1:length(voxelList)) {
    vox <- voxelList[[i]]   
    weightVol[vox] <- weightVol[vox] + finalWeights[i]
    countVol[vox] <- countVol[vox] + 1
  }
  
  nz <- which(countVol > 0)
  weightVol[nz] <- weightVol[nz]/countVol[nz]
  weightVol
}

greedyCombine <- function(dataset, resultList, trainIndex, testIndex, calibrateProbs=FALSE, pruneFrac=.15) {
  #resultFrame <- .summarizeResults(resultList)    
  stopifnot(pruneFrac > 0 && pruneFrac <= 1)
  
  ## compute AUC for resultSet
  #AUC <- do.call(rbind, lapply(resultList, function(x) performance(x)))[,3]
  AUC <- do.call(rbind, lapply(resultList, function(x) performance(x)))[,2]
  
  ## order from highest to lowest
  perforder <- order(AUC, decreasing=TRUE)
  
  ## results to keep according to 'pruneFrac'
  keep <- sort(perforder[1:(pruneFrac * length(perforder))])
 
  Ytrain <- dataset$Y[trainIndex]
  Ytest <- dataset$Y[testIndex] 

  prunedList <- resultList[keep]
  
  predlist <- lapply(prunedList, "[[", "probs")
  
  baggedWeights <- if (length(levels(Ytrain)) == 2) {
    Pred <- do.call(cbind, lapply(predlist, function(x) x[,2]))
    BaggedTwoClassOptAUC(Pred, ifelse(Ytrain == levels(Ytrain)[2], 1,0))
  } else {
    #BaggedMultiClassOptAUC(predlist, Ytrain)
    BaggedMultiClassOptACC(predlist, Ytrain)
  }
  
  baggedWeights <- baggedWeights/sum(baggedWeights)
  
  posWeights <- which(baggedWeights > 0)
  
  ## weights for final model
  finalWeights <- baggedWeights[posWeights]
  
  ## full set of weights
  allWeights <- numeric(length(resultList))
  allWeights[keep[posWeights]] <- finalWeights
  
  positiveList <- prunedList[posWeights]
  voxelList <- lapply(positiveList, function(el) attr(el, "vox"))
  
  predictorList <- if (calibrateProbs) {
    lapply(positiveList, function(pred) {  
      mvpaPred <- MVPAVoxelPredictor(pred$predictor, attr(pred, "vox"))
      CalibratedPredictor(mvpaPred, pred$probs, Ytrain)
    })
  } else {
    lapply(positiveList, function(pred) {  
      MVPAVoxelPredictor(pred$predictor, attr(pred, "vox"))      
    })
  }
  
  ## weighted predictor based on finalWeights
  baggedPredictor <- WeightedPredictor(predictorList, weights=finalWeights)
  
  ## predict on test set
  testVec <- takeVolume(dataset$trainVec, testIndex, merge=TRUE)
  testPreds <- evaluateModel(baggedPredictor, testVec)
  

  result <- classificationResult(Ytest, testPreds$class, testPreds$prob, baggedPredictor)
  
  #weightVol <- .computeWeightVol(voxelList, finalWeights, dataset$mask)
  #AUCVol <- .computeWeightVol(lapply(resultList, function(el) attr(el, "vox")), resultFrame$AUC, dataset$mask)
  
  list(result=result, allWeights=allWeights)
  
}

innerIteration <- function(dataset, vox, trainInd, testInd, model, tuneGrid, autobalance, bootstrap) {
  Ytrain <- dataset$Y[-testInd]
  blockVar <- dataset$blockVar[-testInd]
  X <- series(dataset$trainVec, vox)
  Xtrain <- X[-testInd,] 
  foldIterator <- MatrixFoldIterator(Xtrain, Ytrain, blockVar, balance=autobalance, bootstrap=bootstrap)
  
  ## do we need to compute final pedictor?
  res <- try(crossval_internal(foldIterator, model, tuneGrid, fast=TRUE, ncores=1, returnPredictor=TRUE, featureSelector=NULL))
  
  structure(classificationResult(Ytrain, res$class, res$probs, res$predictor),
                     vox=vox)
                   
}

.searchEnsembleIteration <- function(searchIter, dataset, trainInd, testInd, model, tuneGrid, autobalance=FALSE, bootstrap=FALSE) {
  #searchIter <- itertools::ihasNext(RandomSearchlight(mask, radius))
  
  Ytrain <- dataset$Y[-testInd]
  blockVar <- dataset$blockVar[-testInd]
  
  res <- lapply(searchIter, function(vox) { 
    if (nrow(vox) > 2) {   
      cvres <- innerIteration(dataset, vox, trainInd, testInd, model, tuneGrid, autobalance, bootstrap)

      param_names=names(tuneGrid)
      if (!inherits(cvres, "try-error")) {    
        cvres <- structure(classificationResult(Ytrain, cvres$class, cvres$probs, cvres$predictor),
          vox=vox,
          nvox=nrow(vox),
          param_names=param_names,
          #radius=radius,
          param=paste(sapply(1:ncol(tuneGrid), function(i) paste(param_names[i], tuneGrid[,i], sep=":")), collapse="#")
        )
      }
      
      cvres
    }
  })
  
  res <- res[!sapply(res, function(x) is.null(x) || inherits(x, "try-error"))]
  res
}



superLearners = .setupModels(list(
  #avNNet=expand.grid(size = c(2,3,4), decay=c(.01, .001, .0001), bag=FALSE),
  pls=data.frame(ncomp=1:5),
  sda=data.frame(lambda=c(.01, .1, .5, .9), diagonal=c(FALSE,FALSE,FALSE, FALSE))
  #spls=expand.grid(K=c(1:5), eta=c(.1,.3,.5,.7), kappa=.5)
))

mergeEnsembles <- function(ensembleSet, testOrder, returnPredictor=FALSE) {
  
  #testOrder <- blockIterator$getTestOrder()
  predicted <- unlist(lapply(ensembleSet, function(x) x$result$predicted))[testOrder]
  observed <- unlist(lapply(ensembleSet, function(x) x$result$observed))[testOrder]
  probs <- do.call(rbind, lapply(ensembleSet, function(x) x$result$probs))[testOrder,]
  
  if (returnPredictor) {
    ## predictor equally weights each ensemble
    predictor <- WeightedPredictor(lapply(ensembleSet, function(x) x$result$predictor))
    classificationResult(observed, predicted, probs, predictor)
  } else {
    classificationResult(observed, predicted, probs)
  }
  
  
}                               

mvpa_roi_ensemble <- function(dataset, regionMask, ncores=1, 
                              pruneFrac=.5, 
                              combiner="greedyAUC", 
                              bootstrapSamples=TRUE,
                              autobalance=FALSE,
                              returnPredictor=FALSE,
                              ensembleSpec=NULL) {
  
  if (length(dataset$blockVar) != length(dataset$Y)) {
    stop(paste("length of 'labels' must equal length of 'cross validation blocks'", length(Y), "!=", length(blockVar)))
  }
  
  cl <- makeCluster(ncores, outfile="",useXDR=FALSE, type="FORK")
  registerDoParallel(cl)
  
  ## Get the set of unique ROIs (all unique integers > 0 in provided mask)
  regionSet <- sort(as.integer(unique(regionMask[regionMask > 0])))
  ## block iterator
  blockIterator <- FoldIterator(dataset$Y, blockVar=dataset$blockVar)
  
  ensembleSet <- foreach::foreach(fold = blockIterator, .verbose=FALSE, .errorhandling="pass", .packages=c("rMVPA", "MASS", "neuroim", "caret", dataset$model$library)) %dopar% {   
    ## run model over all rois for all folds
    resultList <- lapply(regionSet, function(roinum) {
      idx <- which(regionMask == roinum)
      ## choke on length == 1
      
      if (length(idx) < 2) {
        stop(paste("ROI number", roinum, "has less than 2 voxels, aborting."))
      }
      
      vox <- indexToGrid(regionMask, idx)
      
      lapply(1:nrow(dataset$tuneGrid), function(j) {
        innerIteration(dataset, vox, fold$trainIndex, fold$testIndex, dataset$model, dataset$tuneGrid[j,,drop=FALSE], autobalance=autobalance, bootstrap=bootstrapSamples)
      })
      
    })
    
    resultList <- unlist(resultList, recursive=FALSE)
    ## ensemble for current fold
    ens <- greedyCombine(dataset, resultList, fold$trainIndex, fold$testIndex, calibrateProbs=FALSE, pruneFrac=pruneFrac)
  }
  
  finalResult <- mergeEnsembles(ensembleSet, blockIterator$getTestOrder(), returnPredictor) 
 
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
#' 
mvpa_searchlight_ensemble <- function(modelSet=superLearners, dataset, mask, radius=12, ncores=1, pruneFrac=.15, 
                                      combiner=c("greedyAUC"), bootstrapSamples=TRUE, autobalance=autobalance, searchMethod=c("replacement", "exhaustive"), 
                                      nsamples=10, calibrateProbs=FALSE, returnPredictor=FALSE) {
  
  if (length(dataset$blockVar) != length(dataset$Y)) {
    stop(paste("length of 'labels' must equal length of 'cross validation blocks'", length(Y), "!=", length(blockVar)))
  }
 
  ## iterator over blocks
  blockIterator <- FoldIterator(dataset$Y, blockVar=dataset$blockVar)
  
  ## generate an ensemble classifier for each training block
  allres <- lapply(blockIterator, function(fold) { 
    resultList <- unlist(lapply(radius, function(rad) {
      unlist(parallel::mclapply(modelSet, function(model) {
        if (searchMethod == "exhaustive") {
          searchIter <- itertools::ihasNext(RandomSearchlight(mask, rad))
        } else if (searchMethod == "replacement") {
          searchIter <- itertools::ihasNext(BootstrapSearchlight(mask, radius, nsamples))          
        }
        
        searchResults <- .searchEnsembleIteration(searchIter, dataset, fold$trainIndex, fold$testIndex, model$model, model$tuneGrid, autobalance=autobalance, bootstrap=bootstrapSamples)     
        
        searchResults <- lapply(searchResults, function(sr) { 
          attr(sr, "modelname") <- model$name 
          sr
        })
      }), recursive=FALSE)
    }), recursive=FALSE)
    
    if (combiner == "greedyAUC") {  
      ens <- greedyCombine(dataset, resultList, fold$trainIndex, fold$testIndex, calibrateProbs=calibrateProbs, pruneFrac=pruneFrac)
    } else if (combiner == "metaLearner") {
      ens <- metaCombine(dataset, resultList, fold, pruneFrac=pruneFrac)
    }
        
  })
  
  weightVol <- Reduce("+", lapply(allres, function(x) x$weightVol))/length(allres)
  AUCVol <- Reduce("+", lapply(allres, function(x) x$AUCVol))/length(allres)
  
  res <- mergeEnsembles(allres, returnPredictor)
  
  attr(res, "AUCVol") <- AUCVol 
  attr(res, "weightVol") <- weightVol
  res
}

#' @export
runAnalysis.EnsembleSearchlightModel <- function(object, dataset, vox, returnPredictor=FALSE, autobalance=FALSE, 
                                                 searchMethod="replacement", nsamples=35, radius=12, pruneFrac=.1, bootstrap=FALSE) {
  mask <- array(0, dim(dataset$mask))
  mask[vox] <- 1
  mask <- LogicalBrainVolume(mask, space(dataset$mask))
  modelSet <- .setupModels(object$baseLearners)
  
  mvpa_searchlight_ensemble(modelSet, dataset, mask, radius=radius, pruneFrac=pruneFrac, combiner="greedyAUC", bootstrapSamples=bootstrap, 
                            autobalance=autobalance, searchMethod=searchMethod, nsamples=nsamples, calibrateProbs=FALSE, returnPredictor)
  
}


## could have an ensemble of ensembles
## foreach block, find greedy ensemble of ensembles...
## AUC per voxel, per ROI, over all ROIs


