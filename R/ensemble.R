
.setupModels <- function(learnerSet) {
  unlist(lapply(names(learnerSet), function(mname) {
    model <- loadModel(mname)
    params <- learnerSet[[mname]]
    lapply(1:nrow(params), function(i) {
      list(name=mname, model=model, tuneGrid=params[i,,drop=FALSE])
    }) 
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

metaCombine <- function(dataset, resultList, fold, blockNum, pruneFrac=1, metaLearner="spls", tuneGrid=expand.grid(K=c(1,2,3,4,5), eta=c(.2, .7, .9), kappa=.5)) {
  
  Ytrain <- dataset$Y[fold$trainIndex]
  #Ytest <- dataset$Y[fold$testIndex] 
  
  resultFrame <- .summarizeResults(resultList)
  aucorder <- order(resultFrame$AUC, decreasing=TRUE)
  nkeep <- pruneFrac * length(aucorder)
  keep <- sort(aucorder[1:nkeep])
  
  predmat <- do.call(cbind, lapply(resultList[keep], "[[", "probs"))
  
  foldIterator <- MatrixFoldIterator(X=predmat, Y=Ytrain, blockVar=dataset$blockVar[dataset$blockVar != blockNum])
  
  tuneControl <- caret::trainControl("cv", verboseIter=TRUE, classProbs=TRUE, index=foldIterator$getTrainSets(), indexOut=foldIterator$getTestSets()) 
  #modelFit <- trainModel(loadModel("spls"), predmat, Ytrain,  NULL, NULL, tuneGrid=expand.grid(K=c(1,2,3,4,5), eta=c(.2, .7, .9), kappa=.5), fast=FALSE, tuneControl)
  modelFit <- trainModel(loadModel(metaLearner), predmat, Ytrain,  NULL, NULL, tuneGrid=tuneGrid, fast=FALSE, tuneControl)
  #modelFit <- trainModel(loadModel("rf"), predmat, Ytrain,  NULL, NULL, tuneGrid=expand.grid(ncomp=1:5), fast=FALSE, tuneControl)
  
  ensemblePredictor <- ListPredictor(lapply(resultList[keep], function(res) MVPAVoxelPredictor(res$predictor, attr(res, "vox"))))
  
  voxlist <- lapply(resultList[keep], function(el) attr(el, "vox"))
  testlist <- lapply(voxlist, function(vox) series(dataset$trainVec,vox)[fold$testIndex,])
  
  allpreds <- do.call(cbind, lapply(seq_along(resultList[keep]), function(j) {
    evaluateModel(resultList[keep][[j]]$predictor, testlist[[j]])$probs  
  }))
  
  pfinal <- evaluateModel(modelFit, as.matrix(allpreds))
  #list(predictor=CaretPredictor(modelFit), 
  
}

greedyCombine <- function(dataset, resultList, trainIndex, testIndex, calibrateProbs=FALSE, pruneFrac=.33) {
  resultFrame <- .summarizeResults(resultList)    
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
    BaggedMultiClassOptAUC(predlist, Ytrain)
  }
  
  baggedWeights <- baggedWeights/sum(baggedWeights)
  
  posWeights <- which(baggedWeights > 0)
  finalWeights <- baggedWeights[posWeights]
  
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
  
  baggedPredictor <- WeightedPredictor(predictorList, weights=finalWeights)
  testVec <- takeVolume(dataset$trainVec, testIndex,merge=TRUE)
  testPreds <- evaluateModel(baggedPredictor, testVec)
  result <- classificationResult(Ytest, testPreds$class, testPreds$prob, baggedPredictor)
  
  weightVol <- BrainVolume(array(0, dim(dataset$mask)), space(mask))
  countVol <- BrainVolume(array(0, dim(dataset$mask)), space(mask))
  
  for (i in 1:length(voxelList)) {
    vox <- voxelList[[i]]   
    weightVol[vox] <- weightVol[vox] + finalWeights[i]
    countVol[vox] <- countVol[vox] + 1
  }
  
  nz <- which(countVol > 0)
  weightVol[nz] <- weightVol[nz]/countVol[nz]
  
  list(result=result, weightVol=weightVol)
  
}


.searchEnsembleIteration <- function(searchIter, dataset, trainInd, testInd, model, tuneGrid, autobalance=FALSE, bootstrap=FALSE) {
  #searchIter <- itertools::ihasNext(RandomSearchlight(mask, radius))
  
  Ytrain <- dataset$Y[-testInd]
  blockVar <- dataset$blockVar[-testInd]
  
  res <- lapply(searchIter, function(vox) { 
    print(nrow(vox))
    if (nrow(vox) > 2) {   
      X <- series(dataset$trainVec, vox)
      Xtrain <- X[-testInd,]    
      foldIterator <- MatrixFoldIterator(Xtrain, Ytrain, blockVar, balance=autobalance, bootstrap=bootstrap)
      cvres <- try(crossval_internal(foldIterator, model, tuneGrid, fast=TRUE, ncores=1, returnPredictor=TRUE))
      
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
  #avNNet=expand.grid(size = c(2,4,8), decay=c(.01, .001, .0001), bag=FALSE),
  #pls=data.frame(ncomp=1:5),
  sda=data.frame(lambda=c(.01, .1, .5, .9), diagonal=c(FALSE,FALSE,FALSE, FALSE)),
  spls=expand.grid(K=c(1:5), eta=c(.1,.3,.5,.7), kappa=.5)
))

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
### TODO add argument sampleMethod - "exhaustive", "replacement"//sampleIter=1000
mvpa_searchlight_ensemble <- function(modelSet=superLearners, dataset, mask, radius=12, ncores=1, pruneFrac=.2, 
                                      combiner=c("optAUC"), bootstrapSamples=TRUE, searchMethod=c("replacement", "exhaustive"), nsamples=100, calibrateProbs=FALSE) {
  
  if (length(dataset$blockVar) != length(dataset$Y)) {
    stop(paste("length of 'labels' must equal length of 'cross validation blocks'", length(Y), "!=", length(blockVar)))
  }
 
  blockIterator <- FoldIterator(dataset$Y, blockVar=dataset$blockVar)
  
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
          print(model$name)
          attr(sr, "modelname") <- model$name 
          sr
        })
      }), recursive=FALSE)
    }), recursive=FALSE)
    
    
    ens <- greedyCombine(dataset, resultList, fold$trainIndex, fold$testIndex, calibrateProbs=calibrateProbs, pruneFrac=pruneFrac)
        
  })
  
  weightVol <- Reduce("+", lapply(allres, function(x) x$weightVol))/length(allres)
}

runAnalysis.EnsembleSearchlightModel <- function(object, dataset, vox, returnPredictor=FALSE, autobalance=FALSE, 
                                                 searchMehod="replacement", nsamples=50, radius=12, pruneFrac=.2) {
  browser()
  mask <- array(0, dim(dataset$mask))
  mask[vox] <- 1
  mask <- LogicalBrainVolume(mask, space(dataset$mask))
  modelSet <- .setupModels(object$baseLearners)
  
  res <- mvpa_searchlight_ensemble(modelSet, dataset, mask, radius=radius, pruneFrac=pruneFrac, combiner="optAUC", bootstrapSamples=TRUE, searchMethod=c("replacement", "exhaustive"), nsamples=100)
  
}



