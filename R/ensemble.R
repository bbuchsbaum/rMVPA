
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


.searchEnsembleIteration <- function(dataset, mask, trainInd, testInd, radius, model, tuneGrid) {
  searchIter <- itertools::ihasNext(RandomSearchlight(mask, radius))
  
  Ytrain <- dataset$Y[-testInd]
  blockVar <- dataset$blockVar[-testInd]
  
  res <- lapply(searchIter, function(vox) { 
    if (nrow(vox) > 2) {   
      X <- series(dataset$trainVec, vox)
      Xtrain <- X[-testInd,]    
      foldIterator <- MatrixFoldIterator(Xtrain, Ytrain, blockVar)
      cvres <- try(crossval_internal(foldIterator, model, tuneGrid, fast=TRUE, ncores=1, returnPredictor=TRUE))
      
      if (!inherits(cvres, "try-error")) {    
        cvres <- structure(classificationResult(Ytrain, cvres$class, cvres$probs, cvres$predictor),
          vox=vox,
          nvox=nrow(vox),
          param_names=names(tuneGrid),
          radius=radius,
          param=paste(sapply(1:ncol(tuneGrid), function(i) paste(param_names[i], tuneGrid[,i], sep=":")), collapse="#")
        )
  
      }
      
      cvres
    }
  })
  
  res <- res[!sapply(res, function(x) is.null(x) || inherits(x, "try-error"))]
  res
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
mvpa_searchlight_ensemble <- function(modelSet, dataset, mask, radius=c(14,10,6), ncores=1, pruneFrac=.2, combiner=c("optAUC")) {
  
  if (length(dataset$blockVar) != length(dataset$Y)) {
    stop(paste("length of 'labels' must equal length of 'cross validation blocks'", length(Y), "!=", length(blockVar)))
  }
 
  blockIterator <- FoldIterator(dataset$blockVar)
  
  foreach
  allres <- lapply(blockIterator, function(fold) { 
    resultList <- unlist(lapply(radius, function(rad) {
      unlist(lapply(models, function(model) {
        searchResults <- .searchEnsembleIteration(dataset, mask, fold$trainIndex, fold$testIndex, rad,  model$model, model$tuneGrid)     
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

runAnalysis.EnsembleSearchlightModel <- function(object, dataset, vox, returnPredictor=FALSE) {
  #browser()
  mask <- array(0, dim(dataset$mask))
  mask[vox] <- 1
  mask <- LogicalBrainVolume(mask, space(dataset$mask))
  modelSet <- .setupModels(object$baseLearners)
  
}



