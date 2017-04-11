


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

#' ConsensusLearner
#'
#' construct specification object for a consensus learner
#' @param method the consensus method
#' @param params additional parameters for the method
#' @export
ConsensusLearner <- function(method="glmnet", params=list()) {
  ret <- list(
    method=method,
    params=params
  )

  class(ret) <- c("ConsensusLearner", "list")
  ret
}

#' consensusWeights
#'
#' compute consensusWeights for a set of classification results
#'
#' @param x the result set
#' @param ... addiitonal args
#'
#' @export
consensusWeights <- function(x,...) {
  UseMethod("consensusWeights")
}

innerIteration <- function(dataset, vox, trainInd, testInd, model, tuneGrid, autobalance, bootstrap) {
  Ytrain <- dataset$Y[-testInd]
  blockVar <- dataset$blockVar[-testInd]
  X <- series(dataset$trainVec, vox)
  Xtrain <- X[-testInd,,drop=FALSE]
  foldIterator <- MatrixFoldIterator(Xtrain, Ytrain, blockVar, balance=autobalance, bootstrap=bootstrap)

  ## do we need to compute final pedictor?
  res <- try(crossval_internal(foldIterator, model, tuneGrid, returnPredictor=TRUE, featureSelector=NULL))

  structure(classification_result(Ytrain, res$class, res$probs, res$predictor),
            vox=vox)

}

#' @export
mvpa_regional_consensus <- function(dataset, model, regionMask, autobalance=FALSE, bootstrap=FALSE,
                                    savePredictors=TRUE, classMetrics=TRUE, featureSelector=NULL,
                                    method=c("greedy", "glmnet", "equal_weights", "auc_weights"), ...) {

  if (length(dataset$blockVar) != length(dataset$Y)) {
    stop(paste("length of 'labels' must equal length of 'cross validation blocks'", length(Y), "!=", length(blockVar)))
  }

  method <- match.arg(method)
  crossVal <- blocked_cross_validation(dataset$blockVar, balance=autobalance, bootstrap=bootstrap)

  lookupMetaLearner <- function(method) {
    switch(method,
           "greedy"=greedyWeights,
           "glmnet"=glmnetWeights,
           "auc_weights"=AUCWeights,
           "equal_weights"=equalWeights)

  }

  learner <- lookupMetaLearner(method)

  regionSet <- sort(as.integer(unique(regionMask[regionMask > 0])))

  blockIterator <- FoldIterator(dataset$Y, blockVar=dataset$blockVar)

  Xtest <- series(dataset$trainVec, which(regionMask>0))

  ## loop over blocks
  ensembleSet <- foreach::foreach(fold = blockIterator, .verbose=FALSE, .packages=c("rMVPA", "MASS", "neuroim", "caret", dataset$model$library)) %dopar% {

    ## loop over rois
    resultList <- lapply(regionSet, function(roinum) {
      idx <- which(regionMask == roinum)

      ## choke on length == 1
      if (length(idx) < 2) {
        stop(paste("ROI number", roinum, "has less than 2 voxels, aborting."))
      }

      vox <- indexToGrid(regionMask, idx)
      ## train on sub indices
      result <- model$run(dataset, vox, crossVal, featureSelector, subIndices=fold$trainIndex)

      attr(result, "ROINUM") <- roinum
      attr(result, "vox") <- vox
      result

    })

    ## learn consensus weights
    wts <- learner(resultList, ...)

    if (all(wts==0)) {
      warning("consensus weights are all zero, assigning equal weight to each classifier")
      wts <- rep(1, length(wts))/length(wts)
    }

    combined <- model$combineResults(resultList)

    ## create weighted predictor for this fold
    wpred <- WeightedPredictor(combined$predictor, weights=wts)

    ## predict on held out data
    fullPred <- evaluateModel(wpred, dataset$trainVec, subIndices=fold$testIndex)
    roiPred <- evaluateModel(combined$predictor, dataset$trainVec, subIndices=fold$testIndex)

    list(predictor=wpred, predictions=fullPred, roiPred=roiPred, weights=wts)

  }

  #browser()

  prob <- do.call(rbind, lapply(ensembleSet, function(x) x$predictions$prob))
  #roiProb <- do.call(rbind, lapply(ensembleSet, function(x) x$roiPred$prob))
  testInd <- unlist(blockIterator$getTestSets())
  prob <- prob[testInd,]
  finalWeights <- colMeans(do.call(rbind, lapply(ensembleSet, function(e) unlist(e$weights))))

  maxid <- apply(prob,1,which.max)
  predclass <- colnames(prob)[maxid]

  finalPredictor <- WeightedPredictor(lapply(ensembleSet, "[[", "predictor"))
  result <- classification_result(dataset$Y, predclass, prob=prob, predictor=finalPredictor)

  list(result=result, weights=finalWeights)

}

#' @export
consensusWeights.ClassificationResultSet <- function(x, method=c("greedy", "glmnet", "equal_weights", "auc_weights"), ...) {
  blocks <- sort(unique(x$blockVar))
  observed <- x$resultList[[1]]$observed

  method <- match.arg(method)
  lookupMetaLearner <- function(meth) {
    switch(meth,
           "greedy"=greedyWeights,
           "glmnet"=glmnetWeights,
           "auc_weights"=AUCWeights,
           "equal_weights"=equalWeights)

  }

  learner <- lookupMetaLearner(method)

  res <- lapply(blocks, function(block) {
    heldout <- which(x$blockVar == block)
    trainInd <- which(x$blockVar != block)

    sres <- lapply(x$resultList, function(result) {
      subResult(result, trainInd)
    })

    bvar <- x$blockVar[trainInd]

    wts <- learner(sres, ...)

    if (all(wts==0)) {
      warning("consensus weights are all zero, assigning equal weight to each classifier")
      wts <- rep(1, length(wts))/length(wts)
    }

    probset <- lapply(x$resultList, function(x) x$probs[heldout,])
    wtprob <- Reduce("+", lapply(1:length(probset), function(i) probset[[i]] * wts[i]))
    list(prob=wtprob, weights=wts, testInd=heldout)
  })

  prob <- do.call(rbind, lapply(res, "[[", "prob"))
  ind <- unlist(lapply(res, "[[", "testInd"))
  weights <- do.call(rbind, lapply(res, "[[", "weights"))
  maxid <- apply(prob,1,which.max)
  predclass <- colnames(prob)[maxid]
  finalWeights <- colMeans(weights)
  finalWeights <- finalWeights/sum(finalWeights)

  result <- if (!is.null(x$resultList[[1]]$predictor)) {
    predictor <- WeightedPredictor(lapply(x$resultList, "[[", "pred"), weights=finalWeights)
    classification_result(observed, predclass, prob=prob[ind,], predictor=predictor)
  } else {
    classification_result(observed, predclass, prob=prob[ind,])
    predictor <- NULL
  }

  list(result=result, weights=finalWeights)
}

equalWeights <- function(resultList) {
  wts <- rep(1,length(resultList))/length(resultList)
}

AUCWeights <- function(resultList, pruneFrac=1, power=2) {
  AUC <- unlist(lapply(resultList, function(x) performance(x)[3]))

  aucorder <- order(AUC, decreasing=TRUE)
  finalWeights <- AUC

  if (pruneFrac < 1) {
    nkeep <- pruneFrac * length(aucorder)
    keep <- sort(aucorder[1:nkeep])
    finalWeights[-keep] <- 0
  }

  finalWeights[finalWeights < 0] <- 0

  finalWeights <- finalWeights^power
  finalWeights <- finalWeights/sum(finalWeights)
  finalWeights
}

#' @importFrom glmnet cv.glmnet
glmnetWeights <- function(resultList, alpha=.5) {
  Ytrain <- as.character(resultList[[1]]$observed)
  probset <- lapply(resultList, "[[", "probs")

  posmat <- do.call(rbind, lapply(seq_along(Ytrain), function(i) {
    lab <- Ytrain[i]
    sapply(seq_along(probset), function(j) probset[[j]][i,lab])
  }))

  negmat <- do.call(rbind, lapply(seq_along(Ytrain), function(i) {
    lab <- Ytrain[i]
    idx <- which(colnames(probset[[1]]) != lab)
    colMeans(sapply(seq_along(probset), function(j) probset[[j]][i,idx]))
  }))

  fullmat <- rbind(posmat, negmat)
  y0 <- factor(rep(c("pos", "neg"), each=nrow(posmat)))

  ctrl <- trainControl(method="cv", number=5, repeats=20)
  cvres <- caret::train(fullmat, y0, method="glmnet", tuneLength=8, trControl=ctrl, lower.limits=0)

  #cvres <- glmnet::cv.glmnet(fullmat, y0, family="binomial", lower.limits=0, alpha=alpha)
  #lambda.min <- cvres$lambda.min
  weights <- coef(cvres$finalModel, s=cvres$bestTune$lambda)[-1,1]

  if (all(weights == 0)) {
    weights <- rep(1/length(weights), length(weights))
  } else {
    weights <- weights/sum(weights)
  }
}


AUCCombine <- function(resultList, blockVar, pruneFrac=.5, power=2) {
  finalWeights <- AUCWeights(resultList, blockVar, pruneFrac, power)
  posWeights <- which(finalWeights>0)

  positiveList <- resultList[posWeights]

  predictorList <- lapply(positiveList, function(pred) {
      MVPAVoxelPredictor(pred$predictor, attr(pred, "vox"))
  })

  ## weighted predictor based on finalWeights
  WeightedPredictor(predictorList, weights=finalWeights[posWeights])

}



glmnetCombine <- function(resultList, alpha=.5) {
  weights <- glmnetWeights(resultList, alpha)
  positiveList <- positiveList[weights > 0]

  predictorList <- lapply(positiveList, function(pred) {
    MVPAVoxelPredictor(pred$predictor, attr(pred, "vox"))
  })

  WeightedPredictor(predictorList, weights=weights[weights>0])

}

greedyWeights <- function(resultList, pruneFrac=1) {

  stopifnot(pruneFrac > 0 && pruneFrac <= 1)
  AUC <- unlist(lapply(resultList, function(x) performance(x)[3]))
  Ytrain <- as.character(resultList[[1]]$observed)

  ## order from highest to lowest
  perforder <- order(AUC, decreasing=TRUE)

  ## results to keep according to 'pruneFrac'
  keep <- sort(perforder[1:(pruneFrac * length(perforder))])
  prunedList <- resultList[keep]

  predlist <- lapply(prunedList, "[[", "probs")

  baggedWeights <- if (length(levels(Ytrain)) == 2) {
    Pred <- do.call(cbind, lapply(predlist, function(x) x[,2]))
    BaggedTwoClassOptAUC(Pred, ifelse(Ytrain == levels(Ytrain)[2], 1,0))
  } else {
    BaggedMultiClassOptAUC(predlist, Ytrain)
  }

  baggedWeights/sum(baggedWeights)
}


greedyCombine <- function(resultList, pruneFrac=1) {
  stopifnot(pruneFrac > 0 && pruneFrac <= 1)
  weights <- greedyWeights(resultList, pruneFrac)
  posIndices <- which(weights > 0)

  ## weights for final model
  finalWeights <- weights[posIndices]
  positiveList <- resultList[posIndices]

  predictorList <- lapply(positiveList, function(pred) {
    MVPAVoxelPredictor(pred$predictor, attr(pred, "vox"))
  })

  ## weighted predictor based on finalWeights
  WeightedPredictor(predictorList, weights=finalWeights)
}



metaCombine <- function(dataset, resultList, fold, pruneFrac=1, metaLearner="spls", tuneGrid=expand.grid(K=c(1,2,3,4,5,6), eta=c(.2, .7, .9), kappa=.5)) {
  AUC <- do.call(rbind, lapply(resultList, function(x) performance(x)))[,3]
  Ytrain <- dataset$Y[fold$trainIndex]
  aucorder <- order(AUC, decreasing=TRUE)
  nkeep <- pruneFrac * length(aucorder)
  keep <- sort(aucorder[1:nkeep])

  predmat <- do.call(cbind, lapply(resultList[keep], "[[", "probs"))

  foldIterator <- MatrixFoldIterator(X=predmat, Y=Ytrain, NULL, blockVar=dataset$blockVar[dataset$blockVar != fold$index])

  tuneControl <- caret::trainControl("cv", verboseIter=TRUE, classProbs=TRUE, index=foldIterator$getTrainSets(), indexOut=foldIterator$getTestSets())
  modelFit <- trainModel(loadModel(metaLearner), predmat, Ytrain,  NULL, NULL, tuneGrid=tuneGrid, tuneControl)


  ## need a MetaPredictor class
  ensemblePredictor <- ListPredictor(lapply(resultList[keep], function(res) MVPAVoxelPredictor(res$predictor, attr(res, "vox"))))
  testVec <- takeVolume(dataset$trainVec, fold$testIndex,merge=TRUE)
  allpreds <- do.call(cbind, lapply(evaluateModel(ensemblePredictor, testVec), function(p) p$probs))
  pfinal <- evaluateModel(modelFit, as.matrix(allpreds))
  ## need a MetaPredicotr class


  result <- classification_result(fold$Ytest, pfinal$class, pfinal$probs)

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


computeWeightVol <- function(ensemble, mask) {
  weightVol <- BrainVolume(array(0, dim(mask)), space(mask))
  countVol <- BrainVolume(array(0, dim(mask)), space(mask))
  wts <- attr(ens, "weights")


  for (i in 1:length(ensemble)) {
    vox <- ensemble[[i]]$voxelGrid
    weightVol[vox] <- weightVol[vox] + wts[i]
  }

  weightVol
}





ensembleIteration <- function(searchIter, dataset, fold, modelspec, autobalance=FALSE, bootstrap=FALSE) {

  Ytrain <- dataset$Y[-fold$testIndex]
  blockVar <- dataset$blockVar[-fold$testIndex]

  res <- foreach(vox=searchIter) %do% {
    if (nrow(vox) > 1) {
      cvres <- innerIteration(dataset, vox, fold$trainIndex, fold$testIndex, modelspec$model, modelspec$tuneGrid, autobalance, bootstrap)
      param_names=names(modelspec$tuneGrid)
      if (!inherits(cvres, "try-error")) {
        cvres <- structure(classification_result(Ytrain, cvres$predicted, cvres$probs, cvres$predictor),
          vox=vox,
          nvox=nrow(vox),
          modelname=modelspec$model$name,
          param_names=param_names,
          param=paste(sapply(1:ncol(modelspec$tuneGrid), function(i) paste(param_names[i], modelspec$tuneGrid[,i], sep=":")), collapse="#")
        )
      }

      cvres
    }
  }

  res <- res[!sapply(res, function(x) is.null(x) || inherits(x, "try-error"))]
  res
}



mergeEnsembles <- function(ensembleSet, testOrder, returnPredictor=FALSE) {

  #testOrder <- blockIterator$getTestOrder()
  predicted <- unlist(lapply(ensembleSet, function(x) x$result$predicted))[testOrder]
  observed <- unlist(lapply(ensembleSet, function(x) x$result$observed))[testOrder]
  probs <- do.call(rbind, lapply(ensembleSet, function(x) x$result$probs))[testOrder,]

  if (returnPredictor) {
    ## predictor equally weights each ensemble
    predictor <- WeightedPredictor(lapply(ensembleSet, function(x) x$result$predictor))
    classification_result(observed, predicted, probs, predictor)
  } else {
    classification_result(observed, predicted, probs)
  }

}

mvpa_roi_ensemble <- function(dataset, regionMask,
                              pruneFrac=.5,
                              combiner="greedyAUC",
                              bootstrapSamples=TRUE,
                              autobalance=FALSE,
                              returnPredictor=FALSE,
                              ensembleSpec=NULL) {

  if (length(dataset$blockVar) != length(dataset$Y)) {
    stop(paste("length of 'labels' must equal length of 'cross validation blocks'", length(Y), "!=", length(blockVar)))
  }

  #cl <- makeCluster(ncores, outfile="",useXDR=FALSE, type="FORK")
  #registerDoParallel(cl)

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
    ens <- greedyCombine(dataset, resultList, fold, pruneFrac=pruneFrac)

    preds <- evaluateModel(ens, dataset$trainVec)
    ## get results
  }



  finalResult <- mergeEnsembles(ensembleSet, blockIterator$getTestOrder(), returnPredictor)

}

genSearchlight <- function(mask, radius, type=c("exhaustive", "replacement"), nsamples=NULL) {
  if (type == "exhaustive") {
    itertools::ihasNext(RandomSearchlight(mask, radius))
  } else if (type == "replacement") {
    itertools::ihasNext(neuroim:::BootstrapSearchlight(mask, radius, nsamples))
  } else {
    stop(paste("unrecognized searchlight type: ", type))
  }
}

combineResults <- function(dataset, resultList, fold, pruneFrac, combiner, ...) {
  if (combiner == "greedy") {
    greedyCombine(dataset, resultList, fold, pruneFrac=pruneFrac,...)
  } else if (combiner == "glmnet") {
    glmnetCombine(dataset, resultList, fold, pruneFrac=pruneFrac,...)
  } else if (combiner == "naive") {
    AUCCombine(dataset, resultList, fold, pruneFrac=pruneFrac,...)
  } else {
    stop(paste("unrecognized 'combiner' ", combiner))
  }

}


#' mvpa_searchlight_ensemble
#' @param modelSet a list of classifiers
#' @param dataset a \code{MVPADataset} instance.
#' @param mask a \code{BrainVolume} a mask volume
#' @param pruneFrac the percentage of models to retain for meta-classifier
#' @param ncores the number of cores for parallel processign (default is 1)
#' @return a named list of \code{BrainVolume} objects, where each name indicates the performance metric and label (e.g. accuracy, AUC)
#' @import itertools
#' @import foreach
#' @import doParallel
#' @import parallel
#' @importFrom assertthat assert_that
#' @export
mvpa_searchlight_ensemble <- function(modelSet=superLearners, dataset, mask, radius=12, pruneFrac=.15,
                                      combiner=c("greedy", "glmnet", "naive")[1],
                                      bootstrapSamples=TRUE, autobalance=autobalance,
                                      searchMethod=c("replacement", "exhaustive"),
                                      nsamples=10, calibrateProbs=FALSE, returnPredictor=FALSE) {

  ### returns global performance + local performance + weightvol

  assert_that(length(dataset$blockVar) == length(dataset$Y))

  ## iterator over blocks
  blockIterator <- FoldIterator(dataset$Y, blockVar=dataset$blockVar)

  ## generate an ensemble classifier for each training block
  allres <- lapply(blockIterator, function(fold) {
    resultList <-
      unlist(parallel::mclapply(modelSet, function(model) {
        searchIter <- genSearchlight(mask, radius, searchMethod, nsamples)

        searchResults <- ensembleIteration(searchIter, dataset,
                                                  fold,
                                                  model,
                                                  autobalance=autobalance,
                                                  bootstrap=bootstrapSamples)
      }), recursive=FALSE)


    ens <- combineResults(dataset, resultList, fold, pruneFrac=pruneFrac, combiner)
    preds <- evaluateModel(ens, subVector(dataset$trainVec, fold$testIndex))
    list(preds=preds, ens=ens, weightvol=computeWeightVol(ens, dataset$mask))

  })

  testInd <- unlist(blockIterator$getTestSets())
  pclass <- unlist(lapply(allres, function(x) x$preds$class))[testInd]
  prob <- do.call(rbind, lapply(allres, function(x) x$preds$prob))[testInd,]
  result <- classification_result(dataset$Y, pclass, prob)

  weightvol <- Reduce("+", lapply(allres, function(x) x$weightvol))/length(allres)

  #AUCVol <- Reduce("+", lapply(allres, function(x) x$AUCVol))/length(allres)

  res <- mergeEnsembles(allres, returnPredictor)

  attr(res, "AUCVol") <- AUCVol
  attr(res, "weightVol") <- weightVol
  res
}




## could have an ensemble of ensembles
## foreach block, find greedy ensemble of ensembles...
## AUC per voxel, per ROI, over all ROIs


