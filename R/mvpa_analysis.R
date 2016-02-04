


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
computePerformance <- function(result, vox, splitList=NULL, classMetrics=FALSE) {
  perf <- t(performance(result, splitList, classMetrics))
  out <- cbind(vox, perf[rep(1, nrow(vox)),])   
}



#' @export
mvpa_crossval <- function(dataset, vox, crossVal, model, tuneGrid=NULL, featureSelector = NULL) {
  X <- series(dataset$trainVec, vox)
  valid.idx <- nonzeroVarianceColumns(X)
    
  X <- X[,valid.idx]
    
  if (ncol(X) == 0) {
    stop("no valid columns")
  }
    
  vox <- vox[valid.idx,]
    
  parcels <- if (!is.null(dataset$parcellation)) {
    dataset$parcellation[vox]
  }
    
    
  if (is.null(tuneGrid)) {
    tuneGrid <- model$grid(X, dataset$Y, 1)
  }
  
  foldIterator <- matrixIter(crossVal, dataset$Y, X, vox)
    
  result <- if (is.null(dataset$testVec)) {
      ### train and test on one set.
    crossval_internal(foldIterator, model, tuneGrid, featureSelector = featureSelector, parcels = parcels)
  } else {
    ### train on one set, test on the other.
    Xtest <- series(dataset$testVec, vox)
    crossval_external(foldIterator, Xtest, dataset$testY, model, tuneGrid, featureSelector = featureSelector, parcels =parcels)
  }
    
  result
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


.doStandard <- function(dataset, model, radius, crossVal, classMetrics=FALSE) {
  searchIter <- itertools::ihasNext(Searchlight(dataset$mask, radius)) 
  
  res <- foreach::foreach(vox = searchIter, .verbose=FALSE) %dopar% {   
    if (nrow(vox) > 1) {
      print(nrow(vox))
      computePerformance(model$run(dataset, vox, crossVal), vox, dataset$testSplits, classMetrics)
    }
  }
  
  .convertResultsToVolumeList(res, dataset$mask)
}
  

.doRandomized <- function(dataset, model, radius, crossVal, classMetrics=FALSE) {
  searchIter <- itertools::ihasNext(RandomSearchlight(dataset$mask, radius))
  
  ## tight inner loop should probably avoid "foreach" as it has a lot of overhead, but c'est la vie for now.
  res <- foreach::foreach(vox = searchIter, .verbose=FALSE, .errorhandling="pass", .packages=c("rMVPA", dataset$model$library)) %do% {   
    if (nrow(vox) > 1) {  
      print(nrow(vox))
      computePerformance(model$run(dataset, vox, crossVal), vox, dataset$testSplits, classMetrics)
    }
  }
  
  .convertResultsToVolumeList(res, dataset$mask)
  
}


## TODO mvpa_hiearchical ??


#' mvpa_regional
#' @param dataset a \code{MVPADataset} instance.
#' @param regionMask a \code{BrainVolume} where each region is identified by a unique integer. Every non-zero set of positive integers will be used to define a set of voxels for clasisifcation analysis.
#' @param savePredictors whether to return prediction model in result (default is \code{FALSE} to save memory)
#' @param autobalance whether to subsample training set so that classes are balanced (default is \code{FALSE})
#' @param bootstrapReplications the number of bootstrapped samples per each cross-validation fold
#' @param featureSelector an option \code{FeatureSelector} object that is used to subselect informative voxels in feature space.
#' @param ensemblePredictor whether returned predictor object averages over cross-validation blocks
#' @return a named list of \code{BrainVolume} objects, where each name indicates the performance metric and label (e.g. accuracy, AUC)
#' @import itertools 
#' @import foreach
#' @import doParallel
#' @import parallel
#' @export
mvpa_regional <- function(dataset, model, regionMask, crossVal=KFoldCrossValidation(length(dataset$Y)), 
                          savePredictors=FALSE, featureSelector=NULL, classMetrics=FALSE) {  
  
  if (!is.null(dataset$parcellation)) {
    if (! all(dim(dataset$parcellation) == dim(regionMask))) {
      stop("dimension of 'featureParcellation' must equals dimensions of 'regionMask'")
    }
  }
  
  ## Get the set of unique ROIs (all unique integers > 0 in provided mask)
  regionSet <- sort(as.integer(unique(regionMask[regionMask > 0])))
  
  res <- foreach::foreach(roinum = regionSet, .verbose=FALSE, .errorhandling="pass", .packages=c("rMVPA", "MASS", "neuroim", "caret", dataset$model$library)) %dopar% {   
    idx <- which(regionMask == roinum)
    
    if (length(idx) > 1) {
      ## what if length is less than 1?
      vox <- indexToGrid(regionMask, idx)
      
      result <- model$run(dataset, vox, crossVal, featureSelector)
    
      attr(result, "ROINUM") <- roinum
      attr(result, "vox") <- vox
      
      perf <- c(ROINUM=roinum, t(performance(result, dataset$testSplits, classMetrics))[1,])     
      perf <- structure(perf, result=result)
      perf    
    }
  }
  
  invalid <- sapply(res, function(x) inherits(x, "simpleError") || is.null(x))
  validRes <- res[!invalid]
  
  if (length(validRes) == 0) {
    print(res)
    flog.error("Regional analysis failed for all of %s ROIs", length(regionSet))
    stop("error in mvpa_regional: aborting")
  } 
  
  if (sum(invalid) > 0) {
    flog.warn("Regional analysis failed for %s ROIs", sum(invalid))
    flog.warn("List of ROIs that failed to complete: ", regionSet[invalid], capture=TRUE)
    print(res[invalid])
  }
  
  
  ## extract additonal results
  results <- lapply(validRes, function(res) {
    attr(res, "result")
  })
 
  ## combine performance metrics into matrix
  perfMat <- do.call(rbind, validRes)
  
  ## put performance metrics in volumetric space
  outVols <- lapply(2:ncol(perfMat), function(cnum) {
    fill(regionMask, cbind(perfMat[, 1], perfMat[,cnum]))    
  })
  
  resultSet <- ClassificationResultSet(dataset$blockVar, results)

  extendedResults <- model$combineResults(results)
  names(outVols) <- colnames(perfMat)[2:ncol(perfMat)]
  list(outVols = outVols, performance=perfMat, 
       resultSet=resultSet, 
       extendedResults=extendedResults, 
       invalid=regionSet[invalid])

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
mvpa_searchlight <- function(dataset, model, crossVal, radius=8, method=c("randomized", "standard"),  niter=4, classMetrics=FALSE) {
  stopifnot(niter > 1)
  
  if (radius < 1 || radius > 100) {
    stop(paste("radius", radius, "outside allowable range (1-100)"))
  }
  
  
  method <- match.arg(method)
  
  flog.info("classification model is: %s", model$model_name)
  
  
  res <- if (method == "standard") {
    .doStandard(dataset, model, radius, crossVal, classMetrics=classMetrics)    
  } else {
    res <- foreach::foreach(i=1:niter) %dopar% {
      flog.info("Running randomized searchlight iteration %s", i)   
      do.call(cbind, .doRandomized(dataset, model, radius, crossVal, classMetrics=classMetrics) )
    }
    
   
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


# #' @export
# combineResults.SimilarityModel <- function(model, resultList) {
#   #ret <- list(sWithin=sWithin, sBetween=sBetween, simMat=simMat)
#   rois <- sapply(resultList, function(res) attr(res, "ROINUM"))
#   
#   simMatList <- lapply(resultList, function(res) {
#     res$simMat   
#   })
#   
#   simWithinList <- lapply(resultList, function(res) {
#     res$simWithinTable  
#   })
#   
#   names(simMatList) <- rois
#   names(simWithinList) <- rois
#   
#   ret <- list(simMatList=simMatList, simWithinList=simWithinList)
#   
#   class(ret) <- c("SimilarityResultList", "list")
#   ret
# }
# 
# 




# combineResults.EnsembleSearchlightModel <- function(model, resultList) {
#   rois <- sapply(resultList, function(res) attr(res, "ROINUM"))
#   
#   predictorList <- lapply(resultList, function(res) {
#     MVPAVoxelPredictor(res$predictor, attr(res, "vox"))
#   })
#   
#   predFrame <- as.data.frame(do.call(rbind, lapply(resultList, function(res) {
#     data.frame(ROI=rep(attr(res, "ROINUM"), length(res$observed)), observed=res$observed, pred=res$predicted, correct=as.character(res$observed) == as.character(res$predicted), prob=res$prob)
#   })))
#   
#   weightVol <- Reduce("+", lapply(resultList, function(res) attr(res, "weightVol")))
#   AUCVol <- Reduce("+", lapply(resultList, function(res) attr(res, "AUCVol")))
#   
#   ret <- list(predictor=ListPredictor(predictorList, rois), predictions=predFrame, weightVol=weightVol, AUCVol=AUCVol)
#   class(ret) <- c("EnsembleSearchlightResultList", "list")
#   ret
#   
# }
# 
# #' @export
# saveResults.EnsembleSearchlightResultList <- function(results, folder) {
#   write.table(format(results$predictions,  digits=2, scientific=FALSE, drop0trailing=TRUE), paste0(paste0(folder, "/prediction_table.txt")), row.names=FALSE, quote=FALSE)  
#   if (!is.null(results$predictor)) {
#     saveRDS(results$predictor, paste0(folder, "/predictor.RDS"))
#   }
#   
#   writeVolume(results$AUCVol, paste0(folder, "/AUC_Scores.nii"))
#   writeVolume(results$weightVol, paste0(folder, "/EnsembleWeights.nii"))
# }
# 
# #' @export
# saveResults.ClassificationResultList <- function(results, folder) {
#   write.table(format(results$predictions,  digits=2, scientific=FALSE, drop0trailing=TRUE), paste0(paste0(folder, "/prediction_table.txt")), row.names=FALSE, quote=FALSE)  
#   if (!is.null(results$predictor)) {
#     saveRDS(results$predictor, paste0(folder, "/predictor.RDS"))
#   }
# }
# 
# #' @export
# saveResults.SimilarityResultList <- function(results, folder) {
#   saveRDS(results$simMatList, paste0(folder, "/similarityMatrices.RDS"))
#   omat <- as.data.frame(do.call(rbind, results$simWithinList))
#   write.table(omat, paste0(folder, "/similarityWithinTable.txt"))
# }

