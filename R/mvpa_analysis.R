


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
computePerformance <- function(result, vox, splitList=NULL) {
  perf <- t(performance(result, splitList))
  out <- cbind(vox, perf[rep(1, nrow(vox)),])   
}

#' run an MVPA analysis on ablock of data
#' @param object the model object
#' @param dataset the MVPADataset instance
#' @return a result object
#' @export
runAnalysis <- function(object, dataset,...) {
  UseMethod("runAnalysis")
}

#' @export
runAnalysis.ClassificationModel <- function(object, dataset, vox, returnPredictor=FALSE) {
  mvpa_crossval(dataset, vox, returnPredictor)
}


#' @export
runAnalysis.SimilarityModel <- function(object, dataset, vox, returnPredictor=FALSE) {
  patternSimilarity(dataset, vox, object$simFun)
}


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
      .computePerformance(runAnalysis(dataset$model, dataset, vox), vox, dataset$testSplits)
    }
  }
  
  .convertResultsToVolumeList(res, dataset$mask)
}
  

.doRandomized <- function(dataset, radius, returnPredictor=FALSE) {
  searchIter <- itertools::ihasNext(RandomSearchlight(dataset$mask, radius))
  
  ## tight inner loop should probbaly avoid "foreach" as it has a lot of overhead.
  res <- foreach::foreach(vox = searchIter, .verbose=FALSE, .errorhandling="pass", .packages=c("rMVPA", dataset$model$library)) %do% {   
    if (nrow(vox) > 1) {  
      print(nrow(vox))
      computePerformance(runAnalysis(dataset$model, dataset, vox), vox, dataset$testSplits)
    }
  }
  
  .convertResultsToVolumeList(res, dataset$mask)
  
}



#runSearchLight <- function(searchIter, blockNum, dataset, model, tuneGrid) {
#    
#}




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
  
  cl <- makeCluster(ncores, outfile="",useXDR=FALSE, type="FORK")
  registerDoParallel(cl)
  

  
  ## Get the set of unique ROIs (all unique integers > 0 in provided mask)
  regionSet <- sort(as.integer(unique(regionMask[regionMask > 0])))
  
  res <- foreach::foreach(roinum = regionSet, .verbose=TRUE, .errorhandling="pass", .packages=c("rMVPA", "MASS", "neuroim", "caret", dataset$model$library)) %dopar% {   
    idx <- which(regionMask == roinum)
    
    if (length(idx) > 1) {
      vox <- indexToGrid(regionMask, idx)
      
      result <- runAnalysis(dataset$model, dataset, vox, savePredictors)
      attr(result, "ROINUM") <- roinum
      attr(result, "vox") <- vox
      perf <- c(ROINUM=roinum, t(performance(result, dataset$testSplits))[1,])     
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
  
  ## combine performance metrics into matrix
  perfMat <- do.call(rbind, validRes)
  
  ## put performance metrics in volumetric space
  outVols <- lapply(2:ncol(perfMat), function(cnum) {
    fill(regionMask, cbind(perfMat[, 1], perfMat[,cnum]))    
  })

  ## extract additonal results
  resultList <- lapply(validRes, function(res) {
    attr(res, "result")
  })
  
  
  extendedResults <- combineResults(dataset$model, resultList)
  names(outVols) <- colnames(perfMat)[2:ncol(perfMat)]
  list(outVols = outVols, performance=perfMat, extendedResults=extendedResults)

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
#' @param ncores the number of cores for parallel processing (default is 1)
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
  
 
  #cl <- makeCluster(ncores)
  #registerDoParallel(cl)
  
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

#' @export
combineResults <- function(resultObject, resultList, ...) {
  UseMethod("combineResults")
}

#' @export
saveResults <- function(results, folder) {
  UseMethod("saveResults")
}

#' @export
combineResults.SimilarityModel <- function(model, resultList) {
  #ret <- list(sWithin=sWithin, sBetween=sBetween, simMat=simMat)
  rois <- sapply(resultList, function(res) attr(res, "ROINUM"))
  
  simMatList <- lapply(resultList, function(res) {
    res$simMat   
  })
  
  simWithinList <- lapply(resultList, function(res) {
    res$simWithinTable  
  })
  
  names(simMatList) <- rois
  names(simWithinList) <- rois
  
  ret <- list(simMatList=simMatList, simWithinList=simWithinList)
  
  class(ret) <- c("SimilarityResultList", "list")
  ret
}


#' @export
combineResults.ClassificationModel <- function(model, resultList) {
  
  rois <- sapply(resultList, function(res) attr(res, "ROINUM"))
  
  predictorList <- lapply(resultList, function(res) {
    MVPAVoxelPredictor(res$predictor, attr(res, "vox"))
  })
  
  predFrame <- as.data.frame(do.call(rbind, lapply(resultList, function(res) {
    data.frame(ROI=rep(attr(res, "ROINUM"), length(res$observed)), observed=res$observed, pred=res$predicted, correct=as.character(res$observed) == as.character(res$predicted), prob=res$prob)
  })))
  
  
  ret <- list(predictor=ListPredictor(predictorList, rois), predictions=predFrame)
  class(ret) <- c("ClassificationResultList", "list")
  ret
}


#' @export
combineResults.EnsembleSearchlightModel <- function(model, resultList) {
  
}

#' @export
saveResults.ClassificationResultList <- function(results, folder) {
  write.table(format(results$predictions,  digits=2, scientific=FALSE, drop0trailing=TRUE), paste0(paste0(folder, "/prediction_table.txt")), row.names=FALSE, quote=FALSE)  
  if (!is.null(results$predictor)) {
    saveRDS(results$predictor, paste0(folder, "/predictor.RDS"))
  }
}

#' @export
saveResults.SimilarityResultList <- function(results, folder) {
  saveRDS(results$simMatList, paste0(folder, "/similarityMatrices.RDS"))
  omat <- as.data.frame(do.call(rbind, results$simWithinList))
  write.table(omat, paste0(folder, "/similarityWithinTable.txt"))
}

