


colACC <- function(X, Y) {
  apply(X, 2, function(p) {
    sum(p == Y)/length(Y)
  })
}


matrixToVolumeList <- function(vox, mat, mask, default=NA) {
  lapply(1:ncol(mat), function(i) {
    vol <- array(default, dim(mask))   
    vol[vox] <- mat[,i]
    BrainVolume(vol, space(mask))
  })
} 

# helper function to compute performance metric for a classifcation result
computePerformance <- function(result, vox, splitList=NULL, classMetrics=FALSE, customFun=NULL) {
  perf <- t(performance(result, splitList, classMetrics))
  if (!is.null(customFun)) {
    perf <- cbind(perf, t(customPerformance(result, customFun, splitList)))
  }
  out <- cbind(vox, perf[rep(1, nrow(vox)),])   
}

#' mvpa_crossval
#' 
#' workhorse function for cross_validation in an MVPA analysis
#' 
#' @param dataset an instance of type \code{MVPADataset}
#' @param vox a \code{matrix} of voxel coordinates or \code{vector} of ids defining the subset of the image dataset to use.
#' @param crossVal a cross-validation instance of type \code{CrossValidation}
#' @param model a \code{caret} model object
#' @param tuneGrid an optional caret-formatted \code{data.frame} containing parameters to be tuned
#' @param featureSelector an optional \code{FeatureSelector} instance for selecting relevant features
#' @param subIndices an optional subset vector of row indices (observations) to run cross-validation on.
#' @return a list consisting of the following named elements:
#' prediction
#' predictor  
#' testIndices 
#' featureMask 
#'  
#' @export
mvpa_crossval <- function(dataset, vox, crossVal, model, tuneGrid=NULL, featureSelector = NULL, subIndices=NULL) {
  
  X <- series(dataset$trainVec, vox)
  
  valid.idx <- nonzeroVarianceColumns(X)
  X <- X[,valid.idx]
    
  if (ncol(X) == 0) {
    stop("mvpa_crossval: no valid columns in data matrix")
  }
    
  vox <- if (is.matrix(vox)) {
    vox[valid.idx,]
  } else {
    vox[valid.idx]
  }
  
  parcels <- if (!is.null(dataset$parcellation)) {
    dataset$parcellation[vox]
  }
    
  if (is.null(tuneGrid)) {
    tuneGrid <- model$grid(X, Y, 1)
  }
  
  foldIterator <- matrixIter(crossVal, dataset$Y, X, subIndices)
    
  result <- if (is.null(dataset$testVec)) {
    ### train and test on one set.
    crossval_internal(foldIterator, model, tuneGrid, featureSelector = featureSelector, parcels = parcels)
  } else {
    ### train on one set, test on the other.
    Xtest <- series(dataset$testVec, vox)
    crossval_external(foldIterator, Xtest, dataset$testY, model, tuneGrid, featureSelector = featureSelector, parcels = parcels)
  }
    
  result
}

# helper method to convert result list to a list of BrainVolumes 
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


# standard searchlight
.doStandard <- function(dataset, model, radius, crossVal, classMetrics=FALSE) {
  searchIter <- itertools::ihasNext(Searchlight(dataset$mask, radius)) 
  
  res <- foreach::foreach(vox = searchIter, .verbose=FALSE, .packages=c("rMVPA", "MASS", "neuroim", "caret", model$model$library)) %dopar% {   
    if (nrow(vox) > 1) {
      print(nrow(vox))
      vox <- ROIVolume(space(dataset$mask), vox)
      result <- model$run(dataset, vox, crossVal)
      perf <- computePerformance(result, coords(vox), dataset$testSplits, classMetrics, model$customPerformance)
    }
  }
  
  .convertResultsToVolumeList(res, dataset$mask)
}
  
# randomized searchlight
.doRandomized <- function(dataset, model, radius, crossVal, classMetrics=FALSE) {
  searchIter <- itertools::ihasNext(RandomSearchlight(dataset$mask, radius))
  
  ## tight inner loop should probably avoid "foreach" as it has a lot of overhead, but c'est la vie for now.
  res <- foreach::foreach(vox = searchIter, .verbose=FALSE, .packages=c("rMVPA", "MASS", "caret", "neuroim", model$model$library)) %do% {   
    if (nrow(vox) > 1) {  
      print(nrow(vox))
      vox <- ROIVolume(space(dataset$mask), vox)
      result <- model$run(dataset, vox, crossVal)
      perf <- computePerformance(result, coords(vox), dataset$testSplits, classMetrics, model$customPerformance)
    }
  }
  
  .convertResultsToVolumeList(res, dataset$mask)
  
}


## mvpa_hiearchical ??


#' mvpa_regional
#' 
#' Run a separate MVPA analysis for multiple disjoint regions of interest.
#' 
#' @param dataset a \code{MVPADataset} instance.
#' @param model a \code{BaseModel} instance usually of type \code{CaretModelWrapper}
#' @param regionMask a \code{BrainVolume} where each region is identified by a unique integer. Every non-zero set of positive integers will be used to define a set of voxels for clasisifcation analysis.
#' @param crossVal
#' @param savePredictors whether to return prediction model in result (default is \code{FALSE} to save memory)
#' @param featureSelector an option \code{FeatureSelector} object that is used to subselect informative voxels in feature space.
#' @param classMetrics
#' 
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
      vox <- ROIVolume(space(regionMask), indexToGrid(regionMask, idx))
      
      result <- model$run(dataset, vox, crossVal, featureSelector)
    
      attr(result, "ROINUM") <- roinum
      attr(result, "vox") <- coords(vox)
      
     
      perf <- if (!is.null(model$customPerformance)) {
        standard <- performance(result, dataset$testSplits, classMetrics)
        custom_perf <- customPerformance(result, model$customPerformance, dataset$testSplits)
        c(standard, custom_perf)
      } else {
        t(performance(result, dataset$testSplits, classMetrics))[1,]
  
      }
      
      perf <- c(ROINUM=roinum, perf)     
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
  list(outVols = outVols, 
       performance=perfMat, 
       resultSet=resultSet, 
       extendedResults=extendedResults, 
       invalid=regionSet[invalid])

}

  
  
#' mvpa_searchlight
#' @param dataset a \code{MVPADataset} instance.
#' @param model
#' @param crossVal
#' @param radius the searchlight radus in mm
#' @param method the type of searchlight (randomized, or standard)
#' @param niter the number of searchlight iterations for 'randomized' method
#' @param classMetrics
#' @return a named list of \code{BrainVolume} objects, where each name indicates the performance metric and label (e.g. accuracy, AUC)
#' @import itertools 
#' @import foreach
#' @import doParallel
#' @import parallel
#' @import futile.logger
#' @export
mvpa_searchlight <- function(dataset, model, crossVal, radius=8, method=c("randomized", "standard"),  
                             niter=4, classMetrics=FALSE) {
  stopifnot(niter > 1)
  
  if (radius < 1 || radius > 100) {
    stop(paste("radius", radius, "outside allowable range (1-100)"))
  }
  
  method <- match.arg(method)
  
  flog.info("model is: %s", model$model_name)
  
  res <- if (method == "standard") {
    .doStandard(dataset, model, radius, crossVal, classMetrics=classMetrics)    
  } else {
    
    res <- foreach(i = 1:niter) %dopar% {
      flog.info("Running randomized searchlight iteration %s", i)   
      do.call(cbind, .doRandomized(dataset, model, radius, crossVal, classMetrics=classMetrics) )
    }
    
    ## average over iterations
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





