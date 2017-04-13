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
computePerformance <- function(result, splitList=NULL, class_metrics=FALSE, customFun=NULL) {
  perf <- t(performance(result, splitList, class_metrics))
  
  if (!is.null(customFun)) {
    cbind(perf, t(customPerformance(result, customFun, splitList)))
  } else {
    perf
  }
  
  
  ##out <- cbind(vox, perf[rep(1, nrow(vox)),])   
}

#' mvpa_crossval
#' 
#' workhorse function for cross_validation in an MVPA analysis
#' 
#' @param dataset an instance of type \code{MVPADataset}
#' @param ROI a class of type \code{ROIVolume} or \code{ROISurface} contianing the training data.
#' @param crossval a cross-validation instance of type \code{CrossValidation}
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
mvpa_crossval <- function(dataset, ROI, crossval, model, tuneGrid=NULL, featureSelector = NULL, subIndices=NULL) {
  
  
  ## valid subset
  valid.idx <- nonzeroVarianceColumns(values(ROI))
  
  if (length(valid.idx) < 2) {
    stop("mvpa_crossval: fewer than 2 valid columns in data matrix")
  }
  
  if (length(valid.idx) != length(ROI)) {
    ## subset ROI using valid coordinates
    ROI <- ROI[valid.idx]
  }
  
  
  parcels <- if (!is.null(dataset$parcellation)) {
    stop("parcellation not supported")
    ## subset parcellation with ROI
    dataset$parcellation[ROI]
  }
  
  if (is.null(tuneGrid)) {
    ## should move within crossval_
    tuneGrid <- model$grid(values(ROI), Y, 1)
  }
  
  
  result <- if (is.null(dataset$testVec)) {
    ### train and test on one set.
    crossval_internal(crossval, dataset$Y, ROI, subIndices, model, tuneGrid, 
                      featureSelector = featureSelector, parcels = parcels)
  } else {
    ### train on one set, test on the other.
    testROI <- dataset$testChunk(indices(ROI))
    crossval_external(crossval, dataset$Y, ROI, dataset$testY, testROI, subIndices, model, tuneGrid, 
                      featureSelector = featureSelector, parcels = parcels)
  }
  
  ## attr(result, "valid.idx") <- valid.idx
  
  result
}

# helper method to convert result list to a list of BrainVolumes 
.convertResultsToVolumeList <- function(validRes, mask) {
  vols <- matrixToVolumeList(validRes[,1:3], validRes[,4:ncol(validRes)], mask)
  names(vols) <- colnames(validRes)[4:ncol(validRes)]
  vols
}

.extractValidResults <- function(results) {
  invalid <- sapply(results, function(x) inherits(x, "simpleError") || is.null(x))
  validRes <- results[!invalid]
  
  if (length(validRes) == 0) {
    print(results)
    stop("no valid results, all tasks failed")
  }
  
  validRes
  
}


.doClustered <- function(dataset, model, nclusters, crossval, class_metrics=FALSE, ncores=1) {
  
  iterlist <- lapply(nclusters, function(nc) neuroim:::ClusteredSearchlight(dataset$mask, nc))
  index_mat <- do.call(cbind, lapply(iterlist, function(it) it$clusters))
  
  
  resultSet <- foreach(searchIter = iterlist) %dopar% {
    message("running clustered iterator")
    res <- lapply(searchIter, function(vox) {
      message("num vox: ", nrow(vox))
      roi <- dataset$trainChunk(vox)
      result <- try(model$run(dataset, roi, crossval))
      result$predictor <- NULL
      result
    })
    
    res
  }
  
  
  perfList <- lapply(1:nrow(index_mat), function(i) {
    ind <- index_mat[i,]
    
    rlist <- mapply(function(result, indices) {
      result[indices]
    }, resultSet, ind)
    
    result <- Reduce(merge_results, rlist)
    result$predicted <- predicted_class(result$probs)
    perf <- computePerformance(result, dataset$testSplits, class_metrics, model$customPerformance)
  })
  
  
  validRes <- .extractValidResults(perfList)
  perfmat <- do.call(rbind, validRes)
  
  ## fix only use valid indices
  dataset$convertScores(which(dataset$mask != 0), perfmat)
  
}

#' @importFrom tibble data_frame
#' @importFrom magrittr %>%
#' @importFrom dplyr select mutate rowwise arrange do
.doRandomized2 <- function(dataset, model, radius, crossval, niter, class_metrics = FALSE) {
  iterlist <- replicate(niter, dataset$searchlight(radius, "randomized"), simplify = FALSE)
  
  resultSet <- foreach(searchIter=iterlist) %:% foreach(vox=searchIter, .combine=rbind, .packages=c("dplyr", "tibble")) %do% { 
    roi <- dataset$trainChunk(vox)
    if (length(roi) > 1) {
      result <- try(model$run(dataset, roi, crossval))
      result$predictor <- NULL
      tibble::data_frame(result = list(result),
                         indices = list(attr(vox, "indices")))
    } else {
      tibble::data_frame(result = NA,
                         indices = list(attr(vox, "indices")))
    }
  }
  
  index_mat <- as.matrix(do.call(cbind, lapply(resultSet, function(result) {
    result %>% mutate(rn=row_number()) %>% dplyr::rowwise() %>% 
      do(tibble::data_frame(indices=.$indices, rn=.$rn)) %>% arrange(indices) %>% select(rn)
  })))
  
  
  perfList <- lapply(1:nrow(index_mat), function(i) {
    
    ind <- index_mat[i, ]
    
    rlist <- mapply(function(result, i) {
      result$result[i]
    }, resultSet, ind)
    
    rlist <- rlist[sapply(rlist, length) != 1]
    
    result <- Reduce(merge_results, rlist)
    result$predicted <- predicted_class(result$probs)
    
    perf <-
      computePerformance(result,
                         dataset$testSplits,
                         class_metrics,
                         model$customPerformance)
  })
  
  validRes <- .extractValidResults(perfList)
  perfmat <- do.call(rbind, validRes)
  
  dataset$convertScores(which(dataset$mask != 0), perfmat)
  
}

# standard searchlight
.doStandard <- function(dataset, model, radius, crossval, class_metrics=FALSE) {
  searchIter <- itertools::ihasNext(dataset$searchlight(radius, "standard")) 
  
  res <- foreach::foreach(vox = searchIter, .verbose=FALSE, .packages=c("rMVPA", "MASS", "neuroim", "caret", model$model$library)) %dopar% {   
    roi <- dataset$trainChunk(vox)
    if (length(roi) > 1) {
      result <- try(model$run(dataset, roi, crossval))
      perf <- computePerformance(result, dataset$testSplits, class_metrics, model$customPerformance)
      center <- attr(vox, "center.index")
      list(perf=perf, center=center)
    } 
  }
  
  validRes <- .extractValidResults(res)
  perfmat <- do.call(rbind, lapply(validRes, "[[", "perf"))
  ids <- sapply(validRes, "[[", "center")
  dataset$convertScores(ids, perfmat)
}



# randomized searchlight
.doRandomized <- function(dataset, model, radius, crossval, class_metrics=FALSE) {
  searchIter <- itertools::ihasNext(dataset$searchlight(radius, "randomized"))
  
  ## tight inner loop should probably avoid "foreach" as it has a lot of overhead, but c'est la vie for now.
  res <- foreach::foreach(vox = searchIter, .verbose=FALSE, .errorhandling="pass", .packages=c("rMVPA", "MASS", "caret", "neuroim", model$model$library)) %do% {   
    roi <- dataset$trainChunk(vox)
    if (length(roi) > 1) {  
      print(length(roi))
      result <- try(model$run(dataset, roi, crossval))
      perf <- computePerformance(result, dataset$testSplits, class_metrics, model$customPerformance)
      perf <- perf[rep(1, length(roi)),]
      list(perf=perf, vox=indices(roi))
    } else {
      list(perf=NA, vox=indices(roi))
    }
  }
  
  validRes <- .extractValidResults(res)
  perfmat <- do.call(rbind, lapply(validRes, "[[", "perf"))
  ids <- unlist(lapply(validRes, "[[", "vox"))
  dataset$convertScores(ids, perfmat)
  
}


## mvpa_hiearchical ??


#' mvpa_regional
#' 
#' Run a separate MVPA analysis for multiple disjoint regions of interest.
#' 
#' @param dataset a \code{MVPADataset} instance.
#' @param model a \code{BaseModel} instance usually of type \code{CaretModelWrapper}
#' @param region_mask a \code{BrainVolume} where each region is identified by a unique integer. Every non-zero set of positive integers will be used to define a set of voxels for clasisifcation analysis.
#' @param crossval
#' @param savePredictors whether to return prediction model in result (default is \code{FALSE} to save memory)
#' @param featureSelector an option \code{FeatureSelector} object that is used to subselect informative voxels in feature space.
#' @param class_metrics
#' 
#' @return a named list of \code{BrainVolume} objects, where each name indicates the performance metric and label (e.g. accuracy, AUC)
#' @import itertools 
#' @import foreach
#' @import doParallel
#' @import parallel
#' @export
mvpa_regional <- function(dataset, model, region_mask, crossval=kfold_cross_validation(length(dataset$Y)), 
                          savePredictors=FALSE, featureSelector=NULL, class_metrics=FALSE) {  
  
  if (!is.null(dataset$parcellation)) {
    if (! all(dim(dataset$parcellation) == dim(region_mask))) {
      stop("mvpa_regional: dimension of 'featureParcellation' must equal dimensions of 'region_mask'")
    }
  }
  
  ## Get the set of unique ROIs (all unique integers > 0 in provided mask)
  regionSet <- sort(as.integer(unique(region_mask[region_mask > 0])))
  if (length(regionSet) < 1) {
    stop("mvpa_regional: invalid ROI mask, number of ROIs = 0")
  }
  
  res <- foreach::foreach(roinum = regionSet, .verbose=FALSE, .errorhandling="pass", .packages=c("rMVPA", "MASS", "neuroim", "caret", dataset$model$library)) %dopar% {   
    idx <- which(region_mask == roinum)
    
    if (length(idx) > 1) {
      ## what if length is less than 1?
      
      roi <- dataset$trainChunk(idx)
      #vox <- ROIVolume(space(region_mask), indexToGrid(region_mask, idx))
      
      result <- model$run(dataset, roi, crossval, featureSelector)
      
      attr(result, "ROINUM") <- roinum
      attr(result, "vox") <- idx
      
      
      perf <- if (!is.null(model$customPerformance)) {
        standard <- performance(result, dataset$testSplits, class_metrics)
        custom_perf <- customPerformance(result, model$customPerformance, dataset$testSplits)
        c(standard, custom_perf)
      } else {
        t(performance(result, dataset$testSplits, class_metrics))[1,]
        
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
    fill(region_mask, cbind(perfMat[, 1], perfMat[,cnum]))    
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


#' mvpa_clustered_searchlight
#' 
#' @param dataset a \code{MVPADataset} instance.
#' @param model
#' @param crossval
#' @param nclusters a vector of integers indicating the number of clusters for each searchlight iteration.
#' @param class_metrics
#' @export
mvpa_clustered_searchlight <- function(dataset, model, crossval, nclusters = NULL, class_metrics=FALSE) {
  if (is.null(nclusters)) {
    nvox <- sum(dataset$mask != 0)
    maxc <- min(floor(nvox/10),nvox)
    minc <- min(ceiling(nvox/100), nvox)
    nclusters <- round(seq(minc, maxc, length.out=5))
  }
  
  .doClustered(dataset, model, nclusters, crossval, class_metrics=class_metrics) 
  
  
}


#' mvpa_searchlight
#' 
#' @param dataset a \code{MVPADataset} instance.
#' @param model
#' @param crossval
#' @param radius the searchlight radus in mm
#' @param method the type of searchlight (randomized, randomized2, or standard)
#' @param niter the number of searchlight iterations for 'randomized' method
#' @param class_metrics
#' @return a named list of \code{BrainVolume} objects, where each name indicates the performance metric and label (e.g. accuracy, AUC)
#' @import itertools 
#' @import foreach
#' @import doParallel
#' @import parallel
#' @importFrom futile.logger flog.info
#' @export
mvpa_searchlight <- function(dataset, model, crossval, radius=8, method=c("randomized", "randomized2", "standard"),  
                             niter=4, class_metrics=FALSE) {
  stopifnot(niter > 1)
  
  if (radius < 1 || radius > 100) {
    stop(paste("radius", radius, "outside allowable range (1-100)"))
  }
  
  method <- match.arg(method)
  
  flog.info("model is: %s", model$model_name)
  
  res <- if (method == "standard") {
    .doStandard(dataset, model, radius, crossval, class_metrics=class_metrics)    
  } else if (method == "randomized") {
    
    res <- foreach(i = 1:niter) %dopar% {
      flog.info("Running randomized searchlight iteration %s", i)   
      .doRandomized(dataset, model, radius, crossval, class_metrics=class_metrics)
    }
    
    dataset$averageOverIterations(res)
  } else if (method == "randomized2") {
    .doRandomized2(dataset, model, radius, crossval, niter, class_metrics=class_metrics)  
  }
  
}

# average_over_iterations <- function(dataset, result) {
#   Xall <- lapply(1:ncol(result[[1]]), function(i) {
#     X <- do.call(cbind, lapply(result, function(M) M[,i]))
#     xmean <- rowMeans(X, na.rm=TRUE)
#     xmean[is.na(xmean)] <- 0
#     BrainVolume(xmean, space(dataset$mask))
#   })
#   
#   names(Xall) <- colnames(result[[1]])
#   Xall
#   
# }





