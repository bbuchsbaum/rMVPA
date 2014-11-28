

#' @export  
matrixToVolumeList <- function(vox, mat, mask, default=NA) {
  lapply(1:ncol(mat), function(i) {
    vol <- array(default, dim(mask))   
    vol[vox] <- mat[,i]
    BrainVolume(vol, space(mask))
  })
}  

.doStandard <- function(dataset, radius, ncores) {
  searchIter <- itertools::ihasNext(Searchlight(dataset$mask, radius)) 
  
  res <- foreach::foreach(vox = searchIter, .combine=rbind, .verbose=FALSE) %dopar% {   
    if (nrow(vox) < 3) {
      NA
    } else {
      result <- fitMVPAModel(dataset, vox, fast=TRUE, finalFit=FALSE)    
      cen <- attr(vox, "center")
      c(cen, performance(result))  
    }
  }
  
  vols <- matrixToVolumeList(res[,1:3], res[4:ncol(res)], dataset$mask)
  names(vols) <- colnames(res)[4:ncol(res)]
  vols
}
  

.doRandomized <- function(dataset, radius) {
  searchIter <- itertools::ihasNext(RandomSearchlight(dataset$mask, radius))
  
  res <- foreach::foreach(vox = searchIter, .verbose=FALSE, .combine=rbind, .errorhandling="pass", .packages=c("rMVPA", dataset$model$library)) %do% {   
    if (nrow(vox) < 3) {
      NULL
    } else {     
      fit <- fitMVPAModel(dataset, vox, fast=TRUE, finalFit=FALSE)
      result <- t(performance(fit))
      out <- cbind(vox, result[rep(1, nrow(vox)),])
      #attr(out, "prob") <- fit$probs
      out
    }
  }
  
  #print(res)
   
  vols <- matrixToVolumeList(res[,1:3], res[,4:ncol(res)], dataset$mask)
  names(vols) <- colnames(res)[4:ncol(res)]
  vols
  
  
}

#.doRegional <- function(regionSet, model) {
#  res <- foreach::foreach(roinum = regionSet, .verbose=TRUE, .errorhandling="pass", .packages=c("rMVPA", "MASS", "neuroim", model$library)) %dopar% {   
#    idx <- which(mask == roinum)
#    if (length(idx) < 2) {
#      NULL
#    } else {
#      vox <- indexToGrid(mask, idx)
#      fit <- fitMVPAModel(model, bvec, Y, blockVar, vox, fast=TRUE, finalFit=TRUE, tuneGrid=tuneGrid)
#      result <- c(ROINUM=roinum, t(performance(fit))[1,])    
#      attr(result, "finalFit") <- fit
#      result
#    }
#  }
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
mvpa_regional <- function(dataset, regionMask, ncores=1) {
  if (length(dataset$blockVar) != length(dataset$Y)) {
    stop(paste("length of 'labels' must equal length of 'cross validation blocks'", length(Y), "!=", length(blockVar)))
  }
  
  cl <- makeCluster(ncores)
  registerDoParallel(cl)
  
  regionSet <- sort(as.integer(unique(regionMask[regionMask > 0])))
  
  if (ncores > 1 && length(regionSet) == 1) {
    mc.cores <- ncores
  } else {
    mc.cores <- 1
  }
  
  #allowParallel <- if (length(regionSet) == 1) TRUE else FALSE
  
  res <- foreach::foreach(roinum = regionSet, .verbose=TRUE, .errorhandling="pass", .packages=c("rMVPA", "MASS", "neuroim", dataset$model$library)) %do% {   
    idx <- which(regionMask == roinum)
    if (length(idx) < 2) {
      NULL
    } else {
      vox <- indexToGrid(regionMask, idx)
      fit <- fitMVPAModel(dataset, vox, fast=TRUE, finalFit=TRUE, ncores=mc.cores)
      result <- c(ROINUM=roinum, t(performance(fit))[1,]) 
      
      predictor <- asPredictor(fit$finalFit, vox)
      attr(result, "predictor") <- predictor
      result
    }
  }
  
  invalid <- sapply(res, function(x) inherits(x, "simpleError") || is.null(x))
  validRes <- res[!invalid]
  
  if (length(validRes) == 0) {
    flog.error("Classification failed for all of %s ROIs", length(regionSet))
    stop("error in mvpa_regional: aborting")
  } 
  
  perfMat <- do.call(rbind, validRes)
  
  outVols <- lapply(2:ncol(perfMat), function(cnum) {
    fill(regionMask, cbind(perfMat[, 1], perfMat[,cnum]))    
  })
  
  predictorList <- lapply(validRes, function(res) {
    attr(res, "predictor")
  })
  
  names(predictorList) <- regionSet[!invalid]
  names(outVols) <- colnames(perfMat)[2:ncol(perfMat)]
  list(outVols = outVols, performance=perfMat, predictorList=predictorList)

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