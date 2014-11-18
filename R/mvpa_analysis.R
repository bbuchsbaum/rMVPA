

#' @export  
matrixToVolumeList <- function(vox, mat, mask, default=NA) {
  lapply(1:ncol(mat), function(i) {
    vol <- array(default, dim(mask))   
    vol[vox] <- mat[,i]
    BrainVolume(vol, space(mask))
  })
}  

.doStandard <- function(model, bvec, Y, blockVar, mask, radius, ncores) {
  searchIter <- itertools::ihasNext(Searchlight(mask, radius)) 
  
  res <- foreach::foreach(vox = searchIter, .combine=rbind, .verbose=FALSE) %dopar% {   
    if (nrow(vox) < 3) {
      NA
    } else {
      result <- fitMVPAModel(model, bvec, Y, blockVar, vox, ncores, fast=TRUE, finalFit=FALSE)    
      cen <- attr(vox, "center")
      c(cen, performance(result))  
    }
  }
  
  vols <- matrixToVolumeList(res[,1:3], res[4:ncol(res)], mask)
  names(vols) <- colnames(res)[4:ncol(res)]
  vols
}
  

.doRandomized <- function(model, bvec, Y, blockVar, mask, radius=8, ncores=1, tuneGrid=NULL) {
  searchIter <- itertools::ihasNext(RandomSearchlight(mask, radius))
  
  res <- foreach::foreach(vox = searchIter, .verbose=FALSE, .combine=rbind, .errorhandling="pass", .packages=c("rMVPA", model$library)) %do% {   
    if (nrow(vox) < 3) {
      NULL
    } else {     
      fit <- fitMVPAModel(model, bvec, Y, blockVar, vox, ncores, tuneGrid=tuneGrid, fast=TRUE, finalFit=FALSE)
      result <- t(performance(fit))
      out <- cbind(vox, result[rep(1, nrow(vox)),])
      #attr(out, "prob") <- fit$probs
      out
    }
  }
  
  #print(res)
   
  vols <- matrixToVolumeList(res[,1:3], res[,4:ncol(res)], mask)
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
#' @param bvec a \code{BrainVector} instance, a 4-dimensional image where the first three dimensions are space (x,y,z) and the 4th dimension is the image/scan/condition.
#' @param Y the dependent variable. If it is a factor, then classification analysis is performed. If it is a continuous variable then regression is performed.
#' @param mask a \code{BrainVolume} instance indicating the inclusion mask for voxels entering the searchlight analysis. 
#' @param blockVar an \code{integer} vector indicating the blocks to be used for cross-validation. This is usually a variable indicating the scanning "run". 
#'        Must be same length as \code{Y}
#' @param modelName the name of the classifcation model to be used
#' @param ncores the number of cores for parallel processign (default is 1)
#' @return a named list of \code{BrainVolume} objects, where each name indicates the performance metric and label (e.g. accuracy, AUC)
#' @import itertools 
#' @import foreach
#' @import doParallel
#' @import parallel
#' @export
mvpa_regional <- function(bvec, Y, mask, blockVar, modelName="corsim", ncores=2, tuneGrid=NULL) {
  if (length(blockVar) != length(Y)) {
    stop(paste("length of 'labels' must equal length of 'cross validation blocks'", length(Y), "!=", length(blockVar)))
  }
  
  regionSet <- sort(unique(mask[mask > 0]))
  model <- loadModel(modelName)
  cl <- makeCluster(ncores, outfile="")
  registerDoParallel(cl)
  
  res <- foreach::foreach(roinum = regionSet, .verbose=TRUE, .errorhandling="pass", .packages=c("rMVPA", "MASS", "neuroim", model$library)) %dopar% {   
    idx <- which(mask == roinum)
    if (length(idx) < 2) {
      NULL
    } else {
      vox <- indexToGrid(mask, idx)
      fit <- fitMVPAModel(model, bvec, Y, blockVar, vox, fast=TRUE, finalFit=TRUE, tuneGrid=tuneGrid)
      result <- c(ROINUM=roinum, t(performance(fit))[1,])
      
    }
  }
  
  invalid <- sapply(res, function(x) inherits(x, "simpleError") || is.null(x))
  validRes <- res[!invalid]
  
  perfMat <- do.call(rbind, validRes)
  
  outVols <- lapply(2:ncol(perfMat), function(cnum) {
     fill(mask, cbind(perfMat[, 1], perfMat[,cnum]))    
  })
  
  names(outVols) <- colnames(perfMat)[2:ncol(perfMat)]
  list(outVols = outVols, performance=perfMat)

}
  
  
  
#' mvpa_searchlight
#' @param bvec a \code{BrainVector} instance, a 4-dimensional image where the first three dimensons are (x,y,z) and the 4th dimension is the dependent class/variable
#' @param Y the dependent variable. If it is a factor, then classification analysis is performed. If it is a continuous variable then regression is performed.
#' @param mask a \code{BrainVolume} instance indicating the inclusion mask for voxels entering the searchlight analysis. 
#' @param blockVar an \code{integer} vector indicating the blocks to be used for cross-validation. This is usually a variable indicating the scanning "run". 
#'        Must be same length as \code{Y}
#' @param radius the searchlight radus in mm
#' @param modelName the name of the classifcation model to be used
#' @param ncores the number of cores for parallel processign (default is 1)
#' @return a named list of \code{BrainVolume} objects, where each name indicates the performance metric and label (e.g. accuracy, AUC)
#' @import itertools 
#' @import foreach
#' @import doParallel
#' @import parallel
#' @import futile.logger
#' @export
mvpa_searchlight <- function(bvec, Y, mask, blockVar, radius=8, modelName="svmLinear", ncores=2, method=c("randomized", "standard"), niter=4, tuneGrid=NULL) {
  if (radius < 1 || radius > 100) {
    stop(paste("radius", radius, "outside allowable range (1-100)"))
  }
  
  if (length(blockVar) != length(Y)) {
    stop(paste("length of 'labels' must equal length of 'cross validation blocks'", length(Y), "!=", length(blockVar)))
  }
  
 
  cl <- makeCluster(ncores)
  registerDoParallel(cl)
  
  
  model <- loadModel(modelName)
  method <- match.arg(method)
  
  
  flog.info("classification model is: %s", modelName)
  flog.info("tuning grid is", tuneGrid, capture=TRUE)
  
  res <- if (method == "standard") {
    .doStandard(model, bvec, Y, blockVar, mask, radius, ncores)    
  } else {
    res <- parallel::mclapply(1:niter, function(i) {
      flog.info("Running randomized searchlight iteration %s", i)   
      do.call(cbind, .doRandomized(model,bvec, Y, blockVar, mask, radius, ncores, tuneGrid) )
    }, mc.cores=ncores)
   
    Xall <- lapply(1:ncol(res[[1]]), function(i) {
      X <- do.call(cbind, lapply(res, function(M) M[,i]))
      xmean <- rowMeans(X, na.rm=TRUE)
      xmean[is.na(xmean)] <- 0
      BrainVolume(xmean, space(mask))
    })
    
    names(Xall) <- colnames(res[[1]])
    Xall
    
  }
  
}