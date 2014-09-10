

  
  
.doStandard <- function(model, bvec, Y, blockVar, mask, radius, ncores) {
  searchIter <- itertools::ihasNext(Searchlight(mask, radius)) 
  foreach::foreach(vox = searchIter, .combine=rbind, .verbose=FALSE) %dopar% {   
    if (nrow(vox) < 3) {
      NA
    } else {
      result <- fitModel(model, bvec, Y, blockVar, vox, ncores)    
      cen <- attr(vox, "center")
      c(cen, performance(result))  
    }
  }
  
}
  

.doRandomized <- function(model, bvec, Y, blockVar, mask, radius, ncores) {
  searchIter <- itertools::ihasNext(RandomSearchlight(mask, radius))
  res <- do.call(rbind, foreach::foreach(vox = searchIter, .verbose=FALSE, .errorhandling="stop", .packages=c("rMVPA", model$library) %dopar% {   
    if (nrow(vox) < 3) {
      NULL
    } else {
      result <- t(performance(fitModel(model, bvec, Y, blockVar, vox, ncores)))
      out <- cbind(vox, result[rep(1, nrow(vox)),])
      cen <- attr(out, "center")
      out
    }
  })
  

  vols <- lapply(4:ncol(res), function(i) {
    vol <- array(NA, dim(mask))
    vol[res[,1:3]] <- res[,i]
    vol
  })
  
  names(vols) <- colnames(res)[4:ncol(res)]
  vols
  
  
}
  



#' searchlight
#' @param bvec a \code{BrainVector} instance, a 4-dimensional image where the first three dimensons are apce (x,y,z) and the 4th dimension is the dependent class/variable
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
#' @export
searchlight <- function(bvec, Y, mask, blockVar, radius=8, modelName="svmLinear", ncores=2, method=c("randomized", "standard"), niter=4) {
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
  
  
  print(model)
  res <- if (method == "standard") {
    .doStandard(model, bvec, Y, blockVar, mask, radius, ncores)    
  } else {
    res <- lapply(1:niter, function(i) {
      do.call(cbind, .doRandomized(model,bvec, Y, blockVar, mask, radius, ncores) )
    })
   
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