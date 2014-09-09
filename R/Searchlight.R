
.doStandard <- function(model, dset, mask, radius, ncores) {
  searchIter <- itertools::ihasNext(Searchlight(dset$mask, radius)) 
  foreach(vox = searchIter, .combine=rbind, .verbose=TRUE) %do% {   
    if (nrow(vox) < 3) {
      NA
    } else {
      result <- fitModel(model, dset, vox, ncores)    
      cen <- attr(vox, "center")
      c(cen, performance(result))  
    }
  }
  
}
  

.doRandomized <- function(model, dset, mask, radius, ncores) {
  searchIter <- itertools::ihasNext(RandomSearchlight(mask, radius))
  res <- do.call(rbind, foreach(vox = searchIter, .verbose=TRUE, .errorhandling="pass") %do% {   
    if (nrow(vox) < 3) {
      NULL
    } else {
      result <- t(performance(fitModel(model, dset, vox, ncores)))
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
#' @param dataset
#' @param radius the searchlight radus in mm
#' @param modelName the name of the classifcation model to be used
#' @param ncores the number of cores for parallel processign (default is 1)
#' @return a named list of \code{BrainVolume} objects, where each name indicates the performance metric and label (e.g. accuracy, AUC)
#' @importFrom itertools ihasNext
#' @export
searchlight <- function(dset, radius=8, modelName="svmLinear", ncores=1, method=c("randomized", "standard"), niter=4) {
  if (radius < 1 || radius > 100) {
    stop(paste("radius", radius, "outside allowable range (1-100)"))
  }
  
  if (length(dset$blockVar) != length(dset$Y)) {
    stop(paste("length of 'labels' must equal length of 'cross validation blocks'", length(dset$Y), "!=", length(dset$blockVar)))
  }
  
  model <- loadModel(modelName)
  method <- match.arg(method)
  
  
  print(model)
  res <- if (method == "standard") {
    .doStandard(model, dset, dset$mask, radius, ncores)    
  } else {
    res <- lapply(1:niter, function(i) {
      do.call(cbind, .doRandomized(model, dset, dset$mask, radius, ncores))
    })
   
    Xall <- lapply(1:ncol(res[[1]]), function(i) {
      X <- do.call(cbind, lapply(res, function(M) M[,i]))
      xmean <- rowMeans(X, na.rm=TRUE)
      xmean[is.na(xmean)] <- 0
      BrainVolume(xmean, space(dset$mask))
    })
    
    names(Xall) <- colnames(res[[1]])
    Xall
    
  }
  
  #browser()
  
  #grid <- res[,1:3]
  
  #bvout <- lapply(4:ncol(res), function(i) {  
  #  bv <- BrainVolume(numeric(length(dset$mask)), space(dset$mask))
  #  bv[grid] <- res[,i]
  #  bv
  #})
  
  #names(bvout) <- colnames(res[,4:ncol(res)])
  #bvout   
}