
patSimHelper <- memoise::memoise(function(runs, labels, cdim) {
  lowerTri <- as.vector(lower.tri(matrix(0, cdim[1], cdim[2])))
  runCoords <- do.call(rbind, lapply(runs, function(r1) {
    t(sapply(runs, function(r2) {
      c(r2,r1)
    }))
  }))
  
  labelCoords <- do.call(rbind, lapply(labels, function(l1) {
    t(sapply(labels, function(l2) {
      c(l1,l2)
    }))
  }))
  
  crossFac <- factor(apply(runCoords, 1, paste0, collapse=":")) 
  crossrun <- runCoords[,1] != runCoords[,2] & lowerTri
  cwithin <- labelCoords[,1] == labelCoords[,2] & lowerTri
  cbetween <- labelCoords[,1] != labelCoords[,2] & lowerTri
   
  validWithin <- crossrun & cwithin
  validLevels <- levels(factor(crossFac[crossrun & cwithin]))
  validBetween <- crossFac %in% validLevels & cbetween
  
  list(crossFac=crossFac, crossrun=crossrun, cwithin=cwithin, cbetween=cbetween, validWithin=validWithin, validBetween=validBetween)
  
})

#' @export
#' @import mefa4
patternSimilarity <- function(dataset, vox, simFun=cor) {
  X <- series(dataset$trainVec, vox)
  valid.idx <- nonzeroVarianceColumns(X)
  
  X <- X[,valid.idx]
  vox <- vox[valid.idx,]
  
  foldIterator <- MatrixFoldIterator(X, dataset$Y, dataset$blockVar)
  
  if (ncol(X) == 0) {
    stop("no valid columns")
  } else {
    centroids <- lapply(foldIterator, function(fold) {
      gmeans = mefa4::groupMeans(fold$Xtest, 1, fold$Ytest)
      run <- rep(fold$index, nrow(gmeans))
      list(gmeans=gmeans, run=run, labels=row.names(gmeans))
    })
     
    runs <- unlist(lapply(centroids, "[[", "run"))
    labels <- unlist(lapply(centroids, "[[", "labels"))
    gmeans <- do.call(rbind, lapply(centroids, "[[", "gmeans"))
    
    corMat <- simFun(t(gmeans))   
    #lowerTri <- as.vector(lower.tri(corMat))
    
    ret <- patSimHelper(runs,labels, dim(corMat))
          
    withinMeans <- aggregate(corMat[ret$validWithin] ~ ret$crossFac[ret$validWithin], FUN=mean)[,2]    
    betweenMeans <- aggregate(corMat[ret$validBetween] ~ ret$crossFac[ret$validBetween], FUN=mean)[,2]
   
    cor.between <- mean(betweenMeans)
    cor.within <- mean(withinMeans)
    
    simWithinTable <- data.frame(label=names(which(ret$validWithin)), sim=corMat[ret$validWithin])
    SimilarityResult(cor.within, cor.between, corMat, simWithinTable)
  }
  
}

