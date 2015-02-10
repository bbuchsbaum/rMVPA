
patSimHelper <- function(runs, labels) {
  
}

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
      mefa4::groupMeans(fold$Xtest, 1, fold$Ytest)
      run <- rep(fold$index, nrow(gmeans))
      list(gmeans=gmeans, run=run, labels=row.names(gmeans))
    })
     
    runs <- unlist(lapply(centroids, "[[", "run"))
    labels <- unlist(lapply(centroids, "[[", "labels"))
    gmeans <- do.call(rbind, lapply(centroids, "[[", "gmeans"))
    
    corMat <- simFun(t(gmeans))
    
    
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
    
    lowerTri <- as.vector(lower.tri(corMat))
    crossrun <- runCoords[,1] != runCoords[,2] & lowerTri
    cwithin <- labelCoords[,1] == labelCoords[,2] & lowerTri
    cbetween <- labelCoords[,1] != labelCoords[,2] & lowerTri
    
    
    
    crossFac <- factor(apply(runCoords, 1, paste0, collapse=":")) 
    
    validWithin <- crossrun & cwithin
    validLevels <- levels(factor(crossFac[crossrun & cwithin]))
    validBetween <- crossFac %in% validLevels & cbetween
      
        
    withinMeans <- aggregate(corMat[validWithin] ~ crossFac[validWithin], FUN=mean)[,2]    
    betweenMeans <- aggregate(corMat[validBetween] ~ crossFac[validBetween], FUN=mean)[,2]
   
    cor.between <- mean(betweenMeans)
    cor.within <- mean(withinMeans)
    
    simWithinTable <- data.frame(label=names(which(validWithin)), sim=corMat[validWithin])
    SimilarityResult(cor.within, cor.between, corMat, corTable, simWithinTable)
  }
  
}

