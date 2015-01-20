
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
      gmeans <- mefa4::groupMeans(fold$Xtest, 1, fold$Ytest)
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
    
    crossrun <- apply(runCoords, 1, function(x) x[1] != x[2])
    cwithin <- apply(labelCoords, 1, function(x) x[1] == x[2])
    cbetween <- apply(labelCoords, 1, function(x) x[1] != x[2])
   
    cor.between <- mean(corMat[crossrun & cbetween])
    cor.within <- mean(corMat[crossrun & cwithin])
    
    SimilarityResult(cor.within, cor.between, corMat)
  }
  
}

