


#' patSimHelper <- memoise::memoise(function(runs, labels, cdim) {
#'   lowerTri <- as.vector(lower.tri(matrix(0, cdim[1], cdim[2])))
#'   
#'   runCoords <- do.call(rbind, lapply(runs, function(r1) {
#'     t(sapply(runs, function(r2) {
#'       c(r2,r1)
#'     }))
#'   }))
#'   
#'   labelCoords <- do.call(rbind, lapply(labels, function(l1) {
#'     t(sapply(labels, function(l2) {
#'       c(l1,l2)
#'     }))
#'   }))
#'   
#'   crossFac <- factor(apply(runCoords, 1, paste0, collapse=":")) 
#'   crossrun <- runCoords[,1] != runCoords[,2] & lowerTri
#'   cwithin <- labelCoords[,1] == labelCoords[,2] & lowerTri
#'   cbetween <- labelCoords[,1] != labelCoords[,2] & lowerTri
#'    
#'   validWithin <- crossrun & cwithin
#'   validLevels <- levels(factor(crossFac[crossrun & cwithin]))
#'   validBetween <- crossFac %in% validLevels & cbetween
#'   
#'   list(crossFac=crossFac, crossrun=crossrun, cwithin=cwithin, cbetween=cbetween, validWithin=validWithin, validBetween=validBetween)
#'   
#' })
#' 
#' #' @export
#' patternSimilarity <- function(dataset, vox, simFun=cor, contrastMatrix=NULL) {
#'   
#'   if (is.null(contrastMatrix)) {
#'     contrastMatrix <- diag(length(levels(dataset$Y)))
#'     colnames(contrastMatrix) <- levels(dataset$Y)
#'     utri <- upper.tri(contrastMatrix)
#'     ltri <- lower.tri(contrastMatrix)
#'     contrastMatrix[utri] <- length(levels(dataset$Y))/sum(utri)/-2
#'     contrastMatrix[ltri] <- length(levels(dataset$Y))/sum(ltri)/-2
#'   }
#'   
#'   assertthat::assert_that(sum(contrastMatrix) < 1e-6)
#'   
#'   ## TODO deals with parcels.
#'   
#'   X <- series(dataset$trainVec, vox)
#'   valid.idx <- nonzeroVarianceColumns(X)
#'   
#'   X <- X[,valid.idx,drop=FALSE]
#'   vox <- vox[valid.idx,]
#'   
#' 
#'   foldIterator <- MatrixFoldIterator(X, dataset$Y, vox, dataset$blockVar)
#'   
#'   if (ncol(X) == 0) {
#'     stop("no valid columns")
#'   } else {
#'     
#'     corMats <- lapply(foldIterator, function(fold) {
#'       gmeans_test <- groupMeans(fold$Xtest, 1, fold$Ytest)
#'       gmeans_train <- groupMeans(fold$Xtrain, 1, fold$Ytrain)
#'       cmat <- simFun(t(gmeans_test), t(gmeans_train))
#'     })
#'      
#'     corMat <- Reduce("+", corMats)/length(corMats)
#'     cons <- sapply(corMats, function(cmat) sum(cmat * contrastMatrix))
#'   
#'     SimilarityResult(corMat=corMat, avgContrast=mean(cons), sdContrast=sd(cons))
#'   }
#'   
#' }
