
#' BlockedCrossValidation
#' 
#' construct a cross-validation specification
#' 
#' @param blockVar an integer vector of indicating the cross-vlaidation blocks. Each block is indicating by an unique integer.
#' @param balance logical indicating whether cross-validation blocks should automatically balanced, using undersampling if necessary.
#' @param bootstrap logical indicating whether training samples should be sampled with replacement.
#' @export
BlockedCrossValidation <- function(blockVar, balance=FALSE, bootstrap=FALSE) {
  ret <- list(blockVar=blockVar, balance=balance, bootstrap=bootstrap, nfolds=length(unique(blockVar)))
  class(ret) <- c("BlockedCrossValidation", "CrossValidation", "list")
  ret
}

#' KFoldCrossValidation
#' 
#' @export
#' @param len the number of observations
#' @param balance 
#' @param boostrap
KFoldCrossValidation <- function(len, nfolds=10, balance=FALSE, bootstrap=FALSE) {
  blockVar <- sample(rep(seq(1, nfolds), length.out=len))
  ret <- list(blockVar=blockVar, balance=balance, bootstrap=bootstrap, nfolds=nfolds)
  class(ret) <- c("KFoldCrossValidation", "CrossValidation", "list")
  ret
}

#' foldIter
#' 
#' @export
foldIter <- function(obj, Y, subIndices) {
  UseMethod("foldIter")
}

#' matrixIter
#' 
#' @export
matrixIter <- function(obj, Y, X, vox, subIndices) {
  UseMethod("matrixIter")
}


#' @export
foldIter.CrossValidation <- function(obj, Y, subIndices=NULL) {
  if (!is.null(subIndices)) {
    FoldIterator(Y[subIndices], obj$blockVar[subIndices], balance=obj$balance, bootstrap=obj$bootstrap)
  } else {
    FoldIterator(Y, obj$blockVar, balance=obj$balance, bootstrap=obj$bootstrap)
  }
}

#' @export
matrixIter.CrossValidation <- function(obj, Y, X, vox, subIndices=NULL) {
  if (!is.null(subIndices)) {
    MatrixFoldIterator(X[subIndices,,drop=FALSE], Y[subIndices], vox, obj$blockVar[subIndices], balance=obj$balance, bootstrap=obj$bootstrap )
  } else {
    MatrixFoldIterator(X, Y, vox, obj$blockVar, balance=obj$balance, bootstrap=obj$bootstrap )
  }
}

#' invertFolds
#' @param foldSplit
#' @param allInd
#' @export
invertFolds <- function(foldSplit, allInd) { 
  index <-  relist(rank(unlist(foldSplit)), foldSplit)
  
  lapply(index, function(split) {
    allInd[-split]
  })
}

#' FoldIterator
#' 
#' Create an iterator from a block variable for iterating through the folds during cross-validation
#' @param Y the labels
#' @param blockVar variable denoting the cross-validation folds
#' @param balance try to balance each training sample so that the frequency of labels is equal across groups
#' @param bootstrap use bootstrap resampling of the training set
#' @param bootstrapMin 
#' @importFrom assertthat assert_that
#' @export
FoldIterator <- function(Y, blockVar,  balance=FALSE, bootstrap=FALSE, bootstrapMin=2) {

  assert_that(length(blockVar) == length(Y))
  assert_that(length(unique(blockVar)) > 1)
  
  index <- 0
  ord <- NULL
  
  testSets <- split(1:length(blockVar), blockVar)
  
  trainSets <- invertFolds(testSets, 1:length(blockVar))
  
  .getTrainSets <- function() {
    trainSets
  }
  
  .getTestSets <- function() {
    testSets
  }
  
  .getIndex <- function() {
    index
  }
  
  .getTestOrder <- function() {
    if (is.null(ord))
      ord <<- order(unlist(testSets))
    ord
  }
  
  .reset <- function() {
    index <<- 0
  }
  
  nextEl <- function() {
    
    if (index < length(trainSets)) { 
      index <<- index + 1    
      trainIndex <- trainSets[[index]]
      testIndex <- testSets[[index]]
       
      if (bootstrap) {
        
        ## make sure that bootstrap samples have at least bootstrapMin instance of each class.
        ## bit could be infinite loop....
        for (i in 1:50) {
          ind <- sort(sample(trainIndex, replace=TRUE))
          stab <- table(Y[ind])
          minClass <- min(stab)
          
          if (length(stab) == length(levels(Y)) && minClass >= bootstrapMin) {
            trainIndex <- ind
            break
          }
          
          if (i == 50) {
            stop("error in bootstrap sampling: after 50 attempts could not find bootstrap sample with at least 'bootstrapMin' instances for every class.")
          }
        }
      }
      
      if (balance) {
        trainIndex <- caret::upSample(trainIndex, Y[trainIndex])[,1]
      }
      
      list(trainIndex=trainIndex, testIndex=testIndex, Ytrain=Y[trainIndex], Ytest=Y[testIndex], index=index)
      
    } else {
      stop('StopIteration')
    }
  }
  
  obj <- list(nextElem=nextEl, blockVar=blockVar, index=.getIndex, getTrainSets=.getTrainSets, 
              getTestSets=.getTestSets, getTestOrder=.getTestOrder, reset=.reset, balance=balance, bootstrap=bootstrap, nfolds=length(testSets))
  
  class(obj) <- c("FoldIterator", 'abstractiter', 'iter')
  obj
  
}
  
#' MatrixFoldIterator
#' 
#' Create an iterator of training and test data from a data matrix \code{X}, a factor \code{Y}, and an optional precomputed blocking variable \code{blockVar} used for determining cross-validation folds.
#' @param X the data matrix
#' @param Y the labels
#' @param vox the voxel coordinates
#' @param blockVar variable denoting the cross-validation folds.
#' @param nfolds the number of cross-validation folds: only relevant when blockVar is not supplied
#' @param balance try to balance each training sample so that the frequency of labels is equal across groups
#' @param bootstrap use bootstrap resampling of the training set
#' @param bootstrapMin the minumum number of clases per bootstrap iteration.
#' @export
MatrixFoldIterator <- function(X, Y, vox, blockVar, balance=FALSE, bootstrap=FALSE, bootstrapMin=2) {
  
  if (nrow(X) != length(Y)) {
    stop("X matrix must have same number of rows as Y variable")
  }
  
  
  foldIter = FoldIterator(Y, blockVar, balance=balance, bootstrap=bootstrap, bootstrapMin)

  nextEl <- function() {
    ret <- foldIter$nextElem()
    ret$Xtrain <- X[ret$trainIndex,]
    ret$Xtest <- X[ret$testIndex,]
    ret
  }
  
  fullSample <- function() {
    if (bootstrap) {
      bootIndices <- sample(1:nrow(X), nrow(X), replace=TRUE)
    }
  }
  

  obj <- list(X=X, Y=Y, vox=vox, blockVar=blockVar, nextElem=nextEl, index=foldIter$index, 
              getTrainSets=foldIter$getTrainSets, getTestSets=foldIter$getTestSets, 
              getTestOrder=foldIter$getTestOrder, reset=foldIter$reset, balance=foldIter$balance, bootstrap=foldIter$bootstrap, nfolds=foldIter$nfolds)
  
  class(obj) <- c("MatrixFoldIterator", 'abstractiter', 'iter')
  obj
  
}

#' NestedFoldIterator
#' 
#' @param Y
#' @param blockVar
#' @param balance
#' @param bootstrap
#' @param bootstrapMin
#' @export
NestedFoldIterator <- function(Y, blockVar,  balance=FALSE, bootstrap=FALSE, bootstrapMin=2) {
  blockids <- sort(unique(blockVar))
  
  if (length(Y) != length(blockVar)) {
    stop("Y must have same length as blockVar")
  }
  
  
  index <- 0
  
  blockNumber <- function() { index }
  
  nextEl <- function() {
    if (index < length(blockids)) {
      index <<- index + 1
      curBlock <- blockids[index]
      train.idx <- which(blockVar != curBlock)
      iter <- FoldIterator(Y[train.idx], blockVar=blockVar[train.idx], balance=balance, bootstrap=bootstrap, bootstrapMin=bootstrapMin)
      attr(iter, "train_indices") <- train.idx
      attr(iter, "test_indices") <- seq(1, length(Y))[-train.idx]
      iter
    } else {
      stop('StopIteration')
    }
  }
  
  obj <- list(nextElem=nextEl, blockNumber=blockNumber)
  class(obj) <- c("NestedFoldIterator", 'abstractiter', 'iter')
  obj
  
}

#' NestedMatrixIterator
#' 
#' @param X
#' @param Y
#' @param vox
#' @param blockVar
#' @param balance
#' @param bootstrap
#' @param bootstrapMin
#' @export
NestedMatrixIterator <- function(X, Y, vox, blockVar, balance=FALSE, bootstrap=FALSE, bootstrapMin=2) {
  blockids <- sort(unique(blockVar))
  
  if (nrow(X) != length(Y)) {
    stop("X matrix must have same number of rows as Y variable")
  }
  
  
  index <- 0
  
  blockNumber <- function() { index }
  
  nextEl <- function() {
    if (index < length(blockids)) {
      index <<- index + 1
      curBlock <- blockids[index]
      train.idx <- which(blockVar != curBlock)
      iter <- MatrixFoldIterator(X[train.idx,], Y[train.idx], vox, blockVar=blockVar[train.idx], balance=balance, bootstrap=bootstrap, bootstrapMin=bootstrapMin)
      
      attr(iter, "train_indices") <- train.idx
      attr(iter, "test_indices") <- seq(1, length(Y))[-train.idx]
      iter
    } else {
      stop('StopIteration')
    }
  }
  
  obj <- list(nextElem=nextEl, blockNumber=blockNumber)
  class(obj) <- c("NestedMatrixFoldIterator", 'abstractiter', 'iter')
  obj
  
}
  
  
  