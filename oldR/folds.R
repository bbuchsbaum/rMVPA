#' foldIter
#' 
#' @export
foldIter <- function(obj, Y, subIndices) {
  UseMethod("foldIter")
}

#' matrixIter
#' 
#' @export
matrixIter <- function(obj, Y, X, subIndices) {
  UseMethod("matrixIter")
}


#' @export
foldIter.CrossValidation <- function(obj, Y, subIndices=NULL) {
  if (!is.null(subIndices)) {
    FoldIterator(Y[subIndices], obj$block_var[subIndices], balance=obj$balance, bootstrap=obj$bootstrap)
  } else {
    FoldIterator(Y, obj$block_var, balance=obj$balance, bootstrap=obj$bootstrap)
  }
}

#' @export
matrixIter.CrossValidation <- function(obj, Y, X, subIndices=NULL) {
  if (!is.null(subIndices)) {
    MatrixFoldIterator(X[subIndices,,drop=FALSE], Y[subIndices], obj$block_var[subIndices], balance=obj$balance, bootstrap=obj$bootstrap )
  } else {
    MatrixFoldIterator(X, Y,obj$block_var, balance=obj$balance, bootstrap=obj$bootstrap )
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
#' @param block_var variable denoting the cross-validation folds
#' @param balance try to balance each training sample so that the frequency of labels is equal across groups
#' @param bootstrap use bootstrap resampling of the training set
#' @param bootstrapMin 
#' @importFrom assertthat assert_that
#' @export
FoldIterator <- function(Y, block_var,  balance=FALSE, bootstrap=FALSE, bootstrapMin=2) {
  
  assert_that(length(block_var) == length(Y))
  assert_that(length(unique(block_var)) > 1)
  
  index <- 0
  ord <- NULL
  
  testSets <- split(1:length(block_var), block_var)
  
  trainSets <- invertFolds(testSets, 1:length(block_var))
  
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
  
  obj <- list(Y=Y, nextElem=nextEl, block_var=block_var, index=.getIndex, getTrainSets=.getTrainSets, 
              getTestSets=.getTestSets, getTestOrder=.getTestOrder, reset=.reset, balance=balance, bootstrap=bootstrap, nfolds=length(testSets))
  
  class(obj) <- c("FoldIterator", 'abstractiter', 'iter')
  obj
  
}


#' Create an iterator of train/test splits from a data matrix \code{X}, a vector \code{Y}, and an optional 
#' blocking variable \code{block_var} used for determining cross-validation folds.
#' 
#' @param X the data matrix
#' @param Y the label vector
#' @param block_var an index variable indicating the cross-validation folds.
#' @param nfolds the number of cross-validation folds: only relevant when block_var is not supplied
#' @param balance try to balance each training sample so that the frequency of labels is equal across groups
#' @param bootstrap use bootstrap resampling of the training set
#' @param bootstrapMin the minumum number of clases per bootstrap iteration.
#' @export
MatrixFoldIterator <- function(X, Y, block_var, balance=FALSE, bootstrap=FALSE, bootstrapMin=2) {
  
  if (nrow(X) != length(Y)) {
    stop("X matrix must have same number of rows as Y variable")
  }
  
  
  foldIter = FoldIterator(Y, block_var, balance=balance, bootstrap=bootstrap, bootstrapMin)
  
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
  
  
  obj <- list(X=X, Y=Y, block_var=block_var, nextElem=nextEl, index=foldIter$index, 
              getTrainSets=foldIter$getTrainSets, getTestSets=foldIter$getTestSets, 
              getTestOrder=foldIter$getTestOrder, reset=foldIter$reset, balance=foldIter$balance, bootstrap=foldIter$bootstrap, nfolds=foldIter$nfolds)
  
  class(obj) <- c("MatrixFoldIterator", 'abstractiter', 'iter')
  obj
  
}

#' NestedFoldIterator
#' 
#' @param Y a \code{vector} of labels
#' @param block_var an index variable indicating the cross-validation folds.
#' @param balance try to balance each training sample so that the frequency of labels is equal across groups
#' @param bootstrap
#' @param bootstrapMin
#' @export
NestedFoldIterator <- function(Y, block_var,  balance=FALSE, bootstrap=FALSE, bootstrapMin=2) {
  blockids <- sort(unique(block_var))
  
  if (length(Y) != length(block_var)) {
    stop("Y must have same length as block_var")
  }
  
  
  index <- 0
  
  blockNumber <- function() { index }
  
  nextEl <- function() {
    if (index < length(blockids)) {
      index <<- index + 1
      curBlock <- blockids[index]
      train.idx <- which(block_var != curBlock)
      iter <- FoldIterator(Y[train.idx], block_var=block_var[train.idx], balance=balance, bootstrap=bootstrap, bootstrapMin=bootstrapMin)
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
#' @param X the data \code{matrix}
#' @param Y a vector of labels
#' @param block_var
#' @param balance
#' @param bootstrap
#' @param bootstrapMin
#' @export
NestedMatrixIterator <- function(X, Y, block_var, balance=FALSE, bootstrap=FALSE, bootstrapMin=2) {
  blockids <- sort(unique(block_var))
  
  if (nrow(X) != length(Y)) {
    stop("X matrix must have same number of rows as Y variable")
  }
  
  
  index <- 0
  
  blockNumber <- function() { index }
  
  nextEl <- function() {
    if (index < length(blockids)) {
      index <<- index + 1
      curBlock <- blockids[index]
      train.idx <- which(block_var != curBlock)
      iter <- MatrixFoldIterator(X[train.idx,], Y[train.idx],  block_var=block_var[train.idx], balance=balance, bootstrap=bootstrap, bootstrapMin=bootstrapMin)
      
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


