#' @export
invertFolds <- function(foldSplit, allInd) { 
  index <-  relist(rank(unlist(foldSplit)), foldSplit)
  
  lapply(index, function(split) {
    allInd[-split]
  })
}

#' Create an iterator froma block variable for iterating through the folds during cross-validation
#' @param X data matrix
#' @param Y the labels
#' @param blockVar variable denoting the cross-validation folds
#' @export
FoldIterator <- function(blockVar) {
   index <- 0
   
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
  
  .reset <- function() {
    index <<- 0
  }
  
  nextEl <- function() {
    if (index < length(trainSets)) { 
      index <<- index + 1
      list(trainIndex=trainSets[[index]], testIndex=testSets[[index]])
      
    } else {
      stop('StopIteration')
    }
  }
  
  obj <- list(nextElem=nextEl, blockVar=blockVar, index=.getIndex, getTrainSets=.getTrainSets, getTestSets=.getTestSets, reset=.reset)
  class(obj) <- c("FoldIterator", 'abstractiter', 'iter')
  obj
  
}
  

#' Create an iterator from a matrix (X), a dependent variable (Y) and a splitting variable (blockVar)
#' @param X data matrix
#' @param Y the labels
#' @param blockVar variable denoting the cross-validation folds
#' @export
MatrixFoldIterator <- function(X, Y, blockVar, trainSets=NULL, testSets=NULL) {
  
  if (is.null(trainSets) & is.null(testSets)) {
    testSets <- split(1:length(blockVar), blockVar)
    trainSets <- invertFolds(testSets, 1:length(blockVar))
  }
  
  if (is.null(trainSets) && !is.null(testSets)) {
    stop("must supply both trainSets and testSets")
  }
  
  if (is.null(testSets) && !is.null(trainSets)) {
    stop("must supply both trainSets and testSets")
  }
  
  if (nrow(X) != length(Y)) {
    stop("X matrix must have same number of rows as Y variable")
  }
  
  
  index <- 0
  ord <- NULL
  
  .getTrainSets <- function() {
    trainSets
  }
  .getTestSets <- function() {
    testSets
  }
  
  .getIndex <- function() {
    index
  }
  
  .reset <- function() {
    index <<- 0
  }
  
  .getTestOrder <- function() {
    if (is.null(ord))
      ord <<- order(unlist(testSets))
    ord
  }
  
  nextEl <- function() {
    if (index < length(trainSets)) { 
      index <<- index + 1
      list(Xtrain=X[trainSets[[index]], ], Ytrain=Y[trainSets[[index]]], Xtest=X[testSets[[index]],], Ytest=Y[testSets[[index]]], index=index)
      
    } else {
      stop('StopIteration')
    }
  }
  
  obj <- list(X=X, Y=Y, blockVar=blockVar, nextElem=nextEl, index=.getIndex, getTrainSets=.getTrainSets, getTestSets=.getTestSets, getTestOrder=.getTestOrder, reset=.reset)
  class(obj) <- c("MatrixFoldIterator", 'abstractiter', 'iter')
  obj
  
}


#' Create an iterator from a matrix (X), a dependent variable (Y) and a splitting variable (blockVar)
#' @param X data matrix
#' @param Y the labels
#' @param blockVar variable denoting the cross-validation folds
#' @export
MatrixSplitIterator <- function(X, Y, blockVar) {
  
  
  blocks <- split(1:length(blockVar), blockVar)

  if (nrow(X) != length(Y)) {
    stop("X matrix must have same number of rows as Y variable")
  }
  
  
  index <- 0
  
  .getBlocks <- function() {
    blocks
  }
  
  .getIndex <- function() {
    index
  }
  
  .reset <- function() {
    index <<- 0
  }
  
  
  nextEl <- function() {
    if (index < length(blocks)) { 
      index <<- index + 1
      list(X=X[blocks[[index]], ], Y=Y[blocks[[index]]])
      
    } else {
      stop('StopIteration')
    }
  }
  
  obj <- list(X=X, Y=Y, blockVar=blockVar, nextElem=nextEl, index=.getIndex, reset=.reset)
  class(obj) <- c("MatrixSplitIterator", 'abstractiter', 'iter')
  obj
  
}