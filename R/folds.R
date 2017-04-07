#' @export
#' @param k Number of folds (an integer).
crossv_k <- function(data, y, k = 5, id = ".id") {
  if (!is.numeric(k) || length(k) != 1) {
    stop("`k` must be a single integer.", call. = FALSE)
  }
  
  n <- nrow(data)
  folds <- sample(rep(1:k, length.out = n))
  
  idx <- seq_len(n)
  fold_idx <- split(idx, folds)
  
  fold <- function(test) {
    tidx <- setdiff(idx, test)
    list(
      ytrain = y[tidx],
      ytest = y[test],
      train = resample(data, setdiff(idx, test)),
      test = resample(data, test)
    )
  }
  
  cols <- purrr::transpose(purrr::map(fold_idx, fold))
  cols[[id]] <- modelr:::id(k)
  
  tibble::as_data_frame(cols)
}


#' @export
#' @param block_var the blocking variable (an integer vector)
#' @importFrom modelr resample
#' @rdname crossv_block
crossv_block <- function(data, y, block_var, id = ".id", exclude_block=NULL) {
 
  if (!length(block_var) == length(y)) {
    stop("length of `block_var` must be equal to length(y)", call. = FALSE)
  }
  

  if (!is.null(exclude_block)) {
    idx <- seq_len(nrow(data))
    keep <- block_var != exclude_block
    idx <- idx[keep]
    fold_idx <- split(idx, block_var[keep])
  } else {
    idx <- seq_len(nrow(data))
    fold_idx <- split(idx, block_var)
  }
  
  
  fold <- function(test) {
    tidx <- setdiff(idx, test)
    list(
      ytrain = y[tidx],
      ytest = y[test],
      train = modelr::resample(data, tidx),
      test = modelr::resample(data, test)
    )
  }
  
  cols <- purrr::transpose(purrr::map(fold_idx, fold))
  cols[[id]] <- modelr:::id(length(fold_idx))
  
  tibble::as_data_frame(cols)
}



#' blocked_cross_validation
#' 
#' construct a cross-validation specification using a predefined blocking variable
#' 
#' @param block_var an integer vector of indicating the cross-validation blocks. Each block is indicating by a unique integer.
#' @param balance logical indicating whether cross-validation blocks should automatically balanced, using undersampling if necessary.
#' @param bootstrap logical indicating whether training samples should be sampled with replacement.
#' @export
blocked_cross_validation <- function(block_var, balance=FALSE, bootstrap=FALSE) {
  ret <- list(block_var=block_var, balance=balance, bootstrap=bootstrap, nfolds=length(unique(block_var)))
  class(ret) <- c("blocked_cross_validation", "CrossValidation", "list")
  ret
}

#' KFoldCrossValidation
#' 
#' @export
#' @param len the number of observations
#' @param balance 
#' @param boostrap
KFoldCrossValidation <- function(len, nfolds=10, balance=FALSE, bootstrap=FALSE) {
  block_var <- sample(rep(seq(1, nfolds), length.out=len))
  ret <- list(block_var=block_var, balance=balance, bootstrap=bootstrap, nfolds=nfolds)
  class(ret) <- c("KFoldCrossValidation", "CrossValidation", "list")
  ret
}

#' crossval_samples
#' 
#' @export
crossval_samples <- function(obj, data, y) { UseMethod("crossval_samples") }


## todo need to implement local version which stores 'y' variable in data.frame (train, test, y_train, y_test)
crossval_samples.KFoldCrossValidation <- function(obj, data,y) { 
  crossv_k(data, y, obj$nfolds)
}

crossval_samples.blocked_cross_validation <- function(obj, data, y, exclude_block=NULL) { 
  crossv_block(data, y, obj$block_var, exclude_block=exclude_block)
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
#' @param X the data \cdoe{matrix}
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
  
  
  