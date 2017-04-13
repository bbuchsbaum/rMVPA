
#' crossv_k
#' 
#' @param data
#' @param y
#' @param k Number of folds (an integer).
#' @param id
#' 
#' @export
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

#' crossv_block
#' 
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
blocked_cross_validation <- function(block_var) {
  ret <- list(block_var=block_var, nfolds=length(unique(block_var)))
  class(ret) <- c("blocked_cross_validation", "cross_validation", "list")
  ret
}

#' kfold_cross_validation
#' 
#' @export
#' @param len the number of observations
#' @param balance 
#' @param boostrap
kfold_cross_validation <- function(len, nfolds=10) {
  block_var <- sample(rep(seq(1, nfolds), length.out=len))
  ret <- list(block_var=block_var, nfolds=nfolds)
  class(ret) <- c("kfold_cross_validation", "cross_validation", "list")
  ret
}

#' crossval_samples
#' 
#' @export
crossval_samples <- function(obj, data, y) { UseMethod("crossval_samples") }


## todo need to implement local version which stores 'y' variable in data.frame (train, test, y_train, y_test)
#' @export
crossval_samples.kfold_cross_validation <- function(obj, data,y) { 
  crossv_k(data, y, obj$nfolds)
}

#' @export
crossval_samples.blocked_cross_validation <- function(obj, data, y, exclude_block=NULL) { 
  crossv_block(data, y, obj$block_var, exclude_block=exclude_block)
}

