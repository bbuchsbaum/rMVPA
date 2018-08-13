
#' crossv_k
#' 
#' sample dataset using k-fold cross-validation 
#' 
#' @param data the training data
#' @param y the response vector
#' @param k Number of folds (an integer).
#' @param id a character id
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

#' crossv_twofold
#' 
#' sample dataset via repeated two-fold cross-validation
#' @import utils
#' @inheritParams crossv_k
#' @param block_var an \code{integer} \code{vector} defining the cross-validation blocks
#' @param block_ind the ordered \code{integer} ids of the blocks
#' @param nreps the number of resamples
#' @export
#' @examples 
#' 
#' X <- data.frame(x1=rnorm(100), x2=rnorm(100))
#' y <- rep(letters[1:4], 25)
#' block_var <- rep(1:4, each=25)
#' cv <- crossv_twofold(X,y,block_var, nreps=10)
crossv_twofold <- function(data, y, block_var, block_ind=NULL, id = ".id", nreps=15) {
  if (!length(block_var) == length(y)) {
    stop("length of `block_var` must be equal to length(y)", call. = FALSE)
  }
  
  if (is.null(block_ind)) {
    block_ind <- seq(1, length(sort(unique(block_var))))
  }

  nhalf <- floor(length(block_ind)/2)
  assert_that(nhalf > 0)
  
  fold_sets <- utils::combn(block_ind, nhalf)
  nreps <- min(nreps, ncol(fold_sets))
  cols <- sample(1:ncol(fold_sets), nreps)
  
  fold_idx <- lapply(1:nreps, function(i) {
    bind <- fold_sets[, cols[i]]
    which(block_var %in% bind)
  })
  
  idx <- seq_len(nrow(data))
  
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


#' crossv_block
#' 
#' @inheritParams crossv_twofold
#' @importFrom modelr resample
#' @export
#' @examples 

#' X <- data.frame(x1=rnorm(100), x2=rnorm(100))
#' y <- rep(letters[1:4], 25)
#' block_var <- rep(1:4, each=25)
#' cv <- crossv_block(X,y,block_var)
crossv_block <- function(data, y, block_var, id = ".id") {
 
  if (!length(block_var) == length(y)) {
    stop("length of `block_var` must be equal to length(y)", call. = FALSE)
  }
  
  idx <- seq_len(nrow(data))
  fold_idx <- split(idx, block_var)

  
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

#' crossv_bootstrap_block
#' 
#' @inheritParams crossv_twofold
#' @importFrom modelr resample
#' @keywords internal
crossv_bootstrap_block <- function(data, y, block_var, nreps=5, id = ".id") {
  
  if (!length(block_var) == length(y)) {
    stop("length of `block_var` must be equal to length(y)", call. = FALSE)
  }
  
  idx <- seq_len(nrow(data))
  block_idx <- split(idx, block_var)
  
  assert_that(length(block_idx) > 1, msg="crossv_bootstrap_block: must have more than one block to bootstrap.")
  
  fold_idx <- lapply(1:length(block_idx), function(heldout) {
    replicate(nreps, sample(unlist(block_idx[-heldout]), replace=TRUE), simplify=FALSE)
  })
  
  fold_idx <- lapply(unlist(fold_idx, recursive=FALSE), sort)
  
  
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

#' bootstrap_blocked_cross_validation
#' 
#' Construct a cross-validation specification using a predefined blocking variable with bootstrap resamples
#' 
#' @rdname cross_validation
#' @export
bootstrap_blocked_cross_validation <- function(block_var, nreps=10) {
  ret <- list(block_var=block_var, nfolds=length(unique(block_var)), block_ind=sort(unique(block_var)), nreps=nreps)
  class(ret) <- c("bootstrap_blocked_cross_validation", "cross_validation", "list")
  ret
}

#' blocked_cross_validation
#' 
#' Construct a cross-validation specification using a predefined blocking variable
#' 
#' @param block_var an integer vector of indicating the cross-validation blocks. Each block is indicating by a unique integer.
#' @rdname cross_validation
#' @export
blocked_cross_validation <- function(block_var) {
  ret <- list(block_var=block_var, nfolds=length(unique(block_var)), block_ind=sort(unique(block_var)))
  class(ret) <- c("blocked_cross_validation", "cross_validation", "list")
  ret
}

#' custom_cross_validation
#' 
#' Construct a cross-validation specification that used a user-supplied set of training and test indices.
#' 
#' @inheritParams blocked_cross_validation
#' @param sample_set a \code{list} of training and test sample indices
#' @rdname cross_validation
#' @export
custom_cross_validation <- function(sample_set) {
  assert_that(is.list(sample_set))
  for (el in sample_set) {
    assert_that(all(names(el) == c("train", "test")))
  }
  
  ret <- list(sample_set=sample_set, nfolds=length(sample_set))
  class(ret) <- c("custom_cross_validation", "cross_validation", "list")
  ret
  
}
  
  
#' twofold_blocked_cross_validation
#' 
#' Construct a cross-validation specification that randomly partitions input set into two sets of blocks.
#' 
#' @export 
#' @inheritParams blocked_cross_validation
#' @param nreps the number of twofold spits
#' @rdname cross_validation
#' @examples 
#' 
#' blockvar <- rep(1:5, each=10)
#' nreps <- 5
#' cval <- twofold_blocked_cross_validation(blockvar, nreps=nreps)
#' samples <- crossval_samples(cval, as.data.frame(matrix(rnorm(50*50),50,50)), y=rep(letters[1:5],10))
#' stopifnot(nrow(samples) == nreps)
twofold_blocked_cross_validation <- function(block_var, nreps=10) {
  block_var <- as.integer(block_var)
  ret <- list(block_var=block_var, nfolds=2, nreps=nreps, block_ind=sort(unique(block_var)))
  class(ret) <- c("twofold_blocked_cross_validation", "cross_validation", "list")
  ret
}

#' kfold_cross_validation
#' 
#' Construct a cross-validation specification that randomly partitions input set into \code{nfolds} folds.
#' 
#' @param len the number of observations.
#' @param nfolds the number of cross-validation folds.
#' @rdname cross_validation
#' @examples 
#' 
#' cval <- kfold_cross_validation(len=100, nfolds=10)
#' samples <- crossval_samples(cval, data=as.data.frame(matrix(rnorm(100*10), 100, 10)), y=rep(letters[1:5],20))
#' stopifnot(nrow(samples) == 10)
#' @export
kfold_cross_validation <- function(len, nfolds=10) {
  block_var <- sample(rep(seq(1, nfolds), length.out=len))
  ret <- list(block_var=block_var, nfolds=nfolds)
  class(ret) <- c("kfold_cross_validation", "cross_validation", "list")
  ret
}


# nest <- function(cval) {
#   clist <- lapply(cval$block_ind, function(i) {
#     blocked_cross_validation(cval$block_var, exclude=i)
#   })
#   
#   class(clist) = c("nested_blocked_cross_validation", "list")
#   clist
# }


#' @export
crossval_samples.kfold_cross_validation <- function(obj, data,y,...) { 
  crossv_k(data, y, obj$nfolds)
}

#' @export
crossval_samples.blocked_cross_validation <- function(obj, data, y,...) { 
  crossv_block(data, y, obj$block_var)
}

#' @export
crossval_samples.bootstrap_blocked_cross_validation <- function(obj, data, y,...) { 
  crossv_bootstrap_block(data, y, block_var=obj$block_var, nreps=obj$nreps)
}

#' @export
crossval_samples.custom_cross_validation <- function(obj, data, y, id = ".id",...) {
  fold <- function(train, test) {
    list(
      ytrain = y[train],
      ytest = y[test],
      train = resample(data, train),
      test = resample(data, test)
    )
  }
  
  cols <- purrr::transpose(purrr::map(obj$sample_set, function(el) fold(el$train, el$test)))
  cols[[id]] <- modelr:::id(length(obj$sample_set))
  
  tibble::as_data_frame(cols)
}

#' @export
crossval_samples.twofold_blocked_cross_validation <- function(obj, data, y,...) { 
  crossv_twofold(data, y, obj$block_var, obj$block_ind, nreps=obj$nreps)
}


#' @export
print.blocked_cross_validation <- function(x,...) {
  cat("cross-validation: blocked \n")
  cat("  nobservations: ", length(x$block_var), "\n")
  cat("  nfolds: ", x$nfolds, "\n")
  cat("  block sizes: ", table(x$block_var), "\n")
}

#' @export
print.twofold_blocked_cross_validation <- function(x,...) {
  cat("twofold cross-validation: blocked \n")
  cat("  nobservations: ", length(x$block_var), "\n")
  cat("  nreps: ", x$nreps, "\n")
  cat("  block sizes: ", table(x$block_var), "\n")
}

#' @export
print.bootstrap_blocked_cross_validation <- function(x,...) {
  cat("cross-validation: bootstrap blocked \n")
  cat("  n observations: ", length(x$block_var))
  cat("  n bootstrap reps: ", x$nreps, "\n")
  cat("  block sizes: ", table(x$block_var), "\n")
}


#' @export
print.twofold_cross_validation <- function(x,...) {
  cat("cross-validation: repeated two-fold \n")
  cat("  nobservations: ", length(x$block_var))
  cat("  nfolds: ", 2, "\n")
  cat("  nreps: ", x$nreps, "\n")
}


#' @export
print.kfold_cross_validation <- function(x,...) {
  cat("cross-validation: k fold \n")
  cat("  nobservations: ", length(x$block_var), "\n")
  cat("  nfolds: ", x$nfolds, "\n")
}
