
#' K-fold Cross-Validation Data Preparation
#'
#' This function prepares the data for k-fold cross-validation by dividing the
#' dataset into k folds. It creates subsets of training and testing data for
#' each fold without performing any analysis or fitting models.
#'
#' @param data A data frame containing the training data.
#' @param y A response vector.
#' @param k An integer specifying the number of folds for cross-validation.
#' @param id A character string specifying the identifier for the output data frame.
#' @return A tibble containing the training and testing data, response vectors, and fold IDs for each fold.
#' @examples
#' data <- iris[,-5]
#' y <- iris$Species
#' result <- crossv_k(data, y, k = 5)
#' #' @importFrom modelr resample
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
      train = modelr::resample(data, setdiff(idx, test)),
      test = modelr::resample(data, test)
    )
  }
  
  cols <- purrr::transpose(purrr::map(fold_idx, fold))
  cols[[id]] <- gen_id(k)
  
  tibble::as_tibble(cols, .name_repair = .name_repair)
}

#' Repeated Two-Fold Cross-Validation Data Preparation
#'
#' This function prepares the data for repeated two-fold cross-validation by
#' dividing the dataset into two folds based on the provided block variable.
#' It creates subsets of training and testing data for each repetition without
#' performing any analysis or fitting models.
#'
#' @param data A data frame containing the training data.
#' @param y A response vector.
#' @param block_var An integer vector defining the cross-validation blocks.
#' @param block_ind A vector containing the ordered integer IDs of the blocks (optional).
#' @param id A character string specifying the identifier for the output data frame.
#' @param nreps An integer specifying the number of repetitions for two-fold cross-validation.
#' @return A tibble containing the training and testing data, response vectors, and fold IDs for each repetition.
#' @examples
#' X <- data.frame(x1 = rnorm(100), x2 = rnorm(100))
#' y <- rep(letters[1:4], 25)
#' block_var <- rep(1:4, each = 25)
#' cv <- crossv_twofold(X, y, block_var, nreps = 10)
crossv_twofold <- function(data, y, block_var, block_ind=NULL, id = ".id", nreps=15) {
  ## every time this is called, it regenerates new indices
  
  if (nreps < 2) {
    stop("'nreps' must be at least 2")
  }
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
  
 
  cols <- as.integer(seq(1, ncol(fold_sets), length.out=nreps))
  #sample(1:ncol(fold_sets), nreps)
 
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
  
  ## this could return a function that when given data, returns a tibble
  ## fun <- function(data) {... for every rows, mutate(train = ..., test = ...)} 
  
  cols <- purrr::transpose(purrr::map(fold_idx, fold))
  cols[[id]] <-gen_id(length(fold_idx))
  
  tibble::as_tibble(cols, .name_repair = .name_repair)
  
  
}


#' Block Cross-Validation Data Preparation
#'
#' This function prepares the data for block cross-validation by dividing the dataset
#' based on the provided block variable. It creates subsets of training and testing
#' data for each block without performing any analysis or fitting models.
#'
#' @param data A data frame containing the training data.
#' @param y A response vector.
#' @param block_var An integer vector defining the cross-validation blocks.
#' @param id A character string specifying the identifier for the output data frame.
#' @return A tibble containing the training and testing data, response vectors, and block IDs for each fold.
#' @examples
#' X <- data.frame(x1 = rnorm(100), x2 = rnorm(100))
#' y <- rep(letters[1:4], 25)
#' block_var <- rep(1:4, each = 25)
#' cv <- crossv_block(X, y, block_var)
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
  cols[[id]] <- gen_id(length(fold_idx))
  
  tibble::as_tibble(cols, .name_repair = .name_repair)
}

#' Block Bootstrap Cross-Validation Data Preparation
#'
#' This function prepares the data for block bootstrap cross-validation by dividing the dataset
#' based on the provided block variable. It creates subsets of training and testing
#' data for each block using bootstrap sampling within the training blocks, without performing any analysis or fitting models.
#'
#' @param data A data frame containing the training data.
#' @param y A response vector.
#' @param block_var An integer vector defining the cross-validation blocks.
#' @param nreps An integer specifying the number of bootstrap repetitions.
#' @param id A character string specifying the identifier for the output data frame.
#' @param weights An optional numeric vector of weights to be used for bootstrap sampling.
#'
#' @details
#' The function first checks if the length of the `block_var` vector matches the length of the response vector `y`.
#' It then creates a list of block indices and ensures there is more than one block to bootstrap. If weights are provided,
#' the function splits the weights according to the block variable.
#'
#' The function performs bootstrap sampling within the training blocks but keeps the test set fixed.
#' For each block, it generates a list of training indices using bootstrap sampling and creates the corresponding
#' training and testing data sets.
#'
#' @return A tibble containing the training and testing data, response vectors, and block IDs for each fold.
#'
#' @keywords internal
crossv_bootstrap_block <- function(data, y, block_var, nreps=5, id = ".id", weights=NULL) {
  
  if (!length(block_var) == length(y)) {
    stop("length of `block_var` must be equal to length(y)", call. = FALSE)
  }
  
  idx <- seq_len(nrow(data))
  block_idx <- split(idx, block_var)
  
  
  assert_that(length(block_idx) > 1, msg="crossv_bootstrap_block: must have more than one block to bootstrap.")
  
 
  ## alter so that you bootstrap within the training blocks but test set is fixed
  fold_idx <- if (!is.null(weights)) {
    block_wts <- split(weights, block_var)
    fold_idx <- lapply(seq_along(block_idx), function(heldout) {
      replicate(nreps, sample(unlist(block_idx[-heldout]), replace=TRUE, prob=unlist(block_wts[-heldout])), simplify=FALSE)
    })
  } else {
    fold_idx <- lapply(seq_along(block_idx), function(heldout) {
      replicate(nreps, sample(unlist(block_idx[-heldout]), replace=TRUE), simplify=FALSE)
    })
  }
  

  
  #fold_idx <- lapply(unlist(fold_idx, recursive=FALSE), sort)
  #fold_idx <- lapply(fold_idx, function(fidx) setdiff(idx, fidx))
  
  
  fold <- function(tidx, block) {
    list(
      ytrain = y[tidx],
      ytest = y[block_idx[[block]]],
      train = modelr::resample(data, tidx),
      test = modelr::resample(data, block_idx[[block]])
    )
  }
  
  cols <- unlist(lapply(seq_along(fold_idx), function(i) {
    lapply(fold_idx[[i]], function(tidx) {
      fold(tidx, i)
    })
  }), recursive=FALSE)
  
  cols <- purrr::transpose(cols)
  cols[[id]] <- gen_id(nreps * length(block_idx))
  
  tibble::as_tibble(cols, .name_repair = .name_repair)
}

#' Sequential Block Cross-Validation Data Preparation
#'
#' This function prepares the data for sequential block cross-validation by dividing the dataset
#' based on the provided block variable. It creates subsets of training and testing
#' data for each block using sequential sampling within the blocks, without performing any analysis or fitting models.
#'
#' @param data A data frame containing the training data.
#' @param y A response vector.
#' @param nfolds An integer specifying the number of folds for cross-validation.
#' @param block_var An integer vector defining the cross-validation blocks.
#' @param nreps An integer specifying the number of repetitions for each fold.
#' @param block_ind An optional integer vector specifying the ordered ids of the blocks.
#' @param id A character string specifying the identifier for the output data frame.
#'
#' @details
#' The function first checks if the length of the `block_var` vector matches the length of the response vector `y`.
#' It then creates a list of block indices and generates a fold sequence using the provided `nfolds` and `nreps` parameters.
#'
#' For each repetition and fold, the function identifies the indices corresponding to the test data and creates the
#' corresponding training and testing data sets.
#'
#' @return A tibble containing the training and testing data, response vectors, and fold IDs for each fold and repetition.
#'
#' @examples
#' X <- data.frame(x1=rnorm(100), x2=rnorm(100))
#' y <- rep(letters[1:4], 25)
#' block_var <- rep(1:4, each=25)
#' cv <- crossv_seq_block(X,y,2, block_var)
crossv_seq_block <- function(data, y, nfolds, block_var, nreps=4, block_ind = NULL, id = ".id") {
  
  if (!length(block_var) == length(y)) {
    stop("length of `block_var` must be equal to length(y)", call. = FALSE)
  }
  
  idx <- seq_len(nrow(data))
  block_idx <- split(idx, block_var)
  
  if (is.null(block_ind)) {
    block_ind <- seq(1, length(sort(unique(block_var))))
  }
  
  foldseq <- replicate(nreps, {
    unlist(lapply(block_idx, function(id) {
      as.integer(as.character(cut(id, nfolds, labels=sample(1:nfolds))))
    }))
    
  }, simplify=FALSE)
  
  fold_idx <- unlist(lapply(1:nreps, function(i) {
    lapply(1:nfolds, function(j) which(foldseq[[i]] == j))
  }), recursive=FALSE)
  
  
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
  cols[[id]] <- gen_id(length(fold_idx))
  
  tibble::as_tibble(cols, .name_repair = .name_repair)

}


#' bootstrap_blocked_cross_validation
#' 
#' Bootstrap Blocked Cross-Validation Specification
#'
#' This function constructs a cross-validation specification using a predefined blocking variable
#' and creates bootstrap resamples within the blocks.
#'
#' @param block_var An integer vector defining the cross-validation blocks.
#' @param nreps An integer specifying the number of repetitions for each fold.
#' @param weights A numeric vector of the same length as `block_var`, representing the weights for each sample.
#'        Higher weights indicate that observations will be sampled more often. If not provided, all samples are treated as equally likely.
#'
#' @details
#' The function first checks if the provided weights are non-negative and normalizes them to sum to 1.
#' It then constructs a list containing the block variable, number of folds, block indices, number of repetitions, and weights.
#' The output list is assigned the class `"bootstrap_blocked_cross_validation"`, `"cross_validation"`, and `"list"`.
#'
#' @return A list containing the cross-validation specification, with class attributes "bootstrap_blocked_cross_validation", "cross_validation", and "list".
#'
#' @examples
#' block_var <- rep(1:5, each=50)
#' weights <- runif(length(block_var))
#' weights[1] = 0
#' cval <- bootstrap_blocked_cross_validation(block_var, weights=weights)
#' X <- matrix(rnorm(length(block_var) * 10), length(block_var), 10)
#' y <- rep(letters[1:5], length.out=length(block_var))
#'
#' sam <- crossval_samples(cval, as.data.frame(X), y)
#' @rdname cross_validation
#' @export
bootstrap_blocked_cross_validation <- function(block_var, nreps=10, weights=NULL) {
  
  if (!is.null(weights)) {
    assert_that(length(weights) == length(block_var))
    assert_that(all(weights >= 0))  
    weights <- weights/sum(weights)
  }

  ret <- list(block_var=block_var, nfolds=length(unique(block_var)), 
              block_ind=sort(unique(block_var)), nreps=nreps, weights=weights)
  class(ret) <- c("bootstrap_blocked_cross_validation", "cross_validation", "list")
  ret
}


#' Blocked Cross-Validation Specification
#'
#' This function constructs a cross-validation specification using a predefined blocking variable.
#'
#' @param block_var An integer vector defining the cross-validation blocks.
#'
#' @details
#' The function constructs a list containing the block variable, number of folds, and block indices.
#' The output list is assigned the class `"blocked_cross_validation"`, `"cross_validation"`, and `"list"`.
#'
#' @return A list containing the cross-validation specification, with class attributes "blocked_cross_validation", "cross_validation", and "list".
#'
#' @examples
#' block_var <- rep(1:5, each=50)
#' cval <- blocked_cross_validation(block_var)
#' X <- matrix(rnorm(length(block_var) * 10), length(block_var), 10)
#' y <- rep(letters[1:5], length.out=length(block_var))
#'
#' sam <- crossval_samples(cval, as.data.frame(X), y)
#' @rdname cross_validation
#' @export
blocked_cross_validation <- function(block_var) {
  ret <- list(block_var=block_var, nfolds=length(unique(block_var)), block_ind=sort(unique(block_var)))
  class(ret) <- c("blocked_cross_validation", "cross_validation", "list")
  ret
}



#' Sequential Blocked Cross-Validation Specification
#'
#' This function constructs a cross-validation specification using a predefined blocking variable, dividing each block into a specified number of folds.
#'
#' @param block_var An integer vector indicating the cross-validation blocks. Each block is indicated by a unique integer.
#' @param nfolds The number of folds to divide each sequence of trials within a block.
#' @param nreps The number of repetitions for the cross-validation procedure.
#'
#' @details
#' The function constructs a list containing the block variable, number of folds, number of repetitions, and block indices.
#' The output list is assigned the class `"sequential_blocked_cross_validation"`, `"cross_validation"`, and `"list"`.
#'
#' @return A list containing the cross-validation specification, with class attributes "sequential_blocked_cross_validation", "cross_validation", and "list".
#'
#' @examples
#' block_var <- rep(1:5, each=50)
#' nfolds <- 2
#' nreps <- 4
#' cval <- sequential_blocked_cross_validation(block_var, nfolds, nreps)
#' X <- matrix(rnorm(length(block_var) * 10), length(block_var), 10)
#' y <- rep(letters[1:5], length.out=length(block_var))
#'
#' sam <- crossval_samples(cval, as.data.frame(X), y)
#' @rdname cross_validation
#' @export
sequential_blocked_cross_validation <- function(block_var, nfolds=2, nreps=4) {
  block_var <- as.integer(as.character(block_var))
  ret <- list(block_var=block_var, nfolds=nfolds, nreps=nreps, block_ind=sort(unique(block_var)))
  class(ret) <- c("sequential_blocked_cross_validation", "cross_validation", "list")
  ret
}



#' Custom Cross-Validation Specification
#'
#'
#' This function constructs a cross-validation specification that uses a user-supplied set of training and test indices.
#'
#' @param sample_set A list of training and test sample indices. Each element of the list must be a named list with two elements: "train" and "test".
#'
#' @details
#' The custom_cross_validation class allows users to define their own cross-validation structure by providing a set of training and test indices. This can be useful in situations where the standard cross-validation methods (e.g., k-fold, leave-one-out) do not adequately represent the desired validation structure.
#'
#' The function constructs a list containing the sample set and the number of folds, derived from the length of the sample set. The output list is assigned the class `"custom_cross_validation"`, `"cross_validation"`, and `"list"`.
#'
#' @return A list containing the custom cross-validation specification, with class attributes "custom_cross_validation", "cross_validation", and "list".
#'
#' @examples
#' sample_set <- list(
#'   list(train = 1:80, test = 81:100),
#'   list(train = 1:60, test = 61:100),
#'   list(train = 1:40, test = 41:100)
#' )
#' cval <- custom_cross_validation(sample_set)
#' X <- matrix(rnorm(100 * 10), 100, 10)
#' y <- rep(letters[1:4], length.out=100)
#'
#' sam <- crossval_samples(cval, as.data.frame(X), y)
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
#' Construct a cross-validation specification that randomly partitions the input set into two sets of blocks.
#'
#' This function creates a cross-validation scheme for cases where data is organized into blocks, and these blocks
#' are divided into two groups for evaluation. This approach can be useful when there is an inherent structure or
#' dependency within the blocks, and separating them can help to avoid biased estimates of model performance.
#' It returns an object of class "twofold_blocked_cross_validation", "cross_validation", and "list".
#'
#' @param block_var An integer vector representing the cross-validation blocks. Each block is indicated by a unique integer.
#' @param nreps An integer specifying the number of repetitions for the twofold split.
#' @return An object of class "twofold_blocked_cross_validation", "cross_validation", and "list" containing the block_var,
#'   nfolds (fixed at 2 for this function), nreps, and block_ind.
#' @export
#' @examples
#' blockvar <- rep(1:5, each=10)
#' nreps <- 5
#' cval <- twofold_blocked_cross_validation(blockvar, nreps=nreps)
#' samples <- crossval_samples(cval, as.data.frame(matrix(rnorm(50*50),50,50)), y=rep(letters[1:5],10))
#' stopifnot(nrow(samples) == nreps)
twofold_blocked_cross_validation <- function(block_var, nreps=10) {
  block_var <- as.integer(as.character(block_var))
  ret <- list(block_var=block_var, nfolds=2, nreps=nreps, block_ind=sort(unique(block_var)))
  class(ret) <- c("twofold_blocked_cross_validation", "cross_validation", "list")
  ret
}

#' kfold_cross_validation
#'
#' Construct a cross-validation specification that randomly partitions the input set into \code{nfolds} folds.
#'
#' This function creates a k-fold cross-validation scheme for cases where data needs to be split into a specified
#' number of folds for evaluation. It returns an object of class "kfold_cross_validation", "cross_validation", and "list".
#'
#' @param len An integer representing the number of observations.
#' @param nfolds An integer specifying the number of cross-validation folds.
#' @return An object of class "kfold_cross_validation", "cross_validation", and "list" containing the block_var and nfolds.
#' @examples
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
crossval_samples.sequential_blocked_cross_validation <- function(obj, data,y,...) { 
  crossv_seq_block(data, y, nfolds=obj$nfolds, block_var=obj$block_var, nreps=obj$nreps, block_ind=obj$block_ind)
}

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
  crossv_bootstrap_block(data, y, block_var=obj$block_var, nreps=obj$nreps, weights=obj$weights)
}

#' @export
#' @importFrom modelr resample
crossval_samples.custom_cross_validation <- function(obj, data, y, id = ".id",...) {
  fold <- function(train, test) {
    list(
      ytrain = y[train],
      ytest = y[test],
      train = modelr::resample(data, train),
      test = modelr::resample(data, test)
    )
  }
  
  cols <- purrr::transpose(purrr::map(obj$sample_set, function(el) fold(el$train, el$test)))
  cols[[id]] <- gen_id(length(obj$sample_set))
  
  tibble::as_tibble(cols, .name_repair = .name_repair)
}

#' @export
crossval_samples.twofold_blocked_cross_validation <- function(obj, data, y,...) { 
  crossv_twofold(data, y, obj$block_var, obj$block_ind, nreps=obj$nreps)
}


#' @export
#' @method print blocked_cross_validation
print.blocked_cross_validation <- function(x,...) {
  cat("cross-validation: blocked \n")
  cat("  nobservations: ", length(x$block_var), "\n")
  cat("  nfolds: ", x$nfolds, "\n")
  cat("  block sizes: ", table(x$block_var), "\n")
}

#' @export
#' @method print twofold_blocked_cross_validation
print.twofold_blocked_cross_validation <- function(x,...) {
  cat("twofold cross-validation: blocked \n")
  cat("  nobservations: ", length(x$block_var), "\n")
  cat("  nreps: ", x$nreps, "\n")
  cat("  block sizes: ", table(x$block_var), "\n")
}

#' @export
#' @method print bootstrap_blocked_cross_validation
print.bootstrap_blocked_cross_validation <- function(x,...) {
  cat("cross-validation: bootstrap blocked \n")
  cat("  n observations: ", length(x$block_var), "\n")
  cat("  n bootstrap reps: ", x$nreps, "\n")
  cat("  block sizes: ", table(x$block_var), "\n")
}


#' @export
#' @method print twofold_cross_validation
print.twofold_cross_validation <- function(x,...) {
  cat("cross-validation: repeated two-fold \n")
  cat("  nobservations: ", length(x$block_var))
  cat("  nfolds: ", 2, "\n")
  cat("  nreps: ", x$nreps, "\n")
}


#' @export
#' @method print kfold_cross_validation
print.kfold_cross_validation <- function(x,...) {
  cat("cross-validation: k fold \n")
  cat("  nobservations: ", length(x$block_var), "\n")
  cat("  nfolds: ", x$nfolds, "\n")
}
