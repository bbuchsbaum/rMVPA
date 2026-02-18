#' @noRd
#' @keywords internal
check_len <- function(y, block_var) {
  futile.logger::flog.debug("Checking length of y and block_var")
  futile.logger::flog.debug("y: %s", paste(dim(y), collapse=" x "))
  futile.logger::flog.debug("block_var: %s", length(block_var))
  if (is.vector(y)) {
    if (!length(block_var) == length(y)) {
      stop("length of `block_var` must be equal to length(y)", call. = FALSE)
    }
  } else if (is.matrix(y)) {
    if (!nrow(y) == length(block_var)) {
      stop("number of rows in `y` must be equal to length(block_var)", call. = FALSE)
    }
  }
}


#' @noRd
#' @keywords internal
subset_y <- function(y, idx) {
  if (is.vector(y) || is.factor(y)) {
    y[idx]
  } else if (is.matrix(y)) {
    y[idx,,drop=FALSE]
  }
}

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
#' @importFrom modelr resample
#' @export
crossv_k <- function(data, y, k = 5, id = ".id") {
  if (!is.numeric(k) || length(k) != 1) {
    stop("`k` must be a single integer.", call. = FALSE)
  }
  
  if (k < 2) {
    stop("`k` must be at least 2 for cross-validation.", call. = FALSE)
  }
  
  n <- nrow(data)
  folds <- sample(rep(1:k, length.out = n))
  check_len(y, folds) # Ensure y and the generated folds are compatible in length
  
  idx <- seq_len(n)
  fold_idx <- split(idx, folds)
  
  fold <- function(test) {
    tidx <- setdiff(idx, test)
    list(
      ytrain = subset_y(y, tidx),
      ytest = subset_y(y, test),
      train = modelr::resample(data, setdiff(idx, test)),
      test = modelr::resample(data, test)
    )
  }
  
  
  cols <- purrr::transpose(purrr::map(fold_idx, fold))
  cols[[id]] <- gen_id(k)
  
  tibble::as_tibble(cols, .name_repair = "unique")
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
#' @noRd
crossv_twofold <- function(data, y, block_var, block_ind=NULL, id = ".id", nreps=15) {
  ## every time this is called, it regenerates new indices
  if (nreps < 2) {
    stop("'nreps' must be at least 2")
  }
  
  check_len(y, block_var)
  
  if (is.null(block_ind)) {
    block_ind <- seq(1, length(sort(unique(block_var))))
  }

  if (length(block_ind) < 2) {
    stop("`block_ind` must contain at least two unique blocks for twofold cross-validation.", call. = FALSE)
  }
  
  nhalf <- floor(length(block_ind)/2)
  assert_that(nhalf > 0)
  
  fold_sets <- utils::combn(block_ind, nhalf)
  nreps <- min(nreps, ncol(fold_sets))
  
 
  cols <- sample(seq_len(ncol(fold_sets)), nreps)
 
  fold_idx <- lapply(1:nreps, function(i) {
    bind <- fold_sets[, cols[i]]
    which(block_var %in% bind)
  })
  
  idx <- seq_len(nrow(data))
  
  fold <- function(test) {
    tidx <- setdiff(idx, test)
    list(
      ytrain = subset_y(y, tidx),
      ytest = subset_y(y, test),
      train = modelr::resample(data, tidx),
      test = modelr::resample(data, test)
    )
  }
  
  ## this could return a function that when given data, returns a tibble
  ## fun <- function(data) {... for every rows, mutate(train = ..., test = ...)} 
  
  cols <- purrr::transpose(purrr::map(fold_idx, fold))
  cols[[id]] <-gen_id(length(fold_idx))
  
  tibble::as_tibble(cols, .name_repair = "unique")
  
  
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
#' @export
crossv_block <- function(data, y, block_var, id = ".id") {
 
  check_len(y, block_var)
  
  idx <- seq_len(nrow(data))
  fold_idx <- split(idx, block_var)

  
  fold <- function(test) {
    tidx <- setdiff(idx, test)
    list(
      ytrain = subset_y(y, tidx),
      ytest = subset_y(y, test),
      train = modelr::resample(data, tidx),
      test = modelr::resample(data, test)
    )
  }
  
  cols <- purrr::transpose(purrr::map(fold_idx, fold))
  cols[[id]] <- gen_id(length(fold_idx))
  
  tibble::as_tibble(cols, .name_repair = "unique")
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
  check_len(y, block_var)
  
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
      ytrain = subset_y(y, tidx),
      ytest = subset_y(y, block_idx[[block]]),
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
  
  tibble::as_tibble(cols, .name_repair = "unique")
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
#' @noRd
crossv_seq_block <- function(data, y, nfolds, block_var, nreps=4, block_ind = NULL, id = ".id") {
  check_len(y, block_var)
  
  idx <- seq_len(nrow(data))
  block_idx <- split(idx, block_var)
  
  if (is.null(block_ind)) {
    block_ind <- seq(1, length(sort(unique(block_var))))
  }
  
  foldseq <- replicate(nreps, {
    unlist(lapply(block_idx, function(id_vec) {
      as.integer(as.character(cut(id_vec, nfolds, labels=sample(1:nfolds))))
    }))
    
  }, simplify=FALSE)
  
  fold_idx <- unlist(lapply(1:nreps, function(i) {
    lapply(1:nfolds, function(j) which(foldseq[[i]] == j))
  }), recursive=FALSE)
  
  
  fold <- function(test) {
    tidx <- setdiff(idx, test)
    list(
      ytrain = subset_y(y, tidx),
      ytest = subset_y(y, test),
      train = modelr::resample(data, tidx),
      test = modelr::resample(data, test)
    )
  }
  
  cols <- purrr::transpose(purrr::map(fold_idx, fold))
  cols[[id]] <- gen_id(length(fold_idx))
  
  tibble::as_tibble(cols, .name_repair = "unique")

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
  if (nfolds < 2) {
    stop("`nfolds` must be at least 2 for sequential blocked cross-validation.", call. = FALSE)
  }
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
  if (nreps < 2) {
    stop("'nreps' must be at least 2")
  }
  block_var <- as.integer(as.character(block_var))
  block_ind <- sort(unique(block_var))
  if (length(block_ind) < 2) {
    stop("`block_var` must contain at least two unique blocks for twofold cross-validation.", call. = FALSE)
  }
  ret <- list(block_var=block_var, nfolds=2, nreps=nreps, block_ind=block_ind)
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
#' sample_data <- as.data.frame(matrix(rnorm(100*10), 100, 10))
#' sample_y <- rep(letters[1:5], 20)
#' samples <- crossval_samples(cval, data = sample_data, y = sample_y)
#' stopifnot(nrow(samples) == 10)
#' @export
kfold_cross_validation <- function(len, nfolds=10) {
  if (nfolds < 2) {
    stop("`nfolds` must be at least 2 for k-fold cross-validation.", call. = FALSE)
  }
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
#' @rdname crossval_samples
crossval_samples.sequential_blocked_cross_validation <- function(obj, data, y,...) {
  crossv_seq_block(data, y, nfolds=obj$nfolds, block_var=obj$block_var, nreps=obj$nreps, block_ind=obj$block_ind)
}

#' @export
#' @rdname crossval_samples
crossval_samples.kfold_cross_validation <- function(obj, data,y,...) {
  crossv_k(data, y, obj$nfolds)
}

#' @export
#' @rdname crossval_samples
crossval_samples.blocked_cross_validation <- function(obj, data, y,...) {
  crossv_block(data, y, obj$block_var)
}

#' @export
#' @rdname crossval_samples
crossval_samples.bootstrap_blocked_cross_validation <- function(obj, data, y,...) {
  crossv_bootstrap_block(data, y, block_var=obj$block_var, nreps=obj$nreps, weights=obj$weights)
}

#' @export
#' @rdname crossval_samples
#' @param id Column name used for the fold identifier column in the returned tibble.
#' @importFrom modelr resample
crossval_samples.custom_cross_validation <- function(obj, data, y, id = ".id",...) {
  fold <- function(train, test) {
    list(
      ytrain = subset_y(y, train),
      ytest = subset_y(y, test),
      train = modelr::resample(data, train),
      test = modelr::resample(data, test)
    )
  }
  
  cols <- purrr::transpose(purrr::map(obj$sample_set, function(el) fold(el$train, el$test)))
  cols[[id]] <- gen_id(length(obj$sample_set))
  
  tibble::as_tibble(cols, .name_repair = "unique")
}

#' @export
#' @rdname crossval_samples
crossval_samples.twofold_blocked_cross_validation <- function(obj, data, y,...) {
  crossv_twofold(data, y, obj$block_var, obj$block_ind, nreps=obj$nreps)
}


#' @export
#' @method print blocked_cross_validation
print.blocked_cross_validation <- function(x, ...) {
  # Ensure crayon is available
  if (!requireNamespace("crayon", quietly = TRUE)) {
    stop("Package 'crayon' is required for pretty printing. Please install it.")
  }
  
  # Define color scheme
  header_style <- crayon::bold$cyan
  section_style <- crayon::yellow
  info_style <- crayon::white
  number_style <- crayon::green
  stat_style <- crayon::italic$blue
  
  # Print header
  cat("\n", header_style("Blocked Cross-Validation"), "\n\n")
  
  # Basic information
  cat(section_style("- Dataset Information"), "\n")
  cat(info_style("  - Observations: "), number_style(format(length(x$block_var), big.mark=",")), "\n")
  cat(info_style("  - Number of Folds: "), number_style(x$nfolds), "\n")
  
  # Block statistics
  block_sizes <- table(x$block_var)
  cat(section_style("- Block Information"), "\n")
  cat(info_style("  - Total Blocks: "), number_style(length(block_sizes)), "\n")
  cat(info_style("  - Mean Block Size: "), 
      number_style(format(mean(block_sizes), digits=2)), 
      stat_style(" (SD: "), 
      number_style(format(sd(block_sizes), digits=2)),
      stat_style(")"), "\n")
  cat(info_style("  - Block Sizes: "), 
      number_style(paste0(names(block_sizes), ": ", block_sizes, collapse=", ")), "\n\n")
}

#' @export
#' @method print twofold_blocked_cross_validation
print.twofold_blocked_cross_validation <- function(x, ...) {
  # Ensure crayon is available
  if (!requireNamespace("crayon", quietly = TRUE)) {
    stop("Package 'crayon' is required for pretty printing. Please install it.")
  }
  
  # Define color scheme
  header_style <- crayon::bold$cyan
  section_style <- crayon::yellow
  info_style <- crayon::white
  number_style <- crayon::green
  stat_style <- crayon::italic$blue
  
  # Print header
  cat("\n", header_style("Two-Fold Blocked Cross-Validation"), "\n\n")
  
  # Basic information
  cat(section_style("- Configuration"), "\n")
  cat(info_style("  - Observations: "), number_style(format(length(x$block_var), big.mark=",")), "\n")
  cat(info_style("  - Number of Folds: "), number_style("2"), "\n")
  cat(info_style("  - Repetitions: "), number_style(x$nreps), "\n")
  
  # Block statistics
  block_sizes <- table(x$block_var)
  cat(section_style("- Block Information"), "\n")
  cat(info_style("  - Total Blocks: "), number_style(length(block_sizes)), "\n")
  cat(info_style("  - Mean Block Size: "), 
      number_style(format(mean(block_sizes), digits=2)), 
      stat_style(" (SD: "), 
      number_style(format(sd(block_sizes), digits=2)),
      stat_style(")"), "\n")
  cat(info_style("  - Block Sizes: "), 
      number_style(paste0(names(block_sizes), ": ", block_sizes, collapse=", ")), "\n\n")
}

#' @export
#' @method print bootstrap_blocked_cross_validation
print.bootstrap_blocked_cross_validation <- function(x, ...) {
  # Ensure crayon is available
  if (!requireNamespace("crayon", quietly = TRUE)) {
    stop("Package 'crayon' is required for pretty printing. Please install it.")
  }
  
  # Define color scheme
  header_style <- crayon::bold$cyan
  section_style <- crayon::yellow
  info_style <- crayon::white
  number_style <- crayon::green
  stat_style <- crayon::italic$blue
  
  # Print header
  cat("\n", header_style("Bootstrap Blocked Cross-Validation"), "\n\n")
  
  # Basic information
  cat(section_style("- Configuration"), "\n")
  cat(info_style("  - Observations: "), number_style(format(length(x$block_var), big.mark=",")), "\n")
  cat(info_style("  - Bootstrap Repetitions: "), number_style(x$nreps), "\n")
  
  # Block statistics
  block_sizes <- table(x$block_var)
  cat(section_style("- Block Information"), "\n")
  cat(info_style("  - Total Blocks: "), number_style(length(block_sizes)), "\n")
  cat(info_style("  - Mean Block Size: "), 
      number_style(format(mean(block_sizes), digits=2)), 
      stat_style(" (SD: "), 
      number_style(format(sd(block_sizes), digits=2)),
      stat_style(")"), "\n")
  cat(info_style("  - Block Sizes: "), 
      number_style(paste0(names(block_sizes), ": ", block_sizes, collapse=", ")), "\n")
  
  # Weight information if present
  cat(section_style("- Sampling Weights"), "\n")
  if (!is.null(x$weights)) {
    cat(info_style("  - Status: "), crayon::green("Present"), "\n")
    cat(info_style("  - Range: "), 
        number_style(sprintf("[%.3f, %.3f]", min(x$weights), max(x$weights))), "\n")
    cat(info_style("  - Non-zero Weights: "), 
        number_style(sum(x$weights > 0)), 
        stat_style(" ("), 
        number_style(sprintf("%.1f%%", 100*mean(x$weights > 0))),
        stat_style(")"), "\n\n")
  } else {
    cat(info_style("  - Status: "), crayon::red("None"), " (uniform sampling)\n\n")
  }
}

#' @export
#' @method print kfold_cross_validation
print.kfold_cross_validation <- function(x, ...) {
  # Ensure crayon is available
  if (!requireNamespace("crayon", quietly = TRUE)) {
    stop("Package 'crayon' is required for pretty printing. Please install it.")
  }
  
  # Define color scheme
  header_style <- crayon::bold$cyan
  section_style <- crayon::yellow
  info_style <- crayon::white
  number_style <- crayon::green
  stat_style <- crayon::italic$blue
  
  # Print header
  cat("\n", header_style("K-Fold Cross-Validation"), "\n\n")
  
  # Basic information
  cat(section_style("- Configuration"), "\n")
  cat(info_style("  - Observations: "), number_style(format(length(x$block_var), big.mark=",")), "\n")
  cat(info_style("  - Number of Folds: "), number_style(x$nfolds), "\n")
  
  # Fold statistics
  fold_sizes <- table(x$block_var)
  cat(section_style("- Fold Information"), "\n")
  cat(info_style("  - Mean Fold Size: "), 
      number_style(format(mean(fold_sizes), digits=2)), 
      stat_style(" (SD: "), 
      number_style(format(sd(fold_sizes), digits=2)),
      stat_style(")"), "\n")
  cat(info_style("  - Fold Sizes: "), 
      number_style(paste0("Fold ", names(fold_sizes), ": ", fold_sizes, collapse=", ")), 
      "\n\n")
}

#' @export
get_nfolds.blocked_cross_validation <- function(obj, ...) {
  obj$nfolds
}

#' @export
train_indices.blocked_cross_validation <- function(obj, fold_num, ...) {
  if (fold_num < 1 || fold_num > obj$nfolds) {
    stop(paste0("Invalid fold_num: ", fold_num, ". Must be between 1 and ", obj$nfolds), call. = FALSE)
  }
  test_block_value <- obj$block_ind[fold_num]
  which(obj$block_var != test_block_value)
}

#' @export
get_nfolds.kfold_cross_validation <- function(obj, ...) {
  obj$nfolds
}

#' @export
train_indices.kfold_cross_validation <- function(obj, fold_num, ...) {
  if (fold_num < 1 || fold_num > obj$nfolds) {
    stop(paste0("Invalid fold_num: ", fold_num, ". Must be between 1 and ", obj$nfolds), call. = FALSE)
  }
  # For kfold_cross_validation, block_var directly contains fold assignments (1 to k)
  which(obj$block_var != fold_num)
}

#' @export
get_nfolds.twofold_blocked_cross_validation <- function(obj, ...) {
  # The number of "folds" here means the number of unique ways to split the blocks in half.
  # This is determined by nreps, which samples from combn(block_ind, nhalf).
  # For compute_cv_means, we consider each repetition as a distinct fold setup.
  obj$nreps 
}

#' @export
train_indices.twofold_blocked_cross_validation <- function(obj, fold_num, ...) {
  if (fold_num < 1 || fold_num > obj$nreps) {
    stop(paste0("Invalid fold_num: ", fold_num, ". Must be between 1 and ", obj$nreps, " (nreps)"), call. = FALSE)
  }
  
  # Reconstruct the specific fold split logic from crossv_twofold
  # This is complex because crossv_twofold samples combinations.
  # For compute_cv_means, we need a deterministic way if we iterate through folds.
  # The current structure of twofold_blocked_cross_validation does not directly store the chosen combinations.
  # This suggests a potential mismatch or a need for cv_spec to store the actual fold definitions if used by compute_cv_means.
  # 
  # *** Interim solution for train_indices.twofold_blocked_cross_validation ***
  # Option 1: Error out, stating it's not directly supported by compute_cv_means due to sampling.
  # Option 2: Re-run the sampling logic with a fixed seed based on fold_num? (complex, not ideal for S3 method)
  # Option 3: The object stores the actual `fold_idx` from a call to `crossval_samples`? (No, spec is created first)
  #
  # Given compute_cv_means iterates 1:n_folds, and n_folds is nreps, 
  # we must assume each "rep" defines a specific train/test split of blocks.
  # However, the object doesn't store which blocks are in the test set for rep `i`.
  #
  # For now, to allow tests to pass if they mock this, but to highlight the issue for real use:
  # Let's assume that `crossval_samples.twofold_blocked_cross_validation` would have been called
  # and the resulting `test_ind` for that specific `fold_num` (rep) is what defines the test set.
  # This method *cannot* know that without `crossval_samples` being run first and its output stored.
  # 
  # This indicates compute_cv_means may not be able to directly use twofold_blocked_cross_validation
  # in its current generic iteration form unless the object is augmented or the iteration logic changes.
  #
  # Let's make it error for now, as providing a meaningful train_indices is not possible
  # from the object's current state for a specific `fold_num` rep.
  stop(paste0("`train_indices.twofold_blocked_cross_validation` cannot deterministically provide training indices for a specific repetition (fold_num) ",
              "based on the current object structure. `compute_crossvalidated_means_sl` might not be directly compatible with this CV type ",
              "without prior generation and storage of explicit fold assignments."), call. = FALSE)
}

#' @export
get_nfolds.bootstrap_blocked_cross_validation <- function(obj, ...) {
  # Each block held out, and nreps bootstraps for each held-out block
  length(unique(obj$block_var)) * obj$nreps
}

#' @export
train_indices.bootstrap_blocked_cross_validation <- function(obj, fold_num, ...) {
  stop("`train_indices.bootstrap_blocked_cross_validation` is not implemented. `compute_crossvalidated_means_sl` is not suitable for bootstrap resampling due to its assumption of partitioning. Use a different CV scheme or a specialized estimator for bootstrap.", call. = FALSE)
}

#' @export
get_nfolds.sequential_blocked_cross_validation <- function(obj, ...) {
  # This CV scheme defines sub-folds within primary blocks.
  # get_nfolds refers to these sub-folds. The nreps is for different random assignments to these sub-folds.
  # compute_cv_means needs a single deterministic set of folds. We use obj$nfolds.
  obj$nfolds * obj$nreps # This assumes each rep * fold combo is a distinct fold for compute_cv_means.
                        # This is consistent with how crossv_seq_block generates fold_idx.
}

#' @export
train_indices.sequential_blocked_cross_validation <- function(obj, fold_num, ...) {
  # `fold_num` here iterates from 1 to (obj$nfolds * obj$nreps)
  # We need to map this back to a specific repetition and sub-fold j within that rep.
  total_sub_folds_per_rep <- obj$nfolds
  rep_num <- ((fold_num - 1) %/% total_sub_folds_per_rep) + 1
  sub_fold_num_in_rep <- ((fold_num - 1) %% total_sub_folds_per_rep) + 1

  if (rep_num > obj$nreps || sub_fold_num_in_rep > obj$nfolds) {
    stop(paste0("Invalid fold_num: ", fold_num, ". Derived rep_num=", rep_num, ", sub_fold_num=", sub_fold_num_in_rep, 
                 " is out of bounds for nreps=", obj$nreps, ", nfolds=", obj$nfolds, "."), call. = FALSE)
  }

  # Reconstruct the fold sequence for the *specific repetition* `rep_num`.
  # This is tricky because crossv_seq_block uses `sample(1:nfolds)` internally per block per rep.
  # To make this deterministic for a given `fold_num` passed to train_indices,
  # we would need to fix the seed OR the object would need to store the fold assignments.
  #
  # For now, let's assume the intent of compute_cv_means with this scheme is that for each of the total_sub_folds_per_rep * nreps
  # iterations, there's a defined test set. The `crossv_seq_block` logic does this by generating
  # `fold_idx <- unlist(lapply(1:nreps, function(i) { lapply(1:nfolds, function(j) which(foldseq[[i]] == j)) }), recursive=FALSE)`
  # The `fold_num` in train_indices would correspond to an element in this unlisted `fold_idx`.
  # This means we can't easily get *training* indices without knowing the *test* indices for that `fold_num`.
  #
  # This is another case where direct compatibility with compute_cv_means' iteration is problematic
  # without the cv_spec object storing the pre-computed fold assignments.
  stop(paste0("`train_indices.sequential_blocked_cross_validation` cannot deterministically provide training indices for a specific combined repetition and sub-fold (fold_num) ",
              "based on the current object structure. `compute_crossvalidated_means_sl` might not be directly compatible with this CV type ",
              "without prior generation and storage of explicit fold assignments."), call. = FALSE)
}

#' @export
get_nfolds.custom_cross_validation <- function(obj, ...) {
  length(obj$sample_set)
}

#' @export
train_indices.custom_cross_validation <- function(obj, fold_num, ...) {
  if (fold_num < 1 || fold_num > length(obj$sample_set)) {
    stop(paste0("Invalid fold_num: ", fold_num, ". Must be between 1 and ", length(obj$sample_set)), call. = FALSE)
  }
  obj$sample_set[[fold_num]]$train
}


#' Generate Cross-Validation Folds
#'
#' Convenience wrapper around \code{\link{crossval_samples}} for use inside
#' \code{\link{fit_roi}} methods.  Separates fold generation from the
#' iteration engine so that models can generate folds directly.
#'
#' @param cv_spec A cross-validation specification object (e.g.,
#'   \code{blocked_cross_validation}).
#' @param data A data.frame or tibble of training data.
#' @param y Response variable (factor, numeric vector, or matrix).
#' @return A tibble with columns \code{ytrain}, \code{ytest}, \code{train},
#'   \code{test}, and \code{.id}.
#' @seealso \code{\link{crossval_samples}}, \code{\link{fit_roi}}
#' @export
generate_folds <- function(cv_spec, data, y) {
  crossval_samples(cv_spec, data, y)
}
