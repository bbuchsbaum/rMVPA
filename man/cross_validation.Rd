% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/crossval.R
\name{bootstrap_blocked_cross_validation}
\alias{bootstrap_blocked_cross_validation}
\alias{blocked_cross_validation}
\alias{sequential_blocked_cross_validation}
\alias{custom_cross_validation}
\title{bootstrap_blocked_cross_validation}
\usage{
bootstrap_blocked_cross_validation(block_var, nreps = 10, weights = NULL)

blocked_cross_validation(block_var)

sequential_blocked_cross_validation(block_var, nfolds = 2, nreps = 4)

custom_cross_validation(sample_set)
}
\arguments{
\item{block_var}{An integer vector indicating the cross-validation blocks. Each block is indicated by a unique integer.}

\item{nreps}{The number of repetitions for the cross-validation procedure.}

\item{weights}{A numeric vector of the same length as `block_var`, representing the weights for each sample.
Higher weights indicate that observations will be sampled more often. If not provided, all samples are treated as equally likely.}

\item{nfolds}{The number of folds to divide each sequence of trials within a block.}

\item{sample_set}{A list of training and test sample indices. Each element of the list must be a named list with two elements: "train" and "test".}
}
\value{
A list containing the cross-validation specification, with class attributes "bootstrap_blocked_cross_validation", "cross_validation", and "list".

A list containing the cross-validation specification, with class attributes "blocked_cross_validation", "cross_validation", and "list".

A list containing the cross-validation specification, with class attributes "sequential_blocked_cross_validation", "cross_validation", and "list".

A list containing the custom cross-validation specification, with class attributes "custom_cross_validation", "cross_validation", and "list".
}
\description{
Bootstrap Blocked Cross-Validation Specification

This function constructs a cross-validation specification using a predefined blocking variable.

This function constructs a cross-validation specification using a predefined blocking variable, dividing each block into a specified number of folds.

This function constructs a cross-validation specification that uses a user-supplied set of training and test indices.
}
\details{
This function constructs a cross-validation specification using a predefined blocking variable
and creates bootstrap resamples within the blocks.


The function first checks if the provided weights are non-negative and normalizes them to sum to 1.
It then constructs a list containing the block variable, number of folds, block indices, number of repetitions, and weights.
The output list is assigned the class `"bootstrap_blocked_cross_validation"`, `"cross_validation"`, and `"list"`.

The function constructs a list containing the block variable, number of folds, and block indices.
The output list is assigned the class `"blocked_cross_validation"`, `"cross_validation"`, and `"list"`.

The function constructs a list containing the block variable, number of folds, number of repetitions, and block indices.
The output list is assigned the class `"sequential_blocked_cross_validation"`, `"cross_validation"`, and `"list"`.

The custom_cross_validation class allows users to define their own cross-validation structure by providing a set of training and test indices. This can be useful in situations where the standard cross-validation methods (e.g., k-fold, leave-one-out) do not adequately represent the desired validation structure.

The function constructs a list containing the sample set and the number of folds, derived from the length of the sample set. The output list is assigned the class `"custom_cross_validation"`, `"cross_validation"`, and `"list"`.
}
\examples{
block_var <- rep(1:5, each=50)
weights <- runif(length(block_var))
weights[1] = 0
cval <- bootstrap_blocked_cross_validation(block_var, weights=weights)
X <- matrix(rnorm(length(block_var) * 10), length(block_var), 10)
y <- rep(letters[1:5], length.out=length(block_var))

sam <- crossval_samples(cval, as.data.frame(X), y)
block_var <- rep(1:5, each=50)
cval <- blocked_cross_validation(block_var)
X <- matrix(rnorm(length(block_var) * 10), length(block_var), 10)
y <- rep(letters[1:5], length.out=length(block_var))

sam <- crossval_samples(cval, as.data.frame(X), y)
block_var <- rep(1:5, each=50)
nfolds <- 2
nreps <- 4
cval <- sequential_blocked_cross_validation(block_var, nfolds, nreps)
X <- matrix(rnorm(length(block_var) * 10), length(block_var), 10)
y <- rep(letters[1:5], length.out=length(block_var))

sam <- crossval_samples(cval, as.data.frame(X), y)
sample_set <- list(
  list(train = 1:80, test = 81:100),
  list(train = 1:60, test = 61:100),
  list(train = 1:40, test = 41:100)
)
cval <- custom_cross_validation(sample_set)
X <- matrix(rnorm(100 * 10), 100, 10)
y <- rep(letters[1:4], length.out=100)

sam <- crossval_samples(cval, as.data.frame(X), y)
}
