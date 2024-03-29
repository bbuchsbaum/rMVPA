% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/crossval.R
\name{crossv_twofold}
\alias{crossv_twofold}
\title{Repeated Two-Fold Cross-Validation Data Preparation}
\usage{
crossv_twofold(data, y, block_var, block_ind = NULL, id = ".id", nreps = 15)
}
\arguments{
\item{data}{A data frame containing the training data.}

\item{y}{A response vector.}

\item{block_var}{An integer vector defining the cross-validation blocks.}

\item{block_ind}{A vector containing the ordered integer IDs of the blocks (optional).}

\item{id}{A character string specifying the identifier for the output data frame.}

\item{nreps}{An integer specifying the number of repetitions for two-fold cross-validation.}
}
\value{
A tibble containing the training and testing data, response vectors, and fold IDs for each repetition.
}
\description{
This function prepares the data for repeated two-fold cross-validation by
dividing the dataset into two folds based on the provided block variable.
It creates subsets of training and testing data for each repetition without
performing any analysis or fitting models.
}
\examples{
X <- data.frame(x1 = rnorm(100), x2 = rnorm(100))
y <- rep(letters[1:4], 25)
block_var <- rep(1:4, each = 25)
cv <- crossv_twofold(X, y, block_var, nreps = 10)
}
