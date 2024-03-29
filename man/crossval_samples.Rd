% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/allgeneric.R
\name{crossval_samples}
\alias{crossval_samples}
\title{crossval_samples}
\usage{
crossval_samples(obj, data, y, ...)
}
\arguments{
\item{obj}{A cross-validation control object.}

\item{data}{A data frame containing the predictor variables.}

\item{y}{A vector containing the response variable.}

\item{...}{Extra arguments passed to the specific cross-validation methods.}
}
\value{
A tibble containing the training and testing sets for each fold, as well as the response variables for both sets.
}
\description{
A generic function that applies a cross-validation scheme to split the data into training and testing sets.
It is used along with cross-validation control objects and S3 implementation functions to perform the cross-validation process.
}
\examples{
# Example with k-fold cross-validation
cval <- kfold_cross_validation(len=100, nfolds=10)
samples <- crossval_samples(cval, data=as.data.frame(matrix(rnorm(100*10), 100, 10)), y=rep(letters[1:5],20))
stopifnot(nrow(samples) == 10)
}
\seealso{
\code{\link{crossval_samples.sequential_blocked_cross_validation}},
  \code{\link{crossval_samples.kfold_cross_validation}},
  \code{\link{crossval_samples.blocked_cross_validation}},
  \code{\link{crossval_samples.bootstrap_blocked_cross_validation}},
  \code{\link{crossval_samples.custom_cross_validation}},
  \code{\link{crossval_samples.twofold_blocked_cross_validation}}
}
