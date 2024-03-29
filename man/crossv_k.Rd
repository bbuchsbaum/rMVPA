% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/crossval.R
\name{crossv_k}
\alias{crossv_k}
\title{K-fold Cross-Validation Data Preparation}
\usage{
crossv_k(data, y, k = 5, id = ".id")
}
\arguments{
\item{data}{A data frame containing the training data.}

\item{y}{A response vector.}

\item{k}{An integer specifying the number of folds for cross-validation.}

\item{id}{A character string specifying the identifier for the output data frame.}
}
\value{
A tibble containing the training and testing data, response vectors, and fold IDs for each fold.
}
\description{
This function prepares the data for k-fold cross-validation by dividing the
dataset into k folds. It creates subsets of training and testing data for
each fold without performing any analysis or fitting models.
}
\examples{
data <- iris[,-5]
y <- iris$Species
result <- crossv_k(data, y, k = 5)
}
