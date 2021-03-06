% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/crossval.R
\name{crossv_block}
\alias{crossv_block}
\title{crossv_block}
\usage{
crossv_block(data, y, block_var, id = ".id")
}
\arguments{
\item{data}{the training data}

\item{y}{the response vector}

\item{block_var}{an \code{integer} \code{vector} defining the cross-validation blocks}

\item{id}{a character id}
}
\description{
crossv_block
}
\examples{
X <- data.frame(x1=rnorm(100), x2=rnorm(100))
y <- rep(letters[1:4], 25)
block_var <- rep(1:4, each=25)
cv <- crossv_block(X,y,block_var)
}
