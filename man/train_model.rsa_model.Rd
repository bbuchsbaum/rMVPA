% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rsa_model.R
\name{train_model.rsa_model}
\alias{train_model.rsa_model}
\title{Train an RSA Model}
\usage{
\method{train_model}{rsa_model}(obj, train_dat, y, indices, ...)
}
\arguments{
\item{obj}{An object of class \code{rsa_model}.}

\item{train_dat}{The training data.}

\item{y}{The response variable.}

\item{indices}{The indices of the training data.}

\item{...}{Additional arguments passed to the training method.}
}
\value{
The trained model.
}
\description{
This function trains an RSA (representational similarity analysis) model using the specified method and distance calculation.
}
