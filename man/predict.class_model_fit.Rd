% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model_fit.R
\name{predict.class_model_fit}
\alias{predict.class_model_fit}
\title{This function predicts class labels and probabilities for new data using a fitted model.}
\usage{
\method{predict}{class_model_fit}(object, newdata, sub_indices = NULL, ...)
}
\arguments{
\item{object}{A fitted model object of class \code{class_model_fit}.}

\item{newdata}{New data to predict on, either as a \code{matrix} or a \code{NeuroVec} or \code{NeuroSurfaceVector} object.}

\item{sub_indices}{The subset of row indices to compute predictions on (optional).}

\item{...}{Additional arguments to be passed to the underlying prediction function.}
}
\value{
A list containing class predictions and probabilities with class attributes "classification_prediction", "prediction", and "list".
}
\description{
This function predicts class labels and probabilities for new data using a fitted model.
}
