% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/resample.R
\name{filter_roi}
\alias{filter_roi}
\title{Filter Region of Interest (ROI)}
\usage{
filter_roi(roi)
}
\arguments{
\item{roi}{A list containing the train and test ROI data.}
}
\value{
A list with filtered train and test ROI data.
}
\description{
This function filters an ROI, keeping only valid columns.
}
\details{
The function filters an ROI by removing columns with missing values (NA) and zero standard deviation.
It returns a list with filtered train and test ROI data.
}
\keyword{internal}
