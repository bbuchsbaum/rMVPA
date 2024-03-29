% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/allgeneric.R
\name{sub_result}
\alias{sub_result}
\title{Extract Row-wise Subset of a Result Object}
\usage{
sub_result(x, indices)
}
\arguments{
\item{x}{The input result object, which should be an instance of a classification or regression result class.}

\item{indices}{A vector of row indices to extract from the result object.}
}
\value{
A new result object containing only the specified row indices from the input result object.
}
\description{
This function extracts a row-wise subset of a classification or regression result object.
}
\examples{
# Create a simple classification result object
observed <- c(1, 0, 1, 1, 0)
predicted <- c(1, 0, 0, 1, 1)
result <- list(observed = observed, predicted = predicted)
class(result) <- c("classification_result", "list")

# Extract a subset of the result object
sub_result(result, c(2, 3, 4))
}
\seealso{
Other sub_result: 
\code{\link{sub_result.binary_classification_result}()},
\code{\link{sub_result.multiway_classification_result}()}
}
\concept{sub_result}
