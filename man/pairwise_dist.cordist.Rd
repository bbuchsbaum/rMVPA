% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/distcalc.R
\name{pairwise_dist.cordist}
\alias{pairwise_dist.cordist}
\title{Compute Pairwise Correlation Distances}
\usage{
\method{pairwise_dist}{cordist}(obj, X)
}
\arguments{
\item{obj}{A distance function object of class \code{c("cordist", "distfun")}.}

\item{X}{A numeric matrix (rows = observations, columns = variables).}
}
\value{
An \strong{N x N numeric matrix} of pairwise distances.
}
\description{
Computes a full NxN matrix of correlation-based distances: \code{1 - cor(t(X), method=obj$method)}.
\strong{No block-based exclusion is performed here.}
}
\details{
This function calculates correlation distances among the rows of \code{X}. 
If you have a block variable and wish to exclude same-block comparisons, 
handle that \emph{after} obtaining this full matrix (e.g., in \code{second_order_similarity}).
}
\examples{
X <- matrix(rnorm(100), 10, 10)
dist_obj <- cordist(method="pearson")
dist_matrix <- pairwise_dist(dist_obj, X)

}
