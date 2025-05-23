% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rsa_model.R
\name{rsa_design}
\alias{rsa_design}
\title{Construct a design for an RSA (Representational Similarity Analysis) model}
\usage{
rsa_design(
  formula,
  data,
  block_var = NULL,
  split_by = NULL,
  keep_intra_run = FALSE
)
}
\arguments{
\item{formula}{A formula expression specifying the dissimilarity-based regression function.}

\item{data}{A named list containing the dissimilarity matrices and any other auxiliary variables.}

\item{block_var}{An optional \code{formula}, \code{character} name or \code{integer} vector designating the block structure.}

\item{split_by}{An optional \code{formula} indicating grouping structure for evaluating test performance.}

\item{keep_intra_run}{A \code{logical} indicating whether to include within-run comparisons (default: FALSE).}
}
\value{
A list with class attributes "rsa_design" and "list", containing:
  \describe{
    \item{formula}{The input formula}
    \item{data}{The input data}
    \item{split_by}{The split_by formula}
    \item{split_groups}{Grouping structure for split_by}
    \item{block_var}{Block structure}
    \item{include}{Logical vector for including/excluding comparisons}
    \item{model_mat}{Model matrix generated by rsa_model_mat}
  }
}
\description{
This function constructs a design for an RSA model using the provided formula, data, and optional parameters.
}
\details{
The function creates an RSA design based on the input parameters. It checks the validity of the input data and
handles splitting conditions for evaluation of test performance. It also processes optional block structures and
within-run comparisons.
}
\examples{
dismat <- dist(matrix(rnorm(100*100), 100, 100))
rdes <- rsa_design(~ dismat, list(dismat=dismat))
}
