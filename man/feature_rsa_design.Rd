% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/feature_rsa_model.R
\name{feature_rsa_design}
\alias{feature_rsa_design}
\title{Create a Feature-Based RSA Design}
\usage{
feature_rsa_design(
  S = NULL,
  F = NULL,
  labels,
  k = 0,
  max_comps = 10,
  block_var = NULL
)
}
\arguments{
\item{S}{A symmetric similarity matrix representing the feature space relationships.
If NULL, you must supply F.}

\item{F}{A feature space matrix (observations by features). If supplied, this overrides S and k.}

\item{labels}{Vector of labels corresponding to the rows/columns of S or observations of F.}

\item{k}{Integer specifying the number of feature dimensions to retain when using S. If 0 (default),
automatically determines dimensions using eigenvalue threshold > 1 (minimum 2 dimensions kept).
This parameter is ignored if F is supplied directly (k becomes ncol(F)).}

\item{max_comps}{Initial upper limit for the number of components to be derived from the
feature space F by subsequent `feature_rsa_model` methods (PCA, PLS, SCCA).
This value is automatically capped by the final feature dimensionality `k`. Default 10.}

\item{block_var}{Optional blocking variable for cross-validation. If provided and
`crossval` is `NULL` in `feature_rsa_model`, a blocked cross-validation
scheme will be generated using this vector.}
}
\value{
A \code{feature_rsa_design} object (S3 class) containing:
  \describe{
    \item{S}{The input similarity matrix (if used)}
    \item{F}{Feature space projection matrix (k dimensions)}
    \item{labels}{Vector of observation labels}
    \item{k}{The final number of feature dimensions used}
    \item{max_comps}{The upper limit on components (<= k)}
    \item{block_var}{Optional blocking variable for cross-validation}
  }
}
\description{
Creates a design for feature-based Representational Similarity Analysis (RSA).
You can either supply a similarity matrix S (and optionally select dimensions)
or directly supply a feature matrix F.
}
\details{
This function defines the feature space representation for the analysis.
If F is supplied directly, it is used as-is, and `k` becomes `ncol(F)`.
If only S is supplied, an eigen decomposition of S is performed.
`k` determines how many eigenvectors form the feature matrix F. If `k=0`,
dimensions with eigenvalues > 1 are kept (minimum 2).
`max_comps` sets an upper bound for the number of components that model-fitting
methods (like PCA, PLS, SCCA in `feature_rsa_model`) can use, and it cannot
exceed the final feature dimensionality `k`.
}
