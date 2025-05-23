% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/common.R
\name{load_model}
\alias{load_model}
\title{Load a Pre-defined MVPA Model}
\usage{
load_model(name)
}
\arguments{
\item{name}{Character string specifying the model to load. Can be either:
\itemize{
  \item A pre-defined MVPA model name:
  \describe{
    \item{corclass}{Correlation-based classifier with template matching}
    \item{sda_notune}{Simple Shrinkage Discriminant Analysis without tuning}
    \item{sda_boot}{SDA with bootstrap resampling}
    \item{glmnet_opt}{Elastic net with EPSGO parameter optimization}
    \item{sparse_sda}{SDA with sparsity constraints}
    \item{sda_ranking}{SDA with automatic feature ranking}
    \item{mgsda}{Multi-Group Sparse Discriminant Analysis}
    \item{lda_thomaz}{Modified LDA for high-dimensional data}
    \item{hdrda}{High-Dimensional Regularized Discriminant Analysis}
  }
  \item Any valid model name from caret's model library (e.g., "rf" for random forest, "svmRadial" for SVM)
}}
}
\value{
A list containing the model specification with the following components:
  \describe{
    \item{type}{Model type: "Classification" or "Regression"}
    \item{library}{Required R package(s) for the model}
    \item{label}{Human-readable model name}
    \item{parameters}{Data frame describing tunable parameters}
    \item{grid}{Function to generate parameter tuning grid}
    \item{fit}{Function to fit the model}
    \item{predict}{Function to generate predictions}
    \item{prob}{Function to generate class probabilities (classification only)}
  }
}
\description{
Retrieves a model specification from either the pre-defined set of MVPA models or from caret's model library.
}
\examples{
# Load custom MVPA model
model <- load_model("sda_notune")

# Load correlation classifier with parameter tuning options
corr_model <- load_model("corclass")
print(corr_model$parameters)  # View tunable parameters

# Load caret's random forest model
rf_model <- load_model("rf")
print(rf_model$parameters)  # View RF parameters

# Load caret's SVM model
svm_model <- load_model("svmRadial")

}
\seealso{
\code{\link{MVPAModels}} for the complete list of available custom MVPA models

\code{\link[caret]{getModelInfo}} for the complete list of available caret models

\code{\link{mvpa_model}} for using these models in MVPA analyses
}
