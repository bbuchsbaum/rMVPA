
#' Create a Feature-Based RSA Design
#'
#' Creates a design for feature-based Representational Similarity Analysis (RSA) that maps neural patterns 
#' to a predefined feature space using dimensionality reduction techniques.
#'
#' @param S A symmetric similarity matrix representing the feature space relationships
#' @param labels Vector of labels corresponding to the rows/columns of S
#' @param k Integer specifying the number of feature dimensions to retain. If 0 (default),
#'          automatically determines dimensions using eigenvalue threshold > 1
#'
#' @return A \code{feature_rsa_design} object (S3 class) containing:
#'   \describe{
#'     \item{S}{The input similarity matrix}
#'     \item{F}{Feature space projection matrix (eigenvectors)}
#'     \item{labels}{Vector of observation labels}
#'   }
#'
#' @details
#' Unlike standard RSA which compares dissimilarity matrices directly, feature RSA projects neural patterns 
#' into a lower-dimensional feature space defined by the eigenvectors of the similarity matrix S. 
#' The number of dimensions k can be specified explicitly or determined automatically based on 
#' eigenvalues > 1. The resulting design is used with methods like sparse CCA or PLS to relate 
#' neural patterns to this feature space.
#'
#' @examples
#' # Create similarity matrix from feature vectors
#' features <- matrix(rnorm(100*10), 100, 10)  # 100 observations, 10 features
#' S <- tcrossprod(scale(features))  # similarity matrix
#' labels <- paste0("obs", 1:100)
#' 
#' # Create design with automatic dimension selection
#' design <- feature_rsa_design(S, labels)
#' 
#' # Create design with fixed number of dimensions
#' design_k5 <- feature_rsa_design(S, labels, k=5)
#'
#' @seealso 
#' \code{\link{feature_rsa_model}} for creating the full model
#' 
#' \code{\link{rsa_design}} for standard RSA designs
#' 
#' \code{\link{vector_rsa_design}} for vector-based RSA designs
#'
#' @importFrom assertthat assert_that
#' @export
feature_rsa_design <- function(S, labels, k=0) {
  assertthat::assert_that(is.matrix(S))
  assertthat::assert_that(nrow(S) == length(labels))
  assertthat::assert_that(isSymmetric(S))
  
  S <- (S + t(S))/2
  
  if (k == 0) {
    eres <- eigen(S)
    k <- max(which(eres$values > 1))
    k <- max(k, 2)
    F <- eres$vectors[,1:k, drop=FALSE]
  } else {
    assertthat::assert_that(k > 0 && k <= nrow(S))
    eres <- eigen(S)
    F <- eres$vectors[, 1:k, drop=FALSE]
  }
  
  ret <- list(
    S=S,
    F=F,
    labels=labels
  )
  
  class(ret) <- "feature_rsa_design"
  ret
}


#' Create a Feature-Based RSA Model
#'
#' Creates a model for feature-based Representational Similarity Analysis (RSA) that relates neural patterns
#' to a predefined feature space using dimensionality reduction and multivariate techniques.
#'
#' @param dataset An \code{mvpa_dataset} object containing the neural data
#' @param design A \code{feature_rsa_design} object specifying the feature space and its structure
#' @param method Character string specifying the analysis method. One of:
#'   \describe{
#'     \item{scca}{Sparse Canonical Correlation Analysis (default)}
#'     \item{pls}{Partial Least Squares}
#'     \item{pca}{Principal Component Analysis}
#'   }
#' @param crossval Optional cross-validation specification. If NULL and design contains block_var,
#'        creates blocked cross-validation using that variable
#'
#' @return A \code{feature_rsa_model} object (S3 class) containing:
#'   \describe{
#'     \item{dataset}{The input \code{mvpa_dataset}}
#'     \item{design}{The input \code{feature_rsa_design}}
#'     \item{crossval}{Cross-validation specification}
#'   }
#'
#' @details
#' Feature RSA models analyze the relationship between neural patterns and a predefined feature space
#' using multivariate techniques. Unlike traditional RSA which compares dissimilarity matrices,
#' this approach directly relates neural patterns to feature space dimensions through methods like
#' sparse CCA or PLS. Cross-validation can be specified explicitly or derived from the design's
#' block structure.
#'
#' @examples
#' # Create feature space and design
#' features <- matrix(rnorm(100*10), 100, 10)
#' S <- tcrossprod(scale(features))
#' labels <- paste0("obs", 1:100)
#' design <- feature_rsa_design(S, labels)
#'
#' # Create dataset (assuming you have neural data)
#' dataset <- mvpa_dataset(train_data, mask=mask_vol)
#'
#' # Create model with sparse CCA
#' model <- feature_rsa_model(dataset, design, method="scca")
#'
#' # Create model with PLS and explicit cross-validation
#' cv <- blocked_cross_validation(block_var)
#' model_pls <- feature_rsa_model(dataset, design, 
#'                               method="pls", 
#'                               crossval=cv)
#'
#' @seealso 
#' \code{\link{feature_rsa_design}} for creating the feature space design
#' 
#' \code{\link{rsa_model}} for traditional RSA models
#' 
#' \code{\link{vector_rsa_model}} for vector-based RSA models
#'
#' @importFrom assertthat assert_that
#' @export
feature_rsa_model <- function(dataset,
                       design,
                       method=c("scca", "pls", "pca"),
                       crossval=NULL) {
  
  assert_that(inherits(dataset, "mvpa_dataset"))
  assert_that(inherits(design, "feature_rsa_design"))
  
  if (is.null(crossval) && !is.null(design$block_var)) {
    crossval <- blocked_cross_validation(design$block_var)
  }
  
  assertthat::assert_that(!is.null(crossval))
  
  ret <- list(method=method,
              dataset=dataset,
              design=design,
              crossval=crossval)
  
  class(ret) <- "feature_rsa_model"
  ret
}


#' Train an RSA Model
#'
#' This function trains an RSA (representational similarity analysis) model using the specified method.
#'
#' @param obj An object of class \code{feature_rsa_model}
#' @param train_dat The training data
#' @param indices The indices of the training data
#' @param ... Additional arguments passed to the training method
#' @return The trained model with updated components
train_model.feature_rsa_model <- function(obj, train_dat, indices, ...) {
  X <- as.matrix(train_dat)
  
  # Train model based on specified method
  trained_model <- switch(obj$method[1],
    "scca" = whitening::scca(X, obj$design$F),
    "pls" = pls::plsr(X ~ obj$design$F),
    "pca" = prcomp(X, scale. = TRUE)
  )
  
  # Store training results
  obj$trained_model <- trained_model
  obj$training_indices <- indices
  
  # Calculate initial performance metrics
  predicted <- predict_model(obj, X)
  obj$performance <- evaluate_model(obj, predicted, obj$design$F[indices,])
  
  return(obj)
}


#' Predict Method for Feature RSA Model
#'
#' Makes predictions using a trained feature RSA model
#'
#' @param object The trained feature RSA model
#' @param newdata New data to make predictions on
#' @param ... Additional arguments passed to prediction method
#' @return Predicted values in the feature space
#' @export
predict_model.feature_rsa_model <- function(object, newdata, ...) {
  X <- as.matrix(newdata)
  # Project new data into the learned feature space
  pred <- predict(object$trained_model, X)
  return(pred)
}

#' Evaluate Method for Feature RSA Model
#'
#' Evaluates performance of a feature RSA model
#'
#' @param object The feature RSA model
#' @param predicted Predicted values
#' @param observed Observed values
#' @param ... Additional arguments
#' @return Performance metrics
#' @export
evaluate_model.feature_rsa_model <- function(object, predicted, observed, ...) {
  # Calculate correlation between predicted and observed feature projections
  cors <- diag(cor(predicted, observed))
  
  # Calculate mean squared error
  mse <- mean((predicted - observed)^2)
  
  list(
    correlations = cors,
    mse = mse
  )
}

#' Print Method for Feature RSA Model
#'
#' @param x The feature RSA model
#' @param ... Additional arguments
#' @export
print.feature_rsa_model <- function(x, ...) {
  cat("Feature RSA Model\n")
  cat("Method:", x$method, "\n")
  cat("Number of features:", ncol(x$design$F), "\n")
  cat("Number of observations:", nrow(x$design$S), "\n")
}

#' Summary Method for Feature RSA Model
#'
#' @param object The feature RSA model
#' @param ... Additional arguments
#' @export
summary.feature_rsa_model <- function(object, ...) {
  print(object)
  if (!is.null(object$trained_model)) {
    cat("\nModel Performance:\n")
    print(object$performance)
  }
}

#' @importFrom neuroim2 indices values
do_feature_rsa <- function(roi, mod_spec, rnum) {
  xtrain <- tibble::as_tibble(neuroim2::values(roi$train_roi), .name_repair=.name_repair)
  ind <- indices(roi$train_roi)
  ret <- try(train_model(mod_spec, xtrain, ind))
  if (inherits(ret, "try-error")) {
    tibble::tibble(result=list(NULL), indices=list(ind), performance=list(ret), id=rnum, error=TRUE, error_message=attr(ret, "condition")$message)
  } else {
    tibble::tibble(result=list(NULL), indices=list(ind), performance=list(ret), id=rnum, error=FALSE, error_message="~")
  }
}

  
#' feature_rsa_iterate
#'
#' Runs feature-based representational similarity analysis (RSA) for each voxel set in a list.
#'
#' @param mod_spec An object of class \code{rsa_model} specifying the RSA model.
#' @param vox_list A \code{list} of voxel indices or coordinates for each voxel set.
#' @param ids A \code{vector} of IDs for each voxel set (defaults to 1:length(vox_list)).
#' @importFrom dplyr do rowwise
feature_rsa_iterate <- function(mod_spec, vox_list, ids=1:length(vox_list)) {
  assert_that(length(ids) == length(vox_list), msg=paste("length(ids) = ", length(ids), "::", "length(vox_list) =", length(vox_list)))
  sframe <- get_samples(mod_spec$dataset, vox_list)
  
  ## iterate over searchlights using parallel futures
  sf <- sframe %>% dplyr::mutate(rnum=ids) 
  fut_feature_rsa(mod_spec,sf, regtype, distmethod)
}


#' @keywords internal
fut_feature_rsa <- function(mod_spec, sf) {
  #mod_spec$dataset <- NULL
  gc()
  
  #browser()
  sf %>% furrr::future_pmap(function(sample, rnum, .id) {
    ## extract_roi?
    do_feature_rsa(as_roi(sample, mod_spec$dataset), mod_spec, rnum)
  }, .options = furrr::furrr_options(seed = T)) %>% dplyr::bind_rows()
  
}

  