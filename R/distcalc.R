
# Default method
pairwise_dist.default <- function(X, dist_obj) {
  stop("pairwise_dist not implemented for objects of class ", class(dist_obj)[1])
}


#' Create Distance Function Object
#'
#' This function constructs an object representing a distance function, 
#' which can be used with generic functions like `pairwise_dist` for computing distances. 
#' The object stores the method of distance calculation, labels associated with data points, 
#' and any additional parameters that may affect the distance computation.
#'
#' @param name A character string specifying the method of distance computation. 
#' This method name is used to dispatch the appropriate distance calculation function. 
#' Common methods might include "euclidean", "manhattan", "mahalanobis", etc.
#'
#' @param labels A vector of labels or identifiers associated with the rows of the data matrix. 
#' These labels are important for reference in distance computations, particularly when adjustments or 
#' restrictions based on groupings or identifiers are needed.
#'
#' @param ... Additional parameters relevant to the specific distance method. 
#' These could include tuning parameters like `lambda` for shrinkage in covariance estimation 
#' or parameters controlling the behavior of the distance computation.
#'
#' @return Returns an object of class `distfun` and the specific method class 
#' (as specified by the `method` parameter). This object encapsulates all information 
#' necessary to compute distances between data points according to the specified method 
#' and additional parameters.
#'
#' @details
#' The `create_dist` function enables the flexible creation of distance function objects. 
#' By specifying a method and associated parameters, users can customize the behavior of 
#' distance calculations. This functionality is especially useful in statistical and 
#' machine learning applications where different distance metrics can have significant 
#' impacts on the results.
#'
#' @examples
#' # Create a Euclidean distance function object
#' dist_obj_euc <- create_dist("euclidean", labels = c("A", "B", "C", "D"))
#'
#' # Create a Mahalanobis distance function object with additional parameters
#' dist_obj_maha <- create_dist("mahalanobis", labels = c("A", "B", "C", "D"))
#'
#' @export
create_dist <- function(name, labels, ...) {
  structure(list(name = name, labels = labels, ...), class = c(name, "distfun"))
}


#' Distance Function Constructors
#'
#' These functions provide convenient constructors for various types of distance functions.
#' Each constructor function initializes a distance function object for use with distance
#' computation functions, specifying the method and any necessary labels.
#'
#' @param labels A vector of labels associated with the data points.
#' @param method The method of distance computation, applicable for `cordist`.
#'
#' @return An object of class `distfun` with a specific method subclass, encapsulating
#'         all information necessary for computing distances according to the specified method.
#'
#' @details
#' The constructors allow for the specification of distance calculation methods and associated labels:
#' - `cordist` creates a correlation distance function.
#' - `mahadist` creates a Mahalanobis distance function.
#' - `eucdist` creates a Euclidean distance function.
#' - `robustmahadist` creates a robust version of the Mahalanobis distance function.
#'
#' @examples
#' dist_obj_1 <- cordist(labels = c("A", "B", "C"), method = "pearson")
#' dist_obj_2 <- mahadist(labels = c("A", "B", "C"))
#' dist_obj_3 <- eucdist(labels = c("A", "B", "C"))
#' dist_obj_4 <- robustmahadist(labels = c("A", "B", "C"))
#'
#' @seealso \code{\link{create_dist}} for the underlying constructor used by these functions.
#'
#' @rdname distance-constructors
#' @export
#' @keywords methods
cordist <- function(labels=NULL, method=c("pearson", "spearman")) {
  method=match.arg(method)
  create_dist(name="cordist", labels=labels, method=method)
}

# Example usage for Mahalanobis distance with labels
#' @rdname distance-constructors
mahadist <- function(labels=NULL) {
  create_dist("mahalanobis", labels)
}

#' @rdname distance-constructors
eucdist <- function(labels=NULL) {
  create_dist("euclidean", labels)
}

#' @rdname distance-constructors
robustmahadist <- function(labels=NULL) {
  create_dist("robustmahadist", labels)
}



#' Compute Pairwise Correlation Distances
#'
#' This method computes the pairwise correlation distances for a matrix `X`, excluding
#' comparisons within the same block as specified by the `dist_obj$block`.
#'
#' @param dist_obj A list containing the method ("correlation") and a block vector to specify
#' which rows in `X` should not be compared to avoid within-block correlation.
#' @param X Numeric matrix where rows represent observations and columns represent variables.
#'
#' @return An object of class `dist` containing the computed correlation distances.
#'
#' @examples
#' X <- matrix(rnorm(100), 10, 10)
#' block <- rep(1:2, each=5)
#' dist_obj <- list(method = "pearson", block = block)
#' dist_matrix <- pairwise_dist.correlation(dist_obj, X)
#' 
#' @export
pairwise_dist.cordist <- function(obj, X) {
  1 - cor(t(X), method=obj$method)
  # block <- dist_obj$block
  # n <- nrow(X)
  # dist_matrix <- matrix(0, n, n)  # initialize with zeros
  # for (i in seq_len(n)) {
  #   valid_indices <- which(block != block[i])
  #   dist_matrix[i, valid_indices] <- 1 - cor(X[i, , drop = FALSE], t(X[valid_indices, , drop = FALSE]), method=dist_obj$method)
  # }
  # as.dist(dist_matrix)  # convert to distance object
}


#' Compute Pairwise Euclidean Distances
#'
#' Computes the pairwise Euclidean distances for a matrix `X`.
#'
#' @param dist_obj A list containing possibly additional parameters, currently unused.
#' @param X Numeric matrix where rows represent observations and columns represent variables.
#'
#' @return An object of class `dist` containing the computed Euclidean distances.
#'
#' @examples
#' X <- matrix(rnorm(100), 10, 10)
#' dist_matrix <- pairwise_dist.euclidean(list(), X)
#'
#' @export
pairwise_dist.euclidean <- function(obj, X) {
  # Estimate the inverse of the shrunken covariance matrix
  as.matrix(dist(X))
}


#' Compute Pairwise Mahalanobis Distances
#'
#' Computes the pairwise Mahalanobis distances using an inverse covariance matrix estimated
#' from the data matrix `X` with shrinkage.
#'
#' @param dist_obj A list that might include additional parameters for distance computation, 
#' currently unused.
#' @param X Numeric matrix where rows represent observations and columns represent variables.
#'
#' @return An object of class `dist` containing the computed Mahalanobis distances.
#'
#' @examples
#' X <- matrix(rnorm(100), 10, 10)
#' dist_matrix <- pairwise_dist.mahalanobis(list(), X)
#'
#' @export
pairwise_dist.mahalanobis <- function(obj, X) {
  # Estimate the inverse of the shrunken covariance matrix
  inv_cov <- invcov.shrink(X)
  
  n <- nrow(X)
  
  # Compute the squared Mahalanobis distances using mahalanobis()
  dist_matrix_sq <- matrix(0, n, n)
  for (i in 1:n) {
    dist_matrix_sq[i, ] <- mahalanobis(X, center = X[i, ], cov = inv_cov, inverted = TRUE)
  }
  
  sqrt(dist_matrix_sq) # Computing the square root of the squared distances
}


#' Compute Pairwise Robust Mahalanobis Distances
#'
#' Computes the pairwise Mahalanobis distances using a robustly estimated covariance matrix,
#' which can be more resistant to outliers.
#'
#' @param dist_obj A list that might include additional parameters for distance computation, 
#' currently unused.
#' @param X Numeric matrix where rows represent observations and columns represent variables.
#'
#' @return An object of class `dist` containing the computed robust Mahalanobis distances.
#'
#' @examples
#' X <- matrix(rnorm(100), 10, 10)
#' dist_matrix <- pairwise_dist.robustmahadist(list(), X)
#'
#' @export
pairwise_dist.robustmahadist <- function(obj, X) {
  # Use robust covariance estimation
  robust_cov <- robustcov::covGK(X)
  inv_cov <- corpcor::invcov.shrink(robust_cov)
  
  n <- nrow(X)
  dist_matrix <- matrix(0, n, n)
  
  for (i in 1:(n-1)) {
    for (j in (i + 1):n) {
      diff <- X[i, ] - X[j, ]
      dist_matrix[i, j] <- sqrt(t(diff) %*% inv_cov %*% diff)
      dist_matrix[j, i] <- dist_matrix[i, j]  # Fill lower triangle
    }
  }
  
  sqrt(dist_matrix)
}


#' Compute Second-Order Similarity Scores
#'
#' This function calculates the second order similarity between two similarity vectors
#' derived from a provided distance function applied to matrix X and a reference
#' similarity matrix S. The calculation takes into account a blocking variable to exclude
#' comparisons within the same block.
#'
#' @param dist_fun A distance function object or a character string specifying the 
#' method used for distance computation. This function should be capable of processing
#' the matrix X to produce a distance matrix.
#' @param X A numeric matrix where each row is an observation and columns are features.
#' Distances will be computed pairwise between rows of this matrix.
#' @param D A numeric matrix, typically a predefined dissimilarity matrix that
#' serves as a reference to compare against the computed distances from X.
#' @param block A vector (numeric or factor) indicating the block or group for each row
#' in X and S. Comparisons are only made between elements of different blocks.
#' @param method The method used for computing correlation between similarity vectors.
#' Defaults to "pearson", but "spearman" or "kendall" could also be used.
#'
#' @return A numeric vector of similarity scores, one for each observation in X, 
#' representing the correlation between distance vectors derived from X and the
#' corresponding vectors in S for non-matching blocks.
#'
#' @details
#' The function computes a distance matrix for X using the specified `dist_fun`. It then
#' compares these distances with the entries in S for each observation, excluding
#' comparisons within the same block as defined by the `block` argument. This is useful
#' for evaluating how well the distances within X align with an external similarity
#' standard, adjusting for within-block dependencies.
#'
#' @examples
#' # Assuming X and S are numeric matrices and block is a factor or numeric vector
#' dist_fun <- "euclidean"  # This should be defined or loaded from your package/environment
#' X <- matrix(rnorm(100), ncol=10)
#' D <- matrix(rnorm(100), ncol=10)
#' block <- rep(1:5, each=20)
#' scores <- second_order_similarity(dist_fun, X, D, block, method = "pearson")
#'
#' @export
second_order_similarity <- function(dist_fun, X, D, block, method = c("pearson", "spearman")) {
  method <- match.arg(method)

  # Compute distances using the provided distance function
  distance_matrix = pairwise_dist(dist_fun, X)
  
  # Initialize scores vector
  scores <- numeric(length(block))
 
  # Calculate trial-wise similarity scores, considering valid blocks
  for (i in seq_along(block)) {
    valid_indices <- which(block != block[i])
    if (length(valid_indices) > 0) {
      sim_vector_x = distance_matrix[i, valid_indices]
      sim_vector_s = D[i, valid_indices]
      scores[i] <- if (length(sim_vector_x) > 0 && length(sim_vector_s) > 0) {
        cor(sim_vector_x, sim_vector_s, method = method)
      } else {
        NA  # Handle cases where no valid comparisons
      }
    } else {
      scores[i] <- NA  # Assign NA if no valid indices
    }
  }
  
  return(scores)
}










