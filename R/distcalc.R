#' @export
#' @noRd
pairwise_dist.default <- function(obj, X, ...) {
  stop("pairwise_dist not implemented for objects of class ", class(X)[1])
}


#' Create a Distance Function Object
#'
#' Constructs a generic distance function object, storing:
#' \itemize{
#'   \item \code{name}: The name (method) of the distance (e.g., "euclidean", "cordist", "mahalanobis").
#'   \item \code{labels}: (Optional) a vector of labels associated with the data rows.
#'   \item \code{...}: Additional parameters relevant to the specific distance method 
#'         (e.g., correlation method for \code{cordist}, number of components for \code{pcadist}, etc.).
#' }
#'
#' This object is used by \code{pairwise_dist()} to compute an \strong{N x N} matrix of pairwise distances 
#' between rows of a data matrix.
#'
#' @param name A character string specifying the distance method (e.g., "euclidean", "cordist").
#' @param labels A vector of row labels (optional), primarily for informational/reference purposes.
#' @param ... Additional parameters for the distance method (e.g. `method="pearson"` for correlation, 
#'            or \code{whiten=TRUE} for PCA-based distances).
#'
#' @return An S3 object with class \code{c(name, "distfun")} that can be passed to \code{pairwise_dist()}.
#'
#' @details
#' The distance function object itself does \emph{not} exclude same-block comparisons or reorder rows by label.
#' Those tasks (if needed) are handled downstream (for example, in \code{second_order_similarity}).
#'
#' @examples
#' # Create a Euclidean distance function object
#' dist_obj_euc <- create_dist("euclidean")
#'
#' # Create a correlation distance function object with a specified correlation method
#' dist_obj_cor <- create_dist("cordist", method="spearman")
#'
#' @export
create_dist <- function(name, labels=NULL, ...) {
  structure(
    list(name = name, labels = labels, ...), 
    class = c(name, "distfun")
  )
}


#' Distance Function Constructors
#'
#' These convenience functions build specific types of distance function objects via \code{\link{create_dist}}.
#' Each returns an S3 object inheriting from \code{c("<method>", "distfun")}.
#'
#' @param labels Optional vector of row labels (not directly used in distance calculation).
#' @param method For \code{cordist}, the correlation method: "pearson" or "spearman".
#' @param ncomp For \code{pcadist}, the number of components (or a function threshold).
#' @param whiten For \code{pcadist}, whether to whiten principal components (logical).
#' @param threshfun For \code{pcadist}, an optional function that determines how many PCs to retain 
#'                  based on \code{pres$sdev^2}.
#' @param dist_method For \code{pcadist}, the base distance measure in PC space ("euclidean", "manhattan", or "cosine").
#'
#' @return An S3 object with class \code{c("<method>", "distfun")}.
#'
#' @details
#' - \code{cordist(labels, method="pearson")} → correlation-based distance.  
#' - \code{mahadist(labels)} → Mahalanobis distance.  
#' - \code{eucdist(labels)} → Euclidean distance.  
#' - \code{robustmahadist(labels)} → Mahalanobis distance using robust covariance.  
#' - \code{pcadist(labels, ...)} → distance in reduced PCA space.
#'
#' @examples
#' dist_obj_1 <- cordist(method="pearson")
#' dist_obj_2 <- mahadist()
#' dist_obj_3 <- eucdist()
#' dist_obj_4 <- robustmahadist()
#' dist_obj_5 <- pcadist(ncomp=2, dist_method="cosine")
#'
#' @seealso \code{\link{create_dist}} for the underlying constructor.
#'
#' @rdname distance-constructors
#' @export
cordist <- function(labels=NULL, method=c("pearson", "spearman")) {
  method <- match.arg(method)
  create_dist(name="cordist", labels=labels, method=method)
}

#' @rdname distance-constructors
#' @export
mahadist <- function(labels=NULL) {
  create_dist("mahalanobis", labels)
}

#' @rdname distance-constructors
#' @export
eucdist <- function(labels=NULL) {
  create_dist("euclidean", labels)
}

#' @rdname distance-constructors
#' @export
robustmahadist <- function(labels=NULL) {
  create_dist("robustmahadist", labels)
}

#' @rdname distance-constructors
#' @export
pcadist <- function(labels=NULL, ncomp=2, whiten=TRUE, threshfun=NULL,
                    dist_method=c("euclidean", "manhattan", "cosine")) {
  dist_method <- match.arg(dist_method)
  if (is.null(threshfun)) {
    # By default, always use the user-specified 'ncomp'
    tfun <- function(x) ncomp
  } else {
    stopifnot(is.function(threshfun))
    tfun <- threshfun
  }
  create_dist("pcadist", labels, whiten=whiten, threshfun=tfun, dist_method=dist_method)
}


#' Compute Pairwise Correlation Distances
#'
#' Computes a full NxN matrix of correlation-based distances: \code{1 - cor(t(X), method=obj$method)}.
#' \strong{No block-based exclusion is performed here.}
#'
#' @param obj A distance function object of class \code{c("cordist", "distfun")}.
#' @param X A numeric matrix (rows = observations, columns = variables).
#'
#' @return An \strong{N x N numeric matrix} of pairwise distances.
#'
#' @details
#' This function calculates correlation distances among the rows of \code{X}. 
#' If you have a block variable and wish to exclude same-block comparisons, 
#' handle that \emph{after} obtaining this full matrix (e.g., in \code{second_order_similarity}).
#'
#' @examples
#' X <- matrix(rnorm(100), 10, 10)
#' dist_obj <- cordist(method="pearson")
#' dist_matrix <- pairwise_dist(dist_obj, X)
#'
#' @export
#' @noRd 
pairwise_dist.cordist <- function(obj, X,...) {
  1 - cor(t(X), method=obj$method)
}


#' Compute Pairwise Euclidean Distances
#'
#' Returns a full NxN matrix of Euclidean distances among the rows of \code{X}.
#' \strong{No block-based exclusion is performed.}
#'
#' @param obj A distance function object of class \code{c("euclidean", "distfun")}.
#' @param X A numeric matrix (rows = observations, columns = variables).
#'
#' @return An \strong{N x N numeric matrix} of pairwise Euclidean distances.
#'
#' @details
#' This function simply calls \code{dist(X)} internally and converts it to a matrix via \code{as.matrix()}.
#'
#' @examples
#' X <- matrix(rnorm(100), 10, 10)
#' dist_obj <- eucdist()
#' dist_matrix <- pairwise_dist(dist_obj, X)
#'
#' @export
#' @noRd
pairwise_dist.euclidean <- function(obj, X,...) {
  as.matrix(dist(X, method="euclidean"))
}


#' Compute Pairwise Mahalanobis Distances
#'
#' Returns a full NxN matrix of Mahalanobis distances among rows of \code{X}, using a shrunken inverse covariance.
#' \strong{No block-based exclusion is performed.}
#'
#' @param obj A distance function object of class \code{c("mahalanobis", "distfun")}. 
#' @param X A numeric matrix (rows = observations, columns = variables).
#'
#' @return An \strong{N x N numeric matrix} of pairwise Mahalanobis distances.
#'
#' @details
#' Uses \code{corpcor::invcov.shrink} on \code{X} to estimate the inverse covariance matrix 
#' and then computes \code{mahalanobis(...)} for each row vs. each other row.
#'
#' @examples
#' X <- matrix(rnorm(100), 10, 10)
#' dist_obj <- mahadist()
#' dist_matrix <- pairwise_dist(dist_obj, X)
#'
#' @export
#' @importFrom corpcor invcov.shrink
#' @importFrom stats mahalanobis
#' @noRd
pairwise_dist.mahalanobis <- function(obj, X,...) {
  inv_cov <- corpcor::invcov.shrink(X)
  n <- nrow(X)
  
  dist_matrix_sq <- matrix(0, n, n)
  for (i in seq_len(n)) {
    dist_matrix_sq[i, ] <- mahalanobis(X, center = X[i, ], cov = inv_cov, inverted = TRUE)
  }
  
  sqrt(dist_matrix_sq)
}


#' Compute Pairwise PCA-Based Distances
#'
#' Returns a full NxN matrix of distances in a PCA-reduced subspace, with an optional whitening step.
#' \strong{No block-based exclusion is performed.}
#'
#' @param obj A distance function object of class \code{c("pcadist", "distfun")}.
#' @param X A numeric matrix.
#'
#' @return An \strong{N x N numeric matrix} of pairwise distances in the reduced PCA space.
#'
#' @details
#' 1. Performs \code{prcomp(X)} (centered, scaled=TRUE).
#' 2. Determines the number of components via \code{obj$threshfun(...)}.
#' 3. Optionally whitens (divide each principal component by its standard deviation).
#' 4. Computes pairwise distances on the reduced data using \code{obj$dist_method}.
#'
#' @examples
#' X <- matrix(rnorm(100), 10, 10)
#' dist_obj <- pcadist(ncomp=3, dist_method="cosine")
#' dist_matrix <- pairwise_dist(dist_obj, X)
#'
#' @export
#' @importFrom stats prcomp dist
#' @noRd
pairwise_dist.pcadist <- function(obj, X,...) {
  pres <- prcomp(X, center = TRUE, scale. = TRUE)
  ncomp <- obj$threshfun(pres$sdev^2)
  if (ncomp < 1) {
    ncomp <- 1
    warning("Number of components set to 1, as threshold function returned < 1.")
  }
  
  # Extract PC space
  pc_space <- pres$x[, seq_len(ncomp), drop=FALSE]
  
  # Optional whitening
  if (obj$whiten) {
    pc_space <- pc_space %*% diag(1 / pres$sdev[seq_len(ncomp)], ncomp, ncomp)
  }
  
  # Distances
  if (obj$dist_method %in% c("euclidean", "manhattan")) {
    as.matrix(dist(pc_space, method = obj$dist_method))
  } else if (obj$dist_method == "cosine") {
    # proxy::dist for 'cosine' distance
    as.matrix(proxy::dist(pc_space, method = "cosine"))
  }
}


#' Compute Pairwise Robust Mahalanobis Distances
#'
#' Returns a full NxN matrix of robust Mahalanobis distances, using a robust covariance estimator.
#' \strong{No block-based exclusion is performed.}
#'
#' @param obj A distance function object of class \code{c("robustmahadist", "distfun")}.
#' @param X A numeric matrix.
#'
#' @return An \strong{N x N numeric matrix} of pairwise robust Mahalanobis distances.
#'
#' @details
#' - Estimates a robust covariance with \code{robustcov::covGK(X)} (make sure the \code{robustcov} package is installed).
#' - Then calls \code{corpcor::invcov.shrink} to get an inverse covariance estimate.
#' - Finally, loops over row pairs to compute \code{(x_i - x_j) * inv_cov * (x_i - x_j)^T}.
#'
#' @examples
#' \dontrun{
#'   library(robustcov)
#'   X <- matrix(rnorm(100), 10, 10)
#'   dist_obj <- robustmahadist()
#'   dist_matrix <- pairwise_dist(dist_obj, X)
#' }
#'
#' @export
#' @noRd
pairwise_dist.robustmahadist <- function(obj, X,...) {
  robust_cov <- robustcov::covGK(X)
  inv_cov <- corpcor::invcov.shrink(robust_cov)
  
  n <- nrow(X)
  dist_matrix <- matrix(0, n, n)
  
  for (i in seq_len(n-1)) {
    for (j in seq((i+1), n)) {
      diff <- X[i, ] - X[j, ]
      dist_val <- sqrt(t(diff) %*% inv_cov %*% diff)
      dist_matrix[i, j] <- dist_val
      dist_matrix[j, i] <- dist_val
    }
  }
  dist_matrix
}






