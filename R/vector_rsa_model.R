#' Construct a design for a vectorized RSA model
#'
#' This function constructs a design for an RSA model using a single distance matrix, labels, and blocks.
#'
#' @param D A representational dissimilarity matrix with row.names indicating labels.
#' @param labels character vector of labels corresponding to rows in another dataset X.
#' @param block_var vector indicating the block (strata) each label belongs to.
#' @return A list containing the elements of the RSA design, class attributes "vector_rsa_design" and "list".
#' @details The function verifies that all labels in the 'labels' are present in the row.names of 'D'.
#' It also sets up for precomputing cross-block pairs if needed.
#' @export
vector_rsa_design <- function(D, labels, block_var) {
  # Verify that all labels are present in row.names of D
  assertthat::assert_that(all(labels %in% rownames(D)), 
                          msg = "All labels must be present in the row.names of D.")
  
  # S <- dist_to_sim(D)
  # Set up data and preprocessing steps if necessary
  # Here, the design is simple and only needs to handle labels and block structure.
  # Further preprocessing for cross-block pairs might be handled in vector_rsa_model_mat
  
  # Process block_var if provided
  block_var <- if (!is.null(block_var)) {
    parse_variable(block_var, data)
  }
  
  assertthat::assert_that(length(labels) == length(block_var), msg = "Length of labels and block_var must match.")
  
  design <- list(
    D = D,
    #S=S,
    labels = labels,
    block = block_var
  )
  
  # Add model matrix to the design list
  mmat <- vector_rsa_model_mat(design)
  design$model_mat <- mmat
  
  class(design) <- c("vector_rsa_design", "list")
  
  
  return(design)
}

#' @noRd
vector_rsa_model_mat <- function(design) {
  # Convert distance matrix to similarity matrix
  #similarity_matrix <- 1 - as.matrix(design$D)
  
  # Ensure the matrix is expanded to match all instances in X
  Dexpanded <- design$D[match(design$labels, rownames(design$D)), match(design$labels, rownames(design$D))]
  
  # Precompute cross-block pairs with indices for fast access
  label_indices <- match(design$labels, rownames(design$D))
  
  cross_block_data <- lapply(sort(unique(design$block)), function(b) {
    inds <- which(design$block != b)
    list(other_labels = design$labels[inds], indices=inds, block=b)
  })
  
  names(cross_block_data) <- unique(design$block)
  
  # Bundle precomputed data
  precomputed <- list(
    Dexpanded=Dexpanded,
    cross_block_data = cross_block_data
  )
  
  return(precomputed)
}

#' Create a vectorized RSA model
#'
#' This function integrates a vector_rsa_design and precomputes data to create a vectorized RSA model.
#' 
#' @param dataset An \code{mvpa_dataset}  object.
#' @param design A vector_rsa_design object.
#' @param distfun an object of \code{distfun} type for computation of pairwise dissimilarities among image rows
#' @param rsa_simfun a character string specifying the similarity function to use for RSA
#' @return An object representing the RSA model, which includes both the design and precomputed data.
#' @details Integrates RSA design and precomputed data to facilitate efficient RSA computations. 
#' It directly incorporates precomputations for cross-block comparisons.
#' @export
vector_rsa_model <- function(dataset, design, distfun=cordist(), rsa_simfun=c("pearson", "spearman")) {
  rsa_simfun = match.arg(rsa_simfun)
  assert_that(inherits(dataset, "mvpa_dataset"))
  # Verify that design is correctly formatted
  assertthat::assert_that(inherits(design, "vector_rsa_design"),
                          msg = "Input must be a 'vector_rsa_design' object.")
  
  # Precompute cross-block pairs
  label_indices <- match(design$labels, rownames(design$D))
  
  cross_block_data <- lapply(unique(design$block), function(b) {
    inds <- which(design$block != b)
    list(other_labels = design$labels[inds], indices=inds) #indices = label_indices[inds])
  })
  
  names(cross_block_data) <- unique(design$block)
  
  # Store precomputed data in the model
  model <- list(
    dataset = dataset,
    design = design,
    distfun=distfun,
    rsa_simfun=rsa_simfun,
    cross_block_data = cross_block_data
  )
  
  create_model_spec("vector_rsa_model", dataset, design,distfun=distfun,
                    rsa_simfun=rsa_simfun,
                    cross_block_data = cross_block_data)
  

}


#' Train a vector RSA model
#'
#' @param obj An object of class \code{vector_rsa_model}.
#' @param X The data matrix where rows correspond to trials.
#' @param train_dat the training data
#' @param indices the spatial indices of the training data
#' @param ... addiitonal arguments
#' @export
train_model.vector_rsa_model <- function(obj, train_dat, indices, ...) {
  # Compute trial scores using the similarity matrix instead of the distance matrix
  scores <- compute_trial_scores(obj, train_dat)
  return(scores)
}




#' @noRd
#' @keywords internal
compute_trial_scores <- function(obj, X) {
  
  # Retrieve precomputed data
  precomputed = obj$design$model_mat
  dissimilarity_matrix = precomputed$Dexpanded  # Unexpanded similarity matrix
  cross_block_data = precomputed$cross_block_data
  X <- as.matrix(X)
  # Calculate similarity vectors using precomputed data
  #browser()
  second_order_similarity(obj$distfun, X, dissimilarity_matrix, obj$design$block, obj$rsa_simfun)
  
}











