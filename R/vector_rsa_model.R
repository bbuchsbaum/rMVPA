#' Construct a design for a vectorized RSA model
#'
#' This function constructs a design for an RSA model using a single distance matrix, labels, and blocks.
#'
#' @param D A distance matrix with row.names indicating labels.
#' @param labels character vector of labels corresponding to rows in another dataset X.
#' @param block_var vector indicating the block (strata) each label belongs to.
#' @return A list containing the elements of the RSA design, class attributes "vector_rsa_design" and "list".
#' @details The function verifies that all labels in the 'labels' are present in the row.names of 'D'.
#' It also sets up for precomputing cross-block pairs if needed.
#' @export
vector_rsa_design <- function(D, labels, block_var, dist_to_sim = function(D) 1-D) {
  # Verify that all labels are present in row.names of D
  assertthat::assert_that(all(labels %in% rownames(D)), 
                          msg = "All labels must be present in the row.names of D.")
  
  S <- dist_to_sim(D)
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
    S=S,
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
  Sexpanded <- design$S[match(design$labels, rownames(design$S)), match(design$labels, rownames(design$S))]
  
  # Precompute cross-block pairs with indices for fast access
  label_indices <- match(design$labels, rownames(design$S))
  
  cross_block_data <- lapply(sort(unique(design$block)), function(b) {
    inds <- which(design$block != b)
    list(other_labels = design$labels[inds], indices=inds, block=b)
  })
  
  names(cross_block_data) <- unique(design$block)
  
  # Bundle precomputed data
  precomputed <- list(
    Sexpanded=Sexpanded,
    cross_block_data = cross_block_data
  )
  
  return(precomputed)
}

#' Create a vectorized RSA model
#'
#' This function integrates a vector_rsa_design and precomputes data to create a vectorized RSA model.
#'
#' @param design A vector_rsa_design object.
#' @return An object representing the RSA model, which includes both the design and precomputed data.
#' @details Integrates RSA design and precomputed data to facilitate efficient RSA computations. 
#' It directly incorporates precomputations for cross-block comparisons.
#' @export
vector_rsa_model <- function(dataset, design) {
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
    cross_block_data = cross_block_data
  )
  
  class(model) <- c("vector_rsa_model", "list")
  return(model)
}


#' Train a vector RSA model
#'
#' @param obj An object of class \code{vector_rsa_model}.
#' @param X The data matrix where rows correspond to trials.
#' @param ... addiitonal arguments
#' @export
train_model.vector_rsa_model <- function(obj, roi, ...) {
  # Compute trial scores using the similarity matrix instead of the distance matrix
  scores <- compute_trial_scores(obj, roi, ...)
  return(scores)
}

#' Perform vector RSA for a subset of data
#'
#' @param roi A subset of data, usually representing one ROI or one trial block.
#' @param mod_spec The RSA model specification.
#' @param method Method for computing similarities.
#' @return A tibble with RSA results and potentially error information.
#' @noRd
do_vector_rsa <- function(roi, mod_spec, rnum, method = "pearson") {
  xtrain <- tibble::as_tibble(neuroim2::values(roi$train_roi), .name_repair=.name_repair)
  ind <- indices(roi$train_roi)
  tryCatch({
    scores <- train_model(mod_spec, xtrain, method = method)
    tibble::tibble(result = list(NULL), indices=list(ind), performance=list(scores), id=rnum, error = FALSE, error_message = "~")
  }, error = function(e) {
    tibble::tibble(result = list(NULL), indices=list(ind), performance=list(NULL), id=rnum, error = TRUE, error_message = e$message)
  })
}



#' Iterate over data sets applying the vector RSA model
#'
#' @param mod_spec The model specification.
#' @param vox_list A list of voxel sets to analyze.
#' @param ids Identifiers for each data set.
#' @param permute Logical indicating whether to permute the labels
#' @noRd
vector_rsa_iterate <- function(mod_spec, vox_list, ids = seq_along(vox_list),  permute=FALSE) {
  # Ensure IDs match the number of data sets
  if (length(ids) != length(vox_list)) {
    stop("Length of ids must match the number of data sets.")
  }
 
  assert_that(length(ids) == length(vox_list), msg=paste("length(ids) = ", length(ids), "::", "length(vox_list) =", length(vox_list)))
  sframe <- get_samples(mod_spec$dataset, vox_list)
  
  ## iterate over searchlights using parallel futures
  sf <- sframe %>% dplyr::mutate(rnum=ids) 
  
  fut_vector_rsa(mod_spec,sf, method)
 
}


#' Apply the RSA model in parallel using futures
#'
#' @param mod_spec The model specification.
#' @param sf A tibble containing the data sets and their identifiers.
#' @param method Method for computing similarities.
#' @return A combined result of all RSA analyses.
fut_vector_rsa <- function(mod_spec, sf, method) {
  gc()
  sf %>% furrr::future_pmap(function(sample, rnum, .id) {
    do_vector_rsa(as_roi(sample, mod_spec$dataset), mod_spec, rnum, method=method)
  }, .options = furrr::furrr_options(seed = T)) %>% dplyr::bind_rows()
   
}


compute_trial_scores <- function(obj, X, method = c("pearson", "spearman")) {
  method <- match.arg(method)
  
  # Retrieve precomputed data
  precomputed = obj$design$model_mat
  similarity_matrix = precomputed$Sexpanded  # Unexpanded similarity matrix
  cross_block_data = precomputed$cross_block_data
  X <- as.matrix(X)
  # Calculate similarity vectors using precomputed data
  
  ## we could return two lists: data_vec, ref_vec
  ## then use proxy::simil(x,y, pairwise=TRUE, by_rows=TRUE)
  scores <- lapply(seq_len(nrow(X)), function(i) {
    block_id = obj$design$block[i]
    valid_indices = cross_block_data[[as.character(block_id)]]$indices
    
    # Compute first-order similarities
    data_vec = cor(X[i,], t(X[valid_indices,]), method=method)[1,]
    
    # Map indices to pull the correct rows from the compact similarity matrix
    #ref_indices = match(obj$design$labels[valid_indices], rownames(similarity_matrix))
    #ref_vec = similarity_matrix[match(obj$design$labels[i], rownames(similarity_matrix)), ref_indices]
    ref_vec <- similarity_matrix[valid_indices, i]
    
    # Compute second-order similarity score
    second_order_sim = cor(data_vec, ref_vec, method = method)
    return(second_order_sim)
  })
  
  return(unlist(scores))
}

#' Compare trial similarity vectors to the reference
#'
#' @param sim_vector The similarity vector for a given trial.
#' @param reference_vector The corresponding reference vector from the expanded D matrix.
#' @return The similarity score.
#' @noRd
compare_to_reference <- function(sim_vector, reference_vector) {
  cor(sim_vector, reference_vector, method = "pearson")
}

#' Expand the D matrix to match trials in X
#'
#' @param D The reference distance matrix.
#' @param labels Labels for each row in X.
#' @return An expanded version of D matching the order of labels in X.
#' @noRd
expand_D <- function(D, labels) {
  indices <- match(labels, rownames(D))
  return(D[indices, indices])
}










