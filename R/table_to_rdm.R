#' Convert Similarity Table to RDM
#'
#' Converts a similarity lookup table to a Representational Dissimilarity Matrix (RDM)
#' for use in RSA analyses. The function maps label pairs through the similarity table
#' and creates a complete RDM, using default values for missing pairs.
#'
#' @param similarity_table A data frame containing pairwise similarity values with columns
#'   for label1, label2, and similarity values
#' @param labels A character vector of labels for which to create the RDM. The output
#'   RDM will have these labels in this order
#' @param label1_col Character string specifying the column name for the first label 
#'   (default: "label1")
#' @param label2_col Character string specifying the column name for the second label
#'   (default: "label2")
#' @param similarity_col Character string specifying the column name for similarity values
#'   (default: "similarity"). Values should be between 0 and 1, where 1 = identical
#' @param default_similarity Numeric value to use for label pairs not found in the table
#'   (default: 0, meaning maximum dissimilarity)
#' @param symmetric Logical indicating whether the similarity table should be treated
#'   as symmetric (default: TRUE). If TRUE, (A,B) = (B,A)
#' @param self_similarity Numeric value for diagonal elements (self-similarity).
#'   Default is 1 (perfect similarity). Set to NA to look up in table
#' @param as_dist Logical; if TRUE return a dist object, otherwise return matrix
#'   (default: TRUE)
#'
#' @return A dist object or symmetric matrix representing dissimilarities (RDM).
#'   Values are computed as 1 - similarity, so that 0 = identical and 1 = maximally different
#'
#' @details
#' This function is useful for hypothesis-driven RSA where you have theoretical
#' predictions about the similarity structure between conditions. The similarity
#' table can be sparse (not all pairs need to be specified), and missing pairs
#' will use the default similarity value.
#' 
#' The function converts similarities to dissimilarities using the formula:
#' dissimilarity = 1 - similarity
#' 
#' @examples
#' # Create a similarity table based on theoretical predictions
#' sim_table <- data.frame(
#'   label1 = c("cat", "cat", "dog", "car"),
#'   label2 = c("dog", "bird", "bird", "plane"),
#'   similarity = c(0.7, 0.5, 0.6, 0.8)
#' )
#' 
#' # Create RDM for specific conditions
#' conditions <- c("cat", "dog", "bird", "car", "plane", "boat")
#' rdm <- table_to_rdm(sim_table, conditions)
#' 
#' # Get as matrix instead of dist
#' rdm_matrix <- table_to_rdm(sim_table, conditions, as_dist = FALSE)
#' 
#' # Use in RSA design
#' \dontrun{
#' rsa_des <- rsa_design(~ theoretical_rdm,
#'                      data = list(theoretical_rdm = rdm))
#' }
#' 
#' @export
#' @importFrom stats as.dist
table_to_rdm <- function(similarity_table,
                        labels,
                        label1_col = "label1",
                        label2_col = "label2",
                        similarity_col = "similarity",
                        default_similarity = 0,
                        symmetric = TRUE,
                        self_similarity = 1,
                        as_dist = TRUE) {
  
  # Input validation
  if (!is.data.frame(similarity_table)) {
    stop("similarity_table must be a data frame")
  }
  
  if (!all(c(label1_col, label2_col, similarity_col) %in% names(similarity_table))) {
    stop(sprintf("similarity_table must contain columns: %s, %s, %s",
                 label1_col, label2_col, similarity_col))
  }
  
  # Extract columns
  label1_vals <- as.character(similarity_table[[label1_col]])
  label2_vals <- as.character(similarity_table[[label2_col]])
  similarity_vals <- as.numeric(similarity_table[[similarity_col]])
  
  # Check similarity values are in valid range
  if (any(!is.na(similarity_vals) & (similarity_vals < 0 | similarity_vals > 1))) {
    warning("Similarity values should be between 0 and 1. Values outside this range will be clipped.")
    similarity_vals <- pmax(0, pmin(1, similarity_vals))
  }
  
  # Check for duplicate entries
  if (symmetric) {
    # For symmetric tables, (A,B) and (B,A) are the same
    pair_keys <- paste(pmin(label1_vals, label2_vals), 
                      pmax(label1_vals, label2_vals), 
                      sep = "___")
  } else {
    # For asymmetric tables, (A,B) and (B,A) are different
    pair_keys <- paste(label1_vals, label2_vals, sep = "___")
  }
  
  if (any(duplicated(pair_keys))) {
    warning("Duplicate label pairs found in similarity_table. Using first occurrence.")
    # Keep only first occurrence
    keep_idx <- !duplicated(pair_keys)
    label1_vals <- label1_vals[keep_idx]
    label2_vals <- label2_vals[keep_idx]
    similarity_vals <- similarity_vals[keep_idx]
  }
  
  # Create lookup structure
  similarity_lookup <- list()
  for (i in seq_along(label1_vals)) {
    key1 <- paste(label1_vals[i], label2_vals[i], sep = "___")
    similarity_lookup[[key1]] <- similarity_vals[i]
    
    if (symmetric && label1_vals[i] != label2_vals[i]) {
      # Add reverse direction for symmetric case
      key2 <- paste(label2_vals[i], label1_vals[i], sep = "___")
      similarity_lookup[[key2]] <- similarity_vals[i]
    }
  }
  
  # Initialize RDM matrix
  n <- length(labels)
  rdm_matrix <- matrix(NA, nrow = n, ncol = n)
  rownames(rdm_matrix) <- labels
  colnames(rdm_matrix) <- labels
  
  # Fill the matrix
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      if (i == j) {
        # Diagonal elements (self-similarity)
        if (is.na(self_similarity)) {
          # Look up self-similarity in table
          key <- paste(labels[i], labels[j], sep = "___")
          sim_val <- similarity_lookup[[key]]
          if (is.null(sim_val)) {
            rdm_matrix[i, j] <- 0  # Default dissimilarity for diagonal
          } else {
            rdm_matrix[i, j] <- 1 - sim_val
          }
        } else {
          # Use specified self-similarity
          rdm_matrix[i, j] <- 1 - self_similarity
        }
      } else {
        # Off-diagonal elements
        key <- paste(labels[i], labels[j], sep = "___")
        sim_val <- similarity_lookup[[key]]
        
        if (is.null(sim_val)) {
          # Use default similarity for missing pairs
          rdm_matrix[i, j] <- 1 - default_similarity
        } else {
          # Convert similarity to dissimilarity
          rdm_matrix[i, j] <- 1 - sim_val
        }
      }
    }
  }
  
  # Ensure matrix is symmetric if specified
  if (symmetric) {
    # Average the upper and lower triangles to ensure perfect symmetry
    rdm_matrix <- (rdm_matrix + t(rdm_matrix)) / 2
  }
  
  # Return in requested format
  if (as_dist) {
    as.dist(rdm_matrix)
  } else {
    rdm_matrix
  }
}


#' Create Hypothesis RDM from Category Structure
#'
#' Creates an RDM based on category membership, where items in the same category
#' have high similarity and items in different categories have low similarity.
#'
#' @param categories A named vector or factor where names/levels are labels and values
#'   are category assignments
#' @param within_category_sim Similarity value for items within the same category
#'   (default: 0.8)
#' @param between_category_sim Similarity value for items in different categories
#'   (default: 0.2)
#' @param as_dist Logical; if TRUE return a dist object, otherwise return matrix
#'   (default: TRUE)
#'
#' @return A dist object or matrix representing the category-based RDM
#'
#' @examples
#' # Create category structure
#' categories <- c(cat = "animal", dog = "animal", bird = "animal",
#'                car = "vehicle", plane = "vehicle", boat = "vehicle")
#' 
#' # Create category-based RDM
#' rdm <- category_rdm(categories)
#' 
#' # Custom similarity values
#' rdm <- category_rdm(categories, 
#'                    within_category_sim = 0.9,
#'                    between_category_sim = 0.1)
#' 
#' @export
category_rdm <- function(categories,
                        within_category_sim = 0.8,
                        between_category_sim = 0.2,
                        as_dist = TRUE) {
  
  if (is.factor(categories)) {
    # For factors, we need to handle this differently
    # The factor levels are the unique values, not the labels
    # We need a mapping of items to categories
    stop("For factors, please provide a named vector where names are items and values are categories")
  } else if (!is.null(names(categories))) {
    # Named vector: names are labels, values are categories
    labels <- names(categories)
    cats <- as.character(categories)
  } else {
    stop("categories must be a named vector where names are item labels and values are category assignments")
  }
  
  n <- length(labels)
  rdm_matrix <- matrix(0, nrow = n, ncol = n)
  rownames(rdm_matrix) <- labels
  colnames(rdm_matrix) <- labels
  
  # Fill matrix based on category membership
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      if (i == j) {
        rdm_matrix[i, j] <- 0  # Self-dissimilarity is 0
      } else if (cats[i] == cats[j]) {
        rdm_matrix[i, j] <- 1 - within_category_sim
      } else {
        rdm_matrix[i, j] <- 1 - between_category_sim
      }
    }
  }
  
  if (as_dist) {
    as.dist(rdm_matrix)
  } else {
    rdm_matrix
  }
}