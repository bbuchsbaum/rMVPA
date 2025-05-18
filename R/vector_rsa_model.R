#' Construct a design for a vectorized RSA model
#'
#' This function constructs a design for an RSA model using a single distance matrix, labels, and blocks.
#'
#' @param D A representational dissimilarity matrix with row.names indicating labels.
#' @param labels character vector of labels corresponding to rows in another dataset X.
#' @param block_var A vector indicating the block (strata) each label belongs to. Must be the same length as `labels`.
#'
#' @return A list containing the elements of the RSA design, class attributes "vector_rsa_design" and "list".
#' @details
#' The function verifies that all `labels` appear in `rownames(D)` and creates an expanded
#' dissimilarity matrix (`Dexpanded`) matching the order of `labels`.
#' 
#' @export
vector_rsa_design <- function(D, labels, block_var) {
  # Verify that all labels are present in row.names of D
  assertthat::assert_that(
    all(labels %in% rownames(D)),
    msg = "All labels must be present in the row.names of D."
  )
  
  # Fail fast if the number of rows in D is not equal to length(labels)
  #assertthat::assert_that(
  #  nrow(D) == length(labels),
  #  msg = "The number of rows in D must match the number of labels."
  #)
  
  # Ensure length consistency
  assertthat::assert_that(
    length(labels) == length(block_var),
    msg = "Length of labels and block_var must match."
  )
  
  # Build the design object
  design <- list(
    D = D,
    labels = labels,
    block = block_var
  )
  
  # Attach precomputed expansions
  design$model_mat <- vector_rsa_model_mat(design)
  
  class(design) <- c("vector_rsa_design", "list")
  return(design)
}


#' @noRd
vector_rsa_model_mat <- function(design) {
  # Reorder D to match 'labels' in the correct order
  row_idx <- match(design$labels, rownames(design$D))
  Dexpanded <- design$D[row_idx, row_idx, drop=FALSE]
  
  # Return expanded matrix
  list(
    Dexpanded = Dexpanded
  )
}


#' Create a vectorized RSA model
#'
#' This function integrates a vector_rsa_design and an mvpa_dataset to create a vectorized RSA model.
#' 
#' @param dataset An \code{mvpa_dataset} object.
#' @param design A \code{vector_rsa_design} object.
#' @param distfun A \code{distfun} (distance function) for computing pairwise dissimilarities among image rows.
#' @param rsa_simfun A character string specifying the similarity function to use for RSA, 
#'                   one of \code{"pearson"} or \code{"spearman"}.
#' @param nperm Integer, number of permutations for statistical testing (default: 0).
#' @param save_distributions Logical, whether to save full permutation distributions (default: FALSE).
#' @param return_predictions Logical, whether to return per-observation similarity scores (default: FALSE).
#'
#' @return A \code{vector_rsa_model} object (S3 class) containing references to the dataset, design, and function parameters.
#'
#' @details
#' If `return_predictions` is TRUE, the output of `run_regional` or `run_searchlight`
#' will include a `prediction_table` tibble containing the observation-level RSA scores.
#' 
#' @export
vector_rsa_model <- function(dataset, design, 
                           distfun = cordist(), 
                           rsa_simfun = c("pearson", "spearman"),
                           nperm=0, 
                           save_distributions=FALSE,
                           return_predictions=FALSE) {
  rsa_simfun <- match.arg(rsa_simfun)
  
  assertthat::assert_that(inherits(dataset, "mvpa_dataset"))
  assertthat::assert_that(inherits(design, "vector_rsa_design"),
                          msg = "Input must be a 'vector_rsa_design' object.")
  
  # Create the model spec, passing permutation and prediction parameters
  create_model_spec(
    "vector_rsa_model",
    dataset = dataset,
    design  = design,
    distfun = distfun,
    rsa_simfun = rsa_simfun,
    nperm = nperm,
    compute_performance = TRUE, # Assume performance is always computed
    save_distributions = save_distributions,
    return_predictions = return_predictions # Pass the new flag
  )
}


#' Train a vector RSA model
#'
#' @param obj An object of class \code{vector_rsa_model}.
#' @rdname train_model
#' @param train_dat A data frame or matrix representing the training subset (e.g., voxel intensities).
#' @param y Not used in vector RSA (here for consistency with other train_model generics).
#' @param indices The spatial indices of the training data (ROI, searchlight, etc.).
#' @param ... Additional arguments.
#'
#' @return A structure containing "scores" or similar second-order similarity results.
#' @method train_model vector_rsa_model
#' @export
train_model.vector_rsa_model <- function(obj, train_dat, y, indices, ...) {
  # "Training" here is effectively computing RSA-based metrics
  # train_dat is already a matrix of the ROI's data
  scores <- compute_trial_scores(obj, train_dat) 

  if (obj$nperm > 0 && !is.null(train_dat)) {
    # Only package roi_train_data if permutations are requested and data is available
    return(list(scores = scores, roi_train_data = as.matrix(train_dat)))
  } else {
    # For nperm = 0, or if train_dat were NULL, just return scores in the list
    # This ensures a consistent return type (list)
    return(list(scores = scores, roi_train_data = NULL))
  }
}


#' Compute trial scores for vector RSA
#'
#' @keywords internal
#' @noRd
compute_trial_scores <- function(obj, X) {

  # Convert X to matrix if needed and check row count
  X <- as.matrix(X)

  # If we have fewer than two observations, there is no meaningful
  # second-order similarity to compute. Return a vector of NA as a
  # sentinel value rather than attempting to compute pairwise distances.
  if (nrow(X) < 2) {
    return(rep(NA_real_, nrow(X)))
  }

  # Retrieve relevant design expansions
  precomputed <- obj$design$model_mat
  dissimilarity_matrix <- precomputed$Dexpanded
  #cross_block_data     <- precomputed$cross_block_data

  # Basic dimensionality checks
  if (nrow(X) != nrow(dissimilarity_matrix)) {
    stop(
      "nrow(X) (", nrow(X), ") must match nrow(precomputed$Dexpanded) (",
      nrow(dissimilarity_matrix), ")"
    )
  }
  if (length(obj$design$block) != nrow(X)) {
    stop(
      "length(obj$design$block) (", length(obj$design$block),
      ") must match nrow(X) (", nrow(X), ")"
    )
  }
  
  # This function computes second-order similarity:
  #   second_order_similarity(distfun, X, D, block_var, rsa_simfun)
  # 'distfun' can compute distances on X, 'D' is the reference matrix
  # 'block_var' is in obj$design$block, and 'rsa_simfun' is the correlation method.
  second_order_similarity(
    distfun   = obj$distfun,
    X         = X,
    Dref      = dissimilarity_matrix,
    block_var = obj$design$block,
    method    = obj$rsa_simfun
  )
}


#' Compute Second-Order Similarity Scores
#'
#' Calculates correlation-based \emph{second order similarity} between:
#' \itemize{
#'   \item A \strong{full NxN distance matrix} computed from \code{X} via \code{distfun}, and 
#'   \item A \code{Dref} matrix (the "reference" dissimilarities).
#' }
#' For each row \code{i}, this excludes same-block comparisons by selecting \code{which(block_var != block_var[i])}.
#'
#' @param distfun An S3 distance object (see \code{\link{create_dist}}) 
#'   specifying how to compute a pairwise distance matrix from \code{X}.
#' @param X A numeric matrix (rows = observations, columns = features).
#' @param Dref A numeric NxN reference matrix of dissimilarities (e.g., from an ROI mask or a prior).
#' @param block_var A vector indicating block/group memberships for each row in \code{X}.
#' @param method Correlation method: "pearson" or "spearman".
#'
#' @return A numeric vector of length \code{nrow(X)}, where each entry is
#' the correlation (using \code{method}) between \code{distance_matrix[i, valid]} and
#' \code{Dref[i, valid]}, with \code{valid = which(block_var != block_var[i])}.
#'
#' @details
#' This function first calls \code{pairwise_dist(distfun, X)}, obtaining an NxN matrix 
#' of \emph{all} pairwise distances. It does not do block-based exclusion internally. 
#' Instead, for each row \code{i}, it excludes same-block rows from the correlation 
#' by subsetting the distances to \code{valid_indices}.
#'
#' @examples
#' # Suppose we have X (10x5), a reference D (10x10), block var, and a correlation distfun:
#' X <- matrix(rnorm(50), 10, 5)
#' D <- matrix(runif(100), 10, 10)
#' block <- rep(1:2, each=5)
#' dist_obj <- cordist(method="pearson")
#' scores <- second_order_similarity(dist_obj, X, D, block, method="spearman")
#'
#' @export
second_order_similarity <- function(distfun, X, Dref, block_var, method=c("pearson", "spearman")) {
  method <- match.arg(method)
  
  # 1) Compute a full NxN distance matrix from X
  distance_matrix <- pairwise_dist(distfun, X)
  
  n <- nrow(X)
  scores <- numeric(n)
  
  # 2) For each row i, exclude same-block comparisons
  for (i in seq_len(n)) {
    valid_indices <- which(block_var != block_var[i])
    if (length(valid_indices) > 0) {
      x_vec <- distance_matrix[i, valid_indices]
      ref_vec <- Dref[i, valid_indices]
      if (length(x_vec) > 1) {
        scores[i] <- cor(x_vec, ref_vec, method = method)
      } else {
        scores[i] <- NA
      }
    } else {
      scores[i] <- NA
    }
  }
  
  scores
}


#' Print Method for vector_rsa_model
#'
#' @param x An object of class \code{vector_rsa_model}.
#' @param ... Additional arguments (ignored).
#' @export
print.vector_rsa_model <- function(x, ...) {
  # Header
  cat(rep("=", 50), "\\n")
  cat("          Vectorized RSA Model          \\n")
  cat(rep("=", 50), "\\n\\n")

  # Dataset information
  if (!is.null(x$dataset)) {
    dims <- dim(x$dataset$train_data)
    dims_str <- if (!is.null(dims)) paste(dims, collapse = " x ") else "Unknown"
    cat("Dataset:\\n")
    cat("  |- Data Dimensions: ", dims_str, "\\n")
    if (!is.null(x$dataset$mask)) {
      cat("  |- Mask Length:     ", length(x$dataset$mask), "\\n")
    } else {
      cat("  |- Mask:            None\\n")
    }
    cat("\\n")
  }

  # Design information
  if (!is.null(x$design)) {
    n_labels <- length(x$design$labels)
    n_blocks <- length(unique(x$design$block))
    Ddim <- dim(x$design$model_mat$Dexpanded)
    cat("Design:\\n")
    cat("  |- Number of Labels:    ", n_labels, "\\n")
    cat("  |- Number of Blocks:    ", n_blocks, "\\n")
    cat("  |- Dissimilarity Matrix:", paste0(Ddim[1], " x ", Ddim[2]), "\\n\\n")
  }

  # Model Specification
  cat("Model Specification:\\n")
  dist_name <- tryCatch({
    if (!is.null(x$distfun$name)) {
      x$distfun$name
    } else {
      class(x$distfun)[1]
    }
  }, error = function(e) class(x$distfun)[1])
  cat("  |- Distance Function:   ", dist_name, "\\n")
  cat("  |- RSA Similarity Func: ", x$rsa_simfun, "\\n\\n")

  # Footer
  cat(rep("=", 50), "\\n")
}

#' Print Method for vector_rsa_design
#'
#' @param x A vector_rsa_design object.
#' @param ... Additional arguments (ignored).
#' @export
print.vector_rsa_design <- function(x, ...) {
  # Create a border line
  border <- paste(rep("=", 50), collapse="")

  # Header
  cat(border, "\\n")
  cat("         Vectorized RSA Design          \\n")
  cat(border, "\\n\\n")

  # Distance Matrix Information
  if (!is.null(x$D)) {
    Ddim <- dim(x$D)
    cat("Distance Matrix:\\n")
    cat("  |- Original Dimensions: ", paste0(Ddim[1], " x ", Ddim[2]), "\\n")
  }

  # Labels Information
  if (!is.null(x$labels)) {
    n_labels <- length(x$labels)
    cat("Labels:\\n")
    cat("  |- Total Labels: ", n_labels, "\\n")
    # Display the first few labels as a sample
    sample_labels <- paste(head(x$labels, 5), collapse = ", ")
    if(n_labels > 5) sample_labels <- paste0(sample_labels, ", ...")
    cat("  |- Sample:         ", sample_labels, "\\n")
  }

  # Block Information
  if (!is.null(x$block)) {
    unique_blocks <- sort(unique(x$block))
    n_blocks <- length(unique_blocks)
    cat("Blocks:\\n")
    cat("  |- Number of Blocks: ", n_blocks, "\\n")
    cat("  |- Block Labels:     ", paste(unique_blocks, collapse = ", "), "\\n")
  }

  # Expanded D Matrix Information
  if (!is.null(x$model_mat) && !is.null(x$model_mat$Dexpanded)) {
    Dexp_dim <- dim(x$model_mat$Dexpanded)
    cat("Expanded D Matrix:\\n")
    cat("  |- Dimensions:       ", paste0(Dexp_dim[1], " x ", Dexp_dim[2]), "\\n")
  }

  # Footer
  cat("\\n", border, "\\n")
}

#' Evaluate model performance for vector RSA
#'
#' Computes the mean second-order similarity score and handles permutation testing.
#'
#' @param object The vector RSA model specification.
#' @param predicted Ignored (vector RSA doesn't predict in the typical sense).
#' @param observed The computed second-order similarity scores (vector from train_model).
#' @param roi_data_for_perm New parameter
#' @param nperm Number of permutations from the model spec.
#' @param save_distributions Logical, whether to save full permutation distributions.
#' @param ... Additional arguments.
#'
#' @return A list containing the mean RSA score (`rsa_score`), raw scores, and 
#'   optional permutation results (`p_values`, `z_scores`, `permutation_distributions`).
#' @importFrom stats sd
#' @export
evaluate_model.vector_rsa_model <- function(object,
                                             predicted, # Ignored
                                             observed,  # These are the scores from train_model
                                             roi_data_for_perm = NULL, # New parameter
                                             nperm = 0,
                                             save_distributions = FALSE,
                                             ...) 
{
  # Primary metric: mean of the second-order similarity scores
  mean_rsa_score <- mean(observed, na.rm = TRUE)

  # If compute_trial_scores returned only NA values (e.g., because the ROI
  # contained fewer than two observations), there is no valid RSA score to
  # compute or permute.  Return NA immediately.
  if (length(observed) < 2 || all(is.na(observed))) {
    return(list(rsa_score = NA_real_, permutation_results = NULL))
  }

  perm_results <- NULL
  if (nperm > 0) {
    if (is.null(roi_data_for_perm) || !is.matrix(roi_data_for_perm)) {
      stop("evaluate_model.vector_rsa_model: roi_data_for_perm is NULL or not a matrix, but nperm > 0. Permutations cannot run without ROI data.")
    }
    
    X_roi_matrix <- roi_data_for_perm # This is now the ROI's data matrix

    perm_means <- numeric(nperm)
    # Permutation loop
    for (p in seq_len(nperm)) {
      if (nrow(X_roi_matrix) < 2) { # Safety for very small ROIs
         # This case should ideally lead to NA/NaN for the perm_mean if it makes sense statistically
         # or be caught earlier if an ROI is too small for any meaningful analysis.
         # For now, setting to NaN to indicate impossibility.
         perm_means[p] <- NaN 
         next
      }
      # Permute the rows of the ROI's data matrix
      perm_indices <- sample(nrow(X_roi_matrix))
      X_perm_matrix <- X_roi_matrix[perm_indices, , drop = FALSE]
      
      # Recompute trial scores under the null
      scores_perm <- compute_trial_scores(object, X_perm_matrix) # 'object' has design, Dref, block_var
      
      # Store the mean RSA score for this permutation
      perm_means[p] <- mean(scores_perm, na.rm = TRUE)
    }
    # Calculate p-value (one-sided: how many null >= observed)
    count_ge <- sum(perm_means >= mean_rsa_score, na.rm = TRUE)
    p_val <- (count_ge + 1) / (nperm + 1)
    # Compute z-score
    perm_mean <- mean(perm_means, na.rm = TRUE)
    perm_sd   <- sd(perm_means, na.rm = TRUE)
    z_score <- if (!is.na(perm_sd) && perm_sd > 0) {
      (mean_rsa_score - perm_mean) / perm_sd
    } else {
      NA_real_
    }
    # Assemble results
    perm_list <- list(
      p_values = setNames(p_val, "rsa_score"),
      z_scores = setNames(z_score, "rsa_score")
    )
    if (save_distributions) {
      perm_list$permutation_distributions <- list(
        rsa_score = perm_means
      )
    }
    perm_results <- perm_list
  }

  list(
    rsa_score           = mean_rsa_score,
    permutation_results = perm_results
  )
}

#' Merge results for vector RSA model
#'
#' Aggregates results (scores) and calls evaluate_model.
#' Vector RSA typically doesn't involve folds in the same way as classifiers,
#' so this mainly formats the output of train_model for the specific ROI/sphere.
#'
#' @param obj The vector RSA model specification (contains nperm etc.).
#' @param result_set A tibble from the processor. Expected to contain the output
#'   of `train_model.vector_rsa_model` (the scores vector) likely within `$result[[1]]`.
#' @param indices Voxel indices for the current ROI/searchlight sphere.
#' @param id Identifier for the current ROI/searchlight center.
#' @param ... Additional arguments.
#'
#' @return A tibble row with the final performance metrics for the ROI/sphere.
#' @importFrom tibble tibble
#' @importFrom futile.logger flog.error flog.warn
#' @method merge_results vector_rsa_model
merge_results.vector_rsa_model <- function(obj, result_set, indices, id, ...) {
  
  # Check for errors from previous steps (processor/train_model)
  if (any(result_set$error)) {
    emessage <- result_set$error_message[which(result_set$error)[1]]
    # Return standard error tibble structure
    return(
      tibble::tibble(
        result       = list(NULL), 
        indices      = list(indices),
        performance  = list(NULL), 
        id           = id,
        error        = TRUE,
        error_message= emessage
      )
    )
  }
  
  # Extract the list returned by train_model. 
  train_model_output <- result_set$result[[1]]
  
  if (is.null(train_model_output) || !is.list(train_model_output) || !"scores" %in% names(train_model_output)) {
     error_msg <- sprintf("merge_results (vector_rsa): train_model output structure invalid for ROI/ID %s. Expected list with 'scores'.", id)
     futile.logger::flog.error("ROI/Sphere ID %s: %s", id, error_msg)
     # Create NA performance matrix
     perf_names <- "rsa_score" # Base metric
     if (!is.null(obj$nperm) && obj$nperm > 0) { # Check obj$nperm directly
         perf_names <- c(perf_names, "p_rsa_score", "z_rsa_score")
     }
     perf_mat <- matrix(NA_real_, nrow=1, ncol=length(perf_names), 
                        dimnames=list(NULL, perf_names))
     return(tibble::tibble(result=list(NULL), indices=list(indices), performance=list(perf_mat), 
                           id=id, error=TRUE, error_message=error_msg))
  }
  
  scores <- train_model_output$scores
  roi_train_data_for_perm <- train_model_output$roi_train_data # This will be NULL if nperm=0 or data was NULL
  
  # Validate the extracted scores
  if (!is.numeric(scores)) {
      error_msg <- sprintf("merge_results (vector_rsa): Extracted scores are not numeric for ROI/Sphere ID %s.", id)
      futile.logger::flog.error(error_msg)
      # Create NA performance matrix
      perf_names <- "rsa_score"
      if (!is.null(obj$nperm) && obj$nperm > 0) {
          perf_names <- c(perf_names, "p_rsa_score", "z_rsa_score")
      }
      perf_mat <- matrix(NA_real_, nrow=1, ncol=length(perf_names), 
                         dimnames=list(NULL, perf_names))
      return(tibble::tibble(result=list(NULL), indices=list(indices), performance=list(perf_mat), 
                            id=id, error=TRUE, error_message=error_msg))
  }
  
  # Call evaluate_model to compute summary performance and permutations
  perf <- evaluate_model.vector_rsa_model(
    object            = obj,           
    predicted         = NULL, # Not used by vector_rsa's evaluate          
    observed          = scores,  
    roi_data_for_perm = roi_train_data_for_perm, # Pass the ROI's data
    nperm             = obj$nperm,     
    save_distributions = obj$save_distributions
  )
  
  # --- Collate performance matrix --- 
  base_metrics <- c(perf$rsa_score)
  base_names <- c("rsa_score")
  
  if (!is.null(perf$permutation_results)) {
      perm_p_values <- perf$permutation_results$p_values
      perm_z_scores <- perf$permutation_results$z_scores
      if (is.null(names(perm_p_values)) || is.null(names(perm_z_scores))){
           p_names <- paste0("p_", base_names)
           z_names <- paste0("z_", base_names)
      } else {
          p_names <- paste0("p_", names(perm_p_values))
          z_names <- paste0("z_", names(perm_z_scores))
      }
      perf_values <- c(base_metrics, perm_p_values, perm_z_scores)
      perf_names <- c(base_names, p_names, z_names)
  } else {
      perf_values <- base_metrics
      perf_names <- base_names
  }
  perf_mat <- matrix(perf_values, nrow = 1, ncol = length(perf_values), dimnames = list(NULL, perf_names))
  perf_mat <- perf_mat[, colSums(is.na(perf_mat)) < nrow(perf_mat), drop = FALSE]

  # --- Prepare results structure based on return_predictions flag --- 
  result_data <- if (isTRUE(obj$return_predictions)) {
      # Return scores structured for later assembly into prediction_table
      # Wrap scores in a list with a standard name
      list(rsa_scores = scores)
  } else {
      NULL # Return NULL if predictions are not requested
  }
  
  # Return the final tibble structure expected by the framework
  tibble::tibble(
    result      = list(result_data), # Store list(rsa_scores=scores) or NULL here
    indices     = list(indices),
    performance = list(perf_mat),
    id          = id,
    error       = FALSE,
    error_message = "~" 
  )
}






