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
#' The function verifies that all `labels` appear in `rownames(D)`. It also creates an expanded version
#' of the dissimilarity matrix (`Dexpanded`) matching the order of `labels`, and precomputes
#' cross-block information for later use.
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
  
  # Precompute cross-block data: for each block, which labels/indices are outside that block
  unique_blocks <- sort(unique(design$block))
  cross_block_data <- lapply(unique_blocks, function(b) {
    # Indices of everything *not* in block b
    inds_not_b <- which(design$block != b)
    list(
      other_labels = design$labels[inds_not_b],
      indices      = inds_not_b,
      block        = b
    )
  })
  names(cross_block_data) <- unique_blocks
  
  # Return as a list
  list(
    Dexpanded        = Dexpanded,
    cross_block_data = cross_block_data
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
#' The model references the already-precomputed cross-block data from the design. 
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
#' @param train_dat A data frame or matrix representing the training subset (e.g., voxel intensities).
#' @param y Not used in vector RSA (here for consistency with other train_model generics).
#' @param indices The spatial indices of the training data (ROI, searchlight, etc.).
#' @param ... Additional arguments.
#' 
#' @return A structure containing "scores" or similar second-order similarity results.
#' @export
train_model.vector_rsa_model <- function(obj, train_dat, y, indices, ...) {
  # "Training" here is effectively computing RSA-based metrics
  scores <- compute_trial_scores(obj, train_dat)
  return(scores)
}


#' Compute trial scores for vector RSA
#'
#' @keywords internal
#' @noRd
compute_trial_scores <- function(obj, X) {
 
  # Retrieve relevant design expansions
  precomputed <- obj$design$model_mat
  dissimilarity_matrix <- precomputed$Dexpanded
  cross_block_data     <- precomputed$cross_block_data
  
  # Convert X to matrix if needed
  X <- as.matrix(X)
  
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
  # Ensure that crayon is available
  if (!requireNamespace("crayon", quietly = TRUE)) {
    stop("Package 'crayon' is required for pretty printing. Please install it.")
  }
  
  # Define styles using crayon
  header_style  <- crayon::bold$cyan
  section_style <- crayon::yellow
  info_style    <- crayon::white
  number_style  <- crayon::green
  method_style  <- crayon::bold$blue
  label_style   <- crayon::magenta
  italic_style  <- crayon::italic$blue
  
  # Create a border line for header/footer
  border <- header_style(strrep("=", 50))
  
  # Header
  cat(border, "\n")
  cat(header_style("          Vectorized RSA Model          \n"))
  cat(border, "\n\n")
  
  # Dataset information
  if (!is.null(x$dataset)) {
    dims <- dim(x$dataset$train_data)
    dims_str <- if (!is.null(dims)) paste(dims, collapse = " x ") else "Unknown"
    cat(section_style("Dataset:\n"))
    cat(info_style("  ├─ Data Dimensions: "), number_style(dims_str), "\n")
    if (!is.null(x$dataset$mask)) {
      cat(info_style("  └─ Mask Length:     "), number_style(length(x$dataset$mask)), "\n")
    } else {
      cat(info_style("  └─ Mask:            "), crayon::red("None"), "\n")
    }
    cat("\n")
  }
  
  # Design information
  if (!is.null(x$design)) {
    n_labels <- length(x$design$labels)
    n_blocks <- length(unique(x$design$block))
    Ddim <- dim(x$design$model_mat$Dexpanded)
    cat(section_style("Design:\n"))
    cat(info_style("  ├─ Number of Labels:    "), number_style(n_labels), "\n")
    cat(info_style("  ├─ Number of Blocks:    "), number_style(n_blocks), "\n")
    cat(info_style("  └─ Dissimilarity Matrix:"), number_style(paste0(Ddim[1], " x ", Ddim[2])), "\n\n")
  }
  
  # Model Specification
  cat(section_style("Model Specification:\n"))
  cat(info_style("  ├─ Distance Function:   "), method_style(deparse(substitute(x$distfun))), "\n")
  cat(info_style("  └─ RSA Similarity Func: "), method_style(x$rsa_simfun), "\n\n")
  
  # Footer
  cat(border, "\n")
}

#' Print Method for vector_rsa_design
#'
#' @param x A vector_rsa_design object.
#' @param ... Additional arguments (ignored).
#' @export
print.vector_rsa_design <- function(x, ...) {
  # Ensure crayon is available
  if (!requireNamespace("crayon", quietly = TRUE)) {
    stop("Package 'crayon' is required for pretty printing. Please install it.")
  }
  
  # Define crayon styles
  header_style  <- crayon::bold$cyan
  section_style <- crayon::yellow
  info_style    <- crayon::white
  number_style  <- crayon::green
  label_style   <- crayon::magenta
  italic_style  <- crayon::italic$blue
  
  # Create a border line
  border <- header_style(strrep("=", 50))
  
  # Header
  cat(border, "\n")
  cat(header_style("         Vectorized RSA Design          \n"))
  cat(border, "\n\n")
  
  # Distance Matrix Information
  if (!is.null(x$D)) {
    Ddim <- dim(x$D)
    cat(section_style("Distance Matrix:\n"))
    cat(info_style("  ├─ Original Dimensions: "), number_style(paste0(Ddim[1], " x ", Ddim[2])), "\n")
  }
  
  # Labels Information
  if (!is.null(x$labels)) {
    n_labels <- length(x$labels)
    cat(section_style("Labels:\n"))
    cat(info_style("  ├─ Total Labels: "), number_style(n_labels), "\n")
    # Display the first few labels as a sample
    sample_labels <- paste(head(x$labels, 5), collapse = ", ")
    if(n_labels > 5) sample_labels <- paste0(sample_labels, ", ...")
    cat(info_style("  └─ Sample:         "), label_style(sample_labels), "\n")
  }
  
  # Block Information
  if (!is.null(x$block)) {
    unique_blocks <- sort(unique(x$block))
    n_blocks <- length(unique_blocks)
    cat(section_style("Blocks:\n"))
    cat(info_style("  ├─ Number of Blocks: "), number_style(n_blocks), "\n")
    cat(info_style("  └─ Block Labels:     "), label_style(paste(unique_blocks, collapse = ", ")), "\n")
  }
  
  # Expanded D Matrix Information
  if (!is.null(x$model_mat) && !is.null(x$model_mat$Dexpanded)) {
    Dexp_dim <- dim(x$model_mat$Dexpanded)
    cat(section_style("Expanded D Matrix:\n"))
    cat(info_style("  └─ Dimensions:       "), number_style(paste0(Dexp_dim[1], " x ", Dexp_dim[2])), "\n")
  }
  
  # Footer
  cat("\n", border, "\n")
}

#' Evaluate model performance for vector RSA
#'
#' Computes the mean second-order similarity score and handles permutation testing.
#'
#' @param object The vector RSA model specification.
#' @param predicted Ignored (vector RSA doesn't predict in the typical sense).
#' @param observed The computed second-order similarity scores (vector from train_model).
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
                                             nperm = 0,
                                             save_distributions = FALSE,
                                             ...) 
{
  # Primary metric: mean of the second-order similarity scores
  mean_rsa_score <- mean(observed, na.rm = TRUE)

  perm_results <- NULL
  if (nperm > 0) {
    # Retrieve original training data X
    X_orig <- object$dataset$train_data
    perm_means <- numeric(nperm)
    # Perform permutations
    for (p in seq_len(nperm)) {
      # Permute the rows of X (shuffles condition assignments)
      X_perm <- X_orig[sample(nrow(X_orig)), , drop = FALSE]
      # Recompute trial scores under the null
      scores_perm <- compute_trial_scores(object, X_perm)
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
#' @export
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
  
  # Extract the scores computed by train_model. 
  if (!"result" %in% names(result_set) || length(result_set$result) == 0 || is.null(result_set$result[[1]])) {
     error_msg <- "merge_results (vector_rsa): result_set missing or has NULL/empty 'result' field where scores were expected."
     futile.logger::flog.error("ROI/Sphere ID %s: %s", id, error_msg)
     # Create NA performance matrix
     perf_names <- "rsa_score"
     if (obj$nperm > 0) {
         perf_names <- c(perf_names, "p_rsa_score", "z_rsa_score")
     }
     perf_mat <- matrix(NA_real_, nrow=1, ncol=length(perf_names), 
                        dimnames=list(NULL, perf_names))
     return(tibble::tibble(result=list(NULL), indices=list(indices), performance=list(perf_mat), 
                           id=id, error=TRUE, error_message=error_msg))
  }
  
  scores <- result_set$result[[1]]
  
  # Validate the extracted scores
  if (!is.numeric(scores)) {
      error_msg <- sprintf("merge_results (vector_rsa): Extracted scores are not numeric for ROI/Sphere ID %s.", id)
      futile.logger::flog.error(error_msg)
      # Create NA performance matrix
      perf_names <- "rsa_score"
      if (obj$nperm > 0) {
          perf_names <- c(perf_names, "p_rsa_score", "z_rsa_score")
      }
      perf_mat <- matrix(NA_real_, nrow=1, ncol=length(perf_names), 
                         dimnames=list(NULL, perf_names))
      return(tibble::tibble(result=list(NULL), indices=list(indices), performance=list(perf_mat), 
                            id=id, error=TRUE, error_message=error_msg))
  }
  
  # Call evaluate_model to compute summary performance and permutations
  perf <- evaluate_model.vector_rsa_model(
    object    = obj,           
    predicted = NULL,          
    observed  = scores,        
    nperm     = obj$nperm,     
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






