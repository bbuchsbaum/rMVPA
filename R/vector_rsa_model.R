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
  assertthat::assert_that(
    nrow(D) == length(labels),
    msg = "The number of rows in D must match the number of labels."
  )
  
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
#'
#' @return A \code{vector_rsa_model} object (S3 class) containing references to the dataset, design, and function parameters.
#'
#' @details
#' The model references the already-precomputed cross-block data from the design. 
#' 
#' @export
vector_rsa_model <- function(dataset, design, distfun = cordist(), rsa_simfun = c("pearson", "spearman")) {
  rsa_simfun <- match.arg(rsa_simfun)
  
  assertthat::assert_that(inherits(dataset, "mvpa_dataset"))
  assertthat::assert_that(inherits(design, "vector_rsa_design"),
                          msg = "Input must be a 'vector_rsa_design' object.")
  
  # Create the model spec
  create_model_spec(
    "vector_rsa_model",
    dataset = dataset,
    design  = design,
    distfun = distfun,
    rsa_simfun = rsa_simfun
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
  
  # This external function presumably computes second-order similarity:
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






