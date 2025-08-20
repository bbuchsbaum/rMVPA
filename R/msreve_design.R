#' Constructor for msreve_design
#'
#' Creates an msreve_design object, which encapsulates the necessary
#' design information for a Multi-Dimensional Signed Representational
#' Voxel Encoding (MS-ReVE) analysis.
#'
#' @param mvpa_design An object of class \\code{mvpa_design}, containing
#'   information about conditions, blocks, and cross-validation.
#' @param contrast_matrix A numeric matrix (\\code{K x Q}) where \\code{K} is
#'   the number of conditions and \\code{Q} is the number of contrasts.
#'   Each column represents a contrast vector. It is highly recommended
#'   that columns are named to identify the contrasts.
#' @param name An optional character string to name the design.
#' @param include_interactions Logical. If TRUE, automatically add pairwise interaction contrasts using \code{\link{add_interaction_contrasts}}.
#' @param nuisance_rdms Optional named list of K x K matrices or \code{dist} objects representing
#'   nuisance RDMs to be included as additional predictors in the MS-ReVE regression.
#'   These are typically temporal or spatial nuisance patterns that should be accounted
#'   for but are not of primary interest.
#'
#' @return An object of class \\code{msreve_design}, which is a list containing:
#'   \\item{mvpa_design}{The input \\code{mvpa_design} object.}
#'   \\item{contrast_matrix}{The input \\code{contrast_matrix}.}
#'   \\item{name}{The name of the design.}
#'
#' @export
#' @examples
#' # Assume 'mvpa_des_obj' is a pre-existing mvpa_design object
#' # e.g. from mvpa_design(data=my_data_frame, formula = ~ condition_labels + run_labels,
#' #                       block_var = "run_labels")
#' # Let\'s say mvpa_des_obj implies 6 conditions based on unique(my_data_frame$condition_labels)
#' K <- 6 # Number of conditions
#' Q <- 2 # Number of contrasts
#'
#' # Example contrast matrix (K x Q)
#' C_mat <- matrix(c(
#'  # C1: Cond 1,2,3 vs 4,5,6
#'   1,  1,  1, -1, -1, -1,
#'  # C2: Cond 1,2 vs 3 (and 0 for 4,5,6 for simplicity here)
#'   1,  1, -2,  0,  0,  0
#' ), nrow = K, ncol = Q, byrow = FALSE)
#' colnames(C_mat) <- c("GroupComparison", "SubComparison")
#'
#' # if (inherits(mvpa_des_obj, "mvpa_design")) {
#' #  design_obj <- msreve_design(mvpa_des_obj, C_mat, name="example_msreve")
#' #  print(design_obj)
#' # }
#' #
#' # # Automatically add pairwise interactions
#' # design_obj_int <- msreve_design(mvpa_des_obj, C_mat,
#' #                                include_interactions = TRUE)
#' # colnames(design_obj_int$contrast_matrix)
msreve_design <- function(mvpa_design, contrast_matrix, name = "msreve_design_01", 
                         include_interactions = FALSE, nuisance_rdms = NULL) {
  # Basic assertions
  if (!inherits(mvpa_design, "mvpa_design")) {
    stop("`mvpa_design` must be an object of class \'mvpa_design\'.")
  }
  if (!is.matrix(contrast_matrix) || !is.numeric(contrast_matrix)) {
    stop("`contrast_matrix` must be a numeric matrix.")
  }
  if (is.null(colnames(contrast_matrix)) && ncol(contrast_matrix) > 0) {
    warning("`contrast_matrix` does not have column names. It is recommended to name your contrasts.")
    colnames(contrast_matrix) <- paste0("Contrast", 1:ncol(contrast_matrix))
  } else if (!is.null(colnames(contrast_matrix))) {
    # Ensure column names are unique
    original_colnames <- colnames(contrast_matrix)
    unique_colnames <- make.unique(original_colnames, sep = "_")
    if (!identical(original_colnames, unique_colnames)) {
        warning(paste("Duplicate contrast names detected and made unique (e.g., name -> name_1). Original:", 
                      paste(original_colnames[duplicated(original_colnames, fromLast=TRUE)], collapse=", ")))
        colnames(contrast_matrix) <- unique_colnames
    }
  }

  # Check: Number of rows in contrast_matrix should match number of conditions
  if (nrow(contrast_matrix) != mvpa_design$ncond) {
    stop(paste0("Number of rows in contrast_matrix (", nrow(contrast_matrix),
                ") must match number of conditions in mvpa_design (", mvpa_design$ncond, ")."))
  }

  # Suggestion: contrasts should be centered (sum of elements in each column is zero)
  col_sums <- colSums(contrast_matrix)
  if (any(abs(col_sums) > 1e-9)) { # Allow for small floating point inaccuracies
    non_centered_contrasts <- colnames(contrast_matrix)[abs(col_sums) > 1e-9]
    warning(paste("The following contrasts may not be centered (column sums are not zero):",
                  paste(non_centered_contrasts, collapse=", ")))
  }

  # Validate nuisance_rdms if provided
  if (!is.null(nuisance_rdms)) {
    if (!is.list(nuisance_rdms)) {
      stop("`nuisance_rdms` must be a named list of matrices or dist objects.")
    }
    if (is.null(names(nuisance_rdms)) || any(names(nuisance_rdms) == "")) {
      stop("`nuisance_rdms` must be a named list with non-empty names.")
    }
    
    # Check each nuisance RDM
    for (nm in names(nuisance_rdms)) {
      rdm <- nuisance_rdms[[nm]]
      if (inherits(rdm, "dist")) {
        # Convert dist to matrix for dimension check
        n_dist <- attr(rdm, "Size")
        if (n_dist != mvpa_design$ncond) {
          stop(paste0("Nuisance RDM '", nm, "' has size ", n_dist, 
                     " but mvpa_design has ", mvpa_design$ncond, " conditions."))
        }
      } else if (is.matrix(rdm)) {
        if (nrow(rdm) != mvpa_design$ncond || ncol(rdm) != mvpa_design$ncond) {
          stop(paste0("Nuisance RDM '", nm, "' must be a ", mvpa_design$ncond, 
                     " x ", mvpa_design$ncond, " matrix."))
        }
        if (!isSymmetric(rdm)) {
          warning(paste0("Nuisance RDM '", nm, "' is not symmetric."))
        }
      } else {
        stop(paste0("Nuisance RDM '", nm, "' must be a matrix or dist object."))
      }
    }
  }
  
  # Check if contrast_matrix is orthonormal
  is_orthonormal <- FALSE
  if (ncol(contrast_matrix) > 0) {
      ctc <- crossprod(contrast_matrix)
      expected_identity <- diag(ncol(contrast_matrix))
      # Use a reasonable tolerance for floating point comparisons
      if (isTRUE(all.equal(ctc, expected_identity, tolerance = 1e-8, check.attributes = FALSE))) {
          is_orthonormal <- TRUE
      }
  }
  
  # Construct the object
  obj <- list(
      mvpa_design = mvpa_design,
      contrast_matrix = contrast_matrix,
      name = name,
      nuisance_rdms = nuisance_rdms
  )
  
  # Derive and store condition-to-block mapping if block_var exists
  if (!is.null(mvpa_design$block_var) && !is.null(mvpa_design$Y)) {
    if (length(mvpa_design$block_var) == length(mvpa_design$Y)) {
      unique_conditions_ordered <- rownames(contrast_matrix) # K conditions
      if (is.null(unique_conditions_ordered) && nrow(contrast_matrix) > 0) {
          # Fallback if contrast_matrix has no rownames, try from mvpa_design if available
          # This assumes mvpa_design$conditions is ordered correctly for the contrast matrix rows
          if (!is.null(mvpa_design$conditions)) {
              unique_conditions_ordered <- mvpa_design$conditions
          } else {
              # Final fallback: unique levels of Y, but order might not match contrast matrix implicitly
              # This path implies contrast_matrix rows must align with sort(unique(Y))
              warning("Contrast matrix has no rownames and mvpa_design$conditions is NULL. Assuming row order matches sorted unique condition labels from Y for block mapping.")
              unique_conditions_ordered <- sort(unique(as.character(mvpa_design$Y)))
          }
      }

      if (!is.null(unique_conditions_ordered) && length(unique_conditions_ordered) == nrow(contrast_matrix)){
          condition_block_list <- lapply(unique_conditions_ordered, function(cond_name) {
            # Get sample-wise block_var for samples belonging to this condition
            blocks_for_cond <- mvpa_design$block_var[mvpa_design$Y == cond_name]
            unique(blocks_for_cond)
          })
          names(condition_block_list) <- unique_conditions_ordered
          obj$condition_block_list <- condition_block_list
      } else {
          warning("Could not determine ordered unique conditions for block mapping; condition_block_list not created.")
      }
    } else {
      warning("mvpa_design$block_var and mvpa_design$Y have different lengths. Cannot create condition_block_list.")
    }
  } else {
      # If no block_var or Y, cannot create the list
      # obj$condition_block_list <- NULL # (implicitly NULL)
  }

  attr(obj, "is_orthonormal") <- is_orthonormal

  obj_final <- structure(
    obj,
    class = c("msreve_design", "list")
  )

  if (include_interactions) {
    obj_final <- add_interaction_contrasts(obj_final)
  }

  obj_final
}

#' @export
#' @method print msreve_design
print.msreve_design <- function(x, ...) {
  has_crayon <- requireNamespace("crayon", quietly = TRUE)

  header_style  <- if (has_crayon) crayon::bold$cyan else function(txt) txt
  section_style <- if (has_crayon) crayon::yellow else function(txt) txt
  info_style    <- if (has_crayon) crayon::white else function(txt) txt
  dim_style     <- if (has_crayon) crayon::green else function(txt) txt
  name_style    <- if (has_crayon) crayon::italic$blue else function(txt) txt

  cat("\\n", header_style("█▀▀ MS-ReVE Design: "), name_style(x$name), header_style(" ▀▀█"), "\\n\\n")

  cat(section_style("├─ MVPA Design Source"), "\\n")
  # This assumes mvpa_design might have a 'name' or can be summarized.
  # For now, just indicate its class.
  cat(info_style("│  └─ Type: "), dim_style(class(x$mvpa_design)[1]), "\\n")

  cat(section_style("└─ Contrast Matrix"), "\\n")
  dims_cm <- dim(x$contrast_matrix)
  cat(info_style("   ├─ Dimensions: "), dim_style(paste0(dims_cm[1], " Conditions × ", dims_cm[2], " Contrasts")), "\\n")

  contrast_names <- colnames(x$contrast_matrix)
  if (!is.null(contrast_names) && length(contrast_names) > 0) {
    cat(info_style("   └─ Contrast Names: "), info_style(paste(contrast_names, collapse=", ")), "\\n")
  } else {
    cat(info_style("   └─ Contrast Names: "), dim_style("Not specified or empty"), "\\n")
  }

  cat("\\n")
  invisible(x)
}

#' Orthogonalize a Contrast Matrix
#'
#' Orthogonalizes the columns of a contrast matrix using QR decomposition.
#' The resulting matrix will have orthonormal columns spanning the same
#' space as the original columns, up to the rank of the input matrix.
#' Sign of the output columns is heuristically aligned with input columns.
#'
#' @param C A numeric matrix (K x Q) where columns represent contrast vectors.
#' @return An orthogonalized matrix. If the input matrix \code{C} is rank-deficient
#'   (rank < Q), the output matrix will have Q columns, but only the first
#'   \code{rank(C)} columns will be non-zero and form an orthonormal basis;
#'   subsequent columns will be zero vectors. Column names from \code{C} are preserved.
#' @export
#' @examples
#' K <- 6 # Number of conditions
#' Q <- 2 # Number of contrasts
#' C_orig <- matrix(c( 1,  1,  1, -1, -1, -1,  # Contrast 1
#'                     1, -1,  0,  1, -1,  0), # Contrast 2 (not orthogonal to C1)
#'                  nrow=K, ncol=Q)
#' colnames(C_orig) <- c("MainEffect", "InteractionLike")
#' C_ortho <- orthogonalize_contrasts(C_orig)
#' # print(round(crossprod(C_ortho), 10)) # Should be close to identity matrix
#'
#' # Example with a rank-deficient matrix (3rd contrast is sum of first two)
#' C_rank_def <- cbind(C_orig, C_orig[,1] + C_orig[,2])
#' colnames(C_rank_def) <- c("C1", "C2", "C3_dependent")
#' C_ortho_def <- orthogonalize_contrasts(C_rank_def)
#' # print(round(crossprod(C_ortho_def), 10))
#' # The 3rd column of C_ortho_def will be zeros.
orthogonalize_contrasts <- function(C) {
  if (!is.matrix(C) || !is.numeric(C)) {
    stop("Input `C` must be a numeric matrix.")
  }
  if (ncol(C) == 0) {
    return(C) # Return empty matrix if no columns
  }

  qr_decomp <- qr(C)
  rank_C <- qr_decomp$rank

  # Economy-size Q from QR decomposition
  # qr.Q by default returns Q with dimensions nrow(C) x rank(C)
  Q_econ <- qr.Q(qr_decomp, complete = FALSE)

  # Economy-size R for sign correction
  R_econ <- qr.R(qr_decomp, complete = FALSE)

  # Correct signs of Q_econ columns to heuristically match original C directions
  # A common convention is to ensure R_econ has positive diagonal elements.
  # If R_econ[j,j] is negative, flip the sign of Q_econ[,j] and R_econ[j,].
  if (rank_C > 0) {
    for (j in 1:rank_C) {
      if (R_econ[j,j] < 0) {
        Q_econ[,j] <- -Q_econ[,j]
        # R_econ[j,] <- -R_econ[j,] # Not strictly needed for Q output
      }
    }
  }
  
  # Initialize output matrix with original dimensions, filled with zeros
  C_ortho <- matrix(0, nrow = nrow(C), ncol = ncol(C))

  # Fill the first rank_C columns with the orthonormal vectors from Q_econ
  if (rank_C > 0) {
    C_ortho[, 1:rank_C] <- Q_econ[, 1:rank_C, drop = FALSE]
  }

  # Preserve column names
  if (!is.null(colnames(C))) {
    colnames(C_ortho) <- colnames(C)
  } else {
    colnames(C_ortho) <- paste0("Ortho_C", 1:ncol(C))
  }

  if (rank_C < ncol(C)) {
    warning(paste0("Input matrix C is rank-deficient (rank ", rank_C, " < ", ncol(C), " columns). ",
                   "The output matrix has ", ncol(C), " columns, but only the first ", rank_C,
                   " are non-zero and form an orthonormal basis. Subsequent columns are zero vectors."))
  }
  
return(C_ortho)
}

#' Add Interaction Contrasts to an msreve_design
#'
#' Creates new contrast columns representing pairwise interactions of existing
#' contrasts in an \code{msreve_design} object. Interactions are computed as
#' element-wise products of the contrast vectors.
#'
#' @details
#' Interaction contrasts are created by element-wise multiplication of pairs
#' of contrast vectors. If the resulting interaction is a zero vector (which
#' occurs when the original contrasts have non-overlapping support, i.e., no
#' conditions where both contrasts are non-zero), the interaction is skipped
#' with an informative message. This commonly happens with contrasts that
#' compare distinct subsets of conditions, such as c(1,-1,0,0) and c(0,0,1,-1).
#'
#' @param design An object of class \code{msreve_design}.
#' @param pairs Optional two-column matrix or list of character vectors
#'   specifying pairs of contrast column names. Default \code{NULL} uses all
#'   pairwise combinations.
#' @param orthogonalize Logical; if \code{TRUE} (default) the expanded contrast
#'   matrix is passed through \code{\link{orthogonalize_contrasts}}.
#'
#' @return The updated \code{msreve_design} object with non-zero interaction 
#'   columns appended. Zero interactions are automatically skipped.
#'   
#' @examples
#' \dontrun{
#' # Example with non-overlapping contrasts (zero interaction)
#' C1 <- matrix(c(1,-1,0,0, 0,0,1,-1), nrow=4, 
#'              dimnames=list(NULL, c("A","B")))
#' # A compares conditions 1 vs 2, B compares 3 vs 4
#' # Their interaction will be zero and skipped
#' 
#' # Example with overlapping contrasts (non-zero interaction)  
#' C2 <- matrix(c(1,1,-1,-1, 1,-1,1,-1), nrow=4,
#'              dimnames=list(NULL, c("Main1","Main2")))
#' # These contrasts overlap and will produce a meaningful interaction
#' }
#' @export
add_interaction_contrasts <- function(design, pairs = NULL, orthogonalize = TRUE) {
  if (!inherits(design, "msreve_design")) {
    stop("`design` must be an object of class 'msreve_design'.")
  }

  C <- design$contrast_matrix
  if (is.null(colnames(C))) {
    stop("contrast_matrix must have column names to define interactions")
  }

  cn <- colnames(C)
  if (is.null(pairs)) {
    if (ncol(C) < 2) return(design)
    pair_list <- combn(cn, 2, simplify = FALSE)
  } else if (is.matrix(pairs) && ncol(pairs) == 2) {
    pair_list <- split(pairs, seq_len(nrow(pairs)))
  } else if (is.list(pairs)) {
    pair_list <- pairs
  } else if (is.character(pairs)) {
    pair_list <- lapply(pairs, function(x) strsplit(x, ":", fixed = TRUE)[[1]])
  } else {
    stop("`pairs` must be NULL, a two-column matrix, a list, or a character vector")
  }

  for (p in pair_list) {
    if (length(p) != 2 || !all(p %in% cn)) {
      stop("Each pair must contain two valid contrast column names")
    }
    new_col <- C[, p[1]] * C[, p[2]]
    
    # Check if interaction is effectively zero (contrasts have non-overlapping support)
    if (all(abs(new_col) < 1e-10)) {
      message(paste0("Interaction ", p[1], "_x_", p[2], 
                     " is zero (contrasts have non-overlapping support) and will be skipped"))
      next
    }
    
    new_name <- paste0(p[1], "_x_", p[2])
    C <- cbind(C, new_col)
    cn <- c(cn, new_name)
  }
  colnames(C) <- cn

  if (orthogonalize) {
    C <- orthogonalize_contrasts(C)
    attr(design, "is_orthonormal") <- TRUE
  } else {
    ctc <- crossprod(C)
    attr(design, "is_orthonormal") <- isTRUE(all.equal(ctc, diag(ncol(C)), tolerance = 1e-8, check.attributes = FALSE))
  }

  design$contrast_matrix <- C
  design
}
