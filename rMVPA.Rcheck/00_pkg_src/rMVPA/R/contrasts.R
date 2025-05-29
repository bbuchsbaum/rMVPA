#' Parse Mini-Formula Specification to Metadata Tibble
#'
#' Internal helper function to convert the mini-formula syntax within
#' `contrasts(spec=...)` into a structured metadata tibble suitable for
#' `model.matrix`.
#'
#' @param labels Character vector of all unique category labels in the desired order.
#' @param spec A formula object containing the mini-DSL (e.g.,
#'   `~ factor1(levelA + levelB ~ levelC + .) + factor2(...)`).
#'
#' @return A tibble with the first column `label` (matching input `labels` order)
#'   and subsequent columns representing the binary factors derived from the spec.
#'   Coding is typically +1 for the left side, -1 for the right side (or `.`), 0 otherwise.
#'
#' @keywords internal
#' @noRd
#' @importFrom stringr str_match str_squish str_split
#' @importFrom tibble tibble
parse_spec_to_metadata <- function(labels, spec) {
  if (!inherits(spec, "formula")) {
    stop("`spec` must be a formula when `metadata` is not provided.", call. = FALSE)
  }
  
  term_labels <- attr(terms(spec), "term.labels")
  if (length(term_labels) == 0) {
      stop("Formula `spec` contains no terms to parse.", call. = FALSE)
  }
  
  metadata_list <- list()
  metadata_list$label <- labels
  n_labels <- length(labels)
  
  # Regex to capture: factor_name( levels_A ~ levels_B )
  # Allows spaces, requires factor name, captures level groups
  pattern <- "^\\s*(\\w+)\\s*\\(\\s*([^~]+?)\\s*~\\s*(.+?)\\s*\\)\\s*$"
  
  parsed_factors <- character() # Keep track of factors created
  
  for (term in term_labels) {
    match <- str_match(term, pattern)
    
    if (is.na(match[1, 1])) {
      # If it doesn't match the factor(...) pattern, it might be an interaction or simple term
      # These are handled later by model.matrix itself using the columns we create here.
      # We just need to make sure we parse the *definitions* first.
      next
    }
    
    factor_name <- match[1, 2]
    levels_A_str <- str_squish(match[1, 3])
    levels_B_str <- str_squish(match[1, 4])
    
    if (factor_name %in% parsed_factors) {
        stop("Factor '", factor_name, "' defined multiple times in the mini-formula spec.", call. = FALSE)
    }
    
    # Split levels by spaces or '+'
    split_levels <- function(level_str) {
        level_str_sq <- str_squish(level_str)
        if (level_str_sq == ".") return(character(0))
        # Split by '+', trim whitespace
        levels <- str_split(level_str_sq, pattern = "\\s*\\+\\s*", simplify = TRUE)[1,]
        return(str_squish(levels[levels != ""])) 
    }
    
    levels_A <- split_levels(levels_A_str)
    levels_B <- split_levels(levels_B_str)
    
    # Check validity of parsed levels against input labels
    if (length(levels_A) == 0 && levels_A_str != ".") {
         stop("Left side of '~' for factor '", factor_name, "' is empty or only whitespace.", call. = FALSE)
    }
    if (length(levels_A) > 0 && !all(levels_A %in% labels)) {
        stop("Level(s) '", paste(setdiff(levels_A, labels), collapse=", "),
             "' in factor '", factor_name, "' not found in provided labels.", call. = FALSE)
    }
    
    if (length(levels_B) > 0 && !all(levels_B %in% labels)) { 
        stop("Level(s) '", paste(setdiff(levels_B, labels), collapse=", "),
             "' in factor '", factor_name, "' not found in provided labels.", call. = FALSE)
    }
    
    # Handle '.' notation for levels_B
    defined_levels <- levels_A
    if (length(levels_B) > 0) {
        defined_levels <- c(defined_levels, levels_B)
        if (any(duplicated(defined_levels))) {
             stop("Label(s) '", paste(defined_levels[duplicated(defined_levels)], collapse=", "),
                 "' appear on both sides of '~' for factor '", factor_name, "'.", call. = FALSE)
        }
    } else {
        levels_B <- setdiff(labels, levels_A)
        if (length(levels_B) == 0 && length(labels) > 0 && length(levels_A) == length(labels)) {
             warning("Using '.' notation for factor '", factor_name, 
                     "' resulted in an empty set for the right side because all labels were on the left.",
                     "\nThis factor column will likely be all zeros or removed.", call. = FALSE)
        } else if (length(levels_B) == 0 && length(labels) > 0) {
            warning("Using '.' notation for factor '", factor_name, "' resulted in an empty set for the right side.",
                    "\nThis factor column will likely be all zeros or removed.", call. = FALSE)
        }
    }
    
    # Create binary column (+1 / -1 coding)
    factor_col <- numeric(n_labels)
    factor_col[match(levels_A, labels)] <- 1
    factor_col[match(levels_B, labels)] <- -1
    
    metadata_list[[factor_name]] <- factor_col
    parsed_factors <- c(parsed_factors, factor_name)
    
  } # End loop over terms
  
  if (length(parsed_factors) == 0) {
      stop("No valid factor definitions found in the mini-formula `spec`.", call. = FALSE)
  }
  
  metadata_list_tibble <- tibble::as_tibble(metadata_list)
  attr(metadata_list_tibble, "parsed_factor_names") <- parsed_factors
  return(metadata_list_tibble)
  
}

#' Generate Contrast Matrices
#'
#' Creates a numeric contrast matrix for use in RSA or encoding models, based on
#' condition labels and a specification.
#'
#' This function provides two main ways to define contrasts:
#' \enumerate{
#'   \item Via a `labels` vector and a `spec` formula using a mini-DSL like
#'         `~ factor1(levelA + levelB ~ levelC + .) + factor2(...)`.
#'   \item Via a `metadata` tibble (containing condition labels and predictor columns)
#'         and a standard R formula `spec` (e.g., `~ pred1 + pred2 + pred1:pred2`).
#' }
#' The function automatically handles centering, scaling, and optional orthogonalization.
#'
#' @param labels Character vector. Required if `metadata` is NULL. Specifies all
#'   unique condition labels in the desired order for the rows of the contrast matrix.
#' @param spec Formula. Defines the contrasts.
#'   If `metadata` is NULL, uses the mini-DSL (see Details).
#'   If `metadata` is provided, uses standard R formula syntax referencing columns
#'   in `metadata` (excluding the `label` column).
#' @param metadata Optional tibble/data.frame. If provided, it must contain a
#'   `label` column matching the conditions, and other columns representing
#'   features or factors used in the `spec` formula. `labels` argument is ignored
#'   if `metadata` is provided.
#' @param data Ignored in this version. Reserved for future extensions allowing
#'   direct input of feature matrices or RDMs for PCA/MDS contrasts.
#' @param centre Logical. If TRUE (default), columns of the resulting matrix are mean-centered.
#' @param scale Character string specifying scaling method after centering (if `orth=FALSE`).
#'   Options: `"none"` (default), `"sd"` (divide by sample standard deviation),
#'   `"l2"` (divide by L2 norm / vector length to get unit vectors).
#'   This argument is *ignored* if `orth = TRUE`.
#' @param orth Logical. If FALSE (default), the matrix columns represent the specified
#'   contrasts directly (after centering/scaling).
#'   If TRUE, an orthonormal basis for the column space is computed via QR decomposition.
#'   Resulting columns will be orthogonal and have unit length (L2 norm = 1).
#' @param keep_attr Logical. If TRUE (default) and `orth = TRUE`, the original column
#'   names (before orthogonalization) are stored in `attr(C, "source")`.
#'
#' @details
#' **Mini-DSL for `spec` (when `metadata` is NULL):**
#' The formula should be of the form `~ name1(levelsA ~ levelsB) + name2(...)`.
#' \itemize{
#'   \item `name1`, `name2`, etc., become the factor/contrast names. These are used
#'         to generate initial binary (+1/-1/0) columns.
#'   \item `levelsA` are condition labels (from `labels` argument) separated by `+`.
#'         These get coded +1 for the named factor.
#'   \item `levelsB` are condition labels separated by `+`, or `.` (period).
#'         These get coded -1 for the named factor. `.` means "all labels not listed in `levelsA`".
#'   \item Labels not mentioned in a factor definition get coded 0 for that factor.
#'   \item Interaction terms (e.g., `factorName1:factorName2`) can be included in `spec`.
#'         These are passed to `model.matrix` which computes them based on the
#'         previously generated factor columns.
#' }
#' If `centre = TRUE` (default), the resulting columns from `model.matrix` are
#' mean-centered. For binary factors created by the DSL (e.g. +1/-1/0 coding),
#' if groups are balanced, they might already be near zero-mean. The explicit
#' centering step ensures this property regardless of input or balance.
#'
#' @section Orthogonalization:
#' If `orth = TRUE`, uses `qr.Q(qr(C))` to find an orthonormal basis. The number of columns
#' in the output will be the rank of the input matrix. Columns are renamed `Orth1`, `Orth2`, etc.
#' Scaling is ignored as the columns already have unit L2 norm.
#' If `keep_attr = TRUE`:
#'   `attr(C_orth, "source")` stores the names of the original columns
#'   that formed the basis for the orthogonalized matrix.
#'   `attr(C_orth, "dropped")` stores the names of original columns that were
#'   linearly dependent and thus not part of the basis, if any.
#' 
#' @section Scaling:
#' Applied *after* centering if `orth=FALSE`.
#' \itemize{
#'   \item `"none"`: No scaling.
#'   \item `"sd"`: `scale(..., center=FALSE, scale=TRUE)`. Uses sample standard deviation (N-1 denominator).
#'         Note that for columns with few unique values (e.g., a centered +/-1 contrast), the SD
#'         can be slightly different depending on whether the number of items is even or odd,
#'         due to the N-1 denominator. This might lead to minor differences in scaled norms.
#'   \item `"l2"`: Divides each column by its L2 norm (`sqrt(sum(x^2))`).
#' }
#'
#' @section Specific Behaviors:
#' \itemize{
#'   \item If `orth = TRUE` and the input matrix has only one column after potential centering,
#'         that column is scaled to unit L2 norm. Centering still depends on the `centre` argument.
#'   \item If `centre = FALSE` and `orth = TRUE`, the QR decomposition is performed on the
#'         *uncentered* columns.
#'   \item If the mini-DSL `. ` notation is used for `levelsB` and `levelsA` already contains
#'         all `labels`, `levelsB` becomes empty, potentially resulting in a constant (zero)
#'         column before centering. A warning is issued in this case.
#' }
#'
#' @section Masking:
#' This function masks the `stats::contrasts` function. To use the base R function,
#' explicitly call `stats::contrasts()`.
#'
#' @return A numeric matrix (K x Q), where K is the number of labels and Q is the
#'   number of contrasts/orthogonal components.
#'   If `orth = TRUE` and `keep_attr = TRUE`, it includes attributes detailing
#'   the source (`"source"`) and any dropped (`"dropped"`) columns due to rank deficiency.
#'
#' @seealso [transform_contrasts()], [make_feature_contrasts()], [stats::contrasts()]
#'
#' @export
#' @importFrom stats terms model.matrix sd
#' @importFrom methods is
#' @importFrom stats reformulate
#' @examples
#' labs <- c("faces","animals","plants","tools",
#'           "vehicles","furniture","buildings","food")
#'
#' # 1) Mini-DSL: 2x2 Factorial (Animacy x Size) + Interaction, Orthonormal
#' C1 <- contrasts(
#'         labels = labs,
#'         spec   = ~ anim( faces + animals + plants + food ~ . )
#'                  + size( faces + animals + tools + furniture ~ . )
#'                  + anim:size,
#'         orth   = TRUE)
#' print(colnames(C1))
#' print(attr(C1, "source"))
#' print(round(crossprod(C1), 5))
#'
#' # 2) Mini-DSL: One-vs-rest, Centered, Unit Length (L2)
#' C2 <- contrasts(labels = labs,
#'                 spec   = ~ faces( faces ~ . ) + tools( tools ~ . ),
#'                 scale = "l2")
#' print(round(colSums(C2^2), 5)) # Should be 1
#'
#' # 3) Metadata + Formula: Centered, Scaled (SD)
#' meta <- tibble::tribble(
#'   ~label,      ~anim, ~size,
#'   "faces",        1,    0,
#'   "animals",      1,    0,
#'   "plants",       1,    1,
#'   "tools",        0,    0,
#'   "vehicles",     0,    1,
#'   "furniture",    0,    0,
#'   "buildings",    0,    1,
#'   "food",         1,    1)
#' # Note: labels argument is ignored here, order comes from meta$label
#' # Also note: This function masks stats::contrasts
#' C3 <- contrasts(metadata = meta,
#'                 spec     = ~ anim + size + anim:size,
#'                 scale    = "sd")
#' print(round(colMeans(C3), 5)) # Should be 0
#' print(round(apply(C3, 2, sd), 5)) # Should be 1
contrasts <- function(labels = NULL,
                      spec,
                      metadata = NULL,
                      data = NULL, # Reserved for future
                      centre = TRUE,
                      scale = c("none", "sd", "l2"),
                      orth = FALSE,
                      keep_attr = TRUE) {

  scale <- match.arg(scale)
  
  # --- Input Mode Handling ---
  rhs_formula <- NULL
  final_labels <- NULL
  
  if (!is.null(metadata)) {
    # Mode 2: Metadata provided
    if (!is.data.frame(metadata)) stop("`metadata` must be a data.frame or tibble.", call. = FALSE)
    if (!"label" %in% colnames(metadata)) stop("`metadata` must contain a 'label' column.", call. = FALSE)
    if (!is.character(metadata$label) && !is.factor(metadata$label)) {
        stop("`metadata$label` column must be character or factor.", call. = FALSE)
    }
    metadata$label <- as.character(metadata$label) # Ensure character
    if (any(duplicated(metadata$label))) {
        stop("Labels in `metadata$label` column must be unique.", call. = FALSE)
    }
    
    final_labels <- metadata$label
    # Check formula type
    if (!inherits(spec, "formula")) stop("`spec` must be a formula when `metadata` is provided.", call. = FALSE)
    rhs_formula <- spec
    
    # Select only predictor columns for model.matrix data
    # predictor_data <- metadata[, setdiff(colnames(metadata), "label"), drop = FALSE]
    
    # AUDIT A-4: Allow numeric or factor columns - REMOVED this explicit check
    # model.matrix will error if necessary columns have wrong type.
    # is_valid_col_type <- sapply(predictor_data, function(x) is.numeric(x) || is.factor(x))
    # if (!all(is_valid_col_type)) {
    #     non_valid_type_cols <- names(is_valid_col_type)[!is_valid_col_type]
    #     stop("Non-numeric/factor columns found in `metadata` predictors: ", 
    #          paste(non_valid_type_cols, collapse=", "), 
    #          ".\nAll predictor columns must be numeric or factor.", call. = FALSE)
    # }
    # meta_for_mm <- as.data.frame(predictor_data) # model.matrix needs data.frame
    
    # model.matrix will select columns from the full metadata based on the formula
    meta_for_mm <- as.data.frame(metadata) 
    
  } else {
    # Mode 1: Mini-DSL spec + labels
    if (is.null(labels)) stop("`labels` must be provided if `metadata` is NULL.", call. = FALSE)
    if (!is.character(labels)) stop("`labels` must be a character vector.", call. = FALSE)
    if (any(duplicated(labels))) stop("`labels` must contain unique values.", call. = FALSE)
    
    final_labels <- labels
    parsed_meta <- parse_spec_to_metadata(labels, spec)
    meta_for_mm <- as.data.frame(parsed_meta[, -1, drop = FALSE]) # Exclude label column
    
    # AUDIT A-3: Update spec to be the RHS formula for model.matrix
    parsed_factor_names <- attr(parsed_meta, "parsed_factor_names")
    all_term_labels <- attr(terms(spec), "term.labels") # Original terms from spec
    # Keep only interaction terms from the original spec, as factors are already defined
    interaction_terms <- all_term_labels[grepl(":", all_term_labels)] 
    # Ensure interaction terms only use known parsed factors, or factors from original spec for robustness
    
    final_terms_for_formula <- parsed_factor_names
    
    # Validate and add interaction terms
    valid_interaction_terms <- character()
    if (length(interaction_terms) > 0) {
        for (it in interaction_terms) {
            # Split interaction into individual components
            components <- strsplit(it, ":")[[1]]
            # Check if all components are among parsed_factors (or potentially other terms if we extend)
            if (all(components %in% parsed_factor_names)) {
                valid_interaction_terms <- c(valid_interaction_terms, it)
            } else {
                warning("Interaction term '", it, "' involves factors not defined in the mini-DSL.",
                        "\nIt will be ignored.", call. = FALSE)
            }
        }
         final_terms_for_formula <- c(final_terms_for_formula, valid_interaction_terms)
    }

    if (length(final_terms_for_formula) > 0) {
        rhs_formula <- reformulate(final_terms_for_formula)
    } else {
        # This case implies original spec might have only factor definitions, no interactions,
        # or interactions were invalid. model.matrix with only intercept would fail later.
        # This should be caught by `ncol(mm) == 0` check later.
        # However, to avoid issues with reformulate(NULL), create a formula that would select nothing or intercept.
        # The parsed_meta already contains the columns, model.matrix(~0 + .) on this would be better
        # but reformulate needs term labels. For safety, if all are defined in parsed_meta:
        rhs_formula <- reformulate(colnames(meta_for_mm)) # Use the generated factor names
        if (length(colnames(meta_for_mm)) == 0) {
             stop("Mini-DSL parsing resulted in no valid factors for the model matrix.", call. = FALSE)
        }
    }
  }
  
  # --- Build Model Matrix ---
  mm <- tryCatch({
      stats::model.matrix(rhs_formula, data = meta_for_mm)
    },
    error = function(e) {
        stop("Failed to build model matrix from spec and metadata. Original error: \n", e$message, call. = FALSE)
    })
    
  # Remove intercept if present
  intercept_col <- which(colnames(mm) == "(Intercept)")
  if (length(intercept_col) > 0) {
    mm <- mm[, -intercept_col, drop = FALSE]
  }
  
  if (ncol(mm) == 0) {
      stop("Specification resulted in an empty contrast matrix (no terms other than intercept?).", call. = FALSE)
  }
  original_colnames <- colnames(mm)
  
  # --- Transformations ---
  
  # 1. Centering
  if (centre) {
    mm <- sweep(mm, 2, colMeans(mm, na.rm = TRUE), FUN = "-", check.margin = FALSE)
  }
  
  # 2. Orthogonalization (takes precedence over scaling)
  if (orth) {
    if (ncol(mm) > 1) {
        qr_obj <- qr(mm)
        rank_mm <- qr_obj$rank
        
        if (rank_mm < ncol(mm)) {
             warning("Original contrast matrix columns are linearly dependent.",
                     "\nOrthogonalized matrix will have fewer columns (", rank_mm,
                     ") than original non-intercept terms (", ncol(mm), ").", call. = FALSE)
        }
        if (rank_mm == 0) { 
             stop("Contrast matrix has rank 0 after centering; cannot orthogonalize.", call. = FALSE)
        }
        
        mm_orth <- qr.Q(qr_obj)[, seq_len(rank_mm), drop = FALSE]
        # Assign Orth* names immediately
        colnames(mm_orth) <- paste0("Orth", seq_len(rank_mm))
        
        if (keep_attr) {
          # AUDIT A-11: Store names of columns forming the basis and dropped columns
          pivot_indices <- qr_obj$pivot
          original_cols_contributing_to_basis <- original_colnames[pivot_indices[seq_len(rank_mm)]]
          attr(mm_orth, "source") <- original_cols_contributing_to_basis 
          
          if (rank_mm < length(pivot_indices)) {
             dropped_column_indices <- pivot_indices[(rank_mm + 1):length(pivot_indices)]
             # Ensure indices are valid
             valid_dropped_indices <- dropped_column_indices[dropped_column_indices <= length(original_colnames)]
             if (length(valid_dropped_indices) > 0) {
                attr(mm_orth, "dropped") <- original_colnames[valid_dropped_indices]
             }
          }
        }
        mm <- mm_orth
    } else if (ncol(mm) == 1) { 
        # AUDIT A-6: Orth=TRUE, 1 col: centre already done if requested above. Just normalize.
        # Use drop=FALSE to ensure mm remains a matrix
        col_values <- mm[, 1, drop = FALSE]
        col_norm <- sqrt(sum(col_values^2, na.rm = TRUE))
        if (col_norm > 1e-10) { 
             mm[, 1] <- col_values / col_norm
        } else {
             # Stop if the single column is zero
             stop("Single contrast column has zero length after centering; cannot orthogonalize/normalize.", call. = FALSE)
        }
        # Ensure consistent naming with multi-column orth
        colnames(mm) <- "Orth1"
        if (keep_attr) {
          attr(mm, "source") <- original_colnames # Use "source"
        }
    }
  } else {
    # 3. Scaling (only if not orthogonalized)
    if (scale != "none") {
      if (scale == "sd") {
          sds <- apply(mm, 2, stats::sd, na.rm = TRUE)
          zero_sd_cols <- sds < 1e-10
          if (any(zero_sd_cols)) {
              warning("One or more columns have near-zero standard deviation;",
                      "\nthese columns will not be scaled by SD.", call. = FALSE)
              sds[zero_sd_cols] <- 1 # Avoid division by zero, effectively no scaling for these
          }
          if (any(!zero_sd_cols)) {
               mm_to_scale <- mm[, !zero_sd_cols, drop = FALSE]
               sds_to_scale <- sds[!zero_sd_cols]
               # Ensure sweep inputs are valid dimensions
               if (ncol(mm_to_scale) > 0 && length(sds_to_scale) == ncol(mm_to_scale)) { 
                   mm[, !zero_sd_cols] <- sweep(mm_to_scale, 2, sds_to_scale, FUN = "/", check.margin = FALSE)
               }
          }
      } else if (scale == "l2") {
          norms <- sqrt(colSums(mm^2, na.rm = TRUE))
          zero_norm_cols <- norms < 1e-10
          if (any(zero_norm_cols)) {
              warning("One or more columns have near-zero L2 norm (length);",
                      "\nthese columns will not be scaled to unit length.", call. = FALSE)
              norms[zero_norm_cols] <- 1 # Avoid division by zero
          }
          if (any(!zero_norm_cols)) {
               mm_to_scale <- mm[, !zero_norm_cols, drop = FALSE]
               norms_to_scale <- norms[!zero_norm_cols]
               # Ensure sweep inputs are valid dimensions
               if (ncol(mm_to_scale) > 0 && length(norms_to_scale) == ncol(mm_to_scale)) { 
                   mm[, !zero_norm_cols] <- sweep(mm_to_scale, 2, norms_to_scale, FUN = "/", check.margin = FALSE)
               }
          }
      }
    } 
  } # End scaling/orth logic
  
  # Add row names based on final label order
  rownames(mm) <- final_labels
  
  # --- Final Naming --- 
  # Ensure column names are set correctly, especially after orthogonalization
  if (orth) {
      colnames(mm) <- paste0("Orth", seq_len(ncol(mm)))
  } else {
      if (!is.null(original_colnames) && ncol(mm) == length(original_colnames)) {
        colnames(mm) <- original_colnames
      }
  }
  # Provide vector-style names attribute for compatibility with testthat::expect_named on matrices
  # names(mm) <- colnames(mm) # REMOVED - incorrect for matrices
  
  return(mm)
}

#' Apply Transformations to an Existing Contrast Matrix
#'
#' Applies centering, scaling, and/or orthogonalization to a pre-existing
#' numeric contrast matrix.
#'
#' This function is useful for post-processing a contrast matrix, especially one
#' that might have been created by combining outputs from different sources
#' (e.g., theory-driven contrasts and data-driven contrasts) or by direct
#' manual construction.
#'
#' @param C A numeric matrix representing the contrasts (conditions x contrasts).
#'   Row names, if present, should correspond to condition labels and will be preserved.
#'   Column names, if present, will be used for the `"source"` attribute if
#'   `orth = TRUE` and `keep_attr = TRUE`.
#' @param centre Logical. If TRUE (default), columns of the matrix are mean-centered.
#' @param scale Character string specifying scaling method after centering (if `orth=FALSE`).
#'   Options: `"none"` (default), `"sd"` (divide by sample standard deviation),
#'   `"l2"` (divide by L2 norm / vector length to get unit vectors).
#'   This argument is *ignored* if `orth = TRUE`.
#' @param orth Logical. If FALSE (default), the matrix columns represent the specified
#'   contrasts directly (after centering/scaling).
#'   If TRUE, an orthonormal basis for the column space is computed via QR decomposition.
#'   Resulting columns will be orthogonal and have unit length (L2 norm = 1).
#' @param keep_attr Logical. If TRUE (default) and `orth = TRUE`, the original column
#'   names (before orthogonalization) are stored in `attr(C_transformed, "source")`,
#'   and names of linearly dependent columns removed during orthogonalization are
#'   stored in `attr(C_transformed, "dropped")`.
#'
#' @inheritSection contrasts Scaling
#' @inheritSection contrasts Orthogonalization
#' @inheritSection contrasts Specific Behaviors
#'
#' @return A transformed numeric matrix. If `orth = TRUE` and `keep_attr = TRUE`,
#'   it includes `"source"` and potentially `"dropped"` attributes.
#'
#' @seealso [contrasts()], [make_feature_contrasts()]
#'
#' @export
#' @examples
#' C_manual <- matrix(c( 1,  1,  0,  0,
#'                      -1,  1,  0,  0,
#'                       0,  0,  1,  1,
#'                       0,  0, -1,  1), nrow = 4, byrow = TRUE,
#'                    dimnames = list(paste0("Cond", 1:4), c("MainA", "MainB")))
#'
#' # Center and make orthonormal
#' C_transformed <- transform_contrasts(C_manual, orth = TRUE)
#' print(C_transformed)
#' print(attr(C_transformed, "source"))
#'
#' # Center and scale to unit L2 norm
#' C_l2 <- transform_contrasts(C_manual, scale = "l2")
#' print(round(colSums(C_l2^2), 5))
#'
transform_contrasts <- function(C,
                                centre = TRUE,
                                scale = c("none", "sd", "l2"),
                                orth = FALSE,
                                keep_attr = TRUE) {
  
  if (!is.matrix(C) || !is.numeric(C)) {
    stop("`C` must be a numeric matrix.", call. = FALSE)
  }
  if (nrow(C) == 0) { # Check only rows first
       stop("`C` must have at least one row.", call.=FALSE)
  }
  if (ncol(C) == 0) { # Check columns
       if (orth || scale != "none") {
          stop("`C` must have columns for the requested operations (orth or scaling).", call. = FALSE)
       }
       # If only centering requested on a 0-col matrix, or no operations, return as is
       return(C) 
  }
  
  scale <- match.arg(scale)
  mm <- C # Use 'mm' internally for consistency with contrasts() function
  original_colnames <- colnames(mm)
  if (is.null(original_colnames) && ncol(mm) > 0) {
      original_colnames <- paste0("Col", seq_len(ncol(mm)))
  }
  original_rownames <- rownames(mm)
  
  # --- Transformations (same logic as in contrasts() function) ---
  
  # 1. Centering
  if (centre) {
    mm <- sweep(mm, 2, colMeans(mm, na.rm = TRUE), FUN = "-", check.margin = FALSE)
  }
  
  # 2. Orthogonalization (takes precedence over scaling)
  if (orth) {
    if (ncol(mm) > 1) {
        qr_obj <- qr(mm)
        rank_mm <- qr_obj$rank
        
        if (rank_mm < ncol(mm)) {
             warning("Input contrast matrix columns are linearly dependent.",
                     "\nOrthogonalized matrix will have fewer columns (", rank_mm,
                     ") than original (", ncol(mm), ").", call. = FALSE)
        }
        if (rank_mm == 0) {
             stop("Contrast matrix has rank 0 after centering; cannot orthogonalize.", call. = FALSE)
        }
        
        mm_orth <- qr.Q(qr_obj)[, seq_len(rank_mm), drop = FALSE]
        # Assign Orth* names immediately
        colnames(mm_orth) <- paste0("Orth", seq_len(rank_mm))
        
        if (keep_attr) {
          # AUDIT A-11: Store names of columns forming the basis and dropped columns
          pivot_indices <- qr_obj$pivot
          original_cols_contributing_to_basis <- original_colnames[pivot_indices[seq_len(rank_mm)]]
          attr(mm_orth, "source") <- original_cols_contributing_to_basis 
          
          if (rank_mm < length(pivot_indices)) {
             dropped_column_indices <- pivot_indices[(rank_mm + 1):length(pivot_indices)]
             # Ensure indices are valid
             valid_dropped_indices <- dropped_column_indices[dropped_column_indices <= length(original_colnames)]
             if (length(valid_dropped_indices) > 0) {
                attr(mm_orth, "dropped") <- original_colnames[valid_dropped_indices]
             }
          }
        }
        mm <- mm_orth
    } else if (ncol(mm) == 1) {
        # Centering already applied if requested above. Just normalize.
        col_values <- mm[, 1, drop = FALSE]
        col_norm <- sqrt(sum(col_values^2, na.rm = TRUE))
        if (col_norm > 1e-10) {
             mm[, 1] <- col_values / col_norm
        } else {
             # Stop if the single column is zero
             stop("Single contrast column has zero length after centering; cannot orthogonalize/normalize.", call. = FALSE)
        }
        # Ensure consistent naming with multi-column orth
        colnames(mm) <- "Orth1"
        if (keep_attr) {
          attr(mm, "source") <- original_colnames # Use "source"
        }
    } # else ncol(mm) == 0, handled by initial check
  } else {
    # 3. Scaling (only if not orthogonalized)
    if (scale != "none") {
      if (scale == "sd") {
          sds <- apply(mm, 2, stats::sd, na.rm = TRUE)
          zero_sd_cols <- sds < 1e-10
          if (any(zero_sd_cols)) {
              warning("One or more columns have near-zero standard deviation;",
                      "\nthese columns will not be scaled by SD.", call. = FALSE)
              sds[zero_sd_cols] <- 1 # Avoid division by zero, effectively no scaling for these
          }
          if (any(!zero_sd_cols)) {
               mm_to_scale <- mm[, !zero_sd_cols, drop = FALSE]
               sds_to_scale <- sds[!zero_sd_cols]
               # Ensure sweep inputs are valid dimensions
               if (ncol(mm_to_scale) > 0 && length(sds_to_scale) == ncol(mm_to_scale)) { 
                   mm[, !zero_sd_cols] <- sweep(mm_to_scale, 2, sds_to_scale, FUN = "/", check.margin = FALSE)
               }
          }
      } else if (scale == "l2") {
          norms <- sqrt(colSums(mm^2, na.rm = TRUE))
          zero_norm_cols <- norms < 1e-10
          if (any(zero_norm_cols)) {
              warning("One or more columns have near-zero L2 norm (length);",
                      "\nthese columns will not be scaled to unit length.", call. = FALSE)
              norms[zero_norm_cols] <- 1 # Avoid division by zero
          }
          if (any(!zero_norm_cols)) {
               mm_to_scale <- mm[, !zero_norm_cols, drop = FALSE]
               norms_to_scale <- norms[!zero_norm_cols]
               # Ensure sweep inputs are valid dimensions
               if (ncol(mm_to_scale) > 0 && length(norms_to_scale) == ncol(mm_to_scale)) { 
                   mm[, !zero_norm_cols] <- sweep(mm_to_scale, 2, norms_to_scale, FUN = "/", check.margin = FALSE)
               }
          }
      }
    } 
  } # End scaling/orth logic
  
  # Restore row names
  if (!is.null(original_rownames)) {
      if (nrow(mm) == length(original_rownames)) {
        rownames(mm) <- original_rownames
      } else {
        warning("Could not restore original row names due to dimension change (e.g., from rank deficiency).")
      }
  }
  
  # --- Final Naming --- 
  # Ensure column names are set correctly, especially after orthogonalization
  if (orth) {
      colnames(mm) <- paste0("Orth", seq_len(ncol(mm)))
  } else {
      if (!is.null(original_colnames) && ncol(mm) == length(original_colnames)) {
        colnames(mm) <- original_colnames
      }
  }
  # Provide vector-style names attribute for compatibility with testthat::expect_named on matrices
  # names(mm) <- colnames(mm) # REMOVED - incorrect for matrices
  
  return(mm)
}

#' Generate Contrasts from a Feature Matrix (Optional PCA)
#'
#' Creates contrasts based on a matrix where rows represent conditions and
#' columns represent features (e.g., neural network embeddings, semantic features).
#' Optionally performs PCA to reduce dimensionality.
#'
#' @param features A numeric matrix (K x P) where K is the number of conditions
#'   and P is the number of features. Row names, if present, should correspond to
#'   condition labels. Column names are recommended.
#' @param labels Optional character vector of condition labels. If provided, rows of
#'   `features` matrix will be reordered to match this order. If NULL, the order
#'   from `rownames(features)` is used (if available).
#' @param use_pca Logical. If TRUE (default), performs Principal Component Analysis
#'   (PCA) on the features. If FALSE, uses the raw features directly.
#' @param centre_pca Logical. If `use_pca = TRUE`, should features be centered
#'   before PCA? (Default: TRUE)
#'   Note: Reordering of rows based on `labels` argument happens *before* PCA.
#' @param scale_pca Logical. If `use_pca = TRUE`, should features be scaled to
#'   unit variance before PCA? (Default: FALSE, as scaling can affect variance explained).
#'   Note: Reordering of rows based on `labels` argument happens *before* PCA.
#' @param pve Numeric (0 to 1). If `use_pca = TRUE`, selects the minimum number
#'   of principal components (PCs) needed to explain at least this proportion of
#'   variance. Ignored if `n_pcs` is specified. (Default: 0.9)
#' @param n_pcs Integer. If `use_pca = TRUE`, selects exactly this number of
#'   principal components. Takes precedence over `pve`. (Default: NULL)
#' @param prefix Character string to prepend to column names of the output matrix
#'   (e.g., "Feat_", "PCA_"). (Default: "Feat_")
#'
#' @return A numeric matrix (K x Q) where K matches the number of conditions/labels
#'   and Q is the number of selected features or principal components.
#'   Rows are ordered according to `labels` or `rownames(features)`.
#'   Columns are named using the `prefix` and either the original feature names
#'   (if `use_pca=FALSE`) or component numbers (e.g., "PCA_PC1", "PCA_PC2").
#'
#' @seealso [contrasts()], [transform_contrasts()]
#'
#' @export
#' @importFrom stats prcomp
#' @examples
#' # Example feature matrix (4 conditions, 5 features)
#' feat_mat <- matrix(rnorm(20), nrow = 4,
#'                    dimnames = list(paste0("Cond", 1:4), paste0("F", 1:5)))
#'
#' # Use raw features (first 3)
#' C_raw <- make_feature_contrasts(feat_mat[, 1:3], use_pca = FALSE, prefix="RawFeat_")
#' print(C_raw)
#'
#' # Use PCA, selecting top 2 PCs
#' C_pca <- make_feature_contrasts(feat_mat, use_pca = TRUE, n_pcs = 2, prefix="PCA_")
#' print(C_pca)
#'
#' # Use PCA, selecting >= 80% variance explained
#' C_pca_pve <- make_feature_contrasts(feat_mat, use_pca = TRUE, pve = 0.8, prefix="PCA_")
#' print(C_pca_pve)
#'
#' # Reorder based on labels
#' C_pca_reorder <- make_feature_contrasts(feat_mat, labels=c("Cond3", "Cond1", "Cond4", "Cond2"),
#'                                       use_pca = TRUE, n_pcs = 2, prefix="PCA_")
#' print(C_pca_reorder)
#'
make_feature_contrasts <- function(features,
                                   labels = NULL,
                                   use_pca = TRUE,
                                   centre_pca = TRUE,
                                   scale_pca = FALSE,
                                   pve = 0.9,
                                   n_pcs = NULL,
                                   prefix = "Feat_") {

  # --- Input Validation ---
  if (!is.matrix(features) || !is.numeric(features)) {
    stop("`features` must be a numeric matrix.", call. = FALSE)
  }
  if (nrow(features) == 0 || ncol(features) == 0) {
    stop("`features` must not be an empty matrix.", call. = FALSE)
  }

  feature_rownames <- rownames(features)
  final_labels <- NULL

  if (!is.null(labels)) {
    if (!is.character(labels)) stop("`labels` must be a character vector.", call. = FALSE)
    if (any(duplicated(labels))) stop("`labels` must contain unique values.", call. = FALSE)
    if (is.null(feature_rownames)) {
      stop("`features` matrix must have row names if `labels` argument is used for reordering.", call. = FALSE)
    }
    if (!all(labels %in% feature_rownames)) {
      stop("Not all values in `labels` are present in `rownames(features)`.", call. = FALSE)
    }
    if (length(labels) != length(unique(feature_rownames))) {
        warning("Length of `labels` does not match the number of unique row names in `features`.", call. = FALSE)
    }
    # Reorder features matrix rows
    features <- features[match(labels, feature_rownames), , drop = FALSE]
    final_labels <- labels
  } else if (!is.null(feature_rownames)) {
    if (any(duplicated(feature_rownames))) {
        stop("`rownames(features)` must be unique if `labels` is not provided.", call. = FALSE)
    }
    final_labels <- feature_rownames
  } else {
      warning("`features` matrix has no row names and `labels` were not provided. Rows will not be named.", call. = FALSE)
  }

  # --- Processing ---
  output_matrix <- NULL

  if (!use_pca) {
    # Return raw features
    output_matrix <- features
    if (!is.null(colnames(output_matrix))) {
      colnames(output_matrix) <- paste0(prefix, colnames(output_matrix))
    } else {
      colnames(output_matrix) <- paste0(prefix, seq_len(ncol(output_matrix)))
    }

  } else {
    # Perform PCA
    if (!is.null(n_pcs) && (!is.numeric(n_pcs) || n_pcs < 1 || n_pcs > ncol(features))) {
      stop("`n_pcs` must be a positive integer less than or equal to the number of features.", call. = FALSE)
    }
    if (!is.numeric(pve) || pve <= 0 || pve > 1) {
       stop("`pve` must be a number between 0 (exclusive) and 1 (inclusive).", call. = FALSE)
    }

    pca_res <- tryCatch({
        stats::prcomp(features, center = centre_pca, scale. = scale_pca)
    }, error = function(e) {
        stop("PCA failed. Original error: ", e$message, call. = FALSE)
    })

    # Determine number of components to keep
    n_comp_to_keep <- ncol(pca_res$x) # Start with max possible

    if (!is.null(n_pcs)) {
      n_comp_to_keep <- min(n_pcs, ncol(pca_res$x))
    } else {
      # Use PVE
      sum_sdev_sq <- sum(pca_res$sdev^2)
      if (sum_sdev_sq < 1e-10) { # If no variance, effectively 0 components
        n_comp_to_keep <- 0
      } else {
        variance_explained <- pca_res$sdev^2 / sum_sdev_sq
        cumulative_variance <- cumsum(variance_explained)
        n_comp_to_keep_pve <- which(cumulative_variance >= pve)[1]
        # If pve is not met even by all components (e.g. sum_sdev_sq > 0 but all var_explained < pve)
        # or if which() returns NA (e.g. for single component not meeting pve)
        # then n_comp_to_keep_pve would be NA. Fallback to all components if positive variance exists.
        if (is.na(n_comp_to_keep_pve)) {
            # This case should ideally not make n_comp_to_keep zero if sum_sdev_sq > 0
            # It means PVE threshold wasn't reached. Max components should be taken.
             n_comp_to_keep <- ncol(pca_res$x) 
        } else {
            n_comp_to_keep <- n_comp_to_keep_pve
        }
      }
    }

    if (n_comp_to_keep == 0) {
        stop("PCA resulted in 0 components selected based on criteria (pve/n_pcs).", call. = FALSE)
    }

    output_matrix <- pca_res$x[, 1:n_comp_to_keep, drop = FALSE]
    colnames(output_matrix) <- paste0(prefix, "PC", 1:n_comp_to_keep)
  }

  # --- Finalize ---
  # Ensure row names are set if possible
  if (!is.null(final_labels) && nrow(output_matrix) == length(final_labels)) {
    rownames(output_matrix) <- final_labels
  }
  # Provide names attribute to satisfy expect_named
  # names(output_matrix) <- colnames(output_matrix) # REMOVED - incorrect for matrices

  return(output_matrix)
} 