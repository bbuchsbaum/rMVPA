#' Average NeuroVec Data by Labels
#'
#' Computes the average brain activation pattern for each unique label/condition,
#' with optional normalization of individual volumes before averaging.
#'
#' @param neurovec A NeuroVec object containing the neuroimaging data
#' @param labels A vector of labels/conditions corresponding to each volume in neurovec
#' @param mask An optional mask (NeuroVol or logical array) to restrict analysis to specific voxels
#' @param normalize Character string specifying normalization method:
#'   \describe{
#'     \item{"none"}{No normalization (default)}
#'     \item{"scale"}{Scale each volume to unit variance (divide by SD)}
#'     \item{"center"}{Center each volume (subtract mean)}
#'     \item{"z"}{Z-score normalization (center and scale)}
#'     \item{"unit"}{Scale to unit norm (L2 normalization)}
#'     \item{"percent"}{Convert to percent signal change from volume mean}
#'   }
#' @param normalize_by Character string specifying normalization scope:
#'   \describe{
#'     \item{"volume"}{Normalize within each volume (default)}
#'     \item{"voxel"}{Normalize each voxel across time}
#'   }
#' @param na.rm Logical, whether to remove NA values during averaging (default TRUE)
#' @param return_matrix Logical, if TRUE returns the averaged data matrix instead of NeuroVec (default FALSE)
#'
#' @return A NeuroVec object with averaged data (one volume per unique label), or a matrix if return_matrix=TRUE.
#'         The returned object has attributes:
#'         \itemize{
#'           \item{condition_labels}{The unique labels in order}
#'           \item{n_averaged}{Number of volumes averaged for each label}
#'           \item{normalization}{The normalization method used}
#'         }
#' 
#' @examples
#' \dontrun{
#' # Basic averaging
#' averaged <- average_labels(scandat, condition_labels, mask)
#' 
#' # With z-score normalization of each volume
#' averaged_norm <- average_labels(scandat, condition_labels, mask, 
#'                                 normalize = "z", normalize_by = "volume")
#' 
#' # Scale to unit norm for RSA
#' averaged_unit <- average_labels(scandat, condition_labels, mask,
#'                                 normalize = "unit")
#'                                 
#' # Get just the data matrix
#' data_mat <- average_labels(scandat, condition_labels, mask,
#'                            return_matrix = TRUE)
#' }
#' 
#' @importFrom neuroim2 series NeuroVec NeuroSpace space spacing
#' @importFrom assertthat assert_that
#' @export
average_labels <- function(neurovec, 
                          labels, 
                          mask = NULL,
                          normalize = c("none", "scale", "center", "z", "unit", "percent"),
                          normalize_by = c("volume", "voxel"),
                          na.rm = TRUE,
                          return_matrix = FALSE) {
  
  # Input validation
  normalize <- match.arg(normalize)
  normalize_by <- match.arg(normalize_by)
  
  if (!inherits(neurovec, c("NeuroVec", "NeuroVecSeq", "SparseNeuroVec"))) {
    stop("neurovec must be a NeuroVec, NeuroVecSeq, or SparseNeuroVec object")
  }
  
  # Get dimensions
  orig_space <- space(neurovec)
  spatial_dims <- dim(orig_space)[1:3]
  n_volumes <- dim(orig_space)[4]
  
  if (length(labels) != n_volumes) {
    stop(sprintf("Length of labels (%d) must match number of volumes (%d)", 
                 length(labels), n_volumes))
  }
  
  # Extract data matrix (volumes x voxels)
  if (!is.null(mask)) {
    # Get mask indices
    if (inherits(mask, "NeuroVol")) {
      mask_indices <- which(mask > 0)
    } else {
      mask_indices <- which(as.logical(mask))
    }
    # Use series for masked extraction
    data_mat <- neuroim2::series(neurovec, mask_indices)
  } else {
    # Extract all voxels using series
    all_indices <- seq_len(prod(spatial_dims))
    data_mat <- neuroim2::series(neurovec, all_indices)
  }
  
  # Apply normalization if requested
  if (normalize != "none") {
    data_mat <- normalize_data_matrix(data_mat, 
                                      method = normalize, 
                                      by = normalize_by)
  }
  
  # Get unique labels preserving order
  if (is.factor(labels)) {
    unique_labels <- levels(labels)
  } else {
    unique_labels <- unique(labels)
  }
  n_conditions <- length(unique_labels)
  
  # Average by group using group_means
  # Note: group_means with margin=1 averages over rows (volumes)
  averaged_mat <- group_means(data_mat, margin = 1, group = labels)
  
  # Ensure correct row ordering (group_means uses factor levels or appearance order)
  if (!is.null(rownames(averaged_mat))) {
    # Reorder to match unique_labels order
    averaged_mat <- averaged_mat[as.character(unique_labels), , drop = FALSE]
  }
  
  # Return matrix if requested
  if (return_matrix) {
    rownames(averaged_mat) <- unique_labels
    attr(averaged_mat, "normalization") <- normalize
    return(averaged_mat)
  }
  
  # Create new NeuroVec with averaged data
  new_dims <- c(spatial_dims, n_conditions)
  new_space <- NeuroSpace(new_dims, neuroim2::spacing(orig_space))
  
  # Reconstruct the full 4D array
  if (!is.null(mask)) {
    # Create full 4D array and insert averaged data at mask locations
    full_array <- array(NA, dim = new_dims)
    
    # Get mask indices - handle both NeuroVol and logical array
    if (inherits(mask, "NeuroVol")) {
      mask_indices <- which(mask > 0)
    } else {
      mask_indices <- which(as.logical(mask))
    }
    
    # For each condition, insert the averaged values
    for (i in seq_len(n_conditions)) {
      vol_data <- array(NA, dim = spatial_dims)
      vol_data[mask_indices] <- averaged_mat[i, ]
      full_array[,,,i] <- vol_data
    }
    
    averaged_vec <- NeuroVec(full_array, new_space)
  } else {
    # Reshape averaged matrix to 4D array
    # Need to transpose because we want volumes in the 4th dimension
    full_array <- array(t(averaged_mat), dim = new_dims)
    averaged_vec <- NeuroVec(full_array, new_space)
  }
  
  # Add metadata as attributes
  attr(averaged_vec, "condition_labels") <- unique_labels
  attr(averaged_vec, "n_averaged") <- table(labels)[unique_labels]
  attr(averaged_vec, "normalization") <- normalize
  
  averaged_vec
}


#' Normalize Data Matrix
#'
#' Internal helper function to normalize a data matrix using various methods.
#' Supports normalization by row (volume) or column (voxel).
#'
#' @param data_mat Matrix to normalize (rows = observations, cols = features)
#' @param method Normalization method ("center", "scale", "z", "unit", "percent")
#' @param by Whether to normalize by row ("volume") or column ("voxel")
#' @return Normalized matrix
#' @keywords internal
#' @noRd
normalize_data_matrix <- function(data_mat, method, by) {
  
  if (by == "volume") {
    # Normalize each row (volume)
    
    if (method == "center") {
      # Center each volume
      row_means <- rowMeans(data_mat, na.rm = TRUE)
      data_mat <- sweep(data_mat, 1, row_means, "-")
      
    } else if (method == "scale") {
      # Scale each volume to unit variance
      row_sds <- apply(data_mat, 1, sd, na.rm = TRUE)
      row_sds[row_sds == 0] <- 1  # Avoid division by zero
      data_mat <- sweep(data_mat, 1, row_sds, "/")
      
    } else if (method == "z") {
      # Z-score each volume
      row_means <- rowMeans(data_mat, na.rm = TRUE)
      row_sds <- apply(data_mat, 1, sd, na.rm = TRUE)
      row_sds[row_sds == 0] <- 1
      data_mat <- sweep(data_mat, 1, row_means, "-")
      data_mat <- sweep(data_mat, 1, row_sds, "/")
      
    } else if (method == "unit") {
      # L2 normalization - scale each volume to unit norm
      row_norms <- sqrt(rowSums(data_mat^2, na.rm = TRUE))
      row_norms[row_norms == 0] <- 1
      data_mat <- sweep(data_mat, 1, row_norms, "/")
      
    } else if (method == "percent") {
      # Percent signal change from volume mean
      row_means <- rowMeans(data_mat, na.rm = TRUE)
      row_means[row_means == 0] <- 1  # Avoid division by zero
      data_mat <- sweep(data_mat, 1, row_means, "-")
      data_mat <- sweep(data_mat, 1, row_means, "/") * 100
    }
    
  } else {  # by == "voxel"
    # Normalize each column (voxel across time)
    
    if (method == "center") {
      # Center each voxel
      col_means <- colMeans(data_mat, na.rm = TRUE)
      data_mat <- sweep(data_mat, 2, col_means, "-")
      
    } else if (method == "scale") {
      # Scale each voxel to unit variance
      col_sds <- apply(data_mat, 2, sd, na.rm = TRUE)
      col_sds[col_sds == 0] <- 1
      data_mat <- sweep(data_mat, 2, col_sds, "/")
      
    } else if (method == "z") {
      # Z-score each voxel
      col_means <- colMeans(data_mat, na.rm = TRUE)
      col_sds <- apply(data_mat, 2, sd, na.rm = TRUE)
      col_sds[col_sds == 0] <- 1
      data_mat <- sweep(data_mat, 2, col_means, "-")
      data_mat <- sweep(data_mat, 2, col_sds, "/")
      
    } else if (method == "unit") {
      # L2 normalization per voxel (less common)
      col_norms <- sqrt(colSums(data_mat^2, na.rm = TRUE))
      col_norms[col_norms == 0] <- 1
      data_mat <- sweep(data_mat, 2, col_norms, "/")
      
    } else if (method == "percent") {
      # Percent signal change from voxel mean
      col_means <- colMeans(data_mat, na.rm = TRUE)
      col_means[col_means == 0] <- 1
      data_mat <- sweep(data_mat, 2, col_means, "-")
      data_mat <- sweep(data_mat, 2, col_means, "/") * 100
    }
  }
  
  data_mat
}