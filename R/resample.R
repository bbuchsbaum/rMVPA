#' @keywords internal
gen_id <- function(n) {
  width <- max(2, nchar(n))  # Minimum width of 2 for consistent formatting
  sprintf(paste0("%0", width, "d"), seq_len(n))
}

#' @keywords internal
.get_samples <- function(obj, voxlist) {
  ret <- lapply(voxlist, function(vox) {
    # Expand searchlight windows to full voxel indices when possible (duck-typed)
    vox_expanded <- tryCatch({
      coords <- tryCatch(vox@coords, error = function(...) NULL)
      sp     <- tryCatch(vox@space,  error = function(...) NULL)
      if (!is.null(coords) && !is.null(sp)) {
        neuroim2::grid_to_index(sp, coords)
      } else {
        vox
      }
    }, error = function(...) vox)

    sam <- data_sample(obj, vox_expanded)
  })
  
  n <- length(ret)
  df <- tibble::tibble(sample = ret)
  df[[".id"]] <- gen_id(n)
  df
}

#' @export
get_samples.mvpa_dataset <- function(obj, vox_list) {
  .get_samples(obj, vox_list)
}

#' @export
get_samples.mvpa_surface_dataset <- function(obj, vox_list) {
  .get_samples(obj, vox_list)
}


#' @rdname data_sample
#' @export
data_sample.mvpa_dataset <- function(obj, vox, ...) {
  structure(
    list(
      #data = obj,
      data=NULL,
      vox=vox
    ),
    class = "data_sample"
  )
}


#' @export
print.data_sample <- function(x, ...) {
  if (is.matrix(x$vox)) {
    cat("data sample with : ", nrow(x$vox), "features")
  } else {
    cat("data sample with : ", length(x$vox), "features")
  }
}

#' @keywords internal
filter_roi.default <- function(roi, preserve = NULL, ...) {
  stop("Unsupported ROI type")
}

#' @keywords internal
#' @importFrom neuroim2 ROIVec space coords values
filter_roi.ROIVec <- function(roi, preserve = NULL, min_voxels = 2, ...) {
  # Extract the train data values
  trdat <- values(roi$train_roi)

  # Find columns with missing values (NA) - vectorized for performance
  # colSums is much faster than apply() for this operation

  nas <- colSums(is.na(trdat)) > 0

  # Find columns with non-zero standard deviation - vectorized
  # Use matrixStats if available for best performance, otherwise use vectorized base R
  if (requireNamespace("matrixStats", quietly = TRUE)) {
    sdnonzero <- matrixStats::colSds(trdat, na.rm = TRUE) > 0
  } else {
    # Vectorized variance calculation: Var(X) = E[X^2] - E[X]^2
    n <- nrow(trdat)
    col_means <- colMeans(trdat, na.rm = TRUE)
    col_sq_means <- colMeans(trdat^2, na.rm = TRUE)
    col_vars <- col_sq_means - col_means^2
    # Adjust for sample variance (n-1) and handle edge cases
    col_vars <- col_vars * n / (n - 1)
    col_vars[is.na(col_vars) | col_vars < 0] <- 0
    sdnonzero <- col_vars > .Machine$double.eps
  }

  # Determine columns to keep
  keep <- !nas & sdnonzero

  # Preserve specified voxel if requested (e.g., center voxel in searchlight)
  if (!is.null(preserve)) {
    gi <- neuroim2::indices(roi$train_roi)
    kp <- match(preserve, gi)
    if (!is.na(kp) && !keep[kp]) {
      futile.logger::flog.debug("Preserving voxel %s that would have been filtered (NA or zero variance)", preserve)
      keep[kp] <- TRUE
    }
  }

  # Need at least 2 valid columns for multivariate analysis
  # (even with preserve, searchlight edges need multiple voxels)
  if (sum(keep) < min_voxels) {
    reasons <- c()
    if (any(nas)) reasons <- c(reasons, sprintf("NA voxels: %d", sum(nas)))
    if (any(!sdnonzero)) reasons <- c(reasons, sprintf("zero-variance voxels: %d", sum(!sdnonzero)))
    futile.logger::flog.debug(
      "filter_roi.ROIVec: %d valid columns (< %d). Reasons: %s",
      sum(keep), min_voxels,
      paste(reasons, collapse = "; ")
    )
    stop(sprintf("filter_roi: roi has %d valid columns (need at least %d)", sum(keep), min_voxels))
  }
  
  # If there's no test ROI data, return filtered train ROI data only
  if (is.null(roi$test_roi)) {
    troi <- ROIVec(space(roi$train_roi), 
                   coords(roi$train_roi)[keep,,drop=FALSE], 
                   data=trdat[,keep,drop=FALSE])
    list(train_roi=troi, test_roi=NULL)
  } else {
    # Filter train ROI data
    troi <- ROIVec(space(roi$train_roi), 
                   coords(roi$train_roi)[keep,,drop=FALSE], 
                   data=trdat[,keep,drop=FALSE])
    
    # Filter test ROI data
    tedat <- values(roi$test_roi)
    teroi <- ROIVec(space(roi$test_roi), 
                    coords(roi$test_roi)[keep,,drop=FALSE], 
                    data=tedat[,keep,drop=FALSE])
    
    list(train_roi=troi, test_roi=teroi)
  }
}

#' @keywords internal
#' @noRd
filter_roi.list <- function(roi, preserve = NULL, min_voxels = 2, ...) {
  # Expect a list with train_roi/test_roi as produced by as_roi()
  if (!is.list(roi) || is.null(roi$train_roi)) {
    stop("filter_roi.list: expected list with train_roi/test_roi")
  }
  filtered <- filter_roi(roi$train_roi, preserve = preserve, min_voxels = min_voxels, ...)
  list(train_roi = filtered$train_roi,
       test_roi  = if (!is.null(roi$test_roi)) {
         # Apply the same voxel filter to the test ROI if present
         # Use the logical mask from filtered train ROI indices
         keep_idx <- neuroim2::indices(filtered$train_roi)
         orig_idx <- neuroim2::indices(roi$train_roi)
         kp <- match(keep_idx, orig_idx)
         # Guard against mismatch
         if (any(is.na(kp))) {
           roi$test_roi
         } else {
           neuroim2::ROIVec(neuroim2::space(roi$test_roi),
                            neuroim2::coords(roi$test_roi)[kp,,drop=FALSE],
                            data = neuroim2::values(roi$test_roi)[, kp, drop = FALSE])
         }
       } else NULL)
}

#' @keywords internal
#' @importFrom neurosurf ROISurfaceVector geometry nodes
filter_roi.ROISurfaceVector <- function(roi, preserve = NULL, min_voxels = 2, ...) {
  # Extract the train data values
  trdat <- roi$train_roi@data

  # Find columns with missing values (NA) - vectorized for performance
  nas <- colSums(is.na(trdat)) > 0

  # Find columns with non-zero standard deviation - vectorized
  if (requireNamespace("matrixStats", quietly = TRUE)) {
    sdnonzero <- matrixStats::colSds(trdat, na.rm = TRUE) > 0
  } else {
    n <- nrow(trdat)
    col_means <- colMeans(trdat, na.rm = TRUE)
    col_sq_means <- colMeans(trdat^2, na.rm = TRUE)
    col_vars <- col_sq_means - col_means^2
    col_vars <- col_vars * n / (n - 1)
    col_vars[is.na(col_vars) | col_vars < 0] <- 0
    sdnonzero <- col_vars > .Machine$double.eps
  }

  # Determine columns to keep
  keep <- !nas & sdnonzero
  
  # Preserve specified voxel if requested (e.g., center voxel in searchlight)
  if (!is.null(preserve)) {
    gi <- roi$train_roi@indices
    kp <- match(preserve, gi)
    if (!is.na(kp) && !keep[kp]) {
      futile.logger::flog.debug("Preserving voxel %s that would have been filtered (NA or zero variance)", preserve)
      keep[kp] <- TRUE
    }
  }

  # Need at least 2 valid columns for multivariate analysis
  # (even with preserve, searchlight edges need multiple voxels)
  if (sum(keep) < min_voxels) {
    stop(sprintf("filter_roi: roi has %d valid columns (need at least %d)", sum(keep), min_voxels))
  }
  
  # If there's no test ROI data, return filtered train ROI data only
  if (is.null(roi$test_roi)) {
    troi <- ROISurfaceVector(geometry=roi$train_roi@geometry,
                            indices=roi$train_roi@indices[keep],
                            data=trdat[,keep,drop=FALSE])
    list(train_roi=troi, test_roi=NULL)
  } else {
    # Filter train ROI data
    troi <- ROISurfaceVector(geometry=roi$train_roi@geometry,
                            indices=roi$train_roi@indices[keep],
                            data=trdat[,keep,drop=FALSE])
    
    # Filter test ROI data
    tedat <- roi$test_roi@data
    teroi <- ROISurfaceVector(geometry=roi$test_roi@geometry,
                             indices=roi$test_roi@indices[keep],
                             data=tedat[,keep,drop=FALSE])
    
    list(train_roi=troi, test_roi=teroi)
  }
}


#' @keywords internal
#' @noRd
compute_mask_indices <- function(mask) {
  if (inherits(mask, c("NeuroVol", "NeuroVec"))) {
    vals <- neuroim2::values(mask)
    if (is.matrix(vals)) {
      vals <- vals[, 1, drop = TRUE]
    }
    which(vals != 0)
  } else {
    which(mask != 0)
  }
}

#' @keywords internal
#' @noRd
#' @importFrom neuroim2 series_roi
#' @method as_roi data_sample
#' @export
as_roi.data_sample <- function(obj, data, ...) {

  # Defensive validation: ensure all voxel indices are within the mask.
  # Use any precomputed mask_indices if present to avoid O(#voxels) work per ROI.
  if (!is.null(data$mask)) {
    mask_indices <- data$mask_indices
    if (is.null(mask_indices)) {
      mask_indices <- compute_mask_indices(data$mask)
    }

    invalid_vox <- setdiff(obj$vox, mask_indices)

    if (length(invalid_vox) > 0) {
      futile.logger::flog.debug("as_roi.data_sample: %d voxel indices are outside the mask. Filtering to mask (normal for edge searchlights).",
                              length(invalid_vox))
      # Filter to only valid mask indices
      obj$vox <- intersect(obj$vox, mask_indices)

      # Don't fail here for <2 voxels - let filter_roi and extract_roi handle that
      # This is important for searchlights where partial spheres at brain edges are expected
      if (length(obj$vox) == 0) {
        futile.logger::flog.debug("as_roi.data_sample: After mask filtering, ROI has 0 voxels. Returning error ROI.")
        # Return a try-error to signal failure
        train_roi <- structure(list(message = "No voxels remaining after mask filtering"),
                              class = "try-error")
        return(list(train_roi = train_roi, test_roi = NULL))
      }
    }
  }

  train_roi <- try(series_roi(data$train_data, obj$vox), silent = TRUE)

  test_roi <- if (has_test_set(data)) {
    try(series_roi(data$test_data, obj$vox), silent = TRUE)
  }

  if (is.null(test_roi)) {
    list(train_roi = train_roi, test_roi = NULL)
  } else {
    list(train_roi = train_roi, test_roi = test_roi)
  }
}

#' @keywords internal
#' @noRd
#' @method as_roi numeric
#' @export
as_roi.numeric <- function(obj, data, ...) {
  ds <- data_sample(data, obj)
  as_roi(ds, data, ...)
}

#' @keywords internal
#' @noRd
#' @method as_roi integer
#' @export
as_roi.integer <- function(obj, data, ...) {
  ds <- data_sample(data, obj)
  as_roi(ds, data, ...)
}

#' @keywords internal
#' @noRd
#' @importFrom neuroim2 space series series_roi
as.data.frame.data_sample <- function(x, row.names = NULL, optional = FALSE, ..., data) {
  train_mat <- neuroim2::series(data$train_data, x$vox)
  
  test_mat <- if (has_test_set(data)) {
    neuroim2::series(data$test_data, x$vox)
  }
  
  cds <- if (is.vector(x$vox)) {
    cds <- neuroim2::index_to_grid(space(data$mask), x$vox)
  } else {
    x$vox
  }
  
  if (!is.null(test_mat)) {
    .type <- rep(c("train", "test"), c(nrow(train_mat), nrow(test_mat)))
    ret <- as.data.frame(rbind(train_mat, test_mat))
    ret$.type <- .type
    attr(ret, "vox") <- cds
    ret
  } else {
    .type <- rep(c("train"), nrow(train_mat))
    ret <- as.data.frame(train_mat)
    ret$.type <- .type
    attr(ret, "vox") <- cds
    ret
  }
  
}
  
 
