
#' @keywords internal
gen_id <- function(n) {
  width <- nchar(n)
  sprintf(paste0("%0", width, "d"), seq_len(n))
}

#' @keywords internal
.get_samples <- function(obj, voxlist) {
  ret <- lapply(voxlist, function(vox) {
    sam <- data_sample(obj, vox)
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


#' @export
data_sample.mvpa_dataset <- function(obj, vox) {
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

#' Filter Region of Interest (ROI)
#'
#' This function filters an ROI, keeping only valid columns.
#'
#' @param roi A list containing the train and test ROI data.
#' @return A list with filtered train and test ROI data.
#' @details
#' The function filters an ROI by removing columns with missing values (NA) and zero standard deviation.
#' It returns a list with filtered train and test ROI data.
#' @keywords internal
filter_roi <- function(roi) {
  # Extract the train data values
  trdat <- values(roi$train_roi)
  
  # Find columns with missing values (NA)
  nas <- apply(trdat, 2, function(v) any(is.na(v)))
  
  # Find columns with non-zero standard deviation
  sdnonzero <- apply(trdat, 2, sd, na.rm=TRUE) > 0
  
  # Determine columns to keep
  keep <- !nas & sdnonzero
  
  # If no valid columns are found, throw an error
  if (sum(keep) == 0) {
    stop("filter_roi: roi has no valid columns")
  }
  
  # If there's no test ROI data, return filtered train ROI data only
  if (is.null(roi$test_roi)) {
    troi <- neuroim2::ROIVec(space(roi$train_roi), coords(roi$train_roi)[keep,,drop=FALSE], data=trdat[,keep,drop=FALSE])
    list(train_roi=troi,
         test_roi=NULL)
  } else  {
    # Filter train ROI data
    troi <- neuroim2::ROIVec(space(roi$train_roi), coords(roi$train_roi)[keep,,drop=FALSE], data=trdat[,keep,drop=FALSE])
    
    # Filter test ROI data
    tedat <- values(roi$test_roi)
    teroi <- neuroim2::ROIVec(space(roi$test_roi), coords(roi$test_roi)[keep,,drop=FALSE], data=tedat[,keep,drop=FALSE])
    
    # Return filtered train and test ROI data
    list(train_roi=troi,
         test_roi=teroi)
  }
}


#' @keywords internal
#' @noRd
#' @importFrom neuroim2 series_roi
as_roi.data_sample <- function(x, data, ...) {
  
  train_roi <- try(series_roi(data$train_data, x$vox))
  
  test_roi <- if (has_test_set(data)) {
    series_roi(data$test_data, x$vox)
  }
  
  #cds <- if (is.vector(x$vox)) {
  #  cds <- indexToGrid(space(x$data$mask), x$vox)
  #} else {
  #  x$vox
  #}

  if (is.null(test_roi)) {
    list(train_roi=train_roi,
         test_roi=NULL)
  } else {
    list(train_roi=train_roi,
         test_roi=test_roi)
  }
  
  
}

#' @keywords internal
#' @noRd
#' @importFrom neuroim2 space series series_roi
as.data.frame.data_sample <- function(x, data, ...) {
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
  
 
