

.get_samples <- function(obj, voxlist) {
  ret <- lapply(voxlist, function(vox) {
    sam <- data_sample(obj, vox)
  })
  
  n <- length(ret)
  df <- tibble::data_frame(sample = ret)
  df[[".id"]] <- modelr:::id(n)
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
      data = obj,
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
#' @noRd
as_roi.data_sample <- function(x, ...) {
  
  train_roi <- try(series_roi(x$data$train_data, x$vox))
  
  test_roi <- if (has_test_set(x$data)) {
    series_roi(x$data$test_data, x$vox)
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
as.data.frame.data_sample <- function(x, ...) {
  train_mat <- series(x$data$train_data, x$vox)
  
  test_mat <- if (has_test_set(x$data)) {
    series(x$data$test_data, x$vox)
  }
  
  cds <- if (is.vector(x$vox)) {
    cds <- indexToGrid(space(x$data$mask), x$vox)
  } else {
    x$vox
  }
  
  if (!is.null(test_mat)) {
    .type <- rep(c("train", "test"), c(nrow(trainmat), nrow(testmat)))
    ret <- as.data.frame(rbind(train_mat, test_mat))
    ret$.type <- .type
    attr(ret, "vox") <- cds
    ret
  } else {
    .type <- rep(c("train"), nrow(trainmat))
    ret <- as.data.frame(train_mat)
    ret$.type <- .type
    attr(ret, "vox") <- cds
    ret
  }
  
}
  
 
