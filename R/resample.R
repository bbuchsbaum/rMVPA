#' @export
#' @param block_var the blocking variable (an integer vector)
#' @rdname crossv_block
crossv_block <- function(data, block_var, id = ".id", exclude_block=NULL) {
  
  if (!length(block_var) == nrow(data)) {
    stop("length of `block_var` must be equal to row(data).", call. = FALSE)
  }
  
  if (!is.null(exclude_block)) {
    idx <- seq_len(nrow(data))
    keep <- block_var != exclude_block
    idx <- idx[keep]
    fold_idx <- split(idx, block_var[keep])
  } else {
    idx <- seq_len(nrow(data))
    fold_idx <- split(idx, block_var)
  }
  
  
  fold <- function(test) {
    list(
      train = resample(data, setdiff(idx, test)),
      test = resample(data, test)
    )
  }
  
  cols <- purrr::transpose(purrr::map(fold_idx, fold))
  cols[[id]] <- modelr:::id(length(fold_idx))
  
  tibble::as_data_frame(cols)
}





#' @export
get_samples.mvpa_dataset <- function(obj, voxlist) {
  ret <- lapply(voxlist, function(vox) {
    sam <- data_sample(obj, vox)
  })
  
  n <- length(ret)
  df <- tibble::data_frame(sample = ret)
  df[[id]] <- modelr:::id(n)
  df
  
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

#' @export
as_roi.data_sample <- function(x, ...) {
  
  train_mat <- series(x$data$train_data, x$vox)
  
  test_mat <- if (has_test_set(x$data)) {
    series(x$data$test_data, x$vox)
  }
  
  cds <- if (is.vector(x$vox)) {
    cds <- indexToGrid(space(x$data$mask), x$vox)
  } else {
    x$vox
  }

  if (is.null(test_mat)) {
    list(train_roi=neuroim::ROIVolume(space(x$data$mask), cds, data=train_mat),
         test_roi=NULL)
  } else {
    list(train_roi=neuroim::ROIVolume(space(x$data$mask), cds, data=train_mat),
        test_roi=neuroim::ROIVolume(space(x$data$mask), cds, data=test_mat))
  }
  
  
}

#' @export
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
  
 
