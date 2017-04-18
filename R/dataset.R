
roi_volume_matrix <- function(mat, refspace, indices, coords) {
  structure(mat,
            refspace=refspace,
            indices=indices,
            coords=coords,
            class=c("roi_volume_matrix", "matrix"))
            
}

roi_surface_matrix <- function(mat, refspace, indices, coords) {
  structure(mat,
            refspace=refspace,
            indices=indices,
            coords=coords,
            class=c("roi_surface_matrix", "matrix"))
  
}



#' mvpa_dataset
#' 
#' @param train_data the training data set: a \code{BrainVector} instance
#' @param test_data the test data set: a \code{BrainVector} instance
#' @param mask the set of voxels to include: a \code{BrainVolume} instance
#' @param design the design: a \code{mvpa_design} instance
#' @importFrom assertthat assert_that
mvpa_dataset <- function(train_data,test_data=NULL, mask, design) {
  assert_that(inherits(design, "mvpa_design"))
  
  ret <- list(
    train_data=train_data,
    test_data=test_data,
    mask=mask,
    design=design
  )
  
  class(ret) <- c("mvpa_image_dataset", "mvpa_dataset", "list")
  ret
    
}


#' mvpa_surface_dataset
#' 
#' @param train_data
#' @param test_data
#' @param mask
#' @param design
#' @param hemisphere
#' @importFrom assertthat assert_that
#' @export
mvpa_surface_dataset <- function(train_data, test_data=NULL, mask=NULL, hemisphere=c("lh", "rh")) {
  
  hemisphere <- match.arg(hemisphere)
  
  if (is.null(mask)) {
    mask <- numeric(length(nodes(train_data@geometry)))
    mask[indices(train_data)] <- 1
  }
  
  ret <- list(
    train_data=train_data,
    test_data=test_data,
    mask=mask,
    hemisphere=hemisphere
  )
  
  class(ret) <- c("mvpa_surface_dataset", "mvpa_dataset", "list")
  ret
  
}

print.mvpa_surface_dataset <- function(x) {
  cat("mvpa_surface_dataset:", "\n")
  cat("  train_data: ")
  print(x$train_data)
  if (is.null(x$test_data)) {
    cat("  test_data: none. \n")
  } else {
    cat("\n")
    cat("  test_data: ")
    print(x$test_data)
  }
  
  mids <- table(x$mask[x$mask!=0])
  midstr <- paste0(names(mids), "/", mids)
  
  cat("  hemisphere: ", x$hemisphere, "\n")
  cat("  mask areas: ", midstr, "\n")
  cat("  mask cardinality: ", sum(x$mask>0), "\n")
}





#' @export
get_searchlight.mvpa_dataset <- function(x, type=c("standard", "randomized"), radius=8) {
  type <- match.arg(type)
  if (type == "standard") {
    neuroim::Searchlight(x$mask, radius=radius)
  } else {
    neuroim::RandomSearchlight(x$mask, radius=radius)
  }
}



#' @export
get_searchlight.mvpa_surface_dataset <- function(x, type=c("standard", "randomized"), radius=8) {
  type <- match.arg(type)
  if (type == "standard") {
    neurosurf::SurfaceSearchlight(geometry(x$train_data), radius, nodeset=which(x$mask>0))
  } else {
    neurosurf::RandomSurfaceSearchlight(geometry(x$train_data), radius, nodeset=which(x$mask>0))
  }
}

wrap_output.mvpa_dataset <- function(obj, vals, indices) {
  BrainVolume(vals, space(obj$mask), indices=indices)
}

wrap_output.mvpa_surface_dataset <- function(obj, vals, indices) {
  BrainSurface(geometry=geometry(obj$train_data), indices=indices, data=vals)
}

#' @export
has_test_set.mvpa_dataset <- function(obj) {
  !is.null(obj$test_data) 
}



