
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

has_test_set.mvpa_dataset <- function(obj) {
  !is.null(obj$design$y_test) 
}

y_train.mvpa_dataset <- function(obj) y_train(obj$design)

y_test.mvpa_dataset <- function(obj) y_test(obj$design)

#' mvpa_dataset
#' @param train_data
#' @param test_data
#' @param mask
#' @param design
#' @importFrom assertthat assert_that
mvpa_dataset <- function(train_data,test_data=NULL, mask, design) {
  assert_that(inherits(design, "mvpa_design"))
  
  ret <- list(
    train_data=train_data,
    test_data=test_data,
    mask=mask,
    design=design
  )
  
  class(ret) <- c("mvpa_dataset", "list")
  ret
    
}


#' mvpa_surface_dataset
#' 
#' @param train_data
#' @param test_data
#' @param mask
#' @param design
#' @importFrom assertthat assert_that
mvpa_surface_dataset <- function(train_data,test_data=NULL, mask, design) {
  assert_that(inherits(design, "mvpa_design"))
  
  ret <- list(
    train_data=train_data,
    test_data=test_data,
    mask=mask,
    design=design
  )
  
  class(ret) <- c("mvpa_surface_dataset", "list")
  ret
  
}


#' @export
get_searchlight.mvpa_dataset <- function(x, type=c("standard", "randomized"), radius=8) {
  type <- match.arg(type)
  if (type == "standard") {
    neuroim::Searchlight(x$mask, radius=8)
  } else {
    neuroim::RandomSearchlight(x$mask, radius=8)
  }
}



#' @export
get_searchlight.mvpa_surface_dataset <- function(x, type=c("standard", "randomized"), radius=8) {
  type <- match.arg(type)
  if (type == "standard") {
    neurosurf::SurfaceSearchlight(geometry(x$train_data), radius)
  } else {
    neurosurf::RandomSurfaceSearchlight(geometry(x$train_data), radius)
  }
}