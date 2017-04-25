
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

#' @param D the data dimensions
#' @param nobs the number of observations
#' @param response_type 'categorical' or 'continuous'
#' @param data_mode 'image' or 'surface'
#' @param spacing
#' @param blocks
#' @param nlevels
gen_sample_dataset <- function(D, nobs, response_type=c("categorical", "continuous"), data_mode=c("image", "surface"),
                              spacing=c(1,1,1), blocks=5, nlevels=5) {
  
  response_type <- match.arg(response_type)
  data_mode <- match.arg(data_mode)
  
  mat <- array(rnorm(prod(D)*nobs), c(D,nobs))
  bspace <- BrainSpace(c(D,nobs), spacing)
  bvec <- BrainVector(mat, bspace)
  
  mask <- as.logical(BrainVolume(array(rep(1, prod(D)), D), BrainSpace(D, spacing)))
  
  Y <- if (response_type == "categorical") {
    sample(factor(rep(letters[1:nlevels], length.out=nobs)))
  } else {
    rnorm(length(obs))
  }
  
  block_var <- rep(1:folds, length.out=nobs)
  des <- mvpa_design(data.frame(Y=Y), block_var=block_var, y_train= ~ Y)
  dset <- mvpa_dataset(bvec, mask=mask)
  list(dataset=dset, design=des)
}




#' mvpa_dataset
#' 
#' @param train_data the training data set: a \code{BrainVector} instance
#' @param test_data the test data set: a \code{BrainVector} instance
#' @param mask the set of voxels to include: a \code{BrainVolume} instance
#' @importFrom assertthat assert_that
mvpa_dataset <- function(train_data,test_data=NULL, mask) {
  assert_that(inherits(design, "mvpa_design"))
  
  ret <- list(
    train_data=train_data,
    test_data=test_data,
    mask=mask
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
mvpa_surface_dataset <- function(train_data, test_data=NULL, mask=NULL, name="") {
  
  
  
  if (is.null(mask)) {
    mask <- numeric(length(nodes(train_data@geometry)))
    mask[indices(train_data)] <- 1
  }
  
  ret <- list(
    train_data=train_data,
    test_data=test_data,
    mask=mask,
    name=name
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
  
  cat("  name: ", x$name, "\n")
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



