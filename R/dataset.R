

#' @keywords @internal
roi_volume_matrix <- function(mat, refspace, indices, coords) {
  structure(mat,
            refspace=refspace,
            indices=indices,
            coords=coords,
            class=c("roi_volume_matrix", "matrix"))
            
}

#' @keywords @internal
roi_surface_matrix <- function(mat, refspace, indices, coords) {
  structure(mat,
            refspace=refspace,
            indices=indices,
            coords=coords,
            class=c("roi_surface_matrix", "matrix"))
  
}

#' gen_sample_dataset
#' 
#' @param D the data dimension(s)
#' @param nobs the number of observations
#' @param response_type 'categorical' or 'continuous'
#' @param data_mode 'image' or 'surface'
#' @param spacing the voxel spacing
#' @param blocks the number of 'blocks' in the data
#' @param nlevels the number of category levels
#' @param external_test is the test set 'external' to the training set
#' @param ntest_obs number of test observations (only relevant if \code{external_test} is true)
#' @param split_by compute performance measures for each level of factor
#' @export
gen_sample_dataset <- function(D, nobs, response_type=c("categorical", "continuous"), data_mode=c("image", "surface"),
                              spacing=c(1,1,1), blocks=5, nlevels=5, external_test=FALSE, ntest_obs=nobs, split_by=NULL) {
  
  response_type <- match.arg(response_type)
  data_mode <- match.arg(data_mode)
  
  
  if (data_mode == "image") {
    mat <- array(rnorm(prod(D)*nobs), c(D,nobs))
    bspace <- neuroim::BrainSpace(c(D,nobs), spacing)
    bvec <- neuroim::BrainVector(mat, bspace)
    mask <- as.logical(neuroim::BrainVolume(array(rep(1, prod(D)), D), neuroim::BrainSpace(D, spacing)))
    
    if (external_test) {
      mat <- array(rnorm(prod(D)*ntest_obs), c(D,ntest_obs))
      bspace <- neuroim::BrainSpace(c(D,ntest_obs), spacing)
      testvec <- neuroim::BrainVector(mat, bspace)
      dset <- mvpa_dataset(train_data=bvec, test_data=testvec, mask=mask)
    } else {
      dset <- mvpa_dataset(train_data=bvec,mask=mask)
    }
  } else {
    fname <- system.file("extdata/std.lh.smoothwm.asc", package="neuroim")
    geom <- neurosurf::loadSurface(fname)
    nvert <- nrow(neurosurf::vertices(geom))
    mat <- matrix(rnorm(nvert*nobs), nvert, nobs)
    bvec <- neurosurf::BrainSurfaceVector(geom, 1:nvert, mat)
    
    if (external_test) {
      test_data <- neurosurf::BrainSurfaceVector(geom, 1:nvert, matrix(rnorm(nvert*ntest_obs), nvert, ntest_obs))
      dset <- mvpa_surface_dataset(train_data=bvec, test_data=test_data)
    } else {
      dset <- mvpa_surface_dataset(train_data=bvec)
    }
  }
  
  Y <- if (response_type == "categorical") {
    sample(factor(rep(letters[1:nlevels], length.out=nobs)))
  } else {
    rnorm(nobs)
  }
  
  Ytest <- if (response_type == "categorical") {
    sample(factor(rep(letters[1:nlevels], length.out=ntest_obs)))
  } else {
    rnorm(ntest_obs)
  }
  
  block_var <- rep(1:blocks, length.out=nobs)
  
  des <- if (external_test) {
    message("external test")
    mvpa_design(data.frame(Y=Y, block_var=block_var), test_design=data.frame(Ytest = Ytest), 
                       block_var= "block_var", y_train= ~ Y, y_test = ~ Ytest, split_by=split_by)
  } else {
    mvpa_design(data.frame(Y=Y, block_var=block_var), block_var="block_var", y_train= ~ Y, split_by=split_by)
  }
  
  list(dataset=dset, design=des)
}




#' mvpa_dataset
#' 
#' A data structure that encapsulate a standard (volumetric) training dataset, an optional test dataset and a voxel 'mask'.
#' 
#' @param train_data the training data set: a \code{BrainVector} instance
#' @param test_data the test data set: a \code{BrainVector} instance
#' @param mask the set of voxels to include: a \code{BrainVolume} instance
#' @importFrom assertthat assert_that
#' @export
mvpa_dataset <- function(train_data, test_data=NULL, mask) {
  assert_that(inherits(train_data, "BrainVector"))
  if (!is.null(test_data)) {
    assert_that(inherits(test_data, "BrainVector"))
  }
  ret <- structure(
    list(
      train_data=train_data,
      test_data=test_data,
      mask=mask
    ),
    class=c("mvpa_image_dataset", "mvpa_dataset", "list")
  )

}


#' mvpa_surface_dataset
#' 
#' Cnstruct a MVPA dataset with a surface-based dataset.
#' 
#' @param train_data the training data, must inherit from \code{BrainSurfaceVector} in \code{neurosurf} package.
#' @param test_data optional test data, must inherit from \code{BrainSurfaceVector} in \code{neurosurf} package.
#' @param mask a binary mask equal to the number of nodes in the training/test data set.
#' @param name label to identify the dataset (e.g. "lh" or "rh" to indicate hemisphere)
#' @importFrom assertthat assert_that
#' @export
mvpa_surface_dataset <- function(train_data, test_data=NULL, mask=NULL, name="") {
  
  assert_that(inherits(train_data, "BrainSurfaceVector"))
  
  if (!is.null(test_data)) {
    assert_that(inherits(test_data, "BrainSurfaceVector"))
  }
  
  if (is.null(mask)) {
    mask <- numeric(length(nodes(train_data@geometry)))
    mask[indices(train_data)] <- 1
  }
  
  structure(
    list(
      train_data=train_data,
      test_data=test_data,
      mask=mask,
      name=name
    ),
    class=c("mvpa_surface_dataset", "mvpa_dataset", "list")
  )

  
}

#' @export
print.mvpa_dataset <- function(x,...) {
  cat("mvpa_dataset:", "\n")
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
  cat("  mask areas: ", midstr, "\n")
  cat("  mask cardinality: ", sum(x$mask>0), "\n")
}

#' @export
print.mvpa_surface_dataset <- function(x,...) {
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
get_searchlight.mvpa_dataset <- function(obj, type=c("standard", "randomized"), radius=8,...) {
  type <- match.arg(type)
  if (type == "standard") {
    neuroim::Searchlight(obj$mask, radius=radius)
  } else {
    neuroim::RandomSearchlight(obj$mask, radius=radius)
  }
}



#' @export

get_searchlight.mvpa_surface_dataset <- function(obj, type=c("standard", "randomized"), radius=8,...) {
  type <- match.arg(type)
  if (type == "standard") {
    neurosurf::SurfaceSearchlight(geometry(obj$train_data), radius, nodeset=which(obj$mask>0))
  } else {
    neurosurf::RandomSurfaceSearchlight(geometry(obj$train_data), radius, nodeset=which(obj$mask>0))
  }
}

#' @keywords internal
#' @noRd
#' @import neuroim
wrap_output.mvpa_dataset <- function(obj, vals, indices) {
  BrainVolume(vals[indices], space(obj$mask), indices=indices)
}


#' @keywords internal
#' @noRd
#' @importFrom neurosurf nodes geometry BrainSurface
wrap_output.mvpa_surface_dataset <- function(obj, vals, indices) {
  
  dvals <- numeric(length(nodes(geometry(obj$train_data))))
  
  #if (length(indices) != length(vals)) {
  #  browser()
  #}
  
  dvals[indices] <- vals[indices]
  ## bit of a hack
  BrainSurface(geometry=geometry(obj$train_data), indices=indices, data=dvals)
}

#' @export
has_test_set.mvpa_dataset <- function(obj) {
  !is.null(obj$test_data) 
}



