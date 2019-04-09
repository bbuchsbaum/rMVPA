

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
    bspace <- neuroim2::NeuroSpace(c(D,nobs), spacing)
    bvec <- neuroim2::NeuroVec(mat, bspace)
    
    mask <- as.logical(neuroim2::NeuroVol(array(rep(0, prod(D)), D), neuroim2::NeuroSpace(D, spacing)))
    roi <- neuroim2::spherical_roi(mask, round((dim(bspace)[1:3])/2), radius=min(dim(bspace)/2))
    mask[coords(roi)] <- 1
    
    if (external_test) {
      mat <- array(rnorm(prod(D)*ntest_obs), c(D,ntest_obs))
      bspace <- neuroim2::NeuroSpace(c(D,ntest_obs), spacing)
      testvec <- neuroim2::NeuroVec(mat, bspace)
      dset <- mvpa_dataset(train_data=bvec, test_data=testvec, mask=mask)
    } else {
      dset <- mvpa_dataset(train_data=bvec,mask=mask)
    }
  } else {
    fname <- system.file("extdata/std.lh.smoothwm.asc", package="neuroim2")
    geom <- neurosurf::read_surf_geometry(fname)
    nvert <- nrow(neurosurf::vertices(geom))
    mat <- matrix(rnorm(nvert*nobs), nvert, nobs)
    bvec <- neurosurf::NeuroSurfaceVector(geom, 1:nvert, mat)
    
    if (external_test) {
      test_data <- neurosurf::NeuroSurfaceVector(geom, 1:nvert, matrix(rnorm(nvert*ntest_obs), nvert, ntest_obs))
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
#' @param train_data the training data set: a \code{NeuroVec} instance
#' @param test_data the optional test data set: a \code{NeuroVec} instance
#' @param mask the set of voxels to include: a \code{NeuroVol} instance
#' @importFrom assertthat assert_that
#' @export
mvpa_dataset <- function(train_data, test_data=NULL, mask) {
  assert_that(inherits(train_data, "NeuroVec"))
  if (!is.null(test_data)) {
    assert_that(inherits(test_data, "NeuroVec"))
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
#' Construct a MVPA dataset with a surface-based dataset.
#' 
#' @param train_data the training data, must inherit from \code{NeuroSurfaceVector} in \code{neurosurf} package.
#' @param test_data optional test data, must inherit from \code{NeuroSurfaceVector} in \code{neurosurf} package.
#' @param mask a binary mask equal to the number of nodes in the training/test data set.
#' @param name label to identify the dataset (e.g. "lh" or "rh" to indicate hemisphere)
#' @importFrom assertthat assert_that
#' @export
mvpa_surface_dataset <- function(train_data, test_data=NULL, mask=NULL, name="") {
  
  assert_that(inherits(train_data, "NeuroSurfaceVector"))
  
  if (!is.null(test_data)) {
    assert_that(inherits(test_data, "NeuroSurfaceVector"))
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
#' @method print mvpa_design
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
#' @method print mvpa_surface_dataset
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
    neuroim2::searchlight(obj$mask, radius=radius,...)
  } else {
    neuroim2::random_searchlight(obj$mask, radius=radius,...)
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
#' @import neuroim2
wrap_output.mvpa_dataset <- function(obj, vals, indices=NULL) {
  if (!is.null(indices)) {
    NeuroVol(vals, space(obj$mask), indices=indices)
  } else {
    NeuroVol(vals, space(obj$mask))
  }
}


#' @keywords internal
#' @noRd
#' @importFrom neurosurf nodes geometry NeuroSurface
wrap_output.mvpa_surface_dataset <- function(obj, vals, indices) {
  
  dvals <- numeric(length(nodes(geometry(obj$train_data))))
  
  #if (length(indices) != length(vals)) {
  #  browser()
  #}
  
  dvals[indices] <- vals[indices]
  ## bit of a hack
  NeuroSurface(geometry=geometry(obj$train_data), indices=indices, data=dvals)
}

#' @export
has_test_set.mvpa_dataset <- function(obj) {
  !is.null(obj$test_data) 
}



