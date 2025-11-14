#' Generate Sample Dataset for MVPA Analysis
#' 
#' Creates a synthetic dataset for testing and demonstration of MVPA analyses.
#'
#' @param D The data dimension(s): vector of length 2 or 3 for image data, or single number for surface data
#' @param nobs The number of observations
#' @param response_type Either 'categorical' or 'continuous'
#' @param data_mode Either 'image' or 'surface'
#' @param spacing The voxel spacing (default: c(1,1,1))
#' @param blocks The number of 'blocks' in the data (for cross-validation)
#' @param nlevels The number of category levels (only used if response_type='categorical')
#' @param external_test Whether to generate an external test set
#' @param split_by Optional factor for splitting analyses
#' @param na_cols The number of columns to randomly set to NA (default: 0)
#' @param ntest_obs The number of test observations (default: nobs)
#'
#' @return A list containing:
#'   \describe{
#'     \item{dataset}{An \code{mvpa_dataset} object containing:
#'       \itemize{
#'         \item \code{train_data}: Training data as \code{NeuroVec} or \code{ROISurface}
#'         \item \code{test_data}: Test data (if external_test=TRUE)
#'         \item \code{mask}: Binary mask indicating valid voxels/vertices
#'       }
#'     }
#'     \item{design}{An \code{mvpa_design} object containing:
#'       \itemize{
#'         \item \code{y_train}: Response variable for training
#'         \item \code{y_test}: Response variable for test set (if external_test=TRUE)
#'         \item \code{block_var}: Block variable for cross-validation
#'         \item \code{split_by}: Optional splitting factor
#'       }
#'     }
#'   }
#'
#' @examples
#' # Generate categorical image dataset
#' dataset <- gen_sample_dataset(
#'   D = c(10,10,10),
#'   nobs = 100,
#'   response_type = "categorical",
#'   data_mode = "image",
#'   blocks = 3,
#'   nlevels = 2
#' )
#'
#' # Generate continuous surface dataset
#' surf_data <- gen_sample_dataset(
#'   D = 1000,  # number of vertices
#'   nobs = 50,
#'   response_type = "continuous",
#'   data_mode = "surface"
#' )
#'
#' # Generate dataset with external test set
#' test_dataset <- gen_sample_dataset(
#'   D = c(8,8,8),
#'   nobs = 80,
#'   response_type = "categorical",
#'   nlevels = 3,
#'   external_test = TRUE
#' )
#'
#' @export
gen_sample_dataset <- function(D, nobs, response_type=c("categorical", "continuous"), 
                             data_mode=c("image", "surface"), spacing=c(1,1,1), 
                             blocks=5, nlevels=5, external_test=FALSE, 
                             ntest_obs=nobs, split_by=NULL, na_cols=0) {
 
  response_type <- match.arg(response_type)
  data_mode <- match.arg(data_mode)
  
  # Ensure na_cols is numeric and has a default
  na_cols <- as.numeric(na_cols)
  if (is.null(na_cols) || is.na(na_cols)) {
    na_cols <- 0
  }
  
  if (data_mode == "image") {
    mat <- array(rnorm(prod(D)*nobs), c(D,nobs))
    if (na_cols > 0) {
      naidx <- sample(dim(mat)[4], na_cols)
      for (naid in naidx) {
        ind <- arrayInd(naid, dim(mat)[1:3])
        mat[ind[1], ind[2], ind[3],] <- NA
      }
    } 
    bspace <- neuroim2::NeuroSpace(c(D,nobs), spacing)
    bvec <- neuroim2::NeuroVec(mat, bspace)
    
    mask <- as.logical(neuroim2::NeuroVol(array(rep(0, prod(D)), D), neuroim2::NeuroSpace(D, spacing)))
    roi <- neuroim2::spherical_roi(mask, ceiling((dim(bspace)[1:3])/2), radius=ceiling(min(dim(bspace)/2)))
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
    fname <- system.file("extdata/std.8_lh.inflated.asc", package="neurosurf")
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
    if (nobs == 1) {
      factor(letters[1], levels = letters[1:max(2, nlevels)])
    } else {
      sample(factor(rep(letters[1:nlevels], length.out=nobs)))
    }
  } else {
    rnorm(nobs)
  }
  
  Ytest <- if (response_type == "categorical") {
    sample(factor(rep(letters[1:nlevels], length.out=ntest_obs)))
  } else {
    rnorm(ntest_obs)
  }
  
  block_var <- if (nobs == 1 && blocks == 1) {
    1L #  block_var should be an integer
  } else {
    as.integer(as.character(cut(1:nobs, blocks, labels=1:blocks)))
  }
  
  if (external_test) {
    message("external test")
    mvdes <- mvpa_design(data.frame(Y=Y, block_var=block_var), test_design=data.frame(Ytest = Ytest), 
                       block_var= "block_var", y_train= ~ Y, y_test = ~ Ytest, split_by=split_by)
  } else {
    mvdes <- mvpa_design(data.frame(Y=Y, block_var=block_var), block_var="block_var", y_train= ~ Y, split_by=split_by)
  }
  
  # Make sure the dataset also has the _has_test_set flag set consistently
  if (is.list(dset) && "dataset" %in% names(dset)) {
    dset$dataset$has_test_set <- external_test
  } else {
    dset$has_test_set <- external_test
  }
  
  list(dataset=dset, design=mvdes)
}


#' Create an MVPA Dataset Object
#'
#' Creates a dataset object for MVPA analysis that encapsulates a training dataset, 
#' an optional test dataset, and a voxel mask.
#'
#' @param train_data The training data set: a \code{NeuroVec} instance
#' @param test_data Optional test data set: a \code{NeuroVec} instance (default: NULL)
#' @param mask The set of voxels to include: a \code{NeuroVol} instance
#'
#' @return An \code{mvpa_dataset} object (S3 class) containing:
#'   \describe{
#'     \item{train_data}{The training data as a \code{NeuroVec} instance}
#'     \item{test_data}{The test data as a \code{NeuroVec} instance (if provided, otherwise NULL)}
#'     \item{mask}{The binary mask defining valid voxels as a \code{NeuroVol} instance}
#'     \item{has_test_set}{Logical flag indicating whether this dataset has a test set}
#'   }
#'
#' @examples
#' # Use gen_sample_dataset helper to create a simple dataset
#' sample_data <- gen_sample_dataset(c(5, 5, 5), nobs = 100, blocks = 4)
#' dataset <- sample_data$dataset
#'
#' # Access components
#' print(dim(dataset$train_data))
#' print(sum(dataset$mask > 0))
#'
#' @seealso 
#' \code{\link{mvpa_surface_dataset}} for creating surface-based MVPA datasets
#' 
#' \code{\link{mvpa_design}} for creating the corresponding design object
#'
#' @importFrom assertthat assert_that
#' @export
mvpa_dataset <- function(train_data, test_data=NULL, mask) {
  assert_that(inherits(train_data, "NeuroVec"))
  if (!is.null(test_data)) {
    assert_that(inherits(test_data, "NeuroVec"))
  }
  assert_that(inherits(mask, "NeuroVol"))
  
  # Check for single-voxel datasets (1,1,1,time)
  mask_dims <- dim(mask)[1:3]
  total_voxels <- prod(mask_dims)
  if (total_voxels <= 1) {
    stop("Invalid dataset: Only 1 voxel detected (dimensions ",
         paste(mask_dims, collapse="x"),
         "). Feature RSA analysis requires multiple voxels.")
  }
  
  # Check for active voxels in mask
  active_voxels <- sum(mask > 0)
  if (active_voxels <= 1) {
    stop("Invalid dataset: Only ", active_voxels, " active voxel(s) in mask. Feature RSA analysis requires multiple active voxels.")
  }
  
  # Store a flag indicating whether this dataset has a test set
  has_test <- !is.null(test_data)
  
  ret <- structure(
    list(
      train_data=train_data,
      test_data=test_data,
      mask=mask,
      has_test_set=has_test  # Add flag for test set presence
    ),
    class=c("mvpa_image_dataset", "mvpa_dataset", "list")
  )
  ret
}


#' Create a Surface-Based MVPA Dataset Object
#'
#' Creates a dataset object for surface-based MVPA analysis that encapsulates a training dataset,
#' an optional test dataset, and a vertex mask.
#'
#' @param train_data The training data set: must inherit from \code{NeuroSurfaceVector}
#' @param test_data Optional test data set: must inherit from \code{NeuroSurfaceVector} (default: NULL)
#' @param mask Optional binary mask for vertices. If NULL, creates mask from training data indices
#' @param name Optional label to identify the dataset (e.g., "lh" or "rh" to indicate hemisphere)
#'
#' @return An \code{mvpa_surface_dataset} object (S3 class) containing:
#'   \describe{
#'     \item{train_data}{The training data as a \code{NeuroSurfaceVector} instance}
#'     \item{test_data}{The test data as a \code{NeuroSurfaceVector} instance (if provided)}
#'     \item{mask}{A numeric vector indicating valid vertices (1) and excluded vertices (0)}
#'     \item{name}{Character string identifier for the dataset}
#'     \item{has_test_set}{Logical flag indicating whether this dataset has a test set}
#'   }
#'
#' @details
#' If no mask is provided, one will be created automatically using the indices from the training data.
#' The mask will be a numeric vector with length equal to the number of nodes in the surface geometry.
#'
#' @examples
#' \dontrun{
#' # Create surface dataset with automatic mask
#' train_surf <- NeuroSurfaceVector(geometry, data)
#' dataset <- mvpa_surface_dataset(train_surf, name="lh")
#'
#' # Create dataset with test data and custom mask
#' test_surf <- NeuroSurfaceVector(geometry, test_data)
#' mask <- numeric(length(nodes(geometry)))
#' mask[roi_indices] <- 1
#' dataset <- mvpa_surface_dataset(train_surf, test_surf, mask, name="rh")
#' }
#'
#' @seealso 
#' \code{\link{mvpa_dataset}} for creating volume-based MVPA datasets
#' 
#' \code{\link{mvpa_design}} for creating the corresponding design object
#'
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
  
  # Store a flag indicating whether this dataset has a test set
  has_test <- !is.null(test_data)
  
  structure(
    list(
      train_data=train_data,
      test_data=test_data,
      mask=mask,
      name=name,
      has_test_set=has_test  # Add flag for test set presence
    ),
    class=c("mvpa_surface_dataset", "mvpa_dataset", "list")
  )
}

#' @export
#' @method print mvpa_dataset
print.mvpa_dataset <- function(x, ...) {
  # Ensure crayon is available
  if (!requireNamespace("crayon", quietly = TRUE)) {
    stop("Package 'crayon' is required for pretty printing. Please install it.")
  }
  
  # Define color scheme
  header_style <- crayon::bold$cyan
  section_style <- crayon::yellow
  info_style <- crayon::white
  number_style <- crayon::green
  dim_style <- crayon::italic$blue
  
  # Print header
  cat("\n", header_style("MVPA Dataset"), "\n\n")
  
  # Training data section
  cat(section_style("- Training Data"), "\n")
  dims <- dim(x$train_data)
  dim_str <- paste0(paste(dims[-length(dims)], collapse=" x "), 
                   " x ", dim_style(dims[length(dims)]), " observations")
  cat(info_style("  - Dimensions: "), number_style(dim_str), "\n")
  cat(info_style("  - Type: "), class(x$train_data)[1], "\n")
  
  # Test data section
  cat(section_style("- Test Data"), "\n")
  if (is.null(x$test_data)) {
    cat(info_style("  - "), crayon::red("None"), "\n")
  } else {
    dims <- dim(x$test_data)
    dim_str <- paste0(paste(dims[-length(dims)], collapse=" x "), 
                     " x ", dim_style(dims[length(dims)]), " observations")
    cat(info_style("  - Dimensions: "), number_style(dim_str), "\n")
    cat(info_style("  - Type: "), class(x$test_data)[1], "\n")
  }
  
  # Mask information
  cat(section_style("- Mask Information"), "\n")
  mids <- table(x$mask[x$mask != 0])
  if (length(mids) > 0) {
    midstr <- paste(names(mids), ":", number_style(mids), collapse = ", ")
    cat(info_style("  - Areas: "), midstr, "\n")
  }
  cat(info_style("  - Active voxels/vertices: "), 
      number_style(format(sum(x$mask > 0), big.mark=",")), "\n\n")
}

#' @export
#' @method print mvpa_surface_dataset
print.mvpa_surface_dataset <- function(x, ...) {
  # Ensure crayon is available
  if (!requireNamespace("crayon", quietly = TRUE)) {
    stop("Package 'crayon' is required for pretty printing. Please install it.")
  }
  
  # Define color scheme
  header_style <- crayon::bold$cyan
  section_style <- crayon::yellow
  info_style <- crayon::white
  number_style <- crayon::green
  dim_style <- crayon::italic$blue
  name_style <- crayon::magenta
  
  # Print header
  cat("\n", header_style("Surface MVPA Dataset"), "\n\n")
  
  # Dataset name
  if (nzchar(x$name)) {
    cat(section_style("- Name: "), name_style(x$name), "\n")
  }
  
  # Training data section
  cat(section_style("- Training Data"), "\n")
  dims <- dim(x$train_data)
  vertices <- length(nodes(geometry(x$train_data)))
  cat(info_style("  - Vertices: "), number_style(format(vertices, big.mark=",")), "\n")
  cat(info_style("  - Observations: "), number_style(dims[length(dims)]), "\n")
  cat(info_style("  - Type: "), class(x$train_data)[1], "\n")
  
  # Test data section
  cat(section_style("- Test Data"), "\n")
  if (is.null(x$test_data)) {
    cat(info_style("  - "), crayon::red("None"), "\n")
  } else {
    dims <- dim(x$test_data)
    cat(info_style("  - Observations: "), number_style(dims[length(dims)]), "\n")
    cat(info_style("  - Type: "), class(x$test_data)[1], "\n")
  }
  
  # Mask information
  cat(section_style("- Mask Information"), "\n")
  mids <- table(x$mask[x$mask != 0])
  if (length(mids) > 0) {
    midstr <- paste(names(mids), ":", number_style(mids), collapse = ", ")
    cat(info_style("  - Areas: "), midstr, "\n")
  }
  cat(info_style("  - Active vertices: "), 
      number_style(format(sum(x$mask > 0), big.mark=",")), "\n\n")
}


#' @export
#' @method get_searchlight mvpa_image_dataset
get_searchlight.mvpa_image_dataset <- function(obj, type=c("standard", "randomized"), radius=8,...) {
  type <- match.arg(type)
  if (type == "standard") {
    neuroim2::searchlight(obj$mask, radius=radius,...)
  } else {
    neuroim2::random_searchlight(obj$mask, radius=radius,...)
  }
}

#' @export
#' @method get_searchlight mvpa_surface_dataset
get_searchlight.mvpa_surface_dataset <- function(obj, type=c("standard", "randomized"), radius=8,...) {
  type <- match.arg(type)
  #browser()
  # Create the iterator once
  slight <- if (type == "standard") {
    neurosurf::SurfaceSearchlight(geometry(obj$train_data), radius, nodeset=which(obj$mask>0), as_deflist=TRUE)
  } else {
    neurosurf::RandomSurfaceSearchlight(geometry(obj$train_data), radius, nodeset=which(obj$mask>0), as_deflist=TRUE)
  }
  
  slight
}



#' @keywords internal
#' @noRd
#' @importFrom neuroim2 NeuroVol
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
  #browser()
  
  dvals <- numeric(length(nodes(geometry(obj$train_data))))
  
  #if (length(indices) != length(vals)) {
  #  browser()
  #}
  
  dvals[indices] <- vals[indices]
  ## bit of a hack
  ## we fill in with non-zero rather than allow indices to be missing
  #NeuroSurface(geometry=geometry(obj$train_data), indices=indices, data=dvals)
  NeuroSurface(geometry=geometry(obj$train_data), indices=seq_len(length(dvals)), data=dvals)
}

#' @export
has_test_set.mvpa_dataset <- function(obj) {
  # Use the stored flag rather than checking for existence of test_data
  isTRUE(obj$has_test_set)
}

#' @export
nobs.mvpa_dataset <- function(x) {
  dims <- dim(x$train_data)
  if (is.null(dims)) {
    length(x$train_data)
  } else {
    dims[length(dims)]
  }
}

