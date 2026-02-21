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
    
    # Use a dense mask by default so small synthetic datasets yield valid
    # searchlights even at the edges (multivariate methods still enforce
    # their own minimum voxel count downstream).
    mask <- as.logical(neuroim2::NeuroVol(array(1, D), neuroim2::NeuroSpace(D, spacing)))
    
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

#' @keywords internal
#' @noRd
validate_image_mask <- function(mask) {
  assert_that(inherits(mask, "NeuroVol"))

  mask_dims <- dim(mask)[1:3]
  total_voxels <- prod(mask_dims)
  if (total_voxels <= 1) {
    stop(
      "Invalid dataset: Only 1 voxel detected (dimensions ",
      paste(mask_dims, collapse = "x"),
      "). Feature RSA analysis requires multiple voxels."
    )
  }

  active_voxels <- sum(mask > 0)
  if (active_voxels <= 1) {
    stop(
      "Invalid dataset: Only ", active_voxels,
      " active voxel(s) in mask. Feature RSA analysis requires multiple active voxels."
    )
  }
}

#' @keywords internal
#' @noRd
split_multibasis_vector <- function(vec, basis_count, ordering = c("event_major", "basis_major"), arg_name = "train_data") {
  ordering <- match.arg(ordering)

  if (is.null(basis_count) || length(basis_count) != 1 || is.na(basis_count)) {
    stop(arg_name, ": `basis_count` must be supplied when using a single concatenated 4D series.")
  }

  basis_count <- as.integer(basis_count)
  if (basis_count < 1) {
    stop(arg_name, ": `basis_count` must be >= 1.")
  }

  dims <- dim(vec)
  if (is.null(dims) || length(dims) < 4) {
    stop(arg_name, ": expected a 4D NeuroVec.")
  }

  nvol <- dims[length(dims)]
  if (nvol %% basis_count != 0) {
    stop(
      arg_name, ": number of volumes (", nvol, ") is not divisible by basis_count (", basis_count, ")."
    )
  }

  n_events <- nvol %/% basis_count
  idx_sets <- if (ordering == "event_major") {
    lapply(seq_len(basis_count), function(b) seq.int(from = b, to = nvol, by = basis_count))
  } else {
    lapply(seq_len(basis_count), function(b) {
      start <- (b - 1L) * n_events + 1L
      seq.int(from = start, length.out = n_events)
    })
  }

  lapply(idx_sets, function(idx) neuroim2::sub_vector(vec, idx))
}

#' @keywords internal
#' @noRd
as_multibasis_series <- function(x, mask, basis_count = NULL, ordering = c("event_major", "basis_major"), arg_name = "train_data") {
  ordering <- match.arg(ordering)
  mask_for_read <- tryCatch(methods::as(mask, "LogicalNeuroVol"), error = function(...) mask)

  read_one <- function(path) {
    if (!file.exists(path)) {
      stop(arg_name, ": file not found: ", path)
    }
    neuroim2::read_vec(path, mask = mask_for_read)
  }

  as_neurovec <- function(item, idx) {
    if (inherits(item, "NeuroVec")) {
      return(item)
    }
    if (is.character(item) && length(item) == 1L) {
      return(read_one(item))
    }
    stop(arg_name, "[[", idx, "]] must be a NeuroVec or a single file path.")
  }

  series_list <- if (inherits(x, "NeuroVec")) {
    split_multibasis_vector(x, basis_count = basis_count, ordering = ordering, arg_name = arg_name)
  } else if (is.character(x) && !is.list(x)) {
    if (length(x) == 1L) {
      vec <- read_one(x)
      split_multibasis_vector(vec, basis_count = basis_count, ordering = ordering, arg_name = arg_name)
    } else {
      lapply(x, read_one)
    }
  } else if (is.list(x)) {
    if (length(x) == 0L) {
      stop(arg_name, " cannot be an empty list.")
    }
    lapply(seq_along(x), function(i) as_neurovec(x[[i]], i))
  } else {
    stop(arg_name, " must be one of: NeuroVec, list of NeuroVec/file paths, character vector of file paths.")
  }

  if (length(series_list) == 0L) {
    stop(arg_name, ": no basis series were resolved.")
  }

  if (!is.null(basis_count)) {
    basis_count <- as.integer(basis_count)
    if (is.na(basis_count) || basis_count < 1L) {
      stop(arg_name, ": `basis_count` must be a positive integer.")
    }
    if (length(series_list) != basis_count) {
      stop(
        arg_name, ": resolved ", length(series_list),
        " basis series but basis_count=", basis_count, "."
      )
    }
  }

  dims0 <- dim(series_list[[1]])
  if (is.null(dims0) || length(dims0) < 4) {
    stop(arg_name, ": resolved basis series must be 4D NeuroVec objects.")
  }

  for (i in seq_along(series_list)) {
    if (!inherits(series_list[[i]], "NeuroVec")) {
      stop(arg_name, "[[", i, "]] is not a NeuroVec.")
    }
    if (!identical(dim(series_list[[i]]), dims0)) {
      stop(arg_name, ": all basis series must have identical dimensions.")
    }
  }

  series_list
}

#' Create a Multibasis MVPA Image Dataset
#'
#' Creates an image dataset where each event has \code{k} basis beta maps. The
#' dataset stores basis images separately and presents event-level observations
#' with basis-concatenated features in downstream extraction methods.
#'
#' @param train_data Training data in one of the following forms:
#'   \itemize{
#'     \item A list of \code{NeuroVec} objects (one per basis).
#'     \item A character vector of 4D image paths (one file per basis).
#'     \item A single 4D \code{NeuroVec} or image path containing concatenated
#'       basis volumes, where \code{basis_count} specifies splitting.
#'   }
#' @param test_data Optional test data in the same format as \code{train_data}.
#' @param mask A \code{NeuroVol} mask.
#' @param basis_count Number of basis functions when using a single concatenated
#'   4D series.
#' @param ordering Ordering of volumes in concatenated series:
#'   \code{"event_major"} means \code{event1(b1..bk), event2(b1..bk), ...};
#'   \code{"basis_major"} means \code{b1(all events), b2(all events), ...}.
#' @param basis_labels Optional character labels for basis functions.
#'
#' @return An object of class \code{mvpa_multibasis_image_dataset}.
#' @examples
#' \donttest{
#'   ds <- gen_sample_dataset(c(5,5,5), 20)
#'   # Create two basis series from same data
#'   mb <- mvpa_multibasis_dataset(
#'     train_data = list(ds$dataset$train_data, ds$dataset$train_data),
#'     mask = ds$dataset$mask
#'   )
#' }
#' @export
mvpa_multibasis_dataset <- function(train_data,
                                    test_data = NULL,
                                    mask,
                                    basis_count = NULL,
                                    ordering = c("event_major", "basis_major"),
                                    basis_labels = NULL) {
  ordering <- match.arg(ordering)
  validate_image_mask(mask)

  train_series <- as_multibasis_series(
    train_data,
    mask = mask,
    basis_count = basis_count,
    ordering = ordering,
    arg_name = "train_data"
  )

  k <- length(train_series)

  test_series <- if (!is.null(test_data)) {
    out <- as_multibasis_series(
      test_data,
      mask = mask,
      basis_count = k,
      ordering = ordering,
      arg_name = "test_data"
    )
    if (length(out) != k) {
      stop("test_data basis count must match train_data basis count.")
    }
    out
  } else {
    NULL
  }

  if (is.null(basis_labels)) {
    basis_labels <- paste0("b", seq_len(k))
  } else {
    if (length(basis_labels) != k) {
      stop("basis_labels length must match the number of basis series (", k, ").")
    }
    basis_labels <- as.character(basis_labels)
  }

  has_test <- !is.null(test_series)

  structure(
    list(
      train_data = train_series,
      test_data = test_series,
      mask = mask,
      basis_count = k,
      basis_labels = basis_labels,
      has_test_set = has_test
    ),
    class = c("mvpa_multibasis_image_dataset", "mvpa_image_dataset", "mvpa_dataset", "list")
  )
}

#' Alias for \code{mvpa_multibasis_dataset}
#'
#' @inheritParams mvpa_multibasis_dataset
#' @return An object of class \code{mvpa_multibasis_image_dataset}.
#' @examples
#' \dontrun{
#'   # See mvpa_multibasis_dataset for examples
#'   ds <- gen_sample_dataset(c(5,5,5), 20)
#'   mb <- mvpa_multibasis_image_dataset(
#'     list(ds$dataset$train_data, ds$dataset$train_data),
#'     mask = ds$dataset$mask
#'   )
#' }
#' @export
mvpa_multibasis_image_dataset <- mvpa_multibasis_dataset


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
#' @method print mvpa_multibasis_image_dataset
print.mvpa_multibasis_image_dataset <- function(x, ...) {
  if (!requireNamespace("crayon", quietly = TRUE)) {
    stop("Package 'crayon' is required for pretty printing. Please install it.")
  }

  header_style <- crayon::bold$cyan
  section_style <- crayon::yellow
  info_style <- crayon::white
  number_style <- crayon::green

  cat("\n", header_style("Multibasis MVPA Dataset"), "\n\n")
  cat(section_style("- Basis Functions"), "\n")
  cat(info_style("  - Count: "), number_style(x$basis_count), "\n")
  cat(info_style("  - Labels: "), paste(x$basis_labels, collapse = ", "), "\n")

  train_dims <- dim(x$train_data[[1]])
  n_events <- train_dims[length(train_dims)]
  cat(section_style("- Training Data"), "\n")
  cat(info_style("  - Event observations: "), number_style(n_events), "\n")
  cat(info_style("  - Basis files/series: "), number_style(length(x$train_data)), "\n")
  cat(info_style("  - Spatial dims: "), paste(train_dims[1:3], collapse = " x "), "\n")

  cat(section_style("- Test Data"), "\n")
  if (is.null(x$test_data)) {
    cat(info_style("  - "), crayon::red("None"), "\n")
  } else {
    test_dims <- dim(x$test_data[[1]])
    cat(info_style("  - Event observations: "), number_style(test_dims[length(test_dims)]), "\n")
    cat(info_style("  - Basis files/series: "), number_style(length(x$test_data)), "\n")
  }

  cat(section_style("- Mask Information"), "\n")
  cat(
    info_style("  - Active voxels: "),
    number_style(format(sum(x$mask > 0), big.mark = ",")),
    "\n\n"
  )
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


#' @keywords internal
#' @noRd
.searchlight_geometry_cache_enabled <- function() {
  TRUE
}

#' @keywords internal
#' @noRd
.searchlight_geometry_cache_max_entries <- function() {
  8L
}

#' @keywords internal
#' @noRd
.searchlight_geometry_cache_env <- local({
  new.env(hash = TRUE, parent = emptyenv())
})

#' @keywords internal
#' @noRd
.searchlight_geometry_cache_state <- local({
  state <- new.env(parent = emptyenv())
  state$order <- character(0)
  state
})

#' @keywords internal
#' @noRd
.searchlight_geometry_cache_clear <- function() {
  keys <- ls(envir = .searchlight_geometry_cache_env, all.names = TRUE)
  if (length(keys) > 0L) {
    rm(list = keys, envir = .searchlight_geometry_cache_env)
  }
  .searchlight_geometry_cache_state$order <- character(0)
  invisible(NULL)
}

#' @keywords internal
#' @noRd
.searchlight_geometry_cache_size <- function() {
  length(.searchlight_geometry_cache_state$order)
}

#' @keywords internal
#' @noRd
.searchlight_geometry_cache_get <- function(key) {
  if (!exists(key, envir = .searchlight_geometry_cache_env, inherits = FALSE)) {
    return(NULL)
  }
  order <- .searchlight_geometry_cache_state$order
  .searchlight_geometry_cache_state$order <- c(order[order != key], key)
  get(key, envir = .searchlight_geometry_cache_env, inherits = FALSE)
}

#' @keywords internal
#' @noRd
.searchlight_geometry_cache_set <- function(key, value) {
  assign(key, value, envir = .searchlight_geometry_cache_env)
  order <- .searchlight_geometry_cache_state$order
  order <- c(order[order != key], key)
  max_entries <- .searchlight_geometry_cache_max_entries()
  if (length(order) > max_entries) {
    drop_n <- length(order) - max_entries
    drop_keys <- order[seq_len(drop_n)]
    for (drop_key in drop_keys) {
      if (exists(drop_key, envir = .searchlight_geometry_cache_env, inherits = FALSE)) {
        rm(list = drop_key, envir = .searchlight_geometry_cache_env)
      }
    }
    order <- order[-seq_len(drop_n)]
  }
  .searchlight_geometry_cache_state$order <- order
  invisible(value)
}

#' @keywords internal
#' @noRd
.searchlight_geometry_raw_signature <- function(raw) {
  x <- as.numeric(as.integer(raw))
  n <- length(x)
  if (n == 0L) {
    return("r0")
  }
  pos <- seq_len(n)
  s1 <- sum(x * ((pos %% 251L) + 1L))
  s2 <- sum(x * ((pos %% 65521L) + 1L))
  paste0(
    "r", n, "_",
    formatC(round(s1 %% 1e12), format = "f", digits = 0), "_",
    formatC(round(s2 %% 1e12), format = "f", digits = 0)
  )
}

#' @keywords internal
#' @noRd
.searchlight_geometry_numeric_signature <- function(x) {
  vals <- as.numeric(x)
  n <- length(vals)
  if (n == 0L) {
    return("n0")
  }
  pos <- seq_len(n)
  s1 <- sum(vals * ((pos %% 8191L) + 1L))
  s2 <- sum(vals * ((pos %% 131071L) + 3L))
  paste0(
    "n", n, "_",
    formatC(round(s1 %% 1e12), format = "f", digits = 0), "_",
    formatC(round(s2 %% 1e12), format = "f", digits = 0)
  )
}

#' @keywords internal
#' @noRd
.searchlight_geometry_dots_signature <- function(...) {
  dots <- list(...)
  if (length(dots) == 0L) {
    return("none")
  }
  raw <- serialize(dots, NULL, ascii = FALSE, version = 2)
  .searchlight_geometry_raw_signature(raw)
}

#' @keywords internal
#' @noRd
.searchlight_geometry_image_signature <- function(obj) {
  mask <- obj$mask
  sp <- neuroim2::space(mask)
  dims <- dim(sp)
  if (length(dims) > 3L) {
    dims <- dims[seq_len(3L)]
  }
  spacing <- neuroim2::spacing(sp)[seq_len(length(dims))]
  origin <- neuroim2::origin(sp)[seq_len(length(dims))]

  idx <- if (!is.null(obj$mask_indices)) {
    as.integer(obj$mask_indices)
  } else {
    vals <- neuroim2::values(mask)
    if (is.matrix(vals)) {
      vals <- vals[, 1, drop = TRUE]
    }
    which(vals != 0)
  }

  paste(
    "img",
    paste(dims, collapse = ","),
    paste(formatC(spacing, digits = 8, format = "fg"), collapse = ","),
    paste(formatC(origin, digits = 8, format = "fg"), collapse = ","),
    .searchlight_geometry_numeric_signature(idx),
    sep = "::"
  )
}

#' @keywords internal
#' @noRd
.searchlight_geometry_surface_signature <- function(obj) {
  geom <- geometry(obj$train_data)
  active_nodes <- which(obj$mask > 0)
  paste(
    "surf",
    length(neurosurf::nodes(geom)),
    .searchlight_geometry_numeric_signature(active_nodes),
    sep = "::"
  )
}

#' @keywords internal
#' @noRd
.searchlight_geometry_key <- function(dataset_class,
                                      type,
                                      radius,
                                      iter,
                                      nonzero,
                                      k,
                                      signature,
                                      dots_signature) {
  paste(
    dataset_class,
    type,
    paste(radius, collapse = ","),
    if (is.null(iter)) "null" else paste(iter, collapse = ","),
    if (is.null(nonzero)) "null" else as.character(nonzero),
    if (is.null(k)) "null" else paste(k, collapse = ","),
    signature,
    dots_signature,
    sep = "||"
  )
}

#' @keywords internal
#' @noRd
.compute_image_searchlight <- function(obj,
                                       type,
                                       radius,
                                       iter,
                                       nonzero,
                                       ...) {
  if (type == "standard") {
    neuroim2::searchlight(obj$mask, radius = radius, nonzero = nonzero, ...)
  } else if (type == "randomized") {
    neuroim2::random_searchlight(obj$mask, radius = radius, ...)
  } else { # resampled
    if (is.null(iter)) iter <- 1L
    neuroim2::resampled_searchlight(obj$mask, radius = radius, iter = iter, ...)
  }
}

#' @keywords internal
#' @noRd
.compute_surface_searchlight <- function(obj,
                                         type,
                                         radius,
                                         ...) {
  active_nodes <- which(obj$mask > 0)
  if (type == "standard") {
    neurosurf::SurfaceSearchlight(geometry(obj$train_data), radius, nodeset = active_nodes, as_deflist = TRUE)
  } else {
    # No direct analogue of resampled searchlight for surfaces; fall back to randomized sampling.
    neurosurf::RandomSurfaceSearchlight(geometry(obj$train_data), radius, nodeset = active_nodes, as_deflist = TRUE)
  }
}

#' @export
#' @method get_searchlight mvpa_image_dataset
#' @importFrom neuroim2 searchlight random_searchlight resampled_searchlight
get_searchlight.mvpa_image_dataset <- function(obj,
                                               type = c("standard", "randomized", "resampled"),
                                               radius = 8,
                                               iter = NULL,
                                               nonzero = TRUE,
                                               k = NULL,
                                               ...) {
  type <- match.arg(type)
  if (!.searchlight_geometry_cache_enabled() || type != "standard") {
    return(.compute_image_searchlight(
      obj = obj,
      type = type,
      radius = radius,
      iter = iter,
      nonzero = nonzero,
      ...
    ))
  }

  dots_signature <- .searchlight_geometry_dots_signature(...)
  key <- .searchlight_geometry_key(
    dataset_class = "mvpa_image_dataset",
    type = type,
    radius = radius,
    iter = iter,
    nonzero = nonzero,
    k = k,
    signature = .searchlight_geometry_image_signature(obj),
    dots_signature = dots_signature
  )
  cached <- .searchlight_geometry_cache_get(key)
  if (!is.null(cached)) {
    futile.logger::flog.debug("searchlight geometry cache hit [image]: radius=%s", paste(radius, collapse = ","))
    return(cached)
  }

  slight <- .compute_image_searchlight(
    obj = obj,
    type = type,
    radius = radius,
    iter = iter,
    nonzero = nonzero,
    ...
  )
  .searchlight_geometry_cache_set(key, slight)
  slight
}

#' @export
#' @method get_searchlight mvpa_surface_dataset
get_searchlight.mvpa_surface_dataset <- function(obj, type = c("standard", "randomized", "resampled"),
                                                 radius = 8, iter = NULL, k = NULL, ...) {
  type <- match.arg(type)
  if (!.searchlight_geometry_cache_enabled() || type != "standard") {
    return(.compute_surface_searchlight(
      obj = obj,
      type = type,
      radius = radius,
      ...
    ))
  }

  dots_signature <- .searchlight_geometry_dots_signature(...)
  key <- .searchlight_geometry_key(
    dataset_class = "mvpa_surface_dataset",
    type = type,
    radius = radius,
    iter = iter,
    nonzero = NULL,
    k = k,
    signature = .searchlight_geometry_surface_signature(obj),
    dots_signature = dots_signature
  )
  cached <- .searchlight_geometry_cache_get(key)
  if (!is.null(cached)) {
    futile.logger::flog.debug("searchlight geometry cache hit [surface]: radius=%s", paste(radius, collapse = ","))
    return(cached)
  }

  slight <- .compute_surface_searchlight(
    obj = obj,
    type = type,
    radius = radius,
    ...
  )
  .searchlight_geometry_cache_set(key, slight)
  slight
}



#' @keywords internal
#' @noRd
#' @importFrom neuroim2 NeuroVol
wrap_output.mvpa_dataset <- function(obj, vals, indices = NULL, ...) {
  if (!is.null(indices)) {
    NeuroVol(vals, space(obj$mask), indices=indices)
  } else {
    NeuroVol(vals, space(obj$mask))
  }
}


#' @keywords internal
#' @noRd
#' @importFrom neurosurf nodes geometry NeuroSurface
wrap_output.mvpa_surface_dataset <- function(obj, vals, indices, ...) {
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

#' @export
nobs.mvpa_multibasis_image_dataset <- function(x) {
  dims <- dim(x$train_data[[1]])
  dims[length(dims)]
}


# ---- Default implementations of the new searchlight generics ----

#' @export
get_center_ids.mvpa_image_dataset <- function(dataset, ...) {
  which(dataset$mask > 0)
}

#' @export
get_center_ids.mvpa_surface_dataset <- function(dataset, ...) {
  which(dataset$mask > 0)
}

#' @export
build_output_map.mvpa_image_dataset <- function(dataset, metric_vector, ids, ...) {
  build_volume_map(dataset, metric_vector, ids)
}

#' @export
build_output_map.mvpa_multibasis_image_dataset <- function(dataset, metric_vector, ids, aggregation = c("mean", "sum", "maxabs"), ...) {
  aggregation <- match.arg(aggregation)

  if (is.null(ids)) {
    return(build_volume_map(dataset, metric_vector, ids))
  }

  if (length(metric_vector) != length(ids)) {
    stop("Length of metric_vector must match length(ids).")
  }

  if (!anyDuplicated(ids)) {
    return(build_volume_map(dataset, metric_vector, ids))
  }

  split_vals <- split(metric_vector, ids)
  agg_fun <- switch(
    aggregation,
    mean = function(x) mean(x, na.rm = TRUE),
    sum = function(x) sum(x, na.rm = TRUE),
    maxabs = function(x) x[which.max(abs(x))]
  )

  ids_u <- as.integer(names(split_vals))
  vals_u <- vapply(split_vals, agg_fun, numeric(1))
  ord <- order(ids_u)

  build_volume_map(dataset, vals_u[ord], ids_u[ord])
}

#' @export
build_output_map.mvpa_surface_dataset <- function(dataset, metric_vector, ids, ...) {
  build_surface_map(dataset, metric_vector, ids)
}

#' @export
searchlight_scope.default <- function(dataset, ...) {
  "searchlight"
}
