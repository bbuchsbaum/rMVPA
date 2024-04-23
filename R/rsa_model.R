
#' @noRd
#' @keywords internal
sanitize <- function(name) {
  name <- gsub(":", ".", name)
  name <- gsub(" ", "", name)
  name <- gsub("[\\(\\)]", ".", name, perl=TRUE)
  name <- gsub(",", "_", name)
  name <- gsub("\\.$", "", name)
  name
}

#' Construct a design for an RSA (Representational Similarity Analysis) model
#'
#' This function constructs a design for an RSA model using the provided formula, data, and optional parameters.
#'
#' @param formula A formula expression specifying the dissimilarity-based regression function.
#' @param data A named list containing the dissimilarity matrices and any other auxiliary variables.
#' @param block_var An optional \code{formula}, \code{character} name or \code{integer} vector designating the block structure.
#' @param split_by An optional \code{formula} indicating grouping structure for evaluating test performance.
#' @param keep_intra_run A \code{logical} indicating whether to exclude within-run comparisons.
#' @return A list containing the elements of the RSA design, with class attributes "rsa_design" and "list".
#' @details
#' The function creates an RSA design based on the input parameters. It checks the validity of the input data and
#' handles splitting conditions for evaluation of test performance. It also processes optional block structures and
#' within-run comparisons.
#' @importFrom assertthat assert_that
#' @export
#' @examples
#' dismat <- dist(matrix(rnorm(100*100), 100, 100))
#' rdes <- rsa_design(~ dismat, list(dismat=dismat))
rsa_design <- function(formula, data, block_var=NULL, split_by=NULL, keep_intra_run=FALSE) {
  assert_that(purrr::is_formula(formula))
  
  # Check that all variables are either matrices, "dist", or vectors
  nr <- sapply(data, function(x) {
    if (is.matrix(x)) {
      nrow(x)
    } else if (inherits(x, "dist")) {
      attr(x, "Size")
    } else if (is.vector(x)) {
      length(x)
    } else {
      stop(paste("illegal variable type", class(x)))
    }
  })
  
  assert_that(all(nr == nr[1]), msg="all elements in 'data' must have the same number of rows")
  
  check_split <- function(split_var) {
    minSplits <- min(table(split_var))
    if (minSplits < 3) {
      stop(paste("error: splitting condition results in fewer than 3 observations in at least one set"))
    }
  }
  
  # Create split groups if split_by is provided
  split_groups <- if (!is.null(split_by)) {
    split_var <- parse_variable(split_by, data)
    split(seq_along(split_var), split_var)
  }
  
  # Process block_var if provided
  block_var <- if (!is.null(block_var)) {
    parse_variable(block_var, data)
  }
  
  # Include/exclude within-run comparisons based on keep_intra_run
  include <- if (!is.null(block_var) && !keep_intra_run) {
    as.vector(dist(block_var)) != 0
  }
  
  # Create the RSA design as a list
  des <- list(
    formula=formula,
    data=data,
    split_by=split_by,
    split_groups=split_groups,
    block_var=block_var,
    include=include
  )
  
  # Add model matrix to the design list
  mmat <- rsa_model_mat(des)
  des$model_mat <- mmat
  
  # Set the class attributes
  class(des) <- c("rsa_design", "list")
  
  # Return the RSA design
  des
}


#' Construct a model matrix for an RSA (Representational Similarity Analysis) design
#'
#' This function constructs a model matrix for the given RSA design by processing distance matrices and other variables.
#'
#' @param rsa_des An RSA design object created by \code{rsa_design()}.
#' @return A named list of vectors, with each vector corresponding to the processed input data of the RSA design.
#' @details
#' The function takes an RSA design object as input and processes the distance matrices and other variables to
#' construct a model matrix. It handles different types of input matrices, including symmetric and asymmetric
#' distance matrices, and can include or exclude within-run comparisons based on the RSA design.
#' @examples
#' dismat <- dist(matrix(rnorm(100*100), 100, 100))
#' rdes <- rsa_design(~ dismat, list(dismat=dismat))
#' rsa_model_mat(rdes)
rsa_model_mat <- function(rsa_des) {
  rvars <- labels(terms(rsa_des$formula))
  denv <- list2env(rsa_des$data)
  vset <- lapply(rvars, function(x) eval(parse(text=x), denv))
  
  # Process input variables to create vectors from distance matrices
  vmatlist <- lapply(vset, function(v) {
    if (inherits(v, "dist")) {
      # An distance matrix of class "dist"
      as.vector(v)
    } else if (isSymmetric(v)) {
      # A full distance matrix
      v[lower.tri(v)]
    } else {
      as.vector(dist(v))
    }
  })
  
  # Include or exclude within-run comparisons based on rsa_des$include
  if (!is.null(rsa_des$include)) {
    vmatlist <- lapply(vmatlist, function(v) v[rsa_des$include])
  }
  
  # Assign sanitized names to the output list
  names(vmatlist) <- sanitize(rvars)
  
  # Return the model matrix as a named list of vectors
  vmatlist
}


#' Construct an RSA (Representational Similarity Analysis) model
#'
#' This function creates an RSA model object by taking an MVPA (Multi-Variate Pattern Analysis) dataset and an RSA design.
#'
#' @param dataset An instance of an \code{mvpa_dataset}.
#' @param design An instance of an \code{rsa_design} created by \code{rsa_design()}.
#' @return A list with two elements: \code{dataset} and \code{design}, with the class attribute set to \code{"rsa_model"} and \code{"list"}.
#' @examples
#' # Create a random MVPA dataset
#' data <- matrix(rnorm(100 * 100), 100, 100)
#' labels <- factor(rep(1:2, each = 50))
#' mvpa_data <- mvpa_dataset(data, labels)
#'
#' # Create an RSA design
#' dismat <- dist(data)
#' rdes <- rsa_design(~ dismat, list(dismat = dismat))
#'
#' # Create an RSA model
#' rsa_mod <- rsa_model(mvpa_data, rdes)
#' @export
rsa_model <- function(dataset, design) {
  assert_that(inherits(dataset, "mvpa_dataset"))
  assert_that(inherits(design, "rsa_design"))
  
  # Create an RSA model object with the dataset and design
  structure(list(dataset = dataset,
                 design = design),
            class = c("rsa_model", "list"))
}




