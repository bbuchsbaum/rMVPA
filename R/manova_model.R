#' Create a MANOVA Design
#'
#' This function creates a MANOVA design object containing a formula expression and a named list of data.
#'
#' @param formula A formula expression specifying the MANOVA regression model.
#' @param data A named list containing the dissimilarity matrices and any other auxiliary variables.
#' @return A MANOVA design object with class attributes "manova_design" and "list".
#' @details
#' The function takes a formula expression and a named list of data as input, and returns a MANOVA design object.
#' The object is a list that contains the formula expression and the named list of data with class attributes "manova_design" and "list".
#' This object can be further used for MANOVA analysis or other related multivariate statistical methods.
#' @importFrom assertthat assert_that
#' @importFrom purrr is_formula
#' @importFrom ffmanova ffmanova
#' @examples
#' # Create a MANOVA design
#' formula <- y ~ x1 + x2
#' data_list <- list(y = dissimilarity_matrix_y, x1 = dissimilarity_matrix_x1, x2 = dissimilarity_matrix_x2)
#' manova_design_obj <- manova_design(formula, data_list)
#' @export
manova_design <- function(formula, data) {
  assert_that(purrr::is_formula(formula))
  
  des <- list(
    formula=formula,
    data=data
  )
  class(des) <- c("manova_design", "list")
  des
}


#' Create a MANOVA Model
#'
#' This function creates a MANOVA model object containing an `mvpa_dataset` instance and a `manova_design` instance.
#'
#' @param dataset An \code{mvpa_dataset} instance.
#' @param design A \code{manova_design} instance.
#' @return A MANOVA model object with class attributes "manova_model" and "list".
#' @details
#' The function takes an `mvpa_dataset` instance and a `manova_design` instance as input, and returns a MANOVA model object.
#' The object is a list that contains the dataset and the design with class attributes "manova_model" and "list".
#' This object can be used for further multivariate statistical analysis using the MANOVA method.
#' @importFrom assertthat assert_that
#' @importFrom purrr is_formula
#' @examples
#' # Create a MANOVA model
#' dataset <- create_mvpa_dataset(data_matrix, labels, subject_ids)
#' formula <- y ~ x1 + x2
#' data_list <- list(y = dissimilarity_matrix_y, x1 = dissimilarity_matrix_x1, x2 = dissimilarity_matrix_x2)
#' design <- manova_design(formula, data_list)
#' manova_model_obj <- manova_model(dataset, design)
#' @export
manova_model <- function(dataset,
                      design) {
  
  assert_that(inherits(dataset, "mvpa_dataset"))
  assert_that(inherits(design, "manova_design"))
  
  create_model_spec("manova_model", dataset, design)
  
  
}


#' Train a MANOVA Model
#'
#' This function trains a multivariate analysis of variance (MANOVA) model using the specified design.
#'
#' @param obj An object of class \code{manova_model}.
#' @param train_dat The training data.
#' @param y the response variable
#' @param indices The indices of the training data.
#' @param ... Additional arguments passed to the training method.
#' @return A named numeric vector of -log(p-values) for each predictor in the MANOVA model.
#' @importFrom stats as.formula
train_model.manova_model <- function(obj, train_dat, y, indices, ...) {
  dframe <- obj$design$data
  dframe$response <- as.matrix(train_dat)
  form <- stats::as.formula(paste("response", paste(as.character(obj$design$formula), collapse='')))
  
  fres=ffmanova(form, data=dframe)
  pvals=fres$pValues
  names(pvals) <- sanitize(names(pvals))   
  lpvals <- -log(pvals)
  lpvals
}


#' @export
#' @method print manova_model
print.manova_model <- function(x, ...) {
  # Ensure crayon is available
  if (!requireNamespace("crayon", quietly = TRUE)) {
    stop("Package 'crayon' is required for pretty printing. Please install it.")
  }
  
  # Define color scheme
  header_style <- crayon::bold$cyan
  section_style <- crayon::yellow
  info_style <- crayon::white
  number_style <- crayon::green
  formula_style <- crayon::italic$blue
  var_style <- crayon::magenta
  
  # Print header
  cat("\n", header_style("█▀▀ MANOVA Model ▀▀█"), "\n\n")
  
  # Formula section
  cat(section_style("├─ Model Specification"), "\n")
  cat(info_style("│  └─ Formula: "), formula_style(deparse(x$design$formula)), "\n")
  
  # Dataset information
  cat(section_style("├─ Dataset"), "\n")
  dims <- dim(x$dataset$train_data)
  dim_str <- paste0(paste(dims[-length(dims)], collapse=" × "), 
                   " × ", number_style(dims[length(dims)]), " observations")
  cat(info_style("│  ├─ Dimensions: "), dim_str, "\n")
  cat(info_style("│  └─ Type: "), class(x$dataset$train_data)[1], "\n")
  
  # Variables section
  cat(section_style("├─ Variables"), "\n")
  predictors <- all.vars(x$design$formula[[3]])  # Get predictor names from RHS of formula
  response <- all.vars(x$design$formula[[2]])    # Get response name from LHS of formula
  cat(info_style("│  ├─ Response: "), var_style(response), "\n")
  cat(info_style("│  └─ Predictors: "), var_style(paste(predictors, collapse=", ")), "\n")
  
  # Data structure
  cat(section_style("└─ Data Structure"), "\n")
  
  # Check if there's a test set
  has_test <- !is.null(x$dataset$test_data)
  cat(info_style("   ├─ Test Set: "), 
      if(has_test) crayon::green("Present") else crayon::red("None"), "\n")
  
  # Check if there's a mask
  if (!is.null(x$dataset$mask)) {
    mask_sum <- sum(x$dataset$mask > 0)
    cat(info_style("   └─ Active Voxels/Vertices: "), 
        number_style(format(mask_sum, big.mark=",")), "\n")
  } else {
    cat(info_style("   └─ Mask: "), crayon::red("None"), "\n")
  }
  
  cat("\n")
}

#' @export
#' @method print manova_design
print.manova_design <- function(x, ...) {
  # Ensure crayon is available
  if (!requireNamespace("crayon", quietly = TRUE)) {
    stop("Package 'crayon' is required for pretty printing. Please install it.")
  }
  
  # Define color scheme
  header_style <- crayon::bold$cyan
  section_style <- crayon::yellow
  info_style <- crayon::white
  formula_style <- crayon::italic$blue
  var_style <- crayon::magenta
  
  # Print header
  cat("\n", header_style("█▀▀ MANOVA Design ▀▀█"), "\n\n")
  
  # Formula section
  cat(section_style("├─ Formula"), "\n")
  cat(info_style("│  └─ "), formula_style(deparse(x$formula)), "\n")
  
  # Data section
  cat(section_style("└─ Variables"), "\n")
  var_names <- names(x$data)
  cat(info_style("   ├─ Total Variables: "), crayon::green(length(var_names)), "\n")
  cat(info_style("   └─ Names: "), var_style(paste(var_names, collapse=", ")), "\n")
  
  cat("\n")
}


