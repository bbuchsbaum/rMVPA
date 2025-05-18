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
#' data_list <- list(
#'   y = dissimilarity_matrix_y,
#'   x1 = dissimilarity_matrix_x1,
#'   x2 = dissimilarity_matrix_x2
#' )
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
#' data_list <- list(
#'   y = dissimilarity_matrix_y,
#'   x1 = dissimilarity_matrix_x1,
#'   x2 = dissimilarity_matrix_x2
#' )
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
#' @rdname train_model
#' @param y the response variable
#' @param indices The indices of the training data.
#' @param ... Additional arguments passed to the training method.
#' @return A named numeric vector of -log(p-values) for each predictor in the MANOVA model.
#' @importFrom stats as.formula
#' @method train_model manova_model
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
  # Print header
  cat("\\n", "MANOVA Model", "\\n")
  cat(rep("-", 20), "\\n\\n")

  # Formula section
  cat("Model Specification:\\n")
  cat("  |- Formula: ", deparse(x$design$formula), "\\n\\n")

  # Dataset information
  cat("Dataset:\\n")
  dims <- dim(x$dataset$train_data)
  dim_str <- paste0(paste(dims[-length(dims)], collapse=" x "),
                   " x ", dims[length(dims)], " observations")
  cat("  |- Dimensions: ", dim_str, "\\n")
  cat("  |- Type: ", class(x$dataset$train_data)[1], "\\n\\n")

  # Variables section
  cat("Variables:\\n")
  predictors <- all.vars(x$design$formula[[3]])  # Get predictor names from RHS of formula
  response <- all.vars(x$design$formula[[2]])    # Get response name from LHS of formula
  cat("  |- Response: ", response, "\\n")
  cat("  |- Predictors: ", paste(predictors, collapse=", "), "\\n\\n")

  # Data structure
  cat("Data Structure:\\n")

  # Check if there's a test set
  has_test <- !is.null(x$dataset$test_data)
  cat("  |- Test Set: ",
      if(has_test) "Present" else "None", "\\n")

  # Check if there's a mask
  if (!is.null(x$dataset$mask)) {
    mask_sum <- sum(x$dataset$mask > 0)
    cat("  |- Active Voxels/Vertices: ",
        format(mask_sum, big.mark=","), "\\n")
  } else {
    cat("  |- Mask: None\\n")
  }

  cat("\\n")
}

#' @export
#' @method print manova_design
print.manova_design <- function(x, ...) {
  # Print header
  cat("\\n", "MANOVA Design", "\\n")
  cat(rep("-", 20), "\\n\\n")

  # Formula section
  cat("Formula:\\n")
  cat("  |- ", deparse(x$formula), "\\n\\n")

  # Data section
  cat("Variables:\\n")
  var_names <- names(x$data)
  cat("  |- Total Variables: ", length(var_names), "\\n")
  cat("  |- Names: ", paste(var_names, collapse=", "), "\\n")

  cat("\\n")
}

#' Merge Results for MANOVA Model
#'
#' This function takes the computed -log(p-values) from `train_model.manova_model` 
#' for a single ROI/searchlight and formats it into the standard output tibble.
#'
#' @param obj The MANOVA model specification.
#' @param result_set A tibble containing the results from the processor function for the current ROI/sphere. 
#'   Expected to have a `$result` column containing the named vector of -log(p-values).
#' @param indices Voxel indices for the current ROI/searchlight sphere.
#' @param id Identifier for the current ROI/searchlight center.
#' @param ... Additional arguments (ignored).
#'
#' @return A tibble row with the formatted performance metrics for the ROI/sphere.
#' @importFrom tibble tibble
#' @importFrom futile.logger flog.error flog.warn
#' @method merge_results manova_model
merge_results.manova_model <- function(obj, result_set, indices, id, ...) {
  
  # Check for errors from previous steps (processor/train_model)
  if (any(result_set$error)) {
    emessage <- result_set$error_message[which(result_set$error)[1]]
    # Return standard error tibble structure
    return(
      tibble::tibble(
        result       = list(NULL), 
        indices      = list(indices),
        performance  = list(NULL), 
        id           = id,
        error        = TRUE,
        error_message= emessage
      )
    )
  }
  
  # Extract the -log(p-values) computed by train_model. 
  if (!"result" %in% names(result_set) || length(result_set$result) == 0 || is.null(result_set$result[[1]])) {
     error_msg <- sprintf("merge_results (manova): result_set missing or has NULL/empty 'result' field where -log(pvals) were expected for ROI/ID %s.", id)
     futile.logger::flog.error(error_msg)
     # Create NA performance matrix based on expected names if possible
     expected_names <- tryCatch(names(train_model.manova_model(obj, matrix(rnorm(2*2),2,2), NULL, 1:2)), 
                                error = function(e) character(0)) # Get expected names via dummy run
     perf_mat <- matrix(NA_real_, nrow=1, ncol=length(expected_names), dimnames=list(NULL, expected_names))
     if (ncol(perf_mat) == 0) perf_mat <- NULL # If we couldn't get names
     
     return(tibble::tibble(result=list(NULL), indices=list(indices), performance=list(perf_mat), 
                           id=id, error=TRUE, error_message=error_msg))
  }
  
  log_pvals <- result_set$result[[1]]
  
  # Validate the extracted results
  if (!is.numeric(log_pvals) || is.null(names(log_pvals))) {
      error_msg <- sprintf("merge_results (manova): Extracted results are not a named numeric vector for ROI/ID %s.", id)
      futile.logger::flog.error(error_msg)
      perf_mat <- matrix(NA_real_, nrow=1, ncol=length(log_pvals), dimnames=list(NULL, names(log_pvals)))
      if (ncol(perf_mat) == 0) perf_mat <- NULL
      return(tibble::tibble(result=list(NULL), indices=list(indices), performance=list(perf_mat), 
                            id=id, error=TRUE, error_message=error_msg))
  }
  
  # Format the results into a 1-row matrix (performance matrix)
  perf_mat <- matrix(log_pvals, nrow = 1, dimnames = list(NULL, names(log_pvals)))
  
  # Return the standard tibble structure
  tibble::tibble(
    result      = list(NULL), # MANOVA model typically doesn't have 'predictions' in the usual sense
    indices     = list(indices),
    performance = list(perf_mat),
    id          = id,
    error       = FALSE,
    error_message = "~" 
  )
}


