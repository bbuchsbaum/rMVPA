


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
  structure(list(dataset=dataset,
                 design=design),
            class=c("manova_model", "list"))
  
  
}
