


#' manova_design
#' 
#' @param formula a formula expression specifying MANOVA regression model
#' @param data a named list containing the dissimilarity matrices and any other auxiliary variables
#' @importFrom assertthat assert_that
manova_design <- function(formula, data) {
  assert_that(purrr::is_formula(formula))
  
  des <- list(
    formula=formula,
    data=data
  )
  class(des) <- c("manova_design", "list")
  des
}


#' manova_model
#' 
#' @param dataset a \code{mvpa_dataset} instance
#' @param design a \code{manova_design} instance
manova_model <- function(dataset,
                      design) {
  
  assert_that(inherits(dataset, "mvpa_dataset"))
  structure(list(dataset=dataset,
                 design=design),
            class=c("manova_model", "list"))
  
  
}
