#' Create Cross-Validation Folds
#'
#' Generates a list of row indices for k-fold cross-validation.
#' Can perform stratified sampling if y is a factor.
#'
#' @param y A vector, typically the response variable.
#' @param k Integer, the number of folds.
#' @param list Logical, if TRUE, return a list of indices for each fold.
#'        If FALSE, return a vector of fold assignments for each observation.
#' @param seed Optional integer for reproducible fold creation.
#' @return If `list=TRUE`, a list of k integer vectors. If `list=FALSE`, an integer
#'         vector of fold assignments.
#' @importFrom rsample vfold_cv assessment
#' @importFrom tidyselect all_of
#' @keywords internal
create_mvpa_folds <- function(y, k = 5, list = TRUE, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  n <- length(y)
  
  # Create a dummy data frame for rsample
  df_for_rsample <- data.frame(.indices = seq_len(n))
  strata_arg <- NULL
  if (is.factor(y) && k < n) { # Stratification possible and meaningful
     df_for_rsample$.response_var_for_stratification <- y
     strata_arg <- ".response_var_for_stratification"
  }
  
  folds_obj <- if (!is.null(strata_arg)) {
    rsample::vfold_cv(df_for_rsample, v = k, strata = all_of(strata_arg), repeats = 1)
  } else {
    rsample::vfold_cv(df_for_rsample, v = k, repeats = 1)
  }
  
  if (list) {
    # Extract assessment (hold-out) indices for each fold
    out_indices <- lapply(folds_obj$splits, function(split) rsample::assessment(split)$.indices)
  } else {
    # Create a vector of fold assignments
    out_indices <- integer(n)
    for (i in seq_along(folds_obj$splits)) {
      out_indices[rsample::assessment(folds_obj$splits[[i]])$.indices] <- i
    }
  }
  return(out_indices)
}