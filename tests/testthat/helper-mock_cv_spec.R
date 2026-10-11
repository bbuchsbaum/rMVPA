library(rMVPA) # Ensure rMVPA generics are in scope

# Mock cv_spec object with S3 methods
mock_cv_spec_s3 <- function(mvpa_design) {
  n_folds_val <- length(unique(mvpa_design$block_var))
  folds_val <- as.integer(mvpa_design$block_var)
  
  obj <- list(
    # Store necessary data for methods
    .folds_val = folds_val,
    .n_folds_val = n_folds_val
  )
  class(obj) <- c("mock_cv_spec", "cross_validation", "list")
  obj
}

# Method for mock_cv_spec - get_nfolds
get_nfolds.mock_cv_spec <- function(obj, ...) {
  obj$.n_folds_val
}

# Method for mock_cv_spec - train_indices
train_indices.mock_cv_spec <- function(obj, fold_num, ...) {
  which(obj$.folds_val != fold_num)
} 