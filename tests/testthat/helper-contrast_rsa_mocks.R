# Mock get_nfolds / train_indices generics for mock_cv_spec objects, used
# with with_mocked_bindings() in the contrast RSA tests.

# Mock functions for get_nfolds
.mock_get_nfolds_contrast <- function(obj, ...) {
  if (inherits(obj, "mock_cv_spec")) {
    return(obj$.n_folds_val)
  }
  stop(".mock_get_nfolds_contrast called with unexpected object type in this test context.")
}

.mock_train_indices_contrast <- function(obj, fold_num, ...) {
  if (inherits(obj, "mock_cv_spec")) {
    return(which(obj$.folds_val != fold_num))
  }
  stop(".mock_train_indices_contrast called with unexpected object type in this test context.")
}
