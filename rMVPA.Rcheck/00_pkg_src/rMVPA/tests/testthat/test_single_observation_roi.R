context("Single observation ROI handling")

library(testthat)

# Test compute_trial_scores when there is only one observation

test_that("compute_trial_scores returns NA for single observation", {
  # For vector_rsa_design with nobs=1:
  D_mat <- matrix(0, 1, 1, dimnames=list("cond1", "cond1"))
  labels_vec <- "cond1"
  block_var_val <- 1
  
  # Create the specific vector_rsa_design
  design <- vector_rsa_design(D = D_mat, labels = labels_vec, block_var = block_var_val)
  
  # Get dataset from gen_sample_dataset, using continuous to avoid factor issues in mvpa_design,
  # as this part is less critical for the specific test of compute_trial_scores with nrow(X)<2.
  dataset <- gen_sample_dataset(c(2,2,2), nobs=1, blocks=1, response_type = "continuous")$dataset
  
  mspec <- vector_rsa_model(dataset, design, distfun=cordist())
  X <- matrix(rnorm(5), nrow=1) # Single observation data matrix
  scores <- compute_trial_scores(mspec, X)
  expect_true(all(is.na(scores)))
})

# Test evaluate_model handling of NA scores

test_that("evaluate_model returns NA score when observed is NA", {
  # For vector_rsa_design with nobs=1:
  D_mat <- matrix(0, 1, 1, dimnames=list("cond1", "cond1"))
  labels_vec <- "cond1"
  block_var_val <- 1
  
  # Create the specific vector_rsa_design
  design <- vector_rsa_design(D = D_mat, labels = labels_vec, block_var = block_var_val)
  
  # Get dataset from gen_sample_dataset
  dataset <- gen_sample_dataset(c(2,2,2), nobs=1, blocks=1, response_type = "continuous")$dataset
  
  mspec <- vector_rsa_model(dataset, design, distfun=cordist())

  perf <- evaluate_model.vector_rsa_model(mspec, predicted=NULL,
                                          observed=rep(NA_real_, 1),
                                          roi_data_for_perm=NULL,
                                          nperm=0)
  expect_true(is.na(perf$rsa_score))
})
