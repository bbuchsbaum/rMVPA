context("Single observation ROI handling")

library(testthat)

# Test compute_trial_scores when there is only one observation

test_that("compute_trial_scores returns NA for single observation", {
  D <- matrix(0, 1, 1)
  rownames(D) <- colnames(D) <- "cond1"
  labels <- "cond1"
  block_var <- 1
  design <- vector_rsa_design(D=D, labels=labels, block_var=block_var)
  dataset <- gen_sample_dataset(c(2,2,2), nobs=1, blocks=1)$dataset
  mspec <- vector_rsa_model(dataset, design, distfun=cordist())
  X <- matrix(rnorm(5), nrow=1)
  scores <- compute_trial_scores(mspec, X)
  expect_true(all(is.na(scores)))
})

# Test evaluate_model handling of NA scores

test_that("evaluate_model returns NA score when observed is NA", {
  D <- matrix(0, 1, 1)
  rownames(D) <- colnames(D) <- "cond1"
  labels <- "cond1"
  block_var <- 1
  design <- vector_rsa_design(D=D, labels=labels, block_var=block_var)
  dataset <- gen_sample_dataset(c(2,2,2), nobs=1, blocks=1)$dataset
  mspec <- vector_rsa_model(dataset, design, distfun=cordist())

  perf <- evaluate_model.vector_rsa_model(mspec, predicted=NULL,
                                          observed=rep(NA_real_, 1),
                                          roi_data_for_perm=NULL,
                                          nperm=0)
  expect_true(is.na(perf$rsa_score))
})
