library(testthat)

context("compute_trial_scores input validation")

# Helper to set up a simple vector RSA model
make_simple_model <- function(n) {
  dset <- gen_sample_dataset(c(2,2,2), n, blocks = 3)
  D <- as.matrix(dist(matrix(rnorm(n * n), n, n)))
  labels <- paste0("L", seq_len(n))
  rownames(D) <- labels
  colnames(D) <- labels
  rdes <- vector_rsa_design(D = D,
                            labels = sample(labels, length(dset$design$block_var), replace = TRUE),
                            block_var = dset$design$block_var)
  vector_rsa_model(dset$dataset, rdes, distfun = cordist())
}

# nrow(X) vs nrow(precomputed$Dexpanded)

test_that("compute_trial_scores errors when rows do not match Dexpanded", {
  model <- make_simple_model(20)
  X_wrong <- matrix(rnorm(25), 5, 5) # 5 rows, but Dexpanded is 20
  expect_error(
    compute_trial_scores(model, X_wrong),
    "nrow\\(X\\).*nrow\\(precomputed\\$Dexpanded\\)"
  )
})

# length(block) vs nrow(X)

test_that("compute_trial_scores errors when block length mismatches", {
  model <- make_simple_model(20)
  X <- matrix(rnorm(100), 20, 5)
  model$design$block <- model$design$block[-1]
  expect_error(
    compute_trial_scores(model, X),
    "length\\(obj\\$design\\$block\\).*nrow\\(X\\)"
  )
})
