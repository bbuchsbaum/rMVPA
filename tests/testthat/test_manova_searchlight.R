context("manova searchlight")

test_that("standard one-way manova_searchlight runs without error", {
  dataset <- gen_sample_dataset(c(4,4,4), 36, blocks=3)
  
  rdes <- manova_design(~ block_var + Y, dataset$design$train_design)
  mspec <- manova_model(dataset$dataset, design=rdes)
  ret1 <- run_searchlight(mspec, radius=2)
 
  expect_true(!is.null(ret1))
  
})

test_that("randomized one-way manova_searchlight runs without error", {
  dataset <- gen_sample_dataset(c(4,4,4), 36, blocks=3)
  
  rdes <- manova_design(~ block_var + Y, dataset$design$train_design)
  mspec <- manova_model(dataset$dataset, design=rdes)
  ret1 <- run_searchlight(mspec, radius=2, method="randomized", niter=1)
  
  expect_true(!is.null(ret1))
  
})

test_that("randomized one-way manova_searchlight runs without error", {
  dataset <- gen_sample_dataset(c(4,4,4), 36, blocks=3)
  
  rdes <- manova_design(~ block_var + Y, dataset$design$train_design)
  mspec <- manova_model(dataset$dataset, design=rdes)
  ret1 <- run_searchlight(mspec, radius=2, method="randomized", niter=1)
  
  expect_true(!is.null(ret1))
  
})

