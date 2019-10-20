context("manova searchlight")

test_that("standard one-way manova_searchlight runs without error", {
  dataset <- gen_sample_dataset(c(5,5,5), 100, blocks=3)
  
  rdes <- manova_design(~ block_var + Y, dataset$design$train_design)
  mspec <- manova_model(dataset$dataset, design=rdes)
  ret1 <- run_searchlight(mspec, radius=4)
 
  expect_true(!is.null(ret1))
  
})

test_that("randomized one-way manova_searchlight runs without error", {
  dataset <- gen_sample_dataset(c(5,5,5), 100, blocks=3)
  
  rdes <- manova_design(~ block_var + Y, dataset$design$train_design)
  mspec <- manova_model(dataset$dataset, design=rdes)
  ret1 <- run_searchlight(mspec, radius=4,method="randomized")
  
  expect_true(!is.null(ret1))
  
})


