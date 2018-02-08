test_that("standard rsa_searchlight runs without error", {
  dataset <- gen_sample_dataset(c(5,5,5), 100, blocks=3)
  
  Dmat <- dist(matrix(rnorm(100*100), 100, 100))
  rdes <- rsa_design(~ Dmat, list(Dmat=Dmat, block=dataset$design$block_var), block_var="block")
  mspec <- rsa_model(dataset$dataset, design=rdes)
  ret <- run_searchlight.rsa_model(mspec, 4, "standard", 4, regtype="lm")
  
})