test_that("standard rsa_searchlight and blocking variable runs without error", {
  dataset <- gen_sample_dataset(c(5,5,5), 100, blocks=3)
  
  Dmat <- dist(matrix(rnorm(100*100), 100, 100))
  rdes <- rsa_design(~ Dmat, list(Dmat=Dmat, block=dataset$design$block_var), block_var="block")
  mspec <- rsa_model(dataset$dataset, design=rdes)
  ret1 <- run_searchlight(mspec, 4, "standard", 4, regtype="lm")
  ret2 <- run_searchlight(mspec, 4, "standard", 4, regtype="rfit")
  ret3 <- run_searchlight(mspec, 4, "standard", 4, regtype="pearson")
  ret4 <- run_searchlight(mspec, 4, "standard", 4, regtype="spearman")
  
})

test_that("standard rsa_searchlight and no blocking variable runs without error", {
  dataset <- gen_sample_dataset(c(5,5,5), 100, blocks=3)
  
  Dmat <- dist(matrix(rnorm(100*100), 100, 100))
  rdes <- rsa_design(~ Dmat, list(Dmat=Dmat))
  mspec <- rsa_model(dataset$dataset, design=rdes)
  ret1 <- run_searchlight(mspec, 4, "standard", 4, regtype="lm")
  ret2 <- run_searchlight(mspec, 4, "standard", 4, regtype="rfit")
  ret3 <- run_searchlight(mspec, 4, "standard", 4, regtype="pearson")
  ret4 <- run_searchlight(mspec, 4, "standard", 4, regtype="spearman")
  
})