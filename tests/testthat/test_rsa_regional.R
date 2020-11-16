test_that("mvpa_regional with 5 ROIS runs without error", {
  dset <- gen_sample_dataset(c(5,5,5), 100, blocks=3)
  
  Dmat <- dist(matrix(rnorm(100*100), 100, 100))
  rdes <- rsa_design(~ Dmat, list(Dmat=Dmat, block=dset$design$block_var), block_var="block")
  mspec <- rsa_model(dataset$dataset, design=rdes)
  region_mask <- NeuroVol(sample(1:5, size=length(dset$dataset$mask), replace=TRUE), space(dset$dataset$mask))
   
  res <- run_regional(mspec, region_mask)
  expect_true(!is.null(res))
  
})

test_that("mvpa_regional with 5 ROIS and multiple distance matrices runs without error", {
  
  dset <- gen_sample_dataset(c(5,5,5), 100, blocks=3)
  
  Dmat1 <- dist(matrix(rnorm(100*100), 100, 100))
  Dmat2 <- dist(matrix(rnorm(100*100), 100, 100))
  rdes <- rsa_design(~ Dmat1 + Dmat2, list(Dmat1=Dmat1, Dmat2=Dmat2, block=dset$design$block_var), block_var="block")
  
  mspec <- rsa_model(dataset$dataset, design=rdes)
  region_mask <- NeuroVol(sample(1:5, size=length(dset$dataset$mask), replace=TRUE), space(dset$dataset$mask))
  
  res <- run_regional(mspec, region_mask)
  expect_true(!is.null(res))

  
})


