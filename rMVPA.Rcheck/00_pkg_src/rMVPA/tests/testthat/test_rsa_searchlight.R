context("rsa searchlight")

test_that("standard rsa_searchlight and blocking variable runs without error", {
  dataset <- gen_sample_dataset(c(5,5,5), 100, blocks=3)
  
  Dmat <- dist(matrix(rnorm(100*100), 100, 100))
  rdes <- rsa_design(~ Dmat, list(Dmat=Dmat, block=dataset$design$block_var), block_var="block")
  mspec <- rsa_model(dataset$dataset, design=rdes, regtype="lm")
  ret1 <- run_searchlight(mspec, radius=4, method="standard")
  mspec <- rsa_model(dataset$dataset, design=rdes, regtype="rfit")
  ret2 <- run_searchlight(mspec, radius=4, method="standard")
  mspec <- rsa_model(dataset$dataset, design=rdes, regtype="pearson")
  ret3 <- run_searchlight(mspec, radius=4, method="standard")
  mspec <- rsa_model(dataset$dataset, design=rdes)
  ret4 <- run_searchlight(mspec, radius=4, method="standard")
  expect_true(!is.null(ret1) && !is.null(ret2) && !is.null(ret3) && !is.null(ret4))
  
})

test_that("standard rsa_searchlight with multiple distance matrices and blocking variable runs without error", {
  dataset <- gen_sample_dataset(c(5,5,5), 100, blocks=3)
  
  Dmat1 <- dist(matrix(rnorm(100*100), 100, 100))
  Dmat2 <- dist(matrix(rnorm(100*100), 100, 100))
  rdes <- rsa_design(~ Dmat1 + Dmat2, list(Dmat1=Dmat1, Dmat2=Dmat2, block=dataset$design$block_var), block_var="block")
  mspec <- rsa_model(dataset$dataset, design=rdes, regtype="lm")
  ret1 <- run_searchlight(mspec, radius=4, method="standard")
  mspec <- rsa_model(dataset$dataset, design=rdes, regtype="rfit")
  ret2 <- run_searchlight(mspec, radius=4, method="standard")
  mspec <- rsa_model(dataset$dataset, design=rdes, regtype="pearson")
  ret3 <- run_searchlight(mspec, radius=4, method="standard")
  mspec <- rsa_model(dataset$dataset, design=rdes)
  ret4 <- run_searchlight(mspec, radius=4, method="standard")
  expect_true(!is.null(ret1) && !is.null(ret2) && !is.null(ret3) && !is.null(ret4))
  
})


test_that("randomized rsa_searchlight and blocking variable runs without error", {
  dataset <- gen_sample_dataset(c(5,5,5), 100, blocks=3)
  
  Dmat <- dist(matrix(rnorm(100*100), 100, 100))
  rdes <- rsa_design(~ Dmat, list(Dmat=Dmat, block=dataset$design$block_var), block_var="block")
  mspec <- rsa_model(dataset$dataset, design=rdes, regtype="lm")
  ret1 <- run_searchlight(mspec, radius=4, "randomized")
  expect_true(!is.null(ret1))
  
})

test_that("standard rsa_searchlight and no blocking variable runs without error", {
  dataset <- gen_sample_dataset(c(5,5,5), 100, blocks=3)
  
  Dmat <- dist(matrix(rnorm(100*100), 100, 100))
  rdes <- rsa_design(~ Dmat, list(Dmat=Dmat))
 
  for (regtype in c("lm", "rfit", "pearson", "spearman")) {
    mspec <- rsa_model(dataset$dataset, design=rdes, regtype=regtype)
    ret <- run_searchlight(mspec, radius=4, method="standard")
    expect_true(!is.null(ret))
  }

  
})

test_that("standard rsa_searchlight and blocking variable runs without error", {
  # Generate a sample MVPA dataset with a design including a blocking variable
  dataset <- gen_sample_dataset(c(5, 5, 5), 100, blocks = 3)
  
  # Create a distance matrix
  Dmat <- dist(matrix(rnorm(100 * 100), 100, 100))
  # Build RSA design including the block variable
  rdes <- rsa_design(~ Dmat, list(Dmat = Dmat, block = dataset$design$block_var), block_var = "block")
  
  # Test different regression types
  mspec <- rsa_model(dataset$dataset, design = rdes, regtype = "lm")
  ret1 <- run_searchlight(mspec, radius = 4, method = "standard")
  
  mspec <- rsa_model(dataset$dataset, design = rdes, regtype = "rfit")
  ret2 <- run_searchlight(mspec, radius = 4, method = "standard")
  
  mspec <- rsa_model(dataset$dataset, design = rdes, regtype = "pearson")
  ret3 <- run_searchlight(mspec, radius = 4, method = "standard")
  
  mspec <- rsa_model(dataset$dataset, design = rdes)
  ret4 <- run_searchlight(mspec, radius = 4, method = "standard")
  
  expect_true(!is.null(ret1) && !is.null(ret2) && !is.null(ret3) && !is.null(ret4))
})

test_that("standard rsa_searchlight with multiple distance matrices and blocking variable runs without error", {
  dataset <- gen_sample_dataset(c(5, 5, 5), 100, blocks = 3)
  
  # Create two distance matrices
  Dmat1 <- dist(matrix(rnorm(100 * 100), 100, 100))
  Dmat2 <- dist(matrix(rnorm(100 * 100), 100, 100))
  
  # Build RSA design with two distance matrices and a blocking variable
  rdes <- rsa_design(~ Dmat1 + Dmat2, list(Dmat1 = Dmat1, Dmat2 = Dmat2, block = dataset$design$block_var), block_var = "block")
  
  mspec <- rsa_model(dataset$dataset, design = rdes, regtype = "lm")
  ret1 <- run_searchlight(mspec, radius = 4, method = "standard")
  
  mspec <- rsa_model(dataset$dataset, design = rdes, regtype = "rfit")
  ret2 <- run_searchlight(mspec, radius = 4, method = "standard")
  
  mspec <- rsa_model(dataset$dataset, design = rdes, regtype = "pearson")
  ret3 <- run_searchlight(mspec, radius = 4, method = "standard")
  
  mspec <- rsa_model(dataset$dataset, design = rdes)
  ret4 <- run_searchlight(mspec, radius = 4, method = "standard")
  
  expect_true(!is.null(ret1) && !is.null(ret2) && !is.null(ret3) && !is.null(ret4))
})

test_that("randomized rsa_searchlight and blocking variable runs without error", {
  dataset <- gen_sample_dataset(c(5, 5, 5), 100, blocks = 3)
  
  Dmat <- dist(matrix(rnorm(100 * 100), 100, 100))
  rdes <- rsa_design(~ Dmat, list(Dmat = Dmat, block = dataset$design$block_var), block_var = "block")
  
  mspec <- rsa_model(dataset$dataset, design = rdes, regtype = "lm")
  ret1 <- run_searchlight(mspec, radius = 4, method = "randomized")
  
  expect_true(!is.null(ret1))
})

test_that("standard rsa_searchlight and no blocking variable runs without error", {
  dataset <- gen_sample_dataset(c(5, 5, 5), 100, blocks = 3)
  
  Dmat <- dist(matrix(rnorm(100 * 100), 100, 100))
  # Build RSA design without a blocking variable
  rdes <- rsa_design(~ Dmat, list(Dmat = Dmat))
  
  for (regtype in c("lm", "rfit", "pearson", "spearman")) {
    mspec <- rsa_model(dataset$dataset, design = rdes, regtype = regtype)
    ret <- run_searchlight(mspec, radius = 4, method = "standard")
    expect_true(!is.null(ret))
  }
})