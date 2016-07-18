library(neuroim)
library(testthat)

gen_dataset <- function(D, nobs, nlevels, spacing=c(1,1,1), folds=5) {
  mat <- array(rnorm(prod(D)*nobs), c(D,nobs))
  bspace <- BrainSpace(c(D,nobs), spacing)
  bvec <- BrainVector(mat, bspace)
  mask <- as.logical(BrainVolume(array(rep(1, prod(D)), D), BrainSpace(D, spacing)))
  Y <- sample(factor(rep(letters[1:nlevels], length.out=nobs)))
  blockVar <- rep(1:folds, length.out=nobs)
  MVPADataset$new(trainVec=bvec, Y=Y, mask=mask, blockVar=blockVar, testVec=NULL, testY=NULL)
}

gen_regression_dataset <- function(D, nobs, spacing=c(1,1,1), folds=5) {
  mat <- array(rnorm(prod(D)*nobs), c(D,nobs))
  bspace <- BrainSpace(c(D,nobs), spacing)
  bvec <- BrainVector(mat, bspace)
  mask <- as.logical(BrainVolume(array(rep(1, prod(D)), D), BrainSpace(D, spacing)))
  Y <- rnorm(nobs)
  blockVar <- rep(1:folds, length.out=nobs)
  MVPADataset$new(trainVec=bvec, Y=Y, mask=mask, blockVar=blockVar, testVec=NULL, testY=NULL)
}



test_that("mvpa_regional with 5 ROIS runs without error", {
  
  dataset <- gen_dataset(c(10,10,2), 100, 3)
  crossVal <- BlockedCrossValidation(dataset$blockVar)
  
  regionMask <- BrainVolume(sample(1:5, size=length(dataset$mask), replace=TRUE), space(dataset$mask))
  model <- loadModel("sda_notune")
  res <- mvpa_regional(dataset, model, regionMask, crossVal)
  
})

test_that("mvpa_regional with 5 ROIS with sda_boot runs without error", {
  
  dataset <- gen_dataset(c(10,10,2), 100, 3)
  crossVal <- BlockedCrossValidation(dataset$blockVar)
  
  regionMask <- BrainVolume(sample(1:5, size=length(dataset$mask), replace=TRUE), space(dataset$mask))
  model <- loadModel("sda_boot")
  res <- mvpa_regional(dataset, model, regionMask, crossVal)
  
})


test_that("mvpa_regional with 5 ROIS with sda_boot and custom_performance runs without error", {
  
  dataset <- gen_dataset(c(10,10,2), 100, 3)
  crossVal <- BlockedCrossValidation(dataset$blockVar)
  
  regionMask <- BrainVolume(sample(1:5, size=length(dataset$mask), replace=TRUE), space(dataset$mask))
  model <- loadModel("sda_boot", list(custom_performance = function(x) {
    print(names(x$prob))
    print(x$testDesign[[]])
  }))
  res <- mvpa_regional(dataset, model, regionMask, crossVal)
  
})

test_that("mvpa_regional with 5 ROIS runs and sparse_sda without error", {
  tuneGrid <- expand.grid(frac=c(.2,.5,.8), lambda=c(.01, .2, .8))
  model <- loadModel("sparse_sda", list(tuneGrid=tuneGrid))
  dataset <- gen_dataset(c(10,10,2), 100, 3)
  crossVal <- BlockedCrossValidation(dataset$blockVar)
 
  regionMask <- BrainVolume(sample(1:5, size=length(dataset$mask), replace=TRUE), space(dataset$mask))
  res <- mvpa_regional(dataset, model, regionMask, crossVal)
  
})

test_that("mvpa_regional with 5 ROIS runs and clusterSVM without error", {
  
  model <- loadModel("clusterSVM", list(tuneGrid=expand.grid(K=c(2,3), lambda=c(.001, .01), cost=1)))
  dataset <- gen_dataset(c(10,10,2), 100, 3)
  crossVal <- BlockedCrossValidation(dataset$blockVar)

  regionMask <- BrainVolume(sample(1:5, size=length(dataset$mask), replace=TRUE), space(dataset$mask))
  
  res <- mvpa_regional(dataset, model, regionMask, crossVal)
  
})

test_that("mvpa_regional with 5 ROIS and consensus learning runs without error", {
  
  dataset <- gen_dataset(c(10,10,2), 100, 3)
  crossVal <- BlockedCrossValidation(dataset$blockVar)
  regionMask <- BrainVolume(sample(1:5, size=length(dataset$mask), replace=TRUE), space(dataset$mask))
  
  model <- loadModel("sda")
  res <- mvpa_regional(dataset, model, regionMask, crossVal)
  
  consResult1 <- consensusWeights(res$resultSet, "glmnet")
  consResult2 <- consensusWeights(res$resultSet, "greedy")
  consResult3 <- consensusWeights(res$resultSet, "auc_weights")
  consResult4 <- consensusWeights(res$resultSet, "equal_weights")
 
})

test_that("mvpa_regional_consensus with 5 ROIS runs without error", {
  
  dataset <- gen_dataset(c(10,10,2), 100, 3)
  regionMask <- BrainVolume(sample(1:5, size=length(dataset$mask), replace=TRUE), space(dataset$mask))
  
  model <- loadModel("sda")
  res <- mvpa_regional_consensus(dataset, model, regionMask)
  
})



test_that("mvpa_regional with 5 ROIs and ANOVA FeatureSelection with topk=10", {
  
  dataset <- gen_dataset(c(10,10,2), 100, 3)
  crossVal <- BlockedCrossValidation(dataset$blockVar)
  regionMask <- BrainVolume(sample(1:5, size=length(dataset$mask), replace=TRUE), space(dataset$mask))
  model <- loadModel("sda")
  fsel <- FeatureSelector("FTest", "topk", 10)
  res <- mvpa_regional(dataset, model, regionMask, crossVal, featureSelector=fsel)
  
})



test_that("mvpa_regional with 5 ROIs and ANOVA FeatureSelection with topp=.4", {
  
  dataset <- gen_dataset(c(10,10,2), 100, 3)
  crossVal <- BlockedCrossValidation(dataset$blockVar)
  regionMask <- BrainVolume(sample(1:5, size=length(dataset$mask), replace=TRUE), space(dataset$mask))
  model <- loadModel("sda")
  
  fsel <- FeatureSelector("FTest", "topp", .4)
  res <- mvpa_regional(dataset, model, regionMask, crossVal, featureSelector=fsel)
})


test_that("mvpa_regional with regression and 5 ROIs runs without error", {
  
  dataset <- gen_regression_dataset(c(10,10,2), 100)
  crossVal <- BlockedCrossValidation(dataset$blockVar)
  
  regionMask <- BrainVolume(sample(1:5, size=length(dataset$mask), replace=TRUE), space(dataset$mask))
  model <- loadModel("glmnet", list(tuneGrid=expand.grid(alpha=c(.1,.5), lambda=c(.001,.2,.5))))
  res <- mvpa_regional(dataset, model, regionMask, crossVal)
  
})



