library(neuroim)

gen_dataset <- function(D, nobs, nlevels, spacing=c(1,1,1), folds=5) {
  mat <- array(rnorm(prod(D)*nobs), c(D,nobs))
  bspace <- BrainSpace(c(D,nobs), spacing)
  bvec <- BrainVector(mat, bspace)
  mask <- as.logical(BrainVolume(array(rep(1, prod(D)), D), BrainSpace(D, spacing)))
  
  
  Y <- sample(factor(rep(letters[1:nlevels], length.out=nobs)))
  blockVar <- rep(1:folds, length.out=nobs)
  MVPADataset(trainVec=bvec, Y=Y, mask=mask, blockVar=blockVar, testVec=NULL, testY=NULL)
}

test_that("mvpa_regional with 5 ROIS runs without error", {
  
  dataset <- gen_dataset(c(10,10,2), 100, 3)
  crossVal <- BlockedCrossValidation(dataset$blockVar)
  
  regionMask <- BrainVolume(sample(1:5, size=length(dataset$mask), replace=TRUE), space(dataset$mask))
  
  res <- mvpa_regional(dataset, regionMask, crossVal)
  
})

test_that("mvpa_regional with 5 ROIS and consensus learning runs without error", {
  
  dataset <- gen_dataset(c(10,10,2), 100, 3)
  crossVal <- BlockedCrossValidation(dataset$blockVar)
  regionMask <- BrainVolume(sample(1:5, size=length(dataset$mask), replace=TRUE), space(dataset$mask))
  res <- mvpa_regional(dataset, regionMask, crossVal)
  consResult1 <- consensusWeights(res$resultSet, "glmnet")
  consResult2 <- consensusWeights(res$resultSet, "greedy")
  consResult3 <- consensusWeights(res$resultSet, "auc_weights")
  consResult4 <- consensusWeights(res$resultSet, "equal_weights")
 
})
test_that("mvpa_regional with 5 ROIs and ANOVA FeatureSelection with topk=10", {
  
  dataset <- gen_dataset(c(10,10,2), 100, 3)
  crossVal <- BlockedCrossValidation(dataset$blockVar)
  regionMask <- BrainVolume(sample(1:5, size=length(dataset$mask), replace=TRUE), space(dataset$mask))
  
  fsel <- FeatureSelector("FTest", "topk", 10)
  res <- mvpa_regional(dataset, regionMask, crossVal, featureSelector=fsel)
  
})

test_that("mvpa_regional with 5 ROIs and ANOVA FeatureSelection with topk=10", {
  
  dataset <- gen_dataset(c(10,10,2), 100, 3)
  crossVal <- BlockedCrossValidation(dataset$blockVar)
  regionMask <- BrainVolume(sample(1:5, size=length(dataset$mask), replace=TRUE), space(dataset$mask))
  
  fsel <- FeatureSelector("FTest", "topk", 10)
  res <- mvpa_regional(dataset, regionMask, crossVal, featureSelector=fsel)
  fsel2 <- FeatureSelector("FTest", "topp", .4)
  res <- mvpa_regional(dataset, regionMask, crossVal, featureSelector=fsel2)
})
