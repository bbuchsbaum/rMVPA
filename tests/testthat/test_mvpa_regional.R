library(neuroim)
library(neurosurf)
library(testthat)



test_that("mvpa_regional with 5 ROIS runs without error", {
  
  dset <- gen_sample_dataset(c(10,10,4), nobs=100, nlevels=3, data_mode="image", response_type="categorical")
  cval <- twofold_blocked_cross_validation(dset$design$block_var)
  
  region_mask <- BrainVolume(sample(1:5, size=length(dset$dataset$mask), replace=TRUE), space(dset$dataset$mask))
  model <- load_model("sda")
  mspec <- mvpa_model(model, dset$dataset, dset$design, model_type="classification", crossval=cval)
  res <- run_regional(mspec, region_mask, return_fits=TRUE)
  
})

test_that("surface_based mvpa_regional with 5 ROIS runs without error", {
  
  dset <- gen_sample_dataset(c(10,10,4), nobs=100, nlevels=3, data_mode="surface", response_type="categorical")
  cval <- blocked_cross_validation(dset$design$block_var)
  
  maskid <- sample(1:5, size=length(dset$dataset$mask), replace=TRUE)
  region_mask <- BrainSurface(dset$dataset$train_data@geometry, indices=nodes(dset$dataset$train_data@geometry), data=maskid)
  
  model <- load_model("sda")
  mspec <- mvpa_model(model, dset$dataset, dset$design, model_type="classification", crossval=cval)
  res <- run_regional(mspec, region_mask, return_fits=TRUE)
  
})

test_that("mvpa_regional with 5 ROIS with sda_boot runs without error", {
  
  dataset <- gen_sample_dataset(c(10,10,2), nobs=100, nlevels=3,response_type="categorical")
  cval <- blocked_cross_validation(dataset$design$block_var)
  
  regionMask <- BrainVolume(sample(1:5, size=length(dataset$dataset$mask), replace=TRUE), neuroim::space(dataset$dataset$mask))
  model <- load_model("sda_boot")
  mspec <- mvpa_model(model, dataset$dataset, dataset$design, model_type="classification", crossval=cval)
  res <- run_regional(mspec, regionMask)
  
})

# 
# test_that("mvpa_regional with 5 ROIS with sda_boot and custom_performance runs without error", {
#   
#   dataset <- gen_dataset(c(10,10,2), 100, 3)
#   crossVal <- blocked_cross_validation(dataset$blockVar)
#   
#   regionMask <- BrainVolume(sample(1:5, size=length(dataset$mask), replace=TRUE), space(dataset$mask))
#   model <- load_model("sda_boot", list(custom_performance = function(x) {
#     print(names(x$prob))
#     #print(x$testDesign[[]])
#   }))
#   res <- mvpa_regional(dataset, model, regionMask, crossVal)
#   
# })

test_that("mvpa_regional with 5 ROIS runs and sparse_sda without error", {
  tuneGrid <- expand.grid(frac=c(.2,.5,.8), lambda=c(.01, .2, .8))
  model <- load_model("sparse_sda")
  
  dataset <- gen_sample_dataset(c(10,10,2), nobs=100, nlevels=4)
  cval <- blocked_cross_validation(dataset$design$block_var)
  
  mspec <- mvpa_model(model, dataset$dataset, dataset$design, model_type="classification", crossval=cval, tune_grid=tuneGrid)
  regionMask <- BrainVolume(sample(1:5, size=length(dataset$dataset$mask), replace=TRUE), space(dataset$dataset$mask))
  res <- run_regional(mspec, regionMask, TRUE)
  
})

test_that("mvpa_regional with 5 ROIS runs and random forest without error", {
  
  model <- load_model("rf")
  dataset <- gen_sample_dataset(c(10,10,2), nobs=100, nlevels=3)
  cval <- blocked_cross_validation(dataset$design$block_var)

  regionMask <- BrainVolume(sample(1:5, size=length(dataset$dataset$mask), replace=TRUE), space(dataset$dataset$mask))
  mspec <- mvpa_model(model, dataset$dataset, dataset$design, model_type="classification", crossval=cval, 
                      tune_grid=data.frame(mtry=c(2,4,6)))
  res <- run_regional(mspec, regionMask)
  
})

# test_that("mvpa_regional with 5 ROIS and consensus learning runs without error", {
#   
#   dataset <- gen_dataset(c(10,10,2), 100, 3)
#   crossVal <- blocked_cross_validation(dataset$blockVar)
#   regionMask <- BrainVolume(sample(1:5, size=length(dataset$mask), replace=TRUE), space(dataset$mask))
#   
#   model <- load_model("sda")
#   res <- mvpa_regional(dataset, model, regionMask, crossVal)
#   
#   consResult1 <- consensusWeights(res$resultSet, "glmnet")
#   consResult2 <- consensusWeights(res$resultSet, "greedy")
#   consResult3 <- consensusWeights(res$resultSet, "auc_weights")
#   consResult4 <- consensusWeights(res$resultSet, "equal_weights")
#  
# })
# 
# test_that("mvpa_regional_consensus with 5 ROIS runs without error", {
#   
#   dataset <- gen_dataset(c(10,10,2), 100, 3)
#   regionMask <- BrainVolume(sample(1:5, size=length(dataset$mask), replace=TRUE), space(dataset$mask))
#   
#   model <- load_model("sda")
#   res <- mvpa_regional_consensus(dataset, model, regionMask)
#   
# })



test_that("mvpa_regional with 5 ROIs and ANOVA FeatureSelection with topk=10", {
  
  dataset <- gen_sample_dataset(c(10,10,5), nobs=100, nlevels=6)
  crossVal <- blocked_cross_validation(dataset$design$block_var)
  regionMask <- BrainVolume(sample(1:10, size=length(dataset$dataset$mask), replace=TRUE), space(dataset$dataset$mask))
  model <- load_model("sda")
  fsel <- feature_selector("FTest", "topk", 10)
  mspec <- mvpa_model(model, dataset$dataset, dataset$design, model_type="classification", crossval=cval, feature_selector=fsel)
  res <- run_regional(mspec, regionMask, return_fits=TRUE)
  
})



test_that("mvpa_regional with 5 ROIs and ANOVA FeatureSelection with topp=.4", {
  dataset <- gen_sample_dataset(c(10,10,10), nobs=100, nlevels=2)
  crossVal <- blocked_cross_validation(dataset$design$block_var)
  regionMask <- BrainVolume(sample(1:5, size=length(dataset$dataset$mask), replace=TRUE), space(dataset$dataset$mask))
  model <- load_model("sda_notune")
  fsel <- feature_selector("FTest", "topp", .4)
  mspec <- mvpa_model(model, dataset$dataset, dataset$design,model_type="classification", crossval=cval, feature_selector=fsel)
  res <- run_regional(mspec, regionMask, return_fits=TRUE)
 
})


test_that("mvpa_regional with regression and 5 ROIs runs without error", {
  
  dataset <- gen_regression_dataset(c(10,10,2), 100)
  cval <- blocked_cross_validation(dataset$design$block_var)
  
  regionMask <- BrainVolume(sample(1:5, size=length(dataset$mask), replace=TRUE), space(dataset$mask))
  tune_grid <- expand.grid(alpha=c(.1,.5), lambda=c(.001,.2,.5))
  model <- load_model("glmnet")
  mspec <- mvpa_model(model, dataset, model_type="regression", crossval=cval, tune_grid=tune_grid)
  res <- run_regional(mspec, regionMask)
  
})



