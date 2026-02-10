# Libraries are loaded in helper-test-logging.R

context("mvpa regional")

test_that("mvpa_regional with 5 ROIS runs without error", {

  dset <- gen_sample_dataset(c(10,10,8), nobs=100, nlevels=3, data_mode="image",
                             response_type="categorical")
  cval <- twofold_blocked_cross_validation(dset$design$block_var)

  region_mask <- NeuroVol(sample(1:5, size=length(dset$dataset$mask), replace=TRUE), space(dset$dataset$mask))
  model <- load_model("sda_notune")
  mspec <- mvpa_model(model, dset$dataset, dset$design, model_type="classification", crossval=cval, return_fits=TRUE)
  res <- run_regional(mspec, region_mask)
  expect_true(!is.null(res))

})

test_that("mvpa_regional works with ClusteredNeuroVol region_mask", {
  dset <- gen_sample_dataset(c(10,10,4), nobs=80, nlevels=2, data_mode="image",
                             response_type="categorical")
  cval <- twofold_blocked_cross_validation(dset$design$block_var)

  # Build a ClusteredNeuroVol with 3 clusters from the dataset mask
  sp <- space(dset$dataset$mask)
  mask_arr <- as.logical(dset$dataset$mask)
  mask_vol <- LogicalNeuroVol(mask_arr, sp)
  n_vox <- sum(mask_arr)
  clusters <- rep(seq_len(3), length.out = n_vox)
  region_mask <- ClusteredNeuroVol(mask_vol, clusters)

  model <- load_model("sda_notune")
  mspec <- mvpa_model(model, dset$dataset, dset$design,
                      model_type = "classification", crossval = cval)
  res <- run_regional(mspec, region_mask)
  expect_true(!is.null(res))
  expect_true(nrow(res$performance_table) >= 1)
})



test_that("mvpa_regional with 5 ROIS runs without error and can access fitted model", {
  
  dset <- gen_sample_dataset(c(10,10,4), nobs=100, nlevels=3, data_mode="image", 
                             response_type="categorical")
  cval <- twofold_blocked_cross_validation(dset$design$block_var)
  
  region_mask <- NeuroVol(sample(1:5, size=length(dset$dataset$mask), replace=TRUE), space(dset$dataset$mask))
  model <- load_model("sda_notune")
  mspec <- mvpa_model(model, dset$dataset, dset$design, model_type="classification", crossval=cval, return_fits=TRUE)
  res <- run_regional(mspec, region_mask, return_fits=TRUE)
  fit1 <- res$fits[[1]]
  expect_true(class(fit1)[1] == "weighted_model")
  
})

test_that("can combine two prediction tables from two regional analyses", {
  
  dset <- gen_sample_dataset(c(10,10,4), nobs=100, nlevels=3, data_mode="image", 
                             response_type="categorical")
  cval <- twofold_blocked_cross_validation(dset$design$block_var)
  
  region_mask <- NeuroVol(sample(1:2, size=length(dset$dataset$mask), replace=TRUE), space(dset$dataset$mask))
  rmask1 <- region_mask
  rmask1[rmask1 == 1] <- 1
  rmask1[rmask1 == 2] <- 0
  
  rmask2 <- region_mask
  rmask2[rmask2 == 1] <- 0
  rmask2[rmask2 == 2] <- 2
  
  
  model <- load_model("sda_notune")
  mspec <- mvpa_model(model, dset$dataset, dset$design, model_type="classification", crossval=cval)
  res1 <- run_regional(mspec, region_mask=rmask1, return_fits=TRUE)
  res2 <- run_regional(mspec, region_mask=rmask2, return_fits=TRUE)
  
  ptab <- combine_prediction_tables(list(res1$prediction_table, res2$prediction_table))
  ptab2 <- combine_prediction_tables(list(res1$prediction_table, res2$prediction_table), collapse_regions=TRUE)
  expect_true(nrow(ptab) == nrow(res1$prediction_table)*2)
  expect_true(nrow(ptab2) == nrow(res1$prediction_table))
})

test_that("combine_prediction_tables supports numeric predictions", {
  set.seed(1)
  observed <- rnorm(8)
  predtab1 <- tibble::tibble(
    .rownum = seq_along(observed),
    roinum = 1L,
    observed = observed,
    predicted = observed + rnorm(length(observed), sd = 0.2)
  )
  predtab2 <- tibble::tibble(
    .rownum = seq_along(observed),
    roinum = 2L,
    observed = observed,
    predicted = observed + rnorm(length(observed), sd = 0.2)
  )

  pooled <- combine_prediction_tables(list(predtab1, predtab2), collapse_regions = TRUE)
  expect_true(is.data.frame(pooled))
  expect_equal(nrow(pooled), length(observed))
  expect_true(all(c(".rownum", "observed", "predicted", "residual", "abs_error", "sq_error") %in% names(pooled)))
  expect_equal(
    pooled$predicted,
    (predtab1$predicted + predtab2$predicted) / 2,
    tolerance = 1e-10
  )
})

test_that("run_regional can return pooled mean and stacked predictions", {
  dset <- gen_sample_dataset(
    c(10, 10, 4),
    nobs = 120,
    nlevels = 2,
    data_mode = "image",
    response_type = "categorical"
  )
  cval <- blocked_cross_validation(dset$design$block_var)

  region_mask <- NeuroVol(
    sample(1:4, size = length(dset$dataset$mask), replace = TRUE),
    space(dset$dataset$mask)
  )
  model <- load_model("sda_notune")
  mspec <- mvpa_model(
    model,
    dset$dataset,
    dset$design,
    model_type = "classification",
    crossval = cval,
    return_predictions = TRUE
  )

  res_mean <- run_regional(mspec, region_mask, pool_predictions = "mean")
  expect_true(is.data.frame(res_mean$pooled_prediction_table))
  expect_true(is.numeric(res_mean$pooled_performance))
  expect_true("Accuracy" %in% names(res_mean$pooled_performance))
  expect_equal(
    nrow(res_mean$pooled_prediction_table),
    dplyr::n_distinct(res_mean$prediction_table$.rownum)
  )

  res_stack <- run_regional(mspec, region_mask, pool_predictions = "stack", stack_folds = 4, stack_seed = 123)
  expect_true(is.data.frame(res_stack$pooled_prediction_table))
  expect_true(is.numeric(res_stack$pooled_performance))
  expect_true("Accuracy" %in% names(res_stack$pooled_performance))
  expect_equal(
    nrow(res_stack$pooled_prediction_table),
    dplyr::n_distinct(res_stack$prediction_table$.rownum)
  )
})

test_that("surface_based mvpa_regional with 5 ROIS runs without error", {
  
  dset <- gen_sample_dataset(c(10,10,4), nobs=100, nlevels=3, data_mode="surface", response_type="categorical")
  cval <- blocked_cross_validation(dset$design$block_var)
  
  maskid <- sample(1:5, size=length(dset$dataset$mask), replace=TRUE)
  region_mask <- NeuroSurface(dset$dataset$train_data@geometry, indices=nodes(dset$dataset$train_data@geometry), 
                              data=maskid)
  
  model <- load_model("sda_notune")
  mspec <- mvpa_model(model, dset$dataset, dset$design, model_type="classification", crossval=cval)
  res <- run_regional(mspec, region_mask, return_fits=TRUE)
  expect_true(!is.null(res))
  
})

test_that("mvpa_regional with 5 ROIS with sda_boot runs without error", {
  
  dataset <- gen_sample_dataset(c(10,10,4), nobs=100, nlevels=3,response_type="categorical")
  cval <- blocked_cross_validation(dataset$design$block_var)
  
  regionMask <- NeuroVol(sample(1:5, size=length(dataset$dataset$mask), replace=TRUE), neuroim2::space(dataset$dataset$mask))
  model <- load_model("sda_boot")
  mspec <- mvpa_model(model, dataset$dataset, dataset$design, model_type="classification", crossval=cval)
  res <- run_regional(mspec, regionMask)
  expect_true(!is.null(res))
})

test_that("mvpa_regional with 5 ROIS with lda_thomaz_boot runs without error", {
  
  dataset <- gen_sample_dataset(c(10,10,4), nobs=100, nlevels=3,response_type="categorical")
  cval <- blocked_cross_validation(dataset$design$block_var)
  
  regionMask <- NeuroVol(sample(1:5, size=length(dataset$dataset$mask), replace=TRUE), neuroim2::space(dataset$dataset$mask))
  model <- load_model("sda_boot")
  mspec <- mvpa_model(model, dataset$dataset, dataset$design, model_type="classification", crossval=cval)
  res <- run_regional(mspec, regionMask)
  expect_true(!is.null(res))
})

 
test_that("mvpa_regional with 5 ROIS with sda_boot and custom_performance runs without error", {

  dataset <- gen_sample_dataset(c(10,10,4), nobs=100, nlevels=3,response_type="categorical")
  cval <- blocked_cross_validation(dataset$design$block_var)
  regionMask <- NeuroVol(sample(1:5, size=length(dataset$dataset$mask), replace=TRUE), neuroim2::space(dataset$dataset$mask))
  model <- load_model("sda_boot")
  
  p <- function(x) { return(list(x=1)) }
  mspec <- mvpa_model(model, dataset$dataset, dataset$design, model_type="classification", crossval=cval, performance=p)
  res <- run_regional(mspec, regionMask)
  expect_true(!is.null(res))
  expect_true("x" %in% names(res$vol_results))
})

test_that("mvpa_regional with 5 ROIS runs and sda without error", {
  tuneGrid <- expand.grid(lambda=c(.01), diagonal=c(TRUE, FALSE))
  model <- load_model("sda")
  
  dataset <- gen_sample_dataset(c(10,10,5), nobs=100, nlevels=4)
  cval <- blocked_cross_validation(dataset$design$block_var)
  
  mspec <- mvpa_model(model, dataset$dataset, dataset$design, model_type="classification", crossval=cval, tune_grid=tuneGrid, tune_reps=10)
  regionMask <- NeuroVol(sample(1:5, size=length(dataset$dataset$mask), replace=TRUE), space(dataset$dataset$mask))
  res <- run_regional(mspec, regionMask, TRUE)
  expect_true(!is.null(res))
})

test_that("mvpa_regional with 5 ROIS and random forest without error", {
  
  model <- load_model("rf")
  dataset <- gen_sample_dataset(c(10,10,5), nobs=100, nlevels=3)
  cval <- blocked_cross_validation(dataset$design$block_var)

  regionMask <- NeuroVol(sample(1:5, size=length(dataset$dataset$mask), replace=TRUE), space(dataset$dataset$mask))
  mspec <- mvpa_model(model, dataset$dataset, dataset$design, model_type="classification", crossval=cval, 
                      tune_grid=data.frame(mtry=c(2,4,6)))
  res <- run_regional(mspec, regionMask)
  expect_true(!is.null(res))
})

test_that("mvpa_regional with 5 ROIS and random forest and k-fold cross-validation without error", {
  
  model <- load_model("rf")
  dataset <- gen_sample_dataset(c(10,10,6), nobs=100, nlevels=3)
  cval <- kfold_cross_validation(length(dataset$design$block_var), nfolds=4)
  
  regionMask <- NeuroVol(sample(1:5, size=length(dataset$dataset$mask), replace=TRUE), space(dataset$dataset$mask))
  mspec <- mvpa_model(model, dataset$dataset, dataset$design, model_type="classification", crossval=cval, 
                      tune_grid=data.frame(mtry=c(2,4,6)))
  res <- run_regional(mspec, regionMask)
  expect_true(!is.null(res))
})

test_that("mvpa_regional with 5 ROIS and corclass and k-fold cross-validation without error", {
  
  model <- load_model("corclass")
  tune_grid <- expand.grid(method=c("pearson", "kendall", "spearman"), robust=c(TRUE,FALSE))
  dataset <- gen_sample_dataset(c(10,10,4), nobs=100, nlevels=3)
  cval <- kfold_cross_validation(length(dataset$design$block_var), nfolds=4)
  
  regionMask <- NeuroVol(sample(1:5, size=length(dataset$dataset$mask), replace=TRUE), space(dataset$dataset$mask))
  mspec <- mvpa_model(model, dataset$dataset, dataset$design, model_type="classification", crossval=cval, 
                      tune_grid=tune_grid)
  res <- run_regional(mspec, regionMask)
  expect_true(!is.null(res))
})

test_that("mvpa_regional with 5 ROIS runs and external test set", {
  
  model <- load_model("rf")
  dataset <- gen_sample_dataset(c(10,10,4), nobs=100, nlevels=3, ntest_obs=200, external_test=TRUE)
  dataset$design$test_design$auxvar = rnorm(nrow(dataset$design$test_design))
  
  cval <- blocked_cross_validation(dataset$design$block_var)
  
  regionMask <- NeuroVol(sample(1:5, size=length(dataset$dataset$mask), replace=TRUE), space(dataset$dataset$mask))
  mspec <- mvpa_model(model, dataset$dataset, dataset$design, 
                      model_type="classification", crossval=cval, 
                      tune_grid=data.frame(mtry=c(2,4,6)))
  
  expect_true(has_test_set(dataset$design))
  res <- run_regional(mspec, regionMask, coalesce_design_vars = TRUE)
  expect_true(!is.null(res))
})



# test_that("mvpa_regional with 5 ROIS and consensus learning runs without error", {
#   
#   dataset <- gen_dataset(c(10,10,2), 100, 3)
#   crossVal <- blocked_cross_validation(dataset$blockVar)
#   regionMask <- NeuroVol(sample(1:5, size=length(dataset$mask), replace=TRUE), space(dataset$mask))
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
#   regionMask <- NeuroVol(sample(1:5, size=length(dataset$mask), replace=TRUE), space(dataset$mask))
#   
#   model <- load_model("sda")
#   res <- mvpa_regional_consensus(dataset, model, regionMask)
#   
# })



test_that("mvpa_regional with 5 ROIs and ANOVA FeatureSelection with topk=10", {
  
  dataset <- gen_sample_dataset(c(10,10,5), nobs=100, nlevels=6)
  cval <- blocked_cross_validation(dataset$design$block_var)
  regionMask <- NeuroVol(sample(1:10, size=length(dataset$dataset$mask), replace=TRUE), space(dataset$dataset$mask))
  model <- load_model("sda")
  fsel <- feature_selector("FTest", "topk", 10)
  mspec <- mvpa_model(model, dataset$dataset, dataset$design, model_type="classification", crossval=cval, feature_selector=fsel)
  res <- run_regional(mspec, regionMask, return_fits=TRUE)
  expect_true(!is.null(res))
})

test_that("mvpa_regional with 5 ROIs and catscore FeatureSelection with top_p=.1", {
  
  dataset <- gen_sample_dataset(c(10,10,5), nobs=100, nlevels=6)
  cval <- blocked_cross_validation(dataset$design$block_var)
  regionMask <- NeuroVol(sample(1:10, size=length(dataset$dataset$mask), replace=TRUE), space(dataset$dataset$mask))
  model <- load_model("sda")
  fsel <- feature_selector("catscore", "top_p", .1)
  mspec <- mvpa_model(model, dataset$dataset, dataset$design, model_type="classification", crossval=cval, feature_selector=fsel)
  res <- run_regional(mspec, regionMask, return_fits=TRUE)
  expect_true(!is.null(res))
})




test_that("mvpa_regional with 5 ROIs and ANOVA FeatureSelection with topp=.4", {
  dataset <- gen_sample_dataset(c(10,10,10), nobs=100, nlevels=2)
  cval <- blocked_cross_validation(dataset$design$block_var)
  regionMask <- NeuroVol(sample(1:5, size=length(dataset$dataset$mask), replace=TRUE), space(dataset$dataset$mask))
  model <- load_model("sda_notune")
  fsel <- feature_selector("FTest", "topp", .5)
  mspec <- mvpa_model(model, dataset$dataset, dataset$design,model_type="classification", crossval=cval, feature_selector=fsel)
  res <- run_regional(mspec, regionMask, return_fits=TRUE)
  expect_true(!is.null(res))
})


test_that("mvpa_regional with regression and 5 ROIs runs without error", {
  
  dataset <- gen_sample_dataset(c(10,10,6), 100, response_type="continuous")
  cval <- blocked_cross_validation(dataset$design$block_var)
  
  regionMask <- NeuroVol(sample(1:4, size=length(dataset$dataset$mask), replace=TRUE), space(dataset$dataset$mask))
  tune_grid <- expand.grid(alpha=c(.1,.5), lambda=c(.001,.2,.5))
  
  model <- load_model("glmnet")
  mspec <- mvpa_model(model, dataset$dataset, dataset$design, model_type="regression", crossval=cval, tune_grid=tune_grid)
  res <- run_regional(mspec, regionMask)
  expect_equal(nobs(dataset$design), 100)
  expect_true(!is.null(res))
  
})

test_that("standard mvpa_regional and custom cross-validation runs without error", {
  dataset <- gen_sample_dataset(c(5,5,5), 100, blocks=3)
  sample_set <- replicate(5, {
    list(train=sample(1:80), test=sample(1:100))
  }, simplify=FALSE)
  
  regionMask <- NeuroVol(sample(1:5, size=length(dataset$dataset$mask), replace=TRUE), space(dataset$dataset$mask))
  cval <- custom_cross_validation(sample_set)
  
  model <- load_model("sda_notune")
  mspec <- mvpa_model(model, dataset$dataset, design=dataset$design, model_type="classification", crossval=cval)
  res <- run_regional(mspec,regionMask)
  expect_true(!is.null(res))
})

test_that("standard mvpa_regional and custom cross-validation and splitting var runs without error", {
  
  dataset <- gen_sample_dataset(c(5,5,5), 100, blocks=3, split_by=factor(rep(1:4, each=25)))
  sample_set <- replicate(5, {
    list(train=sample(1:80), test=sample(1:100))
  }, simplify=FALSE)
  
  regionMask <- NeuroVol(sample(1:5, size=length(dataset$dataset$mask), replace=TRUE), space(dataset$dataset$mask))
  cval <- custom_cross_validation(sample_set)
  
  model <- load_model("sda_notune")
  mspec <- mvpa_model(model, dataset$dataset, design=dataset$design, model_type="classification", crossval=cval)
  res <- run_regional(mspec,regionMask)
  expect_true(!is.null(res))
})


