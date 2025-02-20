test_that("regional feature_rsa_model with direct F matrix runs without error", {
  # Generate a sample dataset: small volume, say 5x5x5, with 100 observations
  dset <- gen_sample_dataset(c(5,5,5), 100, blocks=3)
  
  # Create a feature matrix F, e.g., 100 observations x 20 features
  Fmat <- matrix(rnorm(100*20), 100, 20)
  labels <- paste0("obs", 1:100)
  
  # Create feature_rsa_design using direct F matrix
  fdes <- feature_rsa_design(F=Fmat, labels=labels)
  
  # Create a feature_rsa_model, for example using 'pls'
  mspec <- feature_rsa_model(dset$dataset, fdes, method="pca", crossval=blocked_cross_validation(dset$design$block_var))
  
  # Create a region mask with 5 ROIs
  region_mask <- NeuroVol(sample(1:5, size=length(dset$dataset$mask), replace=TRUE), space(dset$dataset$mask))
  
  # Run regional analysis
  res <- run_regional(mspec, region_mask)
  
  # Check results
  expect_true(!is.null(res))
  expect_s3_class(res, "regional_mvpa_result")
  expect_true(!is.null(res$performance_table))
})

test_that("regional feature_rsa_model with S-based feature extraction runs without error", {
  # Generate a sample dataset: again 5x5x5 volume, 100 obs
  dset <- gen_sample_dataset(c(6,5,5), 100, blocks=3)
  
  # Create a similarity matrix S: must be symmetric and match number of observations
  obs_features <- matrix(rnorm(100*10), 100, 10)
  S <- tcrossprod(scale(obs_features))  # similarity matrix
  labels <- paste0("obs", 1:100)
  
  # Create feature_rsa_design using S
  fdes <- feature_rsa_design(S=S, labels=labels, k=5) # reduce to 5 dims
  
  # Create a feature_rsa_model using scca this time
  mspec <- feature_rsa_model(dset$dataset, fdes, method="scca", crossval=blocked_cross_validation(dset$design$block_var))
  
  # Create a region mask with 5 ROIs
  region_mask <- NeuroVol(sample(1:5, size=length(dset$dataset$mask), replace=TRUE), space(dset$dataset$mask))
  
  # Run regional analysis
  res <- run_regional(mspec, region_mask)
  
  # Check results
  expect_true(!is.null(res))
  expect_s3_class(res, "regional_mvpa_result")
  expect_true(!is.null(res$performance_table))
  
  # If you want, check that performance_table or prediction_table have expected columns
  # e.g., expect_true("mse" %in% names(res$performance_table))
})

test_that("regional feature_rsa_model with pca method runs without error and returns performance", {
  dset <- gen_sample_dataset(c(4,4,4), 80, blocks=4)
  
  # Create an F matrix directly
  Fmat <- matrix(rnorm(80*15), 80, 15)
  labels <- paste0("obs", 1:80)
  
  fdes <- feature_rsa_design(F=Fmat, labels=labels) # no dimension reduction, just use as is
  mspec <- feature_rsa_model(dset$dataset, fdes, method="pca", crossval=blocked_cross_validation(dset$design$block_var))
  
  region_mask <- NeuroVol(sample(1:5, size=length(dset$dataset$mask), replace=TRUE), space(dset$dataset$mask))
  
  res <- run_regional(mspec, region_mask)
  
  expect_true(!is.null(res))
  expect_s3_class(res, "regional_mvpa_result")
  # Check if performance metrics are available
  expect_true("performance_table" %in% names(res))
  # Possibly check at least one performance metric
  # e.g., expect_true("mse" %in% names(res$performance_table))
})


test_that("PCA method can learn to reconstruct simple linear features", {
  set.seed(42)
  
  # Simulate data: X is (n_obs x n_voxels)
  n_obs <- 50
  n_vox <- 10
  X <- matrix(rnorm(n_obs * n_vox), nrow = n_obs, ncol = n_vox)
  
  # Create a true linear mapping from X to F
  # We'll create 3 feature dimensions from the 10 input dims
  B <- matrix(rnorm(n_vox * 3, sd=1), nrow=n_vox, ncol=3)
  # F = X * B + noise
  Ftrue <- X %*% B + matrix(rnorm(n_obs * 3, sd = 0.1), n_obs, 3)
  
  # Minimal dataset. We'll just store X in a "fake" NeuroVec and mask
  # In real usage, you'd have an actual NeuroVec from your fMRI data, etc.
  mask_vol <- neuroim2::NeuroVol(array(1, c(1,1,1)), neuroim2::NeuroSpace(c(1,1,1)))
  # For demonstration, let's store X in a single-voxel NeuroVec with 50 observations
  # This obviously is a contrived usage, but it triggers the methods we need.
  xspace <- neuroim2::NeuroSpace(c(1,1,1,n_obs))
  xarr <- array(X, dim=c(1,1,1,n_obs)) 
  train_vec <- neuroim2::NeuroVec(xarr, xspace)
  
  ds <- mvpa_dataset(train_data=train_vec, mask=mask_vol)
  
  # The design that holds F
  design <- feature_rsa_design(F = Ftrue, labels = 1:n_obs)
  
  # Create model with PCA method
  # We can skip real cross-validation here and do 'crossval = blocked_cross_validation(...)' or similar
  fake_cv <- blocked_cross_validation(block_var=rep(1, n_obs)) 
  model_spec <- feature_rsa_model(dataset=ds, design=design, method="pca", crossval=fake_cv)
  
  # Train
  trained_fit <- train_model.feature_rsa_model(model_spec, train_dat=X, ytrain=Ftrue, indices=1:n_obs)
  
  # Predict on the same data
  F_pred <- predict_model.feature_rsa_model(model_spec, trained_fit, newdata=X)
  
  # Evaluate
  perf <- evaluate_model.feature_rsa_model(model_spec, F_pred, Ftrue)
  
  # Expect correlations to be quite high
  expect_gt(mean(perf$correlations), 0.8)
  # And MSE to be reasonably small
  expect_lt(perf$mse, 0.2)
})

test_that("PLS method can learn to reconstruct simple linear features", {
  set.seed(42)
  
  n_obs <- 50
  n_vox <- 10
  X <- matrix(rnorm(n_obs * n_vox), n_obs, n_vox)
  
  B <- matrix(rnorm(n_vox * 3, sd=1), n_vox, 3)
  Ftrue <- X %*% B + matrix(rnorm(n_obs * 3, sd = 0.1), n_obs, 3)
  
  mask_vol <- neuroim2::NeuroVol(array(1, c(1,1,1)), neuroim2::NeuroSpace(c(1,1,1)))
  xarr <- array(X, dim=c(1,1,1,n_obs)) 
  train_vec <- neuroim2::NeuroVec(xarr, neuroim2::NeuroSpace(c(1,1,1,n_obs)))
  ds <- mvpa_dataset(train_data=train_vec, mask=mask_vol)
  
  design <- feature_rsa_design(F = Ftrue, labels = 1:n_obs)
  fake_cv <- blocked_cross_validation(block_var=rep(1, n_obs)) 
  model_spec <- feature_rsa_model(dataset=ds, design=design, method="pls", crossval=fake_cv)
  
  trained_fit <- train_model.feature_rsa_model(model_spec, train_dat=X, ytrain=Ftrue, indices=1:n_obs)
  F_pred <- predict_model.feature_rsa_model(model_spec, trained_fit, newdata=X)
  perf <- evaluate_model.feature_rsa_model(model_spec, F_pred, Ftrue)
  
  expect_gt(mean(perf$correlations), 0.8)
  expect_lt(perf$mse, 0.2)
})

test_that("SCCA method can learn to reconstruct correlated features", {
  # SCCA typically is used for correlated sets of variables. 
  # We'll build a scenario with correlated X ~ F.
  set.seed(42)
  n_obs <- 60
  n_vox <- 10
  
  X <- matrix(rnorm(n_obs*n_vox), n_obs, n_vox)
  # We'll create 2 correlated feature dims from X:
  B <- matrix(rnorm(n_vox*2, sd=1), n_vox, 2)
  Ftrue <- X %*% B + matrix(rnorm(n_obs*2, sd = 0.1), n_obs, 2)
  
  mask_vol <- neuroim2::NeuroVol(array(1, c(1,1,1)), neuroim2::NeuroSpace(c(1,1,1)))
  xarr <- array(X, dim=c(1,1,1,n_obs)) 
  train_vec <- neuroim2::NeuroVec(xarr, neuroim2::NeuroSpace(c(1,1,1,n_obs)))
  ds <- mvpa_dataset(train_data=train_vec, mask=mask_vol)
  
  design <- feature_rsa_design(F = Ftrue, labels = 1:n_obs)
  fake_cv <- blocked_cross_validation(block_var=rep(1, n_obs)) 
  model_spec <- feature_rsa_model(dataset=ds, design=design, method="scca", crossval=fake_cv)
  
  trained_fit <- train_model.feature_rsa_model(model_spec, train_dat=X, ytrain=Ftrue, indices=1:n_obs)
  F_pred <- predict_model.feature_rsa_model(model_spec, trained_fit, newdata=X)
  perf <- evaluate_model.feature_rsa_model(model_spec, F_pred, Ftrue)
  
  # Because SCCA can be sensitive to scaling and dimension,
  # we won't set the bar too high, but we still expect decent correlation.
  expect_gt(mean(perf$correlations), 0.7)

})


# test_that("can run feature_rsa with different methods", {
#   data_out <- rMVPA::gen_sample_dataset(D = c(6,6,6), nobs = 50, blocks = 4, nlevels = 2)
#   print(data_out)
#   crossval <- blocked_cross_validation(data_out$design$block_var)
#   # Compare different methods
#   methods <- c("scca", "pca", "pls") 
#   results_list <- lapply(methods, function(method) {
#     model <- feature_rsa_model(
#       dataset = data_out$dset,
#       design = feature_design,
#       method = method,
#       crossval = crossval  # Add cross-validation
#     )
#     run_regional(model, region_mask)
#   })
#   
#   # Compare performance
#   for (i in seq_along(methods)) {
#     cat("\nMethod:", methods[i], "\n")
#     print(results_list[[i]]$performance_table)
#   }
# })