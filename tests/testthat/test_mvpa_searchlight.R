library(neuroim)
library(neurosurf)

gen_regression_dataset <- function(D, nobs, spacing=c(1,1,1), folds=5) {
  mat <- array(rnorm(prod(D)*nobs), c(D,nobs))
  bspace <- BrainSpace(c(D,nobs), spacing)
  bvec <- BrainVector(mat, bspace)
  mask <- as.logical(BrainVolume(array(rep(1, prod(D)), D), BrainSpace(D, spacing)))
  Y <- rnorm(nobs)
  blockVar <- rep(1:folds, length.out=nobs)
  des <- mvpa_design(data.frame(Y=Y), block_var=blockVar, y_train= ~ Y)
  mvpa_dataset(bvec, mask=mask, design=des)
}

gen_dataset_with_test <- function(D, nobs, nlevels, spacing=c(1,1,1), folds=5, splitvar=TRUE) {
  mat <- array(rnorm(prod(D)*nobs), c(D,nobs))
  bspace <- BrainSpace(c(D,nobs), spacing)
  bvec <- BrainVector(mat, bspace)
  mask <- as.logical(BrainVolume(array(rep(1, prod(D)), D), BrainSpace(D, spacing)))
  Y <- sample(factor(rep(letters[1:nlevels], length.out=nobs)))
  Ytest <- rev(Y)
  blockVar <- rep(1:folds, length.out=nobs)
  
  
  if (splitvar) {
    tsplit <- factor(rep(1:5, length.out=length(Y)))
    dframe <- data.frame(Y=Y, Ytest=Ytest, tsplit=tsplit)
    des <- mvpa_design(train_design=dframe, test_design=dframe, block_var=blockVar, y_train= ~ Y, y_test= ~ Ytest, split_by= ~ tsplit)
    mvpa_dataset(train_data=bvec, test_data=bvec, mask=mask, design=des)
  } else {
    des <- mvpa_design(data.frame(Y=Y, Ytest=Ytest), block_var=blockVar, y_train= ~ Y, y_test= ~ Ytest)
    mvpa_dataset(train_data=bvec, test_data=bvec, mask=mask, design=des)
  }
    
}

gen_dataset <- function(D, nobs, nlevels, spacing=c(1,1,1), folds=5) {
  
  mat <- array(rnorm(prod(D)*nobs), c(D,nobs))
  xbad <- array(runif(prod(D)) < .02, D)
  xbad.ind <- which(xbad, arr.ind=TRUE)
  
  for (i in 1:nrow(xbad.ind)) {
    ind <- xbad.ind[i,]
    mat[ind[1], ind[2], ind[3],] <- 0
  }
  
  bspace <- BrainSpace(c(D,nobs), spacing)
  bvec <- BrainVector(mat, bspace)
  mask <- as.logical(BrainVolume(array(rep(1, prod(D)), D), BrainSpace(D, spacing)))
  Y <- sample(factor(rep(letters[1:nlevels], length.out=nobs)))
  blockVar <- rep(1:folds, length.out=nobs)
  
  des <- mvpa_design(data.frame(Y=Y), block_var=blockVar, y_train= ~ Y)
  mvpa_dataset(bvec, mask=mask, design=des)
}

gen_surface_dataset <- function(nobs, nlevels, folds=5) {
  library(neurosurf)
  fname <- system.file("extdata/std.lh.smoothwm.asc", package="neuroim")
  geom <- loadSurface(fname)
  nvert <- nrow(vertices(geom))
  mat <- matrix(rnorm(nvert*nobs), nvert, nobs)
  
  bvec <- BrainSurfaceVector(geom, 1:nvert, mat)
  Y <- sample(factor(rep(letters[1:nlevels], length.out=nobs)))
  blockVar <- rep(1:folds, length.out=nobs)
  
  des <- mvpa_design(train_design=data.frame(Y=Y, blockVar=blockVar), y_train= ~ Y, block_var=~blockVar)
  mask <- rep(1, nvert)
  dataset <- mvpa_surface_dataset(bvec, mask=mask, design=des, hemisphere="lh")
  
  
}




test_that("standard mvpa_searchlight runs without error", {
  
  dataset <- gen_dataset(c(5,5,5), 100, 2)
  cval <- blocked_cross_validation(dataset$design$block_var)
  model <- load_model("sda_notune")
  mspec <- mvpa_model(model, dataset,model_type="classification", crossval=cval)
  res <- run_searchlight(mspec,radius=3, method="standard")
  
})

test_that("standard surface-based mvpa_searchlight runs without error", {
  
  dataset <- gen_surface_dataset(100, 6)
  cval <- blocked_cross_validation(dataset$design$block_var)
  model <- load_model("sda_notune")
  mspec <- mvpa_model(model, dataset,model_type="classification", crossval=cval)
  res <- run_searchlight(mspec, radius=8, method="standard")
  
})

test_that("randomized surface-based mvpa_searchlight runs without error", {
  
  dataset <- gen_surface_dataset(100, 6)
  cval <- blocked_cross_validation(dataset$design$block_var)
  model <- load_model("sda_notune")
  mspec <- mvpa_model(model, dataset,model_type="classification", crossval=cval)
  res <- run_searchlight(mspec, radius=8, method="randomized", niter=5)
  
})

test_that("randomized mvpa_searchlight runs without error", {
  
  dataset <- gen_dataset(c(5,5,5), 100, 2)
  cval <- blocked_cross_validation(dataset$design$block_var)
  model <- load_model("sda_notune")
  mspec <- mvpa_model(model, dataset, model_type="classification", crossval=cval)
  res <- run_searchlight(mspec,radius=3, method="randomized", niter=10)
  
})


# test_that("clustered mvpa_searchlight runs without error", {
#   
#   dataset <- gen_dataset(c(5,5,1), 100, 2)
#   crossVal <- blocked_cross_validation(dataset$blockVar)
#   model <- load_model("sda_notune", list(tuneGrid=NULL))
#   res <- mvpa_clustered_searchlight(dataset, model, crossVal, nclusters=c(2,3,4))
#   
# })

 test_that("randomized mvpa_searchlight runs with custom_performance", {
   
   custom <- function(x) {
     
     cnames <- colnames(x$probs)
     y <- x$observed
     
     p1 <- x$probs[cbind(1:nrow(x$probs), as.integer(y))]
     ret <- c(m1 = mean(p1), m2=max(p1))
     ret
   }
   
   dataset <- gen_dataset(c(5,5,1), 100, 2)
   cval <- blocked_cross_validation(dataset$blockVar)
   model <- load_model("sda_notune")$model
   mspec <- mvpa_model(model, dataset, model_type="classification", crossval=cval, performance=custom)
   res <- run_searchlight(mspec,radius=3, method="randomized")
   

})

test_that("standard mvpa_searchlight and tune_grid runs without error", {
  
  dataset <- gen_dataset(c(3,3,3), 50, 2, folds=3)
  cval <- blocked_cross_validation(dataset$design$block_var)
  tuneGrid <- expand.grid(lambda=c(.1,.8), diagonal=c(TRUE, FALSE))
  model <- load_model("sda")
  mspec <- mvpa_model(model, dataset, model_type="classification", crossval=cval, tune_grid=tuneGrid)
  res <- run_searchlight(mspec, radius=3, method="standard")
  
})

test_that("standard mvpa_searchlight and tune_grid with two-fold cross-validation runs without error", {
  
  dataset <- gen_dataset(c(2,2,6), 100, 2, folds=2)
  cval <- blocked_cross_validation(dataset$design$block_var)
  
  tuneGrid <- expand.grid(lambda=c(.1,.8), diagonal=c(TRUE, FALSE))
  model <- load_model("sda")
  mspec <- mvpa_model(model, dataset, model_type="classification", crossval=cval, tune_grid=tuneGrid)
  res <- run_searchlight(mspec, radius=3, method="standard")
  
})

test_that("randomized mvpa_searchlight and tune_grid runs without error", {
  
  dataset <- gen_dataset(c(2,2,6), 100, 2, folds=2)
  cval <- blocked_cross_validation(dataset$design$block_var)
  
  tuneGrid <- expand.grid(lambda=c(.1,.8), diagonal=c(TRUE, FALSE))
  model <- load_model("sda")$model
  mspec <- mvpa_model(model, dataset, model_type="classification", crossval=cval, tune_grid=tuneGrid)
  res <- run_searchlight(mspec, radius=3, method="randomized")
  
})

test_that("randomized mvpa_searchlight works with regression", {
  
  dataset <- gen_regression_dataset(c(4,4,4), 100, folds=3)
  crossVal <- blocked_cross_validation(dataset$design$block_var)
  tuneGrid <- expand.grid(alpha=.5, lambda=c(.1,.01))
  model <- load_model("glmnet")$model
  mspec <- model_spec(model, dataset, model_type="regression", crossval=cval, tune_grid=tuneGrid)
  res <- run_searchlight(mspec, radius=3, niter=2,method="randomized")
  
})

test_that("mvpa_searchlight works with testset", {
  
  dataset <- gen_dataset_with_test(c(4,4,4), 100, 3, folds=3)
  cval <- blocked_cross_validation(dataset$design$block_var)
  tuneGrid <- expand.grid(alpha=.5, lambda=c(.1,.2))
  model <- load_model("glmnet")$model
  mspec <- mvpa_model(model, dataset, model_type="classification", crossval=cval, tune_grid=tuneGrid)
  res <- run_searchlight( mspec, radius=3, method="standard")
  
})

test_that("mvpa_searchlight works with testset and split_var", {
  
  dataset <- gen_dataset_with_test(c(4,4,4), 100, 3, folds=3, splitvar=TRUE)
  crossVal <- blocked_cross_validation(dataset$blockVar)
  tuneGrid <- expand.grid(alpha=.5, lambda=c(.1))
  model <- load_model("glmnet")$model
  mspec <- mvpa_model(model, dataset, model_type="classification", crossval=cval, tune_grid=tuneGrid, class_metrics=FALSE)
  res <- run_searchlight(mspec, radius=3, niter=2,method="randomized")
  
})





