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
  dataset <- mvpa_surface_dataset(bvec, mask=mask, design=des, name="lh")
  
  
}


test_that("standard mvpa_searchlight runs without error", {
  
  dataset <- gen_sample_dataset(c(5,5,5), 100, blocks=3)
  cval <- blocked_cross_validation(dataset$design$block_var)
  model <- load_model("sda_notune")
  mspec <- mvpa_model(model, dataset$dataset, design=dataset$design, model_type="classification", crossval=cval)
  res <- run_searchlight(mspec,radius=4, method="standard", ncores=2)
  
})

test_that("standard mvpa_searchlight and custom cross-validation runs without error", {
  
  dataset <- gen_sample_dataset(c(5,5,5), 100, blocks=3)
  sample_set <- replicate(5, {
    list(train=sample(1:80), test=sample(81:100))
  }, simplify=FALSE)
  cval <- custom_cross_validation(sample_set)
 
  model <- load_model("sda_notune")
  mspec <- mvpa_model(model, dataset$dataset, design=dataset$design, model_type="classification", crossval=cval)
  res <- run_searchlight(mspec,radius=4, method="standard")
  
})

test_that("standard surface-based mvpa_searchlight runs without error", {
  
  dataset <- gen_sample_dataset(D=0, nobs=100,data_mode="surface")
  cval <- blocked_cross_validation(dataset$design$block_var)
  model <- load_model("sda_notune")
  mspec <- mvpa_model(model, dataset$dataset, dataset$design,model_type="classification", crossval=cval)
  res <- run_searchlight(mspec, radius=8, method="standard")
  
})

test_that("randomized surface-based mvpa_searchlight runs without error", {
  
  dataset <- gen_sample_dataset(D=100, nobs=100, data_mode="surface")
  cval <- blocked_cross_validation(dataset$design$block_var)
  model <- load_model("sda_notune")
  mspec <- mvpa_model(model, dataset$dataset, dataset$design,model_type="classification", crossval=cval)
  res <- run_searchlight(mspec, radius=8, method="randomized", niter=5)
  
  
})

test_that("randomized mvpa_searchlight runs without error", {
  
  dataset <- gen_sample_dataset(c(5,5,5), 100, nlevels=2)
  cval <- blocked_cross_validation(dataset$design$block_var)
  model <- load_model("sda_notune")
  mspec <- mvpa_model(model, dataset$dataset, dataset$design, model_type="classification", crossval=cval)
  res <- run_searchlight(mspec,radius=3, method="randomized", niter=4)


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
   
   dataset <- gen_sample_dataset(c(5,5,1), 100, nlevels=3)
   cval <- blocked_cross_validation(dataset$design$block_var)
   model <- load_model("sda_notune")
   mspec <- mvpa_model(model, dataset$dataset, dataset$design, model_type="classification", crossval=cval, performance=custom)
   res <- run_searchlight(mspec,radius=3, method="randomized")
   

})

test_that("standard mvpa_searchlight and tune_grid runs without error", {
  
  dataset <- gen_sample_dataset(c(3,3,3), 50, nlevels=2, blocks=3)
  cval <- blocked_cross_validation(dataset$design$block_var)
  tuneGrid <- expand.grid(lambda=c(.1,.8), diagonal=c(TRUE, FALSE))
  model <- load_model("sda")
  mspec <- mvpa_model(model, dataset$dataset, dataset$design, model_type="classification", crossval=cval, tune_grid=tuneGrid)
  res <- run_searchlight(mspec, radius=3, method="standard")
  
})

test_that("standard mvpa_searchlight and tune_grid with two-fold cross-validation runs without error", {
  
  dataset <- gen_sample_dataset(c(2,2,6), 100, nlevels=2, blocks=2)
  cval <- blocked_cross_validation(dataset$design$block_var)
  
  tuneGrid <- expand.grid(lambda=c(.1,.8), diagonal=c(TRUE, FALSE))
  model <- load_model("sda")
  mspec <- mvpa_model(model, dataset$dataset, dataset$design, model_type="classification", crossval=cval, tune_grid=tuneGrid)
  res <- run_searchlight(mspec, radius=3, method="standard")
  
})

test_that("randomized mvpa_searchlight and tune_grid runs without error", {
  
  dataset <- gen_sample_dataset(c(2,2,6), 100, nlevels=2, blocks=2)
  cval <- blocked_cross_validation(dataset$design$block_var)
  
  tuneGrid <- expand.grid(lambda=c(.1,.8), diagonal=c(TRUE, FALSE))
  model <- load_model("sda")
  mspec <- mvpa_model(model, dataset$dataset, dataset$design, model_type="classification", crossval=cval, tune_grid=tuneGrid)
  res <- run_searchlight(mspec, radius=3, method="randomized")
  
})

test_that("randomized mvpa_searchlight works with regression", {
  
  dataset <- gen_sample_dataset(c(4,4,4), 100, blocks=3, response_type="continuous")
  cval <- blocked_cross_validation(dataset$design$block_var)
  tuneGrid <- expand.grid(alpha=.5, lambda=c(.1,.01))
  model <- load_model("glmnet")
  mspec <- mvpa_model(model, dataset$dataset, dataset$design, model_type="regression", crossval=cval, tune_grid=tuneGrid)
  res <- run_searchlight(mspec, radius=3, niter=2,method="randomized")
  
})

test_that("mvpa_searchlight works with testset", {
  require("sparsediscrim")
  dataset <- gen_sample_dataset(c(4,4,4), 100, response_type="categorical", data_mode="surface", blocks=5, nlevels=4, external_test=TRUE, nobs=100)
  
  cval <- blocked_cross_validation(dataset$design$block_var)

  model <- load_model("lda_thomaz")
  mspec <- mvpa_model(model, dataset$dataset, dataset$design, model_type="classification", crossval=cval)
  res <- run_searchlight( mspec, radius=6, method="randomized")
  
})


 
test_that("mvpa_searchlight works with testset and split_var", {
  dataset <- gen_sample_dataset(c(5,5,5), 100, blocks=3, split_by=factor(rep(1:4, each=25)))

  crossVal <- blocked_cross_validation(dataset$design$block_var)
  tuneGrid <- expand.grid(alpha=.5, lambda=c(.1))
  model <- load_model("glmnet")
  mspec <- mvpa_model(model, dataset$dataset, dataset$design, model_type="classification", crossval=crossVal, class_metrics=FALSE)
  res <- run_searchlight(mspec, radius=3, niter=2,method="randomized")
 
})





