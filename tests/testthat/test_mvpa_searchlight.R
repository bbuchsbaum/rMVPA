library(neuroim2)
library(neurosurf)

gen_regression_dataset <- function(D, nobs, spacing=c(1,1,1), folds=5) {
  mat <- array(rnorm(prod(D)*nobs), c(D,nobs))
  bspace <- NeuroSpace(c(D,nobs), spacing)
  bvec <- NeuroVec(mat, bspace)
  mask <- as.logical(NeuroVol(array(rep(1, prod(D)), D), NeuroSpace(D, spacing)))
  Y <- rnorm(nobs)
  blockVar <- rep(1:folds, length.out=nobs)
  des <- mvpa_design(data.frame(Y=Y), block_var=blockVar, y_train= ~ Y)
  mvpa_dataset(bvec, mask=mask, design=des)
}

gen_dataset_with_test <- function(D, nobs, nlevels, spacing=c(1,1,1), folds=5, splitvar=TRUE) {
  mat <- array(rnorm(prod(D)*nobs), c(D,nobs))
  bspace <- NeuroSpace(c(D,nobs), spacing)
  bvec <- NeuroVec(mat, bspace)
  mask <- as.logical(NeuroVol(array(rep(1, prod(D)), D), NeuroSpace(D, spacing)))
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
  
  bspace <- NeuroSpace(c(D,nobs), spacing)
  bvec <- NeuroVec(mat, bspace)
  mask <- as.logical(NeuroVol(array(rep(1, prod(D)), D), NeuroSpace(D, spacing)))
  Y <- sample(factor(rep(letters[1:nlevels], length.out=nobs)))
  blockVar <- rep(1:folds, length.out=nobs)
  
  des <- mvpa_design(data.frame(Y=Y), block_var=blockVar, y_train= ~ Y)
  mvpa_dataset(bvec, mask=mask, design=des)
}

gen_surface_dataset <- function(nobs, nlevels, folds=5) {
  library(neurosurf)
  fname <- system.file("extdata/std.lh.smoothwm.asc", package="neuroim2")
  geom <- read_surf_geometry(fname)
  nvert <- nrow(vertices(geom))
  mat <- matrix(rnorm(nvert*nobs), nvert, nobs)
  
  bvec <- NeuroSurfaceVector(geom, 1:nvert, mat)
  Y <- sample(factor(rep(letters[1:nlevels], length.out=nobs)))
  blockVar <- rep(1:folds, length.out=nobs)
  
  des <- mvpa_design(train_design=data.frame(Y=Y, blockVar=blockVar), y_train= ~ Y, block_var=~blockVar)
  mask <- rep(1, nvert)
  dataset <- mvpa_surface_dataset(bvec, mask=mask, design=des, name="lh")
  
  
}


test_that("standard mvpa_searchlight runs without error", {
  
  dataset <- gen_sample_dataset(c(8,8,8), 100, blocks=3)
  cval <- blocked_cross_validation(dataset$design$block_var)
  model <- load_model("sda_notune")
  mspec <- mvpa_model(model, dataset$dataset, design=dataset$design, model_type="classification", crossval=cval)
  res <- run_searchlight(mspec,radius=8, method="standard")
  expect_true(!is.null(res))
  
})

test_that("standard mvpa_searchlight with boot_lda_thomaz runs without error", {
  
  dataset <- gen_sample_dataset(c(8,8,8), 100, blocks=3)
  cval <- blocked_cross_validation(dataset$design$block_var)
  model <- load_model("lda_thomaz")
  grid <- data.frame(nreps=2, frac=.5)
  mspec <- mvpa_model(model, dataset$dataset, design=dataset$design, model_type="classification", crossval=cval, tune_grid=grid)
  res <- run_searchlight(mspec,radius=8, method="standard")
  expect_true(!is.null(res))
  
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
  expect_true(!is.null(res))
  
})

test_that("standard surface-based mvpa_searchlight runs without error", {
  
  dataset <- gen_sample_dataset(D=0, nobs=100,data_mode="surface")
  cval <- blocked_cross_validation(dataset$design$block_var)
  model <- load_model("sda_notune")
  mspec <- mvpa_model(model, dataset$dataset, dataset$design,model_type="classification", crossval=cval)
  res <- run_searchlight(mspec, radius=8, method="standard")
  expect_true(!is.null(res))
  
})

test_that("randomized surface-based mvpa_searchlight runs without error", {
  
  dataset <- gen_sample_dataset(D=100, nobs=100, data_mode="surface")
  cval <- blocked_cross_validation(dataset$design$block_var)
  model <- load_model("sda_notune")
  mspec <- mvpa_model(model, dataset$dataset, dataset$design,model_type="classification", crossval=cval)
  res <- run_searchlight(mspec, radius=8, method="randomized", niter=4)
  expect_true(!is.null(res))
  
  
})

test_that("randomized mvpa_searchlight runs without error", {
  
  dataset <- gen_sample_dataset(c(8,8,8), 100, nlevels=2)
  cval <- blocked_cross_validation(dataset$design$block_var)
  model <- load_model("sda_notune")
  mspec <- mvpa_model(model, dataset$dataset, dataset$design, 
                      model_type="classification", crossval=cval)
  res <- run_searchlight(mspec,radius=5, method="randomized", niter=4)
  expect_true(!is.null(res))


})

test_that("mvpa_searchlight with sda_boot", {
  
  dataset <- gen_sample_dataset(c(8,8,8), 100, nlevels=26, blocks=3)
  cval <- blocked_cross_validation(dataset$design$block_var)
  model <- load_model("sda_boot")
  mspec <- mvpa_model(model, dataset$dataset, dataset$design, model_type="classification", crossval=cval)
  res <- run_searchlight(mspec,radius=6, method="standard", niter=4)
  expect_true(!is.null(res))
  
  
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
   expect_true(!is.null(res))
   

})

test_that("standard mvpa_searchlight and tune_grid runs without error", {
  
  dataset <- gen_sample_dataset(c(3,3,3), 50, nlevels=2, blocks=3)
  cval <- blocked_cross_validation(dataset$design$block_var)
  tuneGrid <- expand.grid(lambda=c(.1,.8), diagonal=c(TRUE, FALSE))
  model <- load_model("sda")
  mspec <- mvpa_model(model, dataset$dataset, dataset$design, model_type="classification", crossval=cval, tune_grid=tuneGrid)
  res <- run_searchlight(mspec, radius=3, method="standard")
  expect_true(!is.null(res))
  
})

test_that("standard mvpa_searchlight and tune_grid with two-fold cross-validation runs without error", {
  
  dataset <- gen_sample_dataset(c(2,2,6), 100, nlevels=2, blocks=2)
  cval <- blocked_cross_validation(dataset$design$block_var)
  
  tuneGrid <- expand.grid(lambda=c(.1,.8), diagonal=c(TRUE, FALSE))
  model <- load_model("sda")
  mspec <- mvpa_model(model, dataset$dataset, dataset$design, model_type="classification", crossval=cval, tune_grid=tuneGrid)
  res <- run_searchlight(mspec, radius=3, method="standard")
  expect_true(!is.null(res))
  
})

test_that("randomized mvpa_searchlight and tune_grid runs without error", {
  
  dataset <- gen_sample_dataset(c(6,6,6), 100, nlevels=2, blocks=2)
  cval <- blocked_cross_validation(dataset$design$block_var)
  
  tuneGrid <- expand.grid(lambda=c(.1,.8), diagonal=c(TRUE, FALSE))
  model <- load_model("sda")
  mspec <- mvpa_model(model, dataset$dataset, dataset$design, model_type="classification", crossval=cval, tune_grid=tuneGrid)
  res <- run_searchlight(mspec, radius=3, method="randomized")
  expect_true(!is.null(res))
  
})

test_that("randomized mvpa_searchlight works with regression", {
  
  dataset <- gen_sample_dataset(c(4,4,4), 100, blocks=3, response_type="continuous")
  cval <- blocked_cross_validation(dataset$design$block_var)
  tuneGrid <- expand.grid(K=3, eta=.5, kappa=.5)
  model <- load_model("spls")
  mspec <- mvpa_model(model, dataset$dataset, dataset$design, model_type="regression", crossval=cval, tune_grid=tuneGrid)
  res <- run_searchlight(mspec, radius=3, niter=2,method="randomized")
  expect_true(!is.null(res))
})

test_that("randomized mvpa_searchlight works with xgboost", {
  
  dataset <- gen_sample_dataset(c(6,6,4), 100, nlevels=5, blocks=5)
  cval <- blocked_cross_validation(dataset$design$block_var)
  tuneGrid <- expand.grid(max_depth=2, eta=.5, nrounds=100,gamma=0,colsample_bytree=.6, min_child_weight=1, subsample=.5)
  model <- load_model("xgbTree")
  mspec <- mvpa_model(model, dataset$dataset, dataset$design, crossval=cval, tune_grid=tuneGrid)
  res <- run_searchlight(mspec, radius=3, niter=2,method="randomized")
  expect_true(!is.null(res))
  
})

test_that("mvpa_searchlight works with testset", {
  require("sda")
  dataset <- gen_sample_dataset(c(4,4,4), 100, response_type="categorical", data_mode="image", blocks=5, nlevels=4, external_test=TRUE, nobs=100)
  
  cval <- blocked_cross_validation(dataset$design$block_var)

  model <- load_model("sda")
  mspec <- mvpa_model(model, dataset$dataset, dataset$design, model_type="classification", crossval=cval)
  res <- run_searchlight( mspec, radius=6, method="randomized")
  expect_true(!is.null(res))
})


 
test_that("mvpa_searchlight works with split_var", {
  dataset <- gen_sample_dataset(c(5,5,5), 100, blocks=3, split_by=factor(rep(1:4, each=25)))

  crossVal <- blocked_cross_validation(dataset$design$block_var)
  tuneGrid <- expand.grid(alpha=.5, lambda=c(.1))
  model <- load_model("glmnet")
  mspec <- mvpa_model(model, dataset$dataset, dataset$design, model_type="classification", 
                      crossval=crossVal, class_metrics=FALSE)
  res <- run_searchlight(mspec, radius=3, niter=2,method="randomized")
  expect_true(!is.null(res))
 
})


test_that("mvpa_searchlight on real data set", {
  tdat <- c(
    system.file("extdata", "sub-1005_task-localizer_run-01_bold_space-MNI152NLin2009cAsym_preproc_betas_small.nii.gz", package="rMVPA"),
    system.file("extdata", "sub-1005_task-localizer_run-02_bold_space-MNI152NLin2009cAsym_preproc_betas_small.nii.gz", package="rMVPA"),
    system.file("extdata", "sub-1005_task-localizer_run-03_bold_space-MNI152NLin2009cAsym_preproc_betas_small.nii.gz", package="rMVPA"),
    system.file("extdata", "sub-1005_task-localizer_run-04_bold_space-MNI152NLin2009cAsym_preproc_betas_small.nii.gz", package="rMVPA")
  )
  
  tdes <- c(
    system.file("extdata", "sub-1005_task-localizer_run-01_events.tsv", package="rMVPA"),
    system.file("extdata", "sub-1005_task-localizer_run-02_events.tsv", package="rMVPA"),
    system.file("extdata", "sub-1005_task-localizer_run-03_events.tsv", package="rMVPA"),
    system.file("extdata", "sub-1005_task-localizer_run-04_events.tsv", package="rMVPA")
  )
  mask <- neuroim2::read_vol(system.file("extdata", "sub-1005-MNI152NLin2009cAsym_small_global_mask.nii", package="rMVPA"))
  
  tvec <- neuroim2::read_vec(tdat)
  des <- do.call(rbind, lapply(tdes,read.table, header=TRUE))
  dset <- mvpa_dataset(tvec, mask=mask)  
  mdes <- mvpa_design(des, y_train = ~ BlockType, block_var=~ Run)
  mod <- mvpa_model(load_model("sda_notune"), dset,mdes,crossval=blocked_cross_validation(des$Run))
  res <- run_searchlight(mod, radius=16, niter=2,method="randomized")
  expect_true(!is.null(res))
  
})

test_that("mvpa_searchlight on real data set", {
  tdat <- c(
    system.file("extdata", "sub-1005_task-localizer_run-01_bold_space-MNI152NLin2009cAsym_preproc_betas_small.nii.gz", package="rMVPA"),
    system.file("extdata", "sub-1005_task-localizer_run-02_bold_space-MNI152NLin2009cAsym_preproc_betas_small.nii.gz", package="rMVPA"),
    system.file("extdata", "sub-1005_task-localizer_run-03_bold_space-MNI152NLin2009cAsym_preproc_betas_small.nii.gz", package="rMVPA"),
    system.file("extdata", "sub-1005_task-localizer_run-04_bold_space-MNI152NLin2009cAsym_preproc_betas_small.nii.gz", package="rMVPA")
  )
  
  tdes <- c(
    system.file("extdata", "sub-1005_task-localizer_run-01_events.tsv", package="rMVPA"),
    system.file("extdata", "sub-1005_task-localizer_run-02_events.tsv", package="rMVPA"),
    system.file("extdata", "sub-1005_task-localizer_run-03_events.tsv", package="rMVPA"),
    system.file("extdata", "sub-1005_task-localizer_run-04_events.tsv", package="rMVPA")
  )
  mask <- neuroim2::read_vol(system.file("extdata", "sub-1005-MNI152NLin2009cAsym_small_global_mask.nii", package="rMVPA"))
  
  tvec <- neuroim2::read_vec(tdat)
  des <- do.call(rbind, lapply(tdes,read.table, header=TRUE))
  dset <- mvpa_dataset(tvec, mask=mask)  
  mdes <- mvpa_design(des, y_train = ~ BlockType, block_var=~ Run)
  mod <- mvpa_model(load_model("sda_notune"), dset,mdes,crossval=blocked_cross_validation(des$Run))
  res <- run_searchlight(mod, radius=16, niter=3,method="randomized")
  expect_true(!is.null(res))
  
})

test_that("mvpa_searchlight on real data set with testset", {
  tdat1 <- c(
    system.file("extdata", "sub-1005_task-localizer_run-01_bold_space-MNI152NLin2009cAsym_preproc_betas_small.nii.gz", package="rMVPA"),
    system.file("extdata", "sub-1005_task-localizer_run-02_bold_space-MNI152NLin2009cAsym_preproc_betas_small.nii.gz", package="rMVPA"),
    system.file("extdata", "sub-1005_task-localizer_run-03_bold_space-MNI152NLin2009cAsym_preproc_betas_small.nii.gz", package="rMVPA"),
    system.file("extdata", "sub-1005_task-localizer_run-04_bold_space-MNI152NLin2009cAsym_preproc_betas_small.nii.gz", package="rMVPA")
  )
  
  tdat2 <- c(
    system.file("extdata", "sub-1005_task-wm_run-01_bold_space-MNI152NLin2009cAsym_preproc_betas_small.nii.gz", package="rMVPA"),
    system.file("extdata", "sub-1005_task-wm_run-02_bold_space-MNI152NLin2009cAsym_preproc_betas_small.nii.gz", package="rMVPA"),
    system.file("extdata", "sub-1005_task-wm_run-03_bold_space-MNI152NLin2009cAsym_preproc_betas_small.nii.gz", package="rMVPA"),
    system.file("extdata", "sub-1005_task-wm_run-04_bold_space-MNI152NLin2009cAsym_preproc_betas_small.nii.gz", package="rMVPA")
  )
  
  tdes1 <- c(
    system.file("extdata", "sub-1005_task-localizer_run-01_events.tsv", package="rMVPA"),
    system.file("extdata", "sub-1005_task-localizer_run-02_events.tsv", package="rMVPA"),
    system.file("extdata", "sub-1005_task-localizer_run-03_events.tsv", package="rMVPA"),
    system.file("extdata", "sub-1005_task-localizer_run-04_events.tsv", package="rMVPA")
  )
  
  tdes2 <- c(
    system.file("extdata", "sub-1005_task-wm_run-01_events.tsv", package="rMVPA"),
    system.file("extdata", "sub-1005_task-wm_run-02_events.tsv", package="rMVPA"),
    system.file("extdata", "sub-1005_task-wm_run-03_events.tsv", package="rMVPA"),
    system.file("extdata", "sub-1005_task-wm_run-04_events.tsv", package="rMVPA")
  )
  
  
  mask <- neuroim2::read_vol(system.file("extdata", "sub-1005-MNI152NLin2009cAsym_small_global_mask.nii", package="rMVPA"))
  
  tvec1 <- neuroim2::read_vec(tdat1)
  tvec2 <- neuroim2::read_vec(tdat2)
  des1 <- do.call(rbind, lapply(tdes1,read.table, header=TRUE))
  des2 <- do.call(rbind, lapply(tdes2,read.table, header=TRUE))
  library(dplyr)
  
  des3 <- des2 %>% mutate(combo = 
                    case_when(
                      (Cue == "Faces" & ToBeIgnored == "Scenes") | (Cue == "Scenes" & ToBeIgnored == "Faces") ~ "Faces-Scenes",
                      (Cue == "Bodies" & ToBeIgnored == "Scenes") | (Cue == "Scenes" & ToBeIgnored == "Bodies") ~ "Bodies-Scenes",
                      (Cue == "Objects" & ToBeIgnored == "Scenes") | (Cue == "Scenes" & ToBeIgnored == "Objects") ~ "Objects-Scenes",
                      (Cue == "Objects" & ToBeIgnored == "Faces") | (Cue == "Faces" & ToBeIgnored == "Objects") ~ "Faces-Objects",
                      (Cue == "Faces" & ToBeIgnored == "Bodies") | (Cue == "Bodies" & ToBeIgnored == "Faces") ~ "Faces-Bodies",
                      (Cue == "Objects" & ToBeIgnored == "Bodies") | (Cue == "Bodies" & ToBeIgnored == "Objects") ~ "Objects-Bodies")
  )
                      
                      
  
  dset <- mvpa_dataset(tvec1, tvec2, mask=mask)  
  mdes <- mvpa_design(train_design=des1, y_train = ~ BlockType, test_design=des3, y_test = ~ Cue, block_var=~ Run, split_by = ~ combo)
  mod <- mvpa_model(load_model("sda_notune"), dset,mdes,
                    crossval=blocked_cross_validation(des1$Run), class_metrics=FALSE)
  res <- run_searchlight(mod, radius=16, niter=3,method="randomized", 
                         combiner=function(mspec, good, bad) {
                           good
                         })
    expect_true(!is.null(res))
})










