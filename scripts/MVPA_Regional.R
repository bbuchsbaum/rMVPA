#! /usr/bin/env Rscript

.suppress <- suppressPackageStartupMessages
.suppress(library(neuroim))
.suppress(library(rMVPA))
.suppress(library(optparse))
.suppress(library(caret))
.suppress(library(futile.logger))
.suppress(library(io))
.suppress(library(caret))
.suppress(library(doParallel))
## add "nestingDimension"
## that means features are nested in  a single 3D block
## series(x, ind) # returns unnested matrix.
## needs then to be unravelled into featureMatrix.


## should be able to use "test_subset" without test_data
## that means: train only on one subset and test only on another subset within single design

## if train_design only, the "test_subset" subsets the the train_design. Check to make sure no overlap in columns.
## if test_design provided, the test design must have matching test_data


## should output importance maps for methods that supply them (e.g. sda)



option_list <- list(
                    make_option(c("-t", "--train_design"), type="character", help="the file name of the training design table"),
                    make_option("--test_design", type="character", help="the file name of the test design table"),
                    make_option("--train_data", type="character", help="the name of the training data file as (4D .nii file)"), 
                    make_option("--test_data", type="character", help="the name of the testing data file as (4D .nii file)"),  
                    make_option(c("-n", "--normalize"), action="store_true", type="logical", help="center and scale each image volume"),
                    make_option(c("-m", "--model"), type="character", help="name of the classifier model"),
                    make_option(c("-a", "--mask"), type="character", help="name of binary image mask file (.nii format)"),
                    make_option(c("--parcellation"), type="character", help="name of file containing integer-based image parcellation for cluster-based dimension reduction (.nii format)"),
                    make_option(c("--autobalance"), action="store_true", type="logical", help="balance training samples by upsampling minority classes"),
                    make_option(c("-p", "--pthreads"), type="numeric", help="the number of parallel threads"),
                    make_option(c("-l", "--label_column"), type="character", help="the name of the column in the design file containing the training labels"),
                    make_option(c("--test_label_column"), type="character", help="the name of the column in the test design file containing the test labels"),
                    make_option(c("-o", "--output"), type="character", help="the name of the output folder where results will be placed"),
                    make_option(c("-b", "--block_column"), type="character", help="the name of the column in the design file indicating the block variable used for cross-validation"),
                    make_option(c("-g", "--tune_grid"), type="character", help="string containing grid parameters in the following form: a=\\(1,2\\), b=\\('one', 'two'\\)"),
                    make_option(c("--tune_length"), type="numeric", help="an integer denoting the number of levels for each model tuning parameter"),
                    make_option(c("--savePredictors"), type="logical", action="store_true", help="save model fits (one per ROI) for predicting new data sets (default is FALSE)"),
                    make_option(c("--skipIfFolderExists"), type="logical", action="store_true", help="skip, if output folder already exists"),
                    ## should be skip_if_folder_exists or "overwrite_folder"
                    make_option(c("--output_class_metrics"), type="logical", help="write out performance metrics for each class in multiclass settings"),
                    make_option(c("--ensemble_predictor"), type="logical", help="predictor is based on average prediction of all cross-validated runs"),
                    make_option(c("--bootstrap_replications"), type="numeric", help="number of bootstrapped models to be computed for each cross-validation fold. Final model is average of bootstrapped models."),
                    make_option(c("-c", "--config"), type="character", help="name of configuration file used to specify program parameters"))


oparser <- OptionParser(usage = "MVPA_Regional.R [options]", option_list=option_list)
opt <- parse_args(oparser, positional_arguments=TRUE)
args <- opt$options

flog.info("command line args are ", args, capture=TRUE)

config <- initializeConfiguration(args)
config <- initializeStandardParameters(config, args, "mvpa_regional")


## Regional Specific Params
setArg("savePredictors", config, args, FALSE)
## Regional Specific Params



config <- initializeTuneGrid(args, config)
configParams <- as.list(config)
config <- initializeDesign(config)

rowIndices <- which(config$train_subset)
config$ROIVolume <- loadMask(config)

if (!is.null(config$roi_subset)) {
  form <- try(eval(parse(text=config$roi_subset)))
  if (inherits(form, "try-error")) {
    flog.error("could not parse roi_subset parameter: %s", config$roi_subset)
    stop()
  }
  
  if (class(form) != "formula") {
    flog.error("roi_subset argument must be an expression that starts with a ~ character")
    stop()
  }
  
  res <- as.logical(eval(form[[2]], list(x=config$ROIVolume)))
  
  config$ROIVolume[!res] <- 0
  flog.info("roi_subset contains %s voxels", sum(config$ROIVolume > 0))
  
  if (sum(config$ROIVolume) <= 1) {
    flog.info("ROI must contain more than one voxel, aborting.")
    stop()
  }
}

if (!is.null(config$roi_grouping)) {
  roilist <- lapply(1:length(config$roi_grouping), function(i) {
    grp <- config$roi_grouping[[i]]
    idx <- which(config$ROIVolume %in% grp)
    vol <- makeVolume(refvol=config$ROIVolume)
    vol[idx] <- i
    vol
  })
  
  if (is.null(names(config$roi_grouping))) {
    names(roilist) <- paste0("roi_group_", seq_along(config$roi_grouping))
  } else {
    names(roilist) <- names(config$roi_grouping)
  }
  
  config$ROIVolume <- roilist
} else {
  config$ROIVolume <- list(config$ROIVolume)
}




config$maskVolume <- as(Reduce("+", lapply(config$ROIVolume, function(roivol) as(roivol, "LogicalBrainVolume"))), "LogicalBrainVolume")

parcellationVolume <- if (!is.null(config$parcellation)) {
  loadVolume(config$parcellation)
}

if (!is.null(parcellationVolume)) {
  if (!all(dim(parcellationVolume) == dim(config$maskVolume))) {
    flog.error("dimensions of parcellation volume must equal dimensions of mask volume")
    stop() 
  }
}


config <- initializeData(config)

flog.info("number of training trials: %s", length(rowIndices))
flog.info("max trial index: %s", max(rowIndices))
flog.info("loading training data: %s", config$train_data)
flog.info("mask contains %s voxels", sum(config$maskVolume))

for (i in seq_along(config$ROIVolume)) {
  rvol <- config$ROIVolume[[i]]
  
  flog.info("Region mask contains: %s ROIs", length(unique(rvol[rvol > 0])))
  flog.info(paste("ROIs are for group ", i, "are:"), rvol, capture=TRUE)
 
}


flog.info("Running regional MVPA with parameters:", configParams, capture=TRUE)
flog.info("With %s roi groups", length(config$ROIVolume))

if (length(config$labels) != dim(config$train_datavec)[4]) {
  flog.error("Number of volumes: %s must equal number of labels: %s", dim(config$train_datavec)[4], length(config$labels))
  stop()
}

featureSelector <- if (!is.null(config$feature_selector)) {
  FeatureSelector(config$feature_selector$method, 
                  config$feature_selector$cutoff_type, 
                  as.numeric(config$feature_selector$cutoff_value))
}

flog.info("feature selector: ", featureSelector, capture=TRUE)
flog.info("bootstrap replications: ", config$bootstrap_replications, capture=TRUE)

consensusLearner <- if (!is.null(config$consensus_learner)) {
  flog.info("Will train a %s consensus classifier on ROI results", config$consensus_learner$method)
  ConsensusLearner(config$consensus_learner$method, config$consensus_learner$params)
}

flog.info("consensus learner: ", consensusLearner, capture=TRUE)

if (!is.null(config$custom_performance)) {
  flog.info("custom performance function provided: ", config$custom_performance)
}


dataset <- MVPADataset$new(config$train_datavec, 
                           config$labels, 
                           config$maskVolume, 
                           config$block, 
                           config$test_datavec, 
                           config$testLabels, 
                           parcellation=parcellationVolume,
                           testSplitVar=config$testSplitVar, 
                           testSplits=config$testSplits,
                           trainDesign=config$train_design,
                           testDesign=config$test_design)

model <- loadModel(config$model, list(tuneGrid=config$tune_grid, custom_performance=config$custom_performance))


for (varname in c("test_subset", "train_subset", "roi_subset", "split_by")) {
  if (!is.null(configParams[[varname]]) && is(configParams[[varname]], "formula")) {
    configParams[[varname]] <- Reduce(paste, deparse(configParams[[varname]]))
  }
}


cl <- makeCluster(config$pthreads, outfile="",useXDR=FALSE, type="FORK")
registerDoParallel(cl)


for (roinum in seq_along(config$ROIVolume)) {
  
  gc()
  
  if (length(config$ROIVolume) > 1) {
    roiname <- names(config$ROIVolume)[roinum]
    outdir <- paste0(config$output, "_roigroup_", roiname)
  } else {
    outdir <- config$output
  }
  
  if (config$skipIfFolderExists) {
    if (file.exists(outdir)) {
      flog.info("output folder %s already exists, skipping.", outdir)
      next
    }
  } else {
    outdir <- makeOutputDir(outdir)
  } 
  
  roivol <- config$ROIVolume[[roinum]]
 
  ## need to handle bootstrap reps
  crossVal <- blocked_cross_validation(dataset$blockVar, balance=config$autobalance)
  
  mvpa_res <- mvpa_regional(dataset, model, roivol, crossVal, config$savePredictors, 
                            featureSelector=featureSelector, classMetrics=config$output_class_metrics)
  
  if (!is.null(consensusLearner)) {
    consResult <- mvpa_regional_consensus(dataset, model, roivol, 
                                          featureSelector=featureSelector, 
                                          classMetrics=config$output_class_metrics,
                                          method=consensusLearner$method)
    
    saveRDS(consResult$result$predictor, paste0(config$output, "/consensusPredictor.RDS"))
    rois <- sapply(mvpa_res$resultSet$resultList, attr, "ROINUM")
    #consResult <- consensusWeights(mvpa_res$resultSet, consensusLearner$method)
    consPerf <- performance(consResult$result, dataset$testSplits, config$output_class_metrics)
    consPerf <- as.data.frame(t(consPerf))
    write.table(format(consPerf,  digits=2, scientific=FALSE, drop0trailing=TRUE), paste0(paste0(outdir, "/consensus_performance.txt")), row.names=FALSE, quote=FALSE)
    consWeights <- data.frame(ROINUM=rois, weights=consResult$weights)
    write.table(format(consWeights,  digits=2, scientific=FALSE, drop0trailing=TRUE), paste0(paste0(outdir, "/consensus_weights.txt")), row.names=FALSE, quote=FALSE)
  }

  write.table(format(mvpa_res$performance,  digits=2, scientific=FALSE, drop0trailing=TRUE), paste0(paste0(outdir, "/performance_table.txt")), row.names=FALSE, quote=FALSE)
  
  model$saveResults(mvpa_res$extendedResults, outdir)
  
  lapply(1:length(mvpa_res$outVols), function(i) {
    out <- paste0(outdir, "/", names(mvpa_res$outVols)[i], ".nii")
    writeVolume(mvpa_res$outVols[[i]], out)  
  })
  

  write.table(config$train_subset, paste0(outdir, "/train_table.txt"), row.names=FALSE, quote=FALSE)

  if (!is.null(config$test_subset)) {
    write.table(config$test_subset, paste0(outdir, "/test_table.txt"), row.names=FALSE, quote=FALSE)
  }

  if (!is.null(config$mask)) {
    file.copy(config$mask, paste0(outdir, "/", basename(config$mask)))
  }

  #write.table(format(mvpa_res$predMat,  digits=2, scientific=FALSE, drop0trailing=TRUE), paste0(paste0(config$output, "/predMat_table.txt")), row.names=FALSE, quote=FALSE)
  #saveRDS(mvpa_res$predictorList, paste0(config$output, "/predictorList.RDS"))

  configParams$currentWorkingDir <- getwd()

  configout <- paste0(outdir, "/config.yaml")
  qwrite(configParams, configout)
}


