#! /usr/bin/env Rscript

.suppress <- suppressPackageStartupMessages
.suppress(library(neuroim))
.suppress(library(rMVPA))
.suppress(library(optparse))
.suppress(library(caret))
.suppress(library(futile.logger))
.suppress(library(io))
.suppress(library(caret))

## add "nestingDimension"
## that means features are nested in  a single 3D block
## series(x, ind) # returns unnested matrix.
## needs then to be unravelled into featureMatrix.


## should be able to use "test_subset" without test_data
## that means: train only on one subset and test only on another subset within single design

## if train_design only, the "test_subset" subsets the the train_design. Check to make sure no overlap in columns.
## if test_design provided, the test design must have matching test_data



option_list <- list(
                    make_option(c("-t", "--train_design"), type="character", help="the file name of the design table"),
                    make_option("--test_design", type="character", help="the file name of the design table"),
                    make_option(c("-s", "--type"), type="character", help="the type of searchlight: standard or randomized"),  
                    make_option("--train_data", type="character", help="the name of the training data file as (4D .nii file)"), 
                    make_option("--test_data", type="character", help="the name of the testing data file as (4D .nii file)"),  
                    make_option(c("-n", "--normalize"), action="store_true", type="logical", help="center and scale each volume vector"),
                    make_option(c("-m", "--model"), type="character", help="name of the classifier model"),
                    make_option(c("-a", "--mask"), type="character", help="name of binary image mask file (.nii format)"),
                    make_option(c("--parcellation"), type="character", help="name of file containing integer-based image parcellation for cluster-based dimension reduction (.nii format)"),
                    make_option(c("--autobalance"), action="store_true", type="logical", help="balance training samples by upsampling minority classes"),
                    make_option(c("-p", "--pthreads"), type="numeric", help="the number of parallel threads"),
                    make_option(c("-l", "--label_column"), type="character", help="the name of the column in the design file containing the training labels"),
                    make_option(c("-o", "--output"), type="character", help="the name of the output folder where results will be placed"),
                    make_option(c("-b", "--block_column"), type="character", help="the name of the column in the design file indicating the block variable used for cross-validation"),
                    make_option(c("-g", "--tune_grid"), type="character", help="string containing grid parameters in the following form: a=\\(1,2\\), b=\\('one', 'two'\\)"),
                    make_option(c("--tune_length"), type="numeric", help="an integer denoting the number of levels for each model tuning parameter"),
                    make_option(c("-i", "--niter"), type="character", help="number of randomized searchlight iterations"),
                    make_option(c("--savePredictors"), type="logical", action="store_true", help="save model fits (one per ROI) for predicting new data sets (default is FALSE)"),
                    make_option(c("--skipIfFolderExists"), type="logical", action="store_true", help="skip, if output folder already exists"),
                    make_option(c("--output_class_metrics"), type="logical", help="write out performance metrics for each class in multiclass settings"),
                    make_option(c("--ensemble_predictor"), type="logical", help="predictor is based on average of all cross-validated runs"),
                    make_option(c("--bootstrap_replications"), type="numeric", help="number of bootstrapped models to be computed for each crossvalidation fold"),
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
}

if (!is.null(config$roi_grouping)) {
  roilist <- lapply(1:length(config$roi_grouping), function(i) {
    grp <- config$roi_grouping[[i]]
    idx <- which(config$ROIVolume %in% grp)
    vol <- makeVolume(refvol=config$ROIVolume)
    vol[idx] <- i
    vol
  })
  
  config$ROIVolume <- roilist
} else {
  config$ROIVolume <- list(config$ROIVolume)
}

parcellationVolume <- if (!is.null(config$parcellation)) {
  loadVolume(config$parcellation)
}


config$maskVolume <- as(Reduce("+", lapply(config$ROIVolume, function(roivol) as(roivol, "LogicalBrainVolume"))), "LogicalBrainVolume")

config <- initializeData(config)

flog.info("number of training trials: %s", length(rowIndices))
flog.info("max trial index: %s", max(rowIndices))
flog.info("loading training data: %s", config$train_data)
flog.info("mask contains %s voxels", sum(config$maskVolume))

for (i in seq_along(config$ROIVolume)) {
  rvol <- config$ROIVolume[[i]]
  flog.info("Region mask contains: %s ROIs", length(unique(rvol[rvol > 0])))
}


flog.info("Running regional MVPA with parameters:", configParams, capture=TRUE)

if (length(config$labels) != dim(config$train_datavec)[4]) {
  flog.error("Number of volumes: %s must equal number of labels: %s", dim(config$train_datavec)[4], length(config$labels))
  stop()
}

featureSelector <- if (!is.null(config$feature_selector)) {
  FeatureSelector(config$feature_selector$method, config$feature_selector$cutoff_type,config$feature_selector$cutoff_value)
}

flog.info("feature selector: ", featureSelector, capture=TRUE)
flog.info("bootstrap replications: ", config$bootstrap_replications, capture=TRUE)

dataset <- MVPADataset(config$train_datavec, config$labels, config$maskVolume, config$block, config$test_datavec, config$testLabels, modelName=config$model, tuneGrid=config$tune_grid,
                       tuneLength=config$tune_length, testSplitVar=config$testSplitVar, testSplits=config$testSplits)



for (varname in c("test_subset", "train_subset", "roi_subset", "split_by")) {
  if (!is.null(configParams[[varname]]) && is(configParams[[varname]], "formula")) {
    configParams[[varname]] <- Reduce(paste, deparse(configParams[[varname]]))
  }
}

for (roinum in seq_along(config$ROIVolume)) {
  
  gc()
  
  if (config$skipIfFolderExists) {
    outdir <- paste0(config$output, "_roigroup_", roinum)
    if (file.exists(outdir)) {
      flog.info("output folder %s already exists, skipping.", outdir)
      next
    }
  }
  
  if (length(config$ROIVolume) > 1) {
    outdir <- paste0(config$output, "_roigroup_", roinum)
    outdir <- makeOutputDir(outdir)
  } else {
    outdir <- makeOutputDir(config$output)
  }
  
  
  roivol <- config$ROIVolume[[roinum]]
  
  
  mvpa_res <- mvpa_regional(dataset, roivol, config$pthreads, config$savePredictors, 
                            autobalance=config$autobalance, 
                            bootstrapReplications=config$bootstrap_replications, 
                            featureSelector=featureSelector, 
                            featureParcellation=parcellationVolume, 
                            classMetrics=config$output_class_metrics,
                            ensemblePredictor=config$ensemble_predictor)
  
  lapply(1:length(mvpa_res$outVols), function(i) {
    out <- paste0(outdir, "/", names(mvpa_res$outVols)[i], ".nii")
    writeVolume(mvpa_res$outVols[[i]], out)  
  })

  write.table(format(mvpa_res$performance,  digits=2, scientific=FALSE, drop0trailing=TRUE), paste0(paste0(outdir, "/performance_table.txt")), row.names=FALSE, quote=FALSE)
  saveResults(mvpa_res$extendedResults, outdir)

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


