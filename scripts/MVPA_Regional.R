#! /usr/bin/env Rscript

.suppress <- suppressPackageStartupMessages
.suppress(library(neuroim))
.suppress(library(rMVPA))
.suppress(library(optparse))
.suppress(library(caret))
.suppress(library(futile.logger))
.suppress(library(io))

option_list <- list(
                    make_option(c("-t", "--train_design"), type="character", help="the file name of the design table"),
                    make_option("--test_design", type="character", help="the file name of the design table"),
                    make_option(c("-s", "--type"), type="character", help="the type of searchlight: standard or randomized"),  
                    make_option("--train_data", type="character", help="the name of the training data file as (4D .nii file)"), 
                    make_option("--test_data", type="character", help="the name of the testing data file as (4D .nii file)"),  
                    make_option(c("-n", "--normalize"), action="store_true", type="logical", help="center and scale each volume vector"),
                    make_option(c("-m", "--model"), type="character", help="name of the classifier model"),
                    make_option(c("-a", "--mask"), type="character", help="name of binary image mask file (.nii format)"),
                    make_option(c("-p", "--pthreads"), type="numeric", help="the number of parallel threads"),
                    make_option(c("-l", "--label_column"), type="character", help="the name of the column in the design file containing the training labels"),
                    make_option(c("-o", "--output"), type="character", help="the name of the output folder where results will be placed"),
                    make_option(c("-b", "--block_column"), type="character", help="the name of the column in the design file indicating the block variable used for cross-validation"),
                    make_option(c("-g", "--tune_grid"), type="character", help="string containing grid parameters in the following form: a=\\(1,2\\), b=\\('one', 'two'\\)"),
                    make_option(c("-i", "--niter"), type="character", help="number of randomized searchlight iterations"),
                    make_option(c("-c", "--config"), type="character", help="name of configuration file used to specify program parameters"))


oparser <- OptionParser(usage = "MVPA_Regional.R [options]", option_list=option_list)
opt <- parse_args(oparser, positional_arguments=TRUE)
args <- opt$options

flog.info("command line args are ", args, capture=TRUE)

config <- initializeConfiguration(args)

config <- initializeStandardParameters(config, args, "mvpa_regional")
config <- initializeTuneGrid(args, config)

configParams <- as.list(config)


config <- initializeDesign(config)

rowIndices <- which(config$train_subset)
config$ROIVolume <- loadMask(config)

if (! is.null(config$roi_subset)) {
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
  
  flog.info("roi_subset contains %s voxels", sum(res))
  config$ROIVolume[!res] <- 0
}


config$maskVolume <- as(config$ROIVolume, "LogicalBrainVolume")

config <- initializeData(config)

flog.info("number of trials: %s", length(rowIndices))
flog.info("max trial index: %s", max(rowIndices))
flog.info("loading training data: %s", config$train_data)
flog.info("mask contains %s voxels", sum(config$maskVolume))
flog.info("Region mask contains: %s ROIs", length(unique(config$ROIVolume[config$ROIVolume > 0])))


flog.info("Running regional MVPA with parameters:", configParams, capture=TRUE)

if (length(config$labels) != dim(config$train_datavec)[4]) {
  flog.error("Number of volumes: %s must equal number of labels: %s", dim(config$train_datavec)[4], length(config$labels))
  stop()
}

dataset <- MVPADataset(config$train_datavec, config$labels, config$maskVolume, config$block, config$test_datavec, config$testLabels, modelName=config$model, tuneGrid=config$tune_grid)
mvpa_res <- mvpa_regional(dataset, config$ROIVolume, config$pthreads)

config$output <- makeOutputDir(config$output)

lapply(1:length(mvpa_res$outVols), function(i) {
  out <- paste0(config$output, "/", names(mvpa_res$outVols)[i], ".nii")
  writeVolume(mvpa_res$outVols[[i]], out)  
})

write.table(format(mvpa_res$performance,  digits=2, scientific=FALSE, drop0trailing=TRUE), paste0(paste0(config$output, "/performance_table.txt")), row.names=FALSE, quote=FALSE)
write.table(format(mvpa_res$predMat,  digits=2, scientific=FALSE, drop0trailing=TRUE), paste0(paste0(config$output, "/predMat_table.txt")), row.names=FALSE, quote=FALSE)

saveRDS(mvpa_res$predictorList, paste0(config$output, "/predictorList.RDS"))


if (!is.null(configParams$test_subset)) {
  configParams$test_subset <- Reduce(paste, deparse(configParams$test_subset))
}

if (!is.null(configParams$train_subset)) {
  configParams$train_subset <- Reduce(paste, deparse(configParams$train_subset))
}
        
configout <- paste0(config$output, "/config.yaml")
qwrite(configParams, configout)


