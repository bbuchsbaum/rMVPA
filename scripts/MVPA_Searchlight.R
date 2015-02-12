#! /usr/bin/env Rscript


.suppress <- suppressPackageStartupMessages
.suppress(library(neuroim))
.suppress(library(rMVPA))
.suppress(library(optparse))
.suppress(library(futile.logger))
.suppress(library(io))

option_list <- list(make_option(c("-r", "--radius"), type="numeric", help="the radius in millimeters of the searchlight"),
                    make_option(c("-t", "--train_design"), type="character", help="the file name of the design table"),
                    make_option("--test_design", type="character", help="the file name of the design table"),
                    make_option(c("-s", "--type"), type="character", help="the type of searchlight: standard or randomized"),  
                    make_option("--train_data", type="character", help="the name of the training data file as (4D .nii file)"), 
                    make_option("--test_data", type="character", help="the name of the testing data file as (4D .nii file)"),  
                    make_option(c("-n", "--normalize"), action="store_true", type="logical", help="center and scale each volume vector"),
                    make_option(c("-b", "--autobalance"), action="store_true", type="logical", help="balance training samples by upsampling minority classes"),
                    #make_option(c("-z", "--zscore"), action="store_true", type="logical", help="z-score all voxels/features within block"),
                    make_option(c("-m", "--model"), type="character", help="name of the classifier model"),
                    make_option(c("-a", "--mask"), type="character", help="name of binary image mask file (.nii format)"),
                    make_option(c("-p", "--pthreads"), type="numeric", help="the number of parallel threads"),
                    make_option(c("-l", "--label_column"), type="character", help="the name of the column in the design file containing the training labels"),
                    make_option(c("-o", "--output"), type="character", help="the name of the output folder where results will be placed"),
                    make_option(c("-b", "--block_column"), type="character", help="the name of the column in the design file indicating the block variable used for cross-validation"),
                    make_option(c("-g", "--tune_grid"), type="character", help="string containing grid parameters in the following form: a=\\(1,2\\), b=\\('one', 'two'\\)"),
                    make_option(c("-i", "--niter"), type="character", help="number of randomized searchlight iterations"),
                    make_option(c("-c", "--config"), type="character", help="name of configuration file used to specify program parameters"))
                  

oparser <- OptionParser(usage = "MVPA_Searchlight.R [options]", option_list=option_list)
opt <- parse_args(oparser, positional_arguments=TRUE)
args <- opt$options

flog.info("command line args are ", args, capture=TRUE)

config <- initializeConfiguration(args)

config <- initializeStandardParameters(config, args, "searchlight")

#flog.appender(appender.file(paste0(config$output, "/rMVPA.log")))
#flog.appender(appender.console(), name='my.logger')
#setDefault("autobalance", config, FALSE)
#setArg("tune_grid", config, args, NULL)
#setDefault("method_params", config, list())


## Searchlight Specific Params
setArg("niter", config, args, 16)
setArg("radius", config, args, 8)
setArg("type", config, args, "randomized")
## Searchlight Specific Params


config <- initializeTuneGrid(args, config)
configParams <- as.list(config)

config <- initializeDesign(config)
config$maskVolume <- as(loadMask(config), "LogicalBrainVolume")


rowIndices <- which(config$train_subset)
flog.info("number of trials: %s", length(rowIndices))
flog.info("max trial index: %s", max(rowIndices))
flog.info("loading training data: %s", config$train_data)
flog.info("mask contains %s voxels", sum(config$maskVolume))

config <- initializeData(config)

flog.info("Running searchlight with parameters:", configParams, capture=TRUE)


dataset <- MVPADataset(config$train_datavec, config$labels, config$maskVolume, config$block, config$test_datavec, config$testLabels, modelName=config$model, tuneGrid=config$tune_grid,
                       testSplitVar=config$testSplitVar, testSplits=config$testSplits)

for (lib in dataset$model$library) {
  library(lib, character.only = TRUE)
}

searchres <- mvpa_searchlight(dataset, config$radius,  config$type, config$niter, config$pthreads)

config$output <- makeOutputDir(config$output)

lapply(1:length(searchres), function(i) {
  out <- paste0(config$output, "/", names(searchres)[i], ".nii")
  writeVolume(searchres[[i]], out)  
})


if (!is.null(configParams$test_subset)) {
  configParams$test_subset <- Reduce(paste, deparse(configParams$test_subset))
}

if (!is.null(configParams$train_subset)) {
  configParams$train_subset <- Reduce(paste, deparse(configParams$train_subset))
}

if (!is.null(configParams$split_by)) {
  configParams$split_by <- Reduce(paste, deparse(configParams$split_by))
}

configout <- paste0(config$output, "/config.yaml")
qwrite(as.list(configParams), configout)

