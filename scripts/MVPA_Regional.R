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
config$output <- makeOutputDir(config$output)
config <- initializeStandardParameters(config, args, "mvpa_regional")

configParams <- as.list(config)

config <- initializeTuneGrid(args, config)
config <- initializeDesign(config)



rowIndices <- which(config$train_subset)
config$ROIVolume <- loadMask(config)
config$maskVolume <- as(config$ROIVolume, "LogicalBrainVolume")

config <- initializeData(config)

flog.info("number of trials: %s", length(rowIndices))
flog.info("max trial index: %s", max(rowIndices))
flog.info("loading training data: %s", config$train_data)
flog.info("mask contains %s voxels", sum(config$maskVolume))
flog.info("Region mask contains: %s ROIs", sum(unique(config$ROIVolume[config$ROIVolume > 0])))


flog.info("Running regional MVPA with parameters:", configParams, capture=TRUE)


dataset <- MVPADataset(config$train_datavec, config$labels, config$maskVolume, config$block, config$test_datavec, config$testLabels, modelName=config$model, tuneGrid=config$tune_grid)
mvpa_res <- mvpa_regional(dataset, config$ROIVolume, config$pthreads)


lapply(1:length(mvpa_res$outVols), function(i) {
  out <- paste0(config$output, "/", names(mvpa_res$outVols)[i], ".nii")
  writeVolume(mvpa_res$outVols[[i]], out)  
})

write.table(mvpa_res$performance, paste0(paste0(config$output, "/performance_table.txt")))

configout <- paste0(config$output, "/config.yaml")
qwrite(as.list(configParams), configout)


