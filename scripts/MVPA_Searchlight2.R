#! /usr/bin/env Rscript


.suppress <- suppressPackageStartupMessages
.suppress(library(neuroim))
.suppress(library(rMVPA))
.suppress(library(optparse))
.suppress(library(futile.logger))
.suppress(library(io))
.suppress(library(doParallel))


option_list <- list(make_option(c("-r", "--radius"), type="numeric", help="the radius in millimeters of the searchlight"),
                    make_option(c("-t", "--train_design"), type="character", help="the file name of the design table"),
                    make_option("--test_design", type="character", help="the file name of the design table"),
                    make_option(c("-s", "--type"), type="character", help="the type of searchlight: standard or randomized"),  
                    make_option("--train_data", type="character", help="the name of the training data file as (4D .nii file)"), 
                    make_option("--test_data", type="character", help="the name of the testing data file as (4D .nii file)"),  
                    make_option(c("-n", "--normalize"), action="store_true", type="logical", help="center and scale each volume vector"),
                    make_option(c("--autobalance"), action="store_true", type="logical", help="balance training samples by upsampling minority classes"),
                    #make_option(c("-z", "--zscore"), action="store_true", type="logical", help="z-score all voxels/features within block"),
                    make_option(c("-m", "--model"), type="character", help="name of the classifier model"),
                    make_option(c("-a", "--mask"), type="character", help="name of binary image mask file (.nii format)"),
                    make_option(c("-p", "--pthreads"), type="numeric", help="the number of parallel threads"),
                    make_option(c("-l", "--label_column"), type="character", help="the name of the column in the design file containing the training labels"),
                    make_option(c("--test_label_column"), type="character", help="the name of the column in the test design file containing the test labels"),
                    make_option(c("-o", "--output"), type="character", help="the name of the output folder where results will be placed"),
                    make_option(c("-b", "--block_column"), type="character", help="the name of the column in the design file indicating the block variable used for cross-validation"),
                    make_option(c("-g", "--tune_grid"), type="character", help="string containing grid parameters in the following form: a=\\(1,2\\), b=\\('one', 'two'\\)"),
                    make_option(c("--tune_length"), type="numeric", help="an integer denoting the number of levels for each model tuning parameter"),
                    make_option(c("-i", "--niter"), type="character", help="number of randomized searchlight iterations"),
                    make_option(c("--output_class_metrics"), type="character", help="write out performance metrics for each class in multiclass settings"),
                    make_option(c("-c", "--config"), type="character", help="name of configuration file used to specify program parameters"))
                  

oparser <- OptionParser(usage = "MVPA_Searchlight.R [options]", option_list=option_list)
opt <- parse_args(oparser, positional_arguments=TRUE)
args <- opt$options

flog.info("command line args are ", args, capture=TRUE)


## set up configuration 
config <- rMVPA:::initialize_configuration(args)

## set default parameters
config <- rMVPA:::initialize_standard_parameters(config, args, "searchlight")

## set Searchlight specific params
rMVPA:::set_arg("niter", config, args, 16)
rMVPA:::set_arg("radius", config, args, 8)
rMVPA:::set_arg("type", config, args, "randomized")


config$tune_grid <- rMVPA:::initialize_tune_grid(args, config)
config_params <- as.list(config)

config$design <- rMVPA:::initialize_design(config)


if (config$data_mode == "image")
  config$mask_volume <- as(load_mask(config), "LogicalBrainVolume")
  flog.info("image mask contains %s voxels", sum(config$mask_volume))
}

row_indices <- which(config$train_subset)

flog.info("number of trials: %s", length(row_indices))
flog.info("max trial index: %s", max(row_indices))
flog.info("loading training data: %s", config$train_data)




config <- initialize_data(config)

flog.info("Running searchlight with parameters:", configParams, capture=TRUE)

if (is.numeric(config$labels)) {
  flog.info("labels are continuous, running a regression analysis.")
} else {
  flog.info("labels are categorical, running a classification analysis.")
  config$labels <- as.factor(config$labels)
}

if (!is.null(config$custom_performance)) {
  flog.info("custom performance function provided: ", config$custom_performance)
}

design <- mvpa_design(train_design=config$train_design,
                      y_train=config$labels,
                      test_design=config$test_design,
                      y_test=config$testLabels,
                      block_var=config$block,
                      split_by=NULL)

dataset <- mvpa_dataset(config$train_datavec,
                        config$test_datavec,
                        mask=config$maskVolume)
                        
model <- load_model(config$model)


if (config$type %in% c("standard", "randomized", "randomized2")) {
  flog.info("searchlight type: ", config$type)
  searchres <- mvpa_searchlight(dataset, model, crossVal, config$radius,  
                              method=config$type, niter=config$niter, 
                              classMetrics=config$output_class_metrics)
} else if (config$type == "clustered") {
  flog.info("clustered searchlight with nclusters: ", config$nclusters)
  searchres <- mvpa_clustered_searchlight(dataset, model, crossVal, nclusters=config$nclusters,
                                classMetrics=config$output_class_metrics)
  
} else {
  stop(paste("unrecognized searchlight type: ", config$type))
}

config$output <- makeOutputDir(config$output)

for (i in 1:length(searchres)) {
  out <- paste0(config$output, "/", names(searchres)[i], ".nii")
  writeVolume(searchres[[i]], out)  
}

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

