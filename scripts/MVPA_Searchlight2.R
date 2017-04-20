#! /usr/bin/env Rscript


.suppress <- suppressPackageStartupMessages
.suppress(library(neuroim))
.suppress(library(neurosurf))
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
                    #make_option(c("--autobalance"), action="store_true", type="logical", help="balance training samples by upsampling minority classes"),
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
                    make_option(c("--class_metrics"), type="character", help="write out performance metrics for each class in multiclass settings"),
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

#flog.info("loading training data: %s", config$train_data)

if (config$data_mode == "image") {
  mask_volume <- as(load_mask(config), "LogicalBrainVolume")
  dataset <- rMVPA:::initialize_image_data(config, mask)
  dataset <- list(dataset)
  names(dataset) <- ""
  flog.info("image mask contains %s voxels", sum(config$mask_volume))
} else if (config$data_mode == "surface") {
  dataset <- rMVPA:::initialize_surface_data(config)
} else {
  flog.error("unrecognized data_mode: %s", config$data_mode)
}

row_indices <- which(config$train_subset)

flog.info("number of trials: %s", length(row_indices))
flog.info("max trial index: %s", max(row_indices))


flog.info("Running searchlight with parameters:", config_params, capture=TRUE)

if (!is.null(config$custom_performance)) {
  flog.info("custom performance function provided: ", config$custom_performance)
}

flog.info("initializing design structure")

design <- mvpa_design(train_design=config$train_design,
                      y_train=config$label_column,
                      test_design=config$test_design,
                      y_test=config$test_label_column,
                      block_var=config$block_column,
                      split_by=config$split_by)

crossval <- rMVPA:::initialize_crossval(config, design)
feature_selector <- rMVPA:::initialize_feature_selection(config)

if (is.numeric(design$y_train)) {
  flog.info("labels are continuous, running a regression analysis.")
} else {
  flog.info("labels are categorical, running a classification analysis.")
}


for (i in 1:length(dataset)) {
  dset <- dataset[[i]]
  if (names(dataset)[i] != "") {
    flog.info("running searchlight for %s dataset: ", names(dataset)[i])
  } else {
    flog.info("running searchlight")
  }
  
  
  mvpa_mod <- rMVPA:::load_mvpa_model(config, dset, design,crossval,feature_selector)
  searchres <- rMVPA:::run_searchlight(mvpa_mod, radius=config$radius, method=config$type, niter=config$niter)

  output <- makeOutputDir(paste0(names(dataset)[i], "_", config$output))
  
  for (i in 1:length(searchres)) {
    out <- paste0(output, "/", names(searchres)[i], ".nii")
    writeVolume(searchres[[i]], out)  
  }
  if (!is.null(configParams$test_subset)) {
    configParams$test_subset <- Reduce(paste, deparse(config_params$test_subset))
  }
  
  if (!is.null(configParams$train_subset)) {
    configParams$train_subset <- Reduce(paste, deparse(config_params$train_subset))
  }
  
  if (!is.null(configParams$split_by)) {
    config_params$split_by <- Reduce(paste, deparse(config_params$split_by))
  }
  
  configout <- paste0(config$output, "/config.yaml")
  qwrite(as.list(configParams), configout)
}



