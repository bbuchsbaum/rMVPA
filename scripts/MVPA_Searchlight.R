#! /usr/bin/env Rscript


.suppress <- suppressPackageStartupMessages
.suppress(library(neuroim))
.suppress(library(neurosurf))
.suppress(library(rMVPA))
.suppress(library(optparse))
.suppress(library(futile.logger))
.suppress(library(io))
.suppress(library(doParallel))
.suppress(library(purrr))

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

## TODO check that 'test_label_column' is valid, if test_data is present

## set up configuration 
config <- rMVPA:::initialize_configuration(args)

## set default parameters
config <- rMVPA:::initialize_standard_parameters(config, args, "searchlight")

## set Searchlight specific params
rMVPA:::set_arg("niter", config, args, 16)
rMVPA:::set_arg("radius", config, args, 8)
rMVPA:::set_arg("type", config, args, "randomized")


## tuning parameters for classiifer optimization
config$tune_grid <- rMVPA:::initialize_tune_grid(args, config)

config_params <- as.list(config)

config$design <- rMVPA:::initialize_design(config)

#flog.info("loading training data: %s", config$train_data)
if (config$data_mode == "image") {
  mask_volume <- as(rMVPA:::load_mask(config), "LogicalBrainVolume")
  dataset <- rMVPA:::initialize_image_data(config, mask_volume)
  dataset <- list(dataset)
  names(dataset) <- ""
  flog.info("image mask contains %s voxels", sum(mask_volume))
} else if (config$data_mode == "surface") {
  dataset <- rMVPA:::initialize_surface_data(config)
} else {
  flog.error("unrecognized data_mode: %s", config$data_mode)
}

## the indices of the included observations
row_indices <- which(config$train_subset)

flog.info("number of included trials: %s", length(row_indices))
flog.info("max trial index: %s", max(row_indices))

if (! is.null(config$test_label_column)) {
  flog.info("test_label_column: %s", config$test_label_column)
  #flog.info("test_design has %s rows", nrow(config$test_design))
}


flog.info("Running searchlight with parameters:", config_params, capture=TRUE)

if (!is.null(config$custom_performance)) {
  flog.info("custom performance function provided: ", config$custom_performance)
}

flog.info("initializing design structure")

design <- config$design

crossval <- rMVPA:::initialize_crossval(config, design)
feature_selector <- rMVPA:::initialize_feature_selection(config)

if (is.numeric(design$y_train)) {
  flog.info("labels are continuous, running a regression analysis.")
} else {
  flog.info("labels are categorical, running a classification analysis.")
}


write_output <- function(searchres, name="", output, data_mode="image") {
  if (data_mode == "image") {
    for (i in 1:length(searchres)) {
      oname <- if (name != "") paste0(output, "/", names(searchres)[i], "_", name, ".nii") else paste0(output, "/", names(searchres)[i], ".nii")
      writeVolume(searchres[[i]], oname)  
    }
  } else if (data_mode == "surface") {
    for (i in 1:length(searchres)) {
      out <- paste0(output, "/", names(searchres)[i])
      neurosurf::writeSurfaceData(searchres[[i]], out, name)  
    }
  } else {
    stop(paste("wrong data_mode:", data_mode))
  }
}


output <- rMVPA:::make_output_dir(config$output)

for (i in 1:length(dataset)) {
  dset <- dataset[[i]]
  if (names(dataset)[i] != "") {
    flog.info("running searchlight for %s dataset: ", names(dataset)[i])
  } else {
    flog.info("running searchlight")
  }

  mvpa_mod <- rMVPA:::load_mvpa_model(config, dset, design,crossval,feature_selector)
  searchres <- rMVPA:::run_searchlight(mvpa_mod, radius=config$radius, method=config$type, niter=config$niter)
  write_output(searchres, name=names(dataset)[i], output, data_mode=config$data_mode)
}

if (!is.null(config_params$test_subset)) {
  config_params$test_subset <- Reduce(paste, deparse(config_params$test_subset))
}

if (!is.null(config_params$train_subset)) {
  config_params$train_subset <- Reduce(paste, deparse(config_params$train_subset))
}

if (!is.null(config_params$split_by)) {
  config_params$split_by <- Reduce(paste, deparse(config_params$split_by))
}

if (purrr::is_formula(config_params$label_column)) {
  config_params$label_column= Reduce(paste, deparse(config_params$label_column))
}

configout <- paste0(config$output, "/config.yaml")
qwrite(as.list(config_params), configout)



