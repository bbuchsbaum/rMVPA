###############################################################################
# MVPA Regional Analysis Script
# 
# This script performs regional multivariate pattern analysis (MVPA) on fMRI data.
# It supports both classification and regression analyses using volumetric (NIfTI) or
# surface-based neuroimaging data. Regional analysis can be performed on a specified 
# set of regions within the brain, and the results include performance metrics and
# prediction tables.
# 
# Key Features:
#   - Configurable via YAML or R config files
#   - Supports test and train subsets (using test_subset with a single design file)
#   - Outputs performance maps and prediction tables
#   - Consistent logging and error checking
###############################################################################

#! /usr/bin/env Rscript

.suppress <- suppressPackageStartupMessages
.suppress(library(neuroim2))
.suppress(library(neurosurf))
.suppress(library(rMVPA))
.suppress(library(optparse))
.suppress(library(futile.logger))
.suppress(library(io))
.suppress(library(doParallel))


###############################################################################
# Command-Line Options
###############################################################################
option_list <- list(
  make_option(c("-t", "--train_design"), type="character", 
              help="File name of the training design table"),
  make_option("--test_design", type="character", 
              help="File name of the test design table"),
  make_option("--train_data", type="character", 
              help="File name of the training data (4D NIfTI file)"),
  make_option("--test_data", type="character", 
              help="File name of the testing data (4D NIfTI file)"),
  make_option(c("-n", "--normalize"), action="store_true", default=FALSE,
              help="Center and scale each image volume"),
  make_option(c("-m", "--model"), type="character", 
              help="Name of the classifier/regression model"),
  make_option(c("-a", "--mask"), type="character", 
              help="Name of binary image mask file (.nii format)"),
  make_option(c("-p", "--pthreads"), type="numeric", 
              help="Number of parallel threads"),
  make_option(c("-l", "--label_column"), type="character", 
              help="Name of the column in the design file containing training labels"),
  make_option(c("--test_label_column"), type="character", 
              help="Name of the column in the test design file containing test labels"),
  make_option(c("-o", "--output"), type="character", 
              help="Output folder where results will be stored"),
  make_option(c("-b", "--block_column"), type="character", 
              help="Name of the column indicating the block variable for cross-validation"),
  make_option(c("-g", "--tune_grid"), type="character", 
              help="String containing grid parameters (e.g., a=\(1,2\), b=\('one','two'\))"),
  make_option(c("--tune_length"), type="numeric", 
              help="Number of levels for model tuning parameters"),
  make_option(c("--save_predictors"), type="logical", action="store_true", default=FALSE,
              help="Save model fits for predicting new data sets (default FALSE)"),
  make_option(c("--skip_if_folder_exists"), type="logical", action="store_true", default=FALSE,
              help="Skip analysis if output folder already exists"),
  make_option(c("--class_metrics"), type="logical", default=FALSE,
              help="Write out performance metrics for each class in multiclass settings"),
  make_option(c("--ensemble_predictor"), type="logical", default=FALSE,
              help="Use ensemble prediction based on average prediction of all cross-validated runs"),
  make_option(c("-c", "--config"), type="character",
              help="Configuration file (YAML or R) specifying program parameters")
)

oparser <- OptionParser(usage = "MVPA_Regional.R [options]", option_list=option_list)
opt <- parse_args(oparser, positional_arguments=TRUE)
args <- opt$options

flog.info("Command-line args:", args, capture=TRUE)

###############################################################################
# Configuration and Parameter Initialization
###############################################################################

# Load configuration from file if provided
config <- rMVPA:::initialize_configuration(args)
# Set default parameters specific to regional analysis
config <- rMVPA:::initialize_standard_parameters(config, args, "mvpa_regional")

# Regional-specific parameter: save_predictors, etc.
rMVPA:::set_arg("save_predictors", config, args, FALSE)

# Initialize tuning grid if provided
config$tune_grid <- rMVPA:::initialize_tune_grid(args, config)

# Convert configuration to a list for logging (if needed)
config_params <- as.list(config)

###############################################################################
# Design Initialization
###############################################################################

flog.info("Initializing design structure...")
config$design <- rMVPA:::initialize_design(config)
design <- config$design

###############################################################################
# Data Loading
###############################################################################

flog.info("Loading training data: %s", config$train_data)

if (config$data_mode == "image") {
  # Load the image mask and image data
  region_mask <- rMVPA:::load_mask(config)
  mask_volume <- as(region_mask, "LogicalNeuroVol")
  dataset <- rMVPA:::initialize_image_data(config, region_mask)
  # Store dataset as a list for consistency
  dataset <- list(dataset)
  names(dataset) <- ""
  flog.info("Image mask contains %s voxels", sum(mask_volume))
  flog.info("Image mask defines %s regions", length(table(region_mask)) - 1)
} else if (config$data_mode == "surface") {
  dataset <- rMVPA:::initialize_surface_data(config)
} else {
  flog.error("Unrecognized data_mode: %s", config$data_mode)
  stop()
}

# Log trial/subset information
rownumbers <- which(config$train_subset)
flog.info("Number of training trials: %s", length(rownumbers))
if (length(rownumbers) > 0) flog.info("Max trial index: %s", max(rownumbers))

###############################################################################
# Cross-Validation and Feature Selection Setup
###############################################################################

if (!is.null(config$custom_performance)) {
  flog.info("Custom performance function provided: %s", config$custom_performance)
}

crossval <- rMVPA:::initialize_crossval(config, design)
feature_selector <- rMVPA:::initialize_feature_selection(config)

if (is.numeric(design$y_train)) {
  flog.info("Labels are continuous - running regression analysis.")
} else {
  flog.info("Labels are categorical - running classification analysis.")
}

###############################################################################
# Model Loading and Regional Analysis Execution
###############################################################################

# Create the output directory
output <- rMVPA:::make_output_dir(config$output)

# Loop over datasets (typically one for regional analysis)
for (i in 1:length(dataset)) {
  dset <- dataset[[i]]
  dname <- names(dataset)[i]
  if (dname != "") {
    flog.info("Running mvpa_regional for dataset: %s", dname)
  } else {
    flog.info("Running mvpa_regional")
  }
  
  # Log number of regions based on mask values
  mvals <- as.vector(dset$mask)
  num_regions <- length(table(mvals[mvals != 0]))
  flog.info("Number of regions: %s", num_regions)
  
  # Load the MVPA model
  mvpa_mod <- rMVPA:::load_mvpa_model(config, dset, design, crossval, feature_selector)
  flog.info("MVPA model loaded.")
  
  # Run regional analysis using the regional function
  regional_res <- rMVPA:::run_regional(mvpa_mod, region_mask, return_fits=TRUE)
  
  # Log performance summary
  flog.info("Regional analysis performance:")
  print(regional_res$performance)
  
  # Write output volumes (e.g., importance or performance maps)
  write_output(regional_res$vol_results, name=dname, output, data_mode=config$data_mode)
  
  # Write out performance and prediction tables
  out_perf <- if (dname != "") {
    paste0(output, "/", dname, "_performance_table.txt")
  } else {
    paste0(output, "/performance_table.txt")
  }
  out_pred <- if (dname != "") {
    paste0(output, "/", dname, "_prediction_table.txt")
  } else {
    paste0(output, "/prediction_table.txt")
  }
  
  write.table(regional_res$performance, out_perf, row.names=FALSE, quote=FALSE, sep="\t")
  write.table(regional_res$prediction_table, out_pred, row.names=FALSE, quote=FALSE, sep="\t")
}

###############################################################################
# Completion Message
###############################################################################

flog.info("Regional MVPA analysis complete! Results saved to: %s", output)


