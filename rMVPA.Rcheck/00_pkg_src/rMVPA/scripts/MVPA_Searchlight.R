#!/usr/bin/env Rscript
###############################################################################
# MVPA Searchlight Script
# 
# This script parses command-line arguments, loads an rMVPA configuration,
# and performs a searchlight analysis (standard or randomized) on either 
# volumetric or surface-based fMRI data.
#
# Major improvements:
#   - Organized into sections (Parsing, Configuration, Execution)
#   - Colorful console messages via crayon
#   - More robust error checking
#   - More readable and maintainable code
###############################################################################

# ----- Library Imports -----
suppressPackageStartupMessages({
  library(neuroim2)
  library(neurosurf)
  library(rMVPA)
  library(optparse)
  library(futile.logger)
  library(io)
  library(doParallel)
  library(purrr)
  library(future)
  library(stringr)
  library(crayon)
})

# ----- Define Command-Line Options -----
option_list <- list(
  make_option(c("-r", "--radius"), type="numeric", 
              help="Radius (mm) of the searchlight sphere."),
  make_option(c("-t", "--train_design"), type="character", 
              help="File path for the training design table."),
  make_option("--test_design", type="character", 
              help="File path for the test design table."),
  make_option(c("-s", "--type"), type="character", default="randomized",
              help="Searchlight type: 'standard' or 'randomized'. [default: %default]"),
  make_option("--train_data", type="character", 
              help="4D .nii file for training data."),
  make_option("--test_data", type="character", 
              help="4D .nii file for test data."),
  make_option(c("-n", "--normalize"), action="store_true", default=FALSE,
              help="Center and scale each volume vector."),
  make_option(c("-m", "--model"), type="character", default="corsim",
              help="Classifier/regression model name. [default: %default]"),
  make_option(c("-a", "--mask"), type="character",
              help="Binary image mask file (.nii)."),
  make_option(c("-p", "--ncores"), type="numeric", default=1,
              help="Number of CPU cores to use in parallel."),
  make_option(c("-l", "--label_column"), type="character", default="labels",
              help="Name of the training label column in the design file."),
  make_option(c("--test_label_column"), type="character",
              help="Name of the test label column in the test design file."),
  make_option(c("-o", "--output"), type="character", default="searchlight_results",
              help="Output folder to store results. [default: %default]"),
  make_option(c("-b", "--block_column"), type="character",
              help="Block (strata) column for cross-validation."),
  make_option(c("-g", "--tune_grid"), type="character",
              help="Parameter grid, e.g. a=\\(1,2\\), b=\\('one','two'\\)."),
  make_option(c("--tune_length"), type="numeric", 
              help="Number of levels for each model tuning parameter."),
  make_option(c("-i", "--niter"), type="character", default="16",
              help="Number of randomized searchlight iterations. [default: %default]"),
  make_option(c("--class_metrics"), type="character", default="TRUE",
              help="Write out performance metrics for each class in multiclass settings."),
  make_option(c("-c", "--config"), type="character",
              help="YAML or R config file specifying program parameters.")
)

# ----- Parse Command-Line Arguments -----
oparser <- OptionParser(usage="MVPA_Searchlight.R [options]", option_list=option_list)
opt     <- parse_args(oparser, positional_arguments=TRUE)
args    <- opt$options

# Increase allowed size for futures if large data
options(future.globals.maxSize = 4024 * 1024^2)

# ----- Pretty Print: Start-Up Info -----
cat(
  bold(cyan("========================================\n")),
  bold(cyan("         MVPA Searchlight Script        \n")),
  bold(cyan("========================================\n"))
)
cat(
  yellow("Command-line arguments:\n"), 
  green(paste(capture.output(str(args)), collapse = "\n")), 
  "\n\n"
)

flog.info("Parsed command-line arguments.")
flog.info("Initializing configuration and parameters...")

# ----- Configuration Initialization -----
# We assume these helper functions exist in your codebase (common.R).
config <- rMVPA:::initialize_configuration(args)
config <- rMVPA:::initialize_standard_parameters(config, args, "searchlight")

# Set specific parameters with default fallbacks
rMVPA:::set_arg("niter", config, args, 16)
rMVPA:::set_arg("radius", config, args, 8)
rMVPA:::set_arg("type", config, args, "randomized")
rMVPA:::set_arg("ncores", config, args, 1)

# Tuning parameters for classifier optimization
config$tune_grid <- rMVPA:::initialize_tune_grid(args, config)

# Build the design object (depending on block columns, etc.)
config$design <- rMVPA:::initialize_design(config)

# ----- Load Data -----
if (config$data_mode == "image") {
  mask_volume <- as(rMVPA:::load_mask(config), "LogicalNeuroVol")
  dataset     <- rMVPA:::initialize_image_data(config, mask_volume)
  dataset     <- list(dataset)
  names(dataset) <- ""
  flog.info("Image mask contains %s voxels", sum(mask_volume))
} else if (config$data_mode == "surface") {
  dataset <- rMVPA:::initialize_surface_data(config)
} else {
  stop(
    red("unrecognized data_mode: "), 
    yellow(config$data_mode)
  )
}

# Indices of included observations (training subset)
row_indices <- which(config$train_subset)
flog.info("Number of included trials: %s", length(row_indices))
flog.info("Max trial index: %s", max(row_indices))

if (!is.null(config$test_label_column)) {
  flog.info("Test label column: %s", config$test_label_column)
}

# ----- Summarize Configuration -----
cat(yellow("Running searchlight with parameters:\n"))
cat(green(paste(capture.output(str(as.list(config))), collapse = "\n")), "\n\n")

if (!is.null(config$custom_performance)) {
  flog.info("Custom performance function provided: %s", config$custom_performance)
}

flog.info("Initializing design structure.")
design    <- config$design
crossval  <- rMVPA:::initialize_crossval(config, design)
feat_sel  <- rMVPA:::initialize_feature_selection(config)

if (is.numeric(design$y_train)) {
  flog.info("Labels are continuous -> regression analysis.")
} else {
  flog.info("Labels are categorical -> classification analysis.")
}

# ----- Output Folder -----
output_folder <- rMVPA:::make_output_dir(config$output)
flog.info("Output will be saved to: %s", output_folder)

# ----- Parallel/Multicore Setup -----
if (as.numeric(config$ncores) > 1) {
  flog.info("Multi-threaded processing with %s cores.", config$ncores)
  future::plan(future::multicore, workers = as.numeric(config$ncores))
}

# ----- Main Searchlight Execution -----
write_output <- function(searchres, name="", output, data_mode="image") {
  # Writes results to disk, either volumetric (NIfTI) or surface data
  if (data_mode == "image") {
    for (i in seq_along(searchres)) {
      outfname <- if (nzchar(name)) {
        file.path(output, paste0(names(searchres)[i], "_", name, ".nii"))
      } else {
        file.path(output, paste0(names(searchres)[i], ".nii"))
      }
      resobj <- searchres[[i]]
      if (inherits(resobj, "NeuroVol")) {
        write_vol(resobj, outfname)  
      } else if (inherits(resobj, "NeuroVec")) {
        write_vec(resobj, outfname) 
      }
    }
  } else if (data_mode == "surface") {
    for (i in seq_along(searchres)) {
      outbase <- file.path(output, names(searchres)[i])
      neurosurf::write_surf_data(searchres[[i]], outbase, name)  
    }
  } else {
    stop("Unknown data_mode: ", data_mode)
  }
}

# We might have multiple dataset entries for surface data 
for (i in seq_along(dataset)) {
  dset_name <- names(dataset)[i]
  if (nzchar(dset_name)) {
    cat(cyan(sprintf("Running searchlight for dataset: %s\n", dset_name)))
  } else {
    cat(cyan("Running searchlight...\n"))
  }
  
  dset    <- dataset[[i]]
  mvpaMod <- rMVPA:::load_mvpa_model(config, dset, design, crossval, feat_sel)
  
  # Actual searchlight call
  searchres <- rMVPA:::run_searchlight(
    mvpaMod, 
    radius = config$radius, 
    method = config$type, 
    niter  = config$niter
  )
  
  write_output(searchres, name=dset_name, output=output_folder, data_mode=config$data_mode)
}

# ----- Final Output (Write Config) -----
# Convert config_params for YAML writing
config_params <- as.list(config)

# Clean up some possibly large formula expressions
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
  config_params$label_column <- Reduce(paste, deparse(config_params$label_column))
}

# Write configuration out to YAML for reproducibility
config_out <- file.path(output_folder, "config.yaml")
io::qwrite(config_params, config_out)

# ----- Completion Message -----
cat(
  bold(cyan("\nSearchlight complete! Results saved to:\n")), 
  green(output_folder), 
  "\n\n"
)

