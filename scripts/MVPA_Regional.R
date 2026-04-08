#!/usr/bin/env Rscript

if (!requireNamespace("optparse", quietly = TRUE)) {
  stop(
    "scripts/MVPA_Regional.R requires the 'optparse' package. Install it with install.packages('optparse').",
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(rMVPA)
})

`%||%` <- function(x, y) if (is.null(x)) y else x

option_list <- list(
  optparse::make_option(c("-t", "--train_design"), type = "character",
              help = "File name of the training design table"),
  optparse::make_option("--test_design", type = "character",
              help = "File name of the test design table"),
  optparse::make_option("--train_data", type = "character",
              help = "File name of the training data (4D NIfTI file)"),
  optparse::make_option("--test_data", type = "character",
              help = "File name of the testing data (4D NIfTI file)"),
  optparse::make_option(c("-n", "--normalize"), action = "store_true", default = FALSE,
              help = "Center and scale each image volume"),
  optparse::make_option(c("-m", "--model"), type = "character",
              help = "Name of the classifier/regression model"),
  optparse::make_option(c("-a", "--mask"), type = "character",
              help = "Name of binary image mask file (.nii format)"),
  optparse::make_option(c("-p", "--pthreads"), type = "numeric",
              help = "Number of parallel threads"),
  optparse::make_option(c("-l", "--label_column"), type = "character",
              help = "Name of the column in the design file containing training labels"),
  optparse::make_option(c("--test_label_column"), type = "character",
              help = "Name of the column in the test design file containing test labels"),
  optparse::make_option(c("-o", "--output"), type = "character",
              help = "Output folder where results will be stored"),
  optparse::make_option(c("-b", "--block_column"), type = "character",
              help = "Name of the column indicating the block variable for cross-validation"),
  optparse::make_option(c("-g", "--tune_grid"), type = "character",
              help = "String containing grid parameters"),
  optparse::make_option(c("--tune_length"), type = "numeric",
              help = "Number of levels for model tuning parameters"),
  optparse::make_option(c("--save_predictors"), action = "store_true", default = FALSE,
              help = "Save model fits for predicting new data sets"),
  optparse::make_option(c("--skip_if_folder_exists"), action = "store_true", default = FALSE,
              help = "Skip analysis if output folder already exists"),
  optparse::make_option(c("--class_metrics"), type = "logical", default = FALSE,
              help = "Write out performance metrics for each class in multiclass settings"),
  optparse::make_option(c("--ensemble_predictor"), type = "logical", default = FALSE,
              help = "Use ensemble prediction based on average prediction of all cross-validated runs"),
  optparse::make_option(c("-c", "--config"), type = "character",
              help = "Configuration file (YAML or R) specifying program parameters")
)

parser <- optparse::OptionParser(usage = "MVPA_Regional.R [options]", option_list = option_list)
parsed <- optparse::parse_args(parser, positional_arguments = TRUE)
args <- parsed$options

cfg <- do.call(rMVPA::mvpa_config, c(list(mode = "regional"), args))
analysis <- rMVPA::build_analysis(cfg)
result <- rMVPA::run_analysis(analysis, preflight = "warn")

out_dir <- cfg$output %||% "mvpa_regional_analysis"
rMVPA::save_results(result, dir = out_dir, level = "complete")

message("Regional analysis complete. Results written to: ", out_dir)
