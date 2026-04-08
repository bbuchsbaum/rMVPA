#!/usr/bin/env Rscript

if (!requireNamespace("optparse", quietly = TRUE)) {
  stop(
    "scripts/MVPA_Searchlight.R requires the 'optparse' package. Install it with install.packages('optparse').",
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(rMVPA)
})

`%||%` <- function(x, y) if (is.null(x)) y else x

option_list <- list(
  optparse::make_option(c("-r", "--radius"), type = "numeric",
              help = "Radius (mm) of the searchlight sphere."),
  optparse::make_option(c("-t", "--train_design"), type = "character",
              help = "File path for the training design table."),
  optparse::make_option("--test_design", type = "character",
              help = "File path for the test design table."),
  optparse::make_option(c("-s", "--type"), type = "character", default = "randomized",
              help = "Searchlight type: 'standard' or 'randomized'. [default: %default]"),
  optparse::make_option("--train_data", type = "character",
              help = "4D .nii file for training data."),
  optparse::make_option("--test_data", type = "character",
              help = "4D .nii file for test data."),
  optparse::make_option(c("-n", "--normalize"), action = "store_true", default = FALSE,
              help = "Center and scale each volume vector."),
  optparse::make_option(c("-m", "--model"), type = "character", default = "corsim",
              help = "Classifier/regression model name. [default: %default]"),
  optparse::make_option(c("-a", "--mask"), type = "character",
              help = "Binary image mask file (.nii)."),
  optparse::make_option(c("-p", "--ncores"), type = "numeric", default = 1,
              help = "Number of CPU cores to use in parallel."),
  optparse::make_option(c("--batch-size"), type = "numeric", default = NULL,
              help = "Number of searchlights per batch."),
  optparse::make_option(c("-l", "--label_column"), type = "character", default = "labels",
              help = "Name of the training label column in the design file."),
  optparse::make_option(c("--test_label_column"), type = "character",
              help = "Name of the test label column in the test design file."),
  optparse::make_option(c("-o", "--output"), type = "character", default = "searchlight_results",
              help = "Output folder to store results. [default: %default]"),
  optparse::make_option(c("-b", "--block_column"), type = "character",
              help = "Block (strata) column for cross-validation."),
  optparse::make_option(c("-g", "--tune_grid"), type = "character",
              help = "Parameter grid, e.g. 'alpha=c(0.1,0.5,1.0)'."),
  optparse::make_option(c("--tune_length"), type = "numeric",
              help = "Number of levels for each model tuning parameter."),
  optparse::make_option(c("-i", "--niter"), type = "numeric", default = 16,
              help = "Number of randomized searchlight iterations. [default: %default]"),
  optparse::make_option(c("--class_metrics"), action = "store_true", default = FALSE,
              help = "Write out per-class metrics in multiclass settings."),
  optparse::make_option(c("-c", "--config"), type = "character",
              help = "YAML or R config file specifying program parameters.")
)

parser <- optparse::OptionParser(usage = "MVPA_Searchlight.R [options]", option_list = option_list)
parsed <- optparse::parse_args(parser, positional_arguments = TRUE)
args <- parsed$options

cfg <- do.call(rMVPA::mvpa_config, c(list(mode = "searchlight"), args))
analysis <- rMVPA::build_analysis(cfg)
result <- rMVPA::run_analysis(analysis, preflight = "warn")

out_dir <- cfg$output %||% "searchlight_results"
rMVPA::save_results(result, dir = out_dir, level = "complete")

message("Searchlight analysis complete. Results written to: ", out_dir)
