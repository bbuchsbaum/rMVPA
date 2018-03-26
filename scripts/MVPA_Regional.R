#! /usr/bin/env Rscript

.suppress <- suppressPackageStartupMessages
.suppress(library(neuroim))
.suppress(library(neurosurf))
.suppress(library(rMVPA))
.suppress(library(optparse))
.suppress(library(caret))
.suppress(library(futile.logger))
.suppress(library(io))
.suppress(library(doParallel))




## should be able to use "test_subset" without test_data
## that means: train only on one subset and test only on another subset within single design

## if train_design only, the "test_subset" subsets the the train_design. Check to make sure no overlap in columns.
## if test_design provided, the test design must have matching test_data

## should output importance maps for methods that supply them (e.g. sda)



option_list <- list(
                    make_option(c("-t", "--train_design"), type="character", help="the file name of the training design table"),
                    make_option("--test_design", type="character", help="the file name of the test design table"),
                    make_option("--train_data", type="character", help="the name of the training data file as (4D .nii file)"), 
                    make_option("--test_data", type="character", help="the name of the testing data file as (4D .nii file)"),  
                    make_option(c("-n", "--normalize"), action="store_true", type="logical", help="center and scale each image volume"),
                    make_option(c("-m", "--model"), type="character", help="name of the classifier model"),
                    make_option(c("-a", "--mask"), type="character", help="name of binary image mask file (.nii format)"),
                    make_option(c("-p", "--pthreads"), type="numeric", help="the number of parallel threads"),
                    make_option(c("-l", "--label_column"), type="character", help="the name of the column in the design file containing the training labels"),
                    make_option(c("--test_label_column"), type="character", help="the name of the column in the test design file containing the test labels"),
                    make_option(c("-o", "--output"), type="character", help="the name of the output folder where results will be placed"),
                    make_option(c("-b", "--block_column"), type="character", help="the name of the column in the design file indicating the block variable used for cross-validation"),
                    make_option(c("-g", "--tune_grid"), type="character", help="string containing grid parameters in the following form: a=\\(1,2\\), b=\\('one', 'two'\\)"),
                    make_option(c("--tune_length"), type="numeric", help="an integer denoting the number of levels for each model tuning parameter"),
                    make_option(c("--save_predictors"), type="logical", action="store_true", help="save model fits (one per ROI) for predicting new data sets (default is FALSE)"),
                    make_option(c("--skip_if_folder_exists"), type="logical", action="store_true", help="skip, if output folder already exists"),
                    ## should be skip_if_folder_exists or "overwrite_folder"
                    make_option(c("--class_metrics"), type="logical", help="write out performance metrics for each class in multiclass settings"),
                    make_option(c("--ensemble_predictor"), type="logical", help="predictor is based on average prediction of all cross-validated runs"),
                    make_option(c("-c", "--config"), type="character", help="name of configuration file used to specify program parameters"))


oparser <- OptionParser(usage = "MVPA_Regional.R [options]", option_list=option_list)
opt <- parse_args(oparser, positional_arguments=TRUE)
args <- opt$options

flog.info("command line args are ", args, capture=TRUE)

## set up configuration 
config <- rMVPA:::initialize_configuration(args)
## set default parameters
config <- rMVPA:::initialize_standard_parameters(config, args, "mvpa_regional")


## Regional Specific Params
rMVPA:::set_arg("save_predictors", config, args, FALSE)


config$tune_grid <- rMVPA:::initialize_tune_grid(args, config)
config_params <- as.list(config)

flog.info("initializing design structure")
config$design <- rMVPA:::initialize_design(config)
design <- config$design


flog.info("loading training data: %s", config$train_data)



if (config$data_mode == "image") {
  region_mask <- rMVPA:::load_mask(config)
  mask_volume <- as(region_mask, "LogicalBrainVolume")
  dataset <- rMVPA:::initialize_image_data(config, region_mask)
  dataset <- list(dataset)
  names(dataset) <- ""
  flog.info("image mask contains %s voxels", sum(mask_volume))
  flog.info("image mask contains %s regions", length(table(region_mask)) - 1)
} else if (config$data_mode == "surface") {
  dataset <- rMVPA:::initialize_surface_data(config)
} else {
  flog.error("unrecognized data_mode: %s", config$data_mode)
}

print(dataset)
row_indices <- which(config$train_subset)

flog.info("number of trials: %s", length(row_indices))
flog.info("max trial index: %s", max(row_indices))


#config$ROIVolume <- loadMask(config)

# if (!is.null(config$roi_subset)) {
#   form <- try(eval(parse(text=config$roi_subset)))
#   if (inherits(form, "try-error")) {
#     flog.error("could not parse roi_subset parameter: %s", config$roi_subset)
#     stop()
#   }
#   
#   if (class(form) != "formula") {
#     flog.error("roi_subset argument must be an expression that starts with a ~ character")
#     stop()
#   }
#   
#   res <- as.logical(eval(form[[2]], list(x=config$ROIVolume)))
#   
#   config$ROIVolume[!res] <- 0
#   flog.info("roi_subset contains %s voxels", sum(config$ROIVolume > 0))
#   
#   if (sum(config$ROIVolume) <= 1) {
#     flog.info("ROI must contain more than one voxel, aborting.")
#     stop()
#   }
# }


#config$maskVolume <- as(Reduce("+", lapply(config$ROIVolume, function(roivol) as(roivol, "LogicalBrainVolume"))), "LogicalBrainVolume")
#config <- initializeData(config)

if (!is.null(config$custom_performance)) {
  flog.info("custom performance function provided: ", config$custom_performance)
}


crossval <- rMVPA:::initialize_crossval(config, design)
feature_selector <- rMVPA:::initialize_feature_selection(config)

if (is.numeric(design$y_train)) {
  flog.info("labels are continuous, running a regression analysis.")
} else {
  flog.info("labels are categorical, running a classification analysis.")
}

write_output <- function(res, name="", output, data_mode="image") {
  if (data_mode == "image") {
    for (i in 1:length(res)) {
      out <- paste0(output, "/", names(res)[i], "_", name, ".nii")
      writeVolume(res[[i]], out)  
    }
  } else if (data_mode == "surface") {
    for (i in 1:length(res)) {
      out <- paste0(output, "/", names(res)[i], "_", name)
      neurosurf::writeSurfaceData(res[[i]], out)  
    }
  } else {
    stop(paste("wrong data_mode:", data_mode))
  }
}

output <- rMVPA:::make_output_dir(config$output)

for (i in 1:length(dataset)) {
  dset <- dataset[[i]]
  dname <- names(dataset)[i]
  
  if (names(dataset)[i] != "") {
    flog.info("running mvpa_regional for %s dataset: ", names(dataset)[i])
  } else {
    flog.info("running mvpa_regional")
  }
  
  mvals <- as.vector(dset$mask)
  
  flog.info("number of regions is %s", length(table(mvals[mvals != 0])))
  
  mvpa_mod <- rMVPA:::load_mvpa_model(config, dset, design,crossval,feature_selector)
  
  flog.info("mvpa model: ", mvpa_mod, capture=TRUE)
  
  regional_res <- rMVPA:::run_regional(mvpa_mod, region_mask, return_fits=TRUE)
  
  print(regional_res$performance)
  
  write_output(regional_res$vol_results, name=names(dataset)[i], output, data_mode=config$data_mode)
  
  out_perf <-
    if (dname != "")
      paste0(paste0(output, paste0("/", dname, "_performance_table.txt")))
  else paste0(output, "/" , "performance_table.txt")
  
  out_pred <-
    if (dname != "")
      paste0(paste0(output, paste0("/", dname, "_prediction_table.txt")))
  else paste0(output, "/" , "prediction_table.txt")
  
  write.table(regional_res$performance, out_perf, row.names=FALSE, quote=FALSE)
  
  write.table(regional_res$prediction_table, out_pred, row.names=FALSE, quote=FALSE)
  
}
 
  


