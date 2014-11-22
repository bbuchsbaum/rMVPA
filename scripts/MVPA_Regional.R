#! /usr/bin/env Rscript

.suppress <- suppressPackageStartupMessages
.suppress(library(neuroim))
.suppress(library(rMVPA))
.suppress(library(optparse))
.suppress(library(futile.logger))
.suppress(library(io))

option_list <- list(make_option(c("-t", "--table"), type="character", help="the file name of the design table"),
                    make_option(c("-d", "--train_data"), type="character", help="the name of the training data file as (4D .nii file)"),  
                    make_option(c("-m", "--model"), type="character", help="name of the classifier model"),
                    make_option(c("-a", "--mask"), type="character", help="name of binary image mask file (.nii format)"),
                    make_option(c("-n", "--normalize"), type="character", help="center and scale each volume vector"),
                    make_option(c("-p", "--pthreads"), type="numeric", help="the number of parallel threads"),
                    make_option(c("-l", "--label_column"), type="character", help="the name of the column in the design file containing the training labels"),
                    make_option(c("-o", "--output"), type="character", help="the name of the output folder where results will be placed"),
                    make_option(c("-b", "--block_column"), type="character", help="the name of the column in the design file indicating the block variable used for cross-validation"),
                    make_option(c("-g", "--tune_grid"), type="character", help="string containing grid parameters in the following form: a=\\(1,2\\), b=\\('one', 'two'\\)"),
                    make_option(c("-c", "--config"), type="character", help="name of configuration file used to specify program parameters"))


oparser <- OptionParser(usage = "MVPA_Regional.R [options]", option_list=option_list)
opt <- parse_args(oparser, positional_arguments=TRUE)
args <- opt$options

config <- new.env()
config$output <- makeOutputDir(config$output)
flog.appender(appender.file(paste0(config$output, "/rMVPA.log")))



if (!is.null(args$config)) {
  if (! file.exists(args$config)) {
    flog.error("cannot find configuration file: %s", args$config)
    stop()
  } else {
    source(args$config, config)
    flog.info("found configuration file with parameters: %s", str(as.list(config)))
  }
}

setArg("table", config, args, "mvpa_design.txt")
setArg("model", config, args, "corsim")
setArg("pthreads", config, args, 1)
setArg("label_column", config, args, "labels")
setArg("output", config, args, paste0("regional_", config$labelColumn))
setArg("block_column", config, args, "block")
setArg("normalize", config, args, FALSE)
#setDefault("autobalance", config, FALSE)
setArg("tune_grid", config, args, NULL)
#setDefault("method_params", config, list())
setArg("mask", config, args, NULL)


if (!is.null(args$tune_grid)) {
  params <- try(expand.grid(eval(parse(text=args$tune_grid))))
  if (inherits(params, "try-error")) {
    stop("could not parse tune_grid expresson: ", args$tune_grid)
  }
  flog.info("tuning grid is", params, capture=TRUE)
  config$tune_grid <- params
}

flog.info("Running regional mvpa with parameters: %s", str(as.list(config)))

configParams <- as.list(config)


config$full_design <- read.table(config$table, header=TRUE, comment.char=";")
config$train_subset <- loadSubset(config$full_design, config)
config$train_design <- config$full_design[config$train_subset,]
config$labels <- loadLabels(config$train_design, config)
config$blockVar <- loadBlockColumn(config, config$train_design)



config$ROIVolume <- loadMask(config)
config$maskVolume <- LogicalBrainVolume(as.logical(config$ROIVolume > 0), space(config$ROIVolume))
config$train_datavec <- loadBrainData(config, indices=which(config$train_subset))

if (config$normalize) {
  flog.info("Normalizing: entering and scaling each volume of training data")
  norm_datavec <- do.call(cbind, eachVolume(config$train_datavec, function(x) scale(x), mask=config$maskVolume))
  config$train_datavec <- SparseBrainVector(norm_datavec, space(config$train_datavec), mask=config$maskVolume)
}

flog.info(paste("subset contains", nrow(config$train_design), "of", nrow(config$full_design), "rows."))

caret_model <- loadModel(config$model)
library(caret_model$library, character.only=TRUE)

#dataset <- MVPADataset(config$train_datavec, config$labels, config$maskVolume, config$blockVar)  
mvpa_res <- mvpa_regional(config$train_datavec, config$labels, config$ROIVolume, config$blockVar, config$method, ncores=config$ncores, tuneGrid=config$tune_grid)


lapply(1:length(mvpa_res$outVols), function(i) {
  out <- paste0(config$output, "/", names(mvpa_res$outVols)[i], ".nii")
  writeVolume(mvpa_res$outVols[[i]], out)  
})

write.table(mvpa_res$performance, paste0(paste0(config$output, "/performance_table.txt")))



