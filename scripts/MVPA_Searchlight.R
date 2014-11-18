#! /usr/bin/env Rscript


.suppress <- suppressPackageStartupMessages
.suppress(library(neuroim))
.suppress(library(rMVPA))
.suppress(library(optparse))
.suppress(library(futile.logger))
.suppress(library(io))

option_list <- list(make_option(c("-r", "--radius"), type="numeric", help="the radius in millimeters of the searchlight"),
                    make_option(c("-t", "--table"), type="character", help="the file name of the design table"),
                    make_option(c("-s", "--type"), type="character", help="the type of searchlight: standard or randomized"),  
                    make_option(c("-d", "--train_data"), type="character", help="the name of the training data file as (4D .nii file)"),  
                    make_option(c("-m", "--model"), type="character", help="name of the classifier model"),
                    make_option(c("-a", "--mask"), type="character", help="name of binary image mask file (.nii format)"),
                    make_option(c("-p", "--pthreads"), type="numeric", help="the number of parallel threads"),
                    make_option(c("-l", "--label_column"), type="character", help="the name of the column in the design file containing the training labels"),
                    make_option(c("-o", "--output"), type="character", help="the name of the output folder where results will be placed"),
                    make_option(c("-b", "--block_column"), type="character", help="the name of the column in the design file indicating the block variable used for cross-validation"),
                    make_option(c("-g", "--tune_grid"), type="character", help="string containing grid parameters in the following form: a=\\(1,2\\), b=\\('one', 'two'\\)"),
                    make_option(c("-i", "--niter"), type="character", help="number of randomized searchlight iterations"),
                    make_option(c("-c", "--config"), type="character", help="name of configuration file used to specify program parameters"))
                  

oparser <- OptionParser(usage = "MVPA_Searchlight.R [options]", option_list=option_list)
opt <- parse_args(oparser, positional_arguments=TRUE)
args <- opt$options

config <- new.env()

if (!is.null(args$config)) {
  if (! file.exists(args$config)) {
    flog.error("cannot find configuration file: %s", args$config)
    stop()
  } else {
    source(args$config, config)
    flog.info("found configuration file with parameters: %s", str(as.list(config)))
  }
}

config$output <- makeOutputDir(config$output)

flog.appender(appender.file(paste0(config$output, "/rMVPA.log")))

flog.info("args are ", str(args), capture=TRUE)

setArg("radius", config, args, 8)
setArg("table", config, args, "mvpa_design.txt")
setArg("type", config, args, "randomized")
setArg("model", config, args, "corsim")
setArg("pthreads", config, args, 1)
setArg("label_column", config, args, "labels")
setArg("output", config, args, paste0("searchlight_", config$labelColumn))
setArg("block_column", config, args, "block")
#setDefault("autobalance", config, FALSE)
setArg("tune_grid", config, args, NULL)
#setDefault("method_params", config, list())
setArg("niter", config, args, 4)
setArg("mask", config, args, NULL)

configParams <- as.list(config)


if (!is.null(args$tune_grid) && !args$tune_grid == "NULL") {
  params <- try(expand.grid(eval(parse(text=args$tune_grid))))
  if (inherits(params, "try-error")) {
    stop("could not parse tune_grid expresson: ", args$tune_grid)
  }
  flog.info("tuning grid is", params, capture=TRUE)
  config$tune_grid <- params
}

flog.info("Running searchlight with parameters: %s", str(as.list(config)))


config$output <- makeOutputDir(config$output)

#config$logFile <- file(paste0(config$output, "/searchlight.log"), "w")


#config$maskVolume[,,c(1:11, 13:26)] <- 0
setArg("train_data", config, args, NULL)


config$full_design <- read.table(config$table, header=TRUE, comment.char=";")
config$train_subset <- loadSubset(config$full_design, config)
config$train_design <- config$full_design[config$train_subset,]
config$labels <- loadLabels(config$train_design, config)
config$blockVar <- loadBlockColumn(config, config$train_design)



config$maskVolume <- loadMask(config)

rowIndices <- which(config$train_subset)
flog.info("number of trials: %s", length(rowIndices))
flog.info("max trial index: %s", max(rowIndices))

config$train_datavec <- loadBrainData(config, indices=which(config$train_subset))
print(paste("subset contains", nrow(config$train_design), "of", nrow(config$full_design), "rows."))

caret_model <- loadModel(config$model)

library(caret_model$library, character.only=TRUE)

dataset <- MVPADataset(config$train_datavec, config$labels, config$maskVolume, config$blockVar)
searchres <- mvpa_searchlight(dataset$trainVec, dataset$Y, dataset$mask,dataset$blockVar, 
                              config$radius, config$model, ncores=config$pthreads, 
                              niter=config$niter, tuneGrid=config$tune_grid)


lapply(1:length(searchres), function(i) {
  out <- paste0(config$output, "/", names(searchres)[i], ".nii")
  writeVolume(searchres[[i]], out)  
})

configout <- paste0(config$output, "/config.yaml")
qwrite(as.list(configParams), configout)

