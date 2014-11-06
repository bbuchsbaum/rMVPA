#! /usr/bin/env Rscript

library(neuroim)
library(rMVPA)

carg <- commandArgs(trailingOnly=TRUE)[[1]]

config <- new.env()

source(carg, config)

setDefault("table", config, "design.txt")
setDefault("method", config, "corsim")
setDefault("ncores", config, 1)
setDefault("labelColumn", config, "labels")
setDefault("output", config, paste0("Out_", config$labelColumn))
setDefault("blockColumn", config, "run")
setDefault("autobalance", config, FALSE)
setDefault("tuneGrid", config, NULL)
setDefault("method_params", config, list())

config$output <- makeOutputDir(config$output)
config$logFile <- file(paste0(config$output, "/rMVPA.log"), "w")


config$ROIVolume <- loadMask(config)
config$maskVolume <- LogicalBrainVolume(as.logical(config$ROIVolume > 0), space(config$ROIVolume))
  
  
#config$maskVolume[,,c(1:11, 13:26)] <- 0

config$full_design <- read.table(config$table, header=TRUE, comment.char=";")


config$train_subset <- loadSubset(config$full_design, config)
config$train_design <- config$full_design[config$train_subset,]
config$labels <- loadLabels(config$train_design, config)
config$blockVar <- loadBlockColumn(config, config$train_design)
config$train_datavec <- loadBrainData(config, indices=which(config$train_subset))
print(paste("subset contains", nrow(config$train_design), "of", nrow(config$full_design), "rows."))

config$model <- loadModel(config$method)
library(config$model$library, character.only=TRUE)

dataset <- MVPADataset(config$train_datavec, config$labels, config$maskVolume, config$blockVar)
mvpa_res <- mvpa_regional(dataset$trainVec, dataset$Y, config$ROIVolume, dataset$blockVar, config$method, ncores=config$ncores, tuneGrid=config$tuneGrid)


lapply(1:length(mvpa_res$outVols), function(i) {
  out <- paste0(config$output, "/", names(mvpa_res$outVols)[i], ".nii")
  writeVolume(mvpa_res$outVols[[i]], out)  
})


