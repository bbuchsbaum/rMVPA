#! /usr/bin/env Rscript

library(neuroim)
library(rMVPA)

carg <- commandArgs(trailingOnly=TRUE)[[1]]

MVPA_CONFIG <- new.env()

source(carg, MVPA_CONFIG)

setDefault("radius", MVPA_CONFIG, 8)
setDefault("table", MVPA_CONFIG, "design.txt")
setDefault("type", MVPA_CONFIG, "standard")
setDefault("method", MVPA_CONFIG, "corsim")
setDefault("ncores", MVPA_CONFIG, 1)
setDefault("labelColumn", MVPA_CONFIG, "labels")
setDefault("output", MVPA_CONFIG, paste0("SearchOut_", config$labelColumn))
setDefault("blockColumn", MVPA_CONFIG, "run")
setDefault("autobalance", MVPA_CONFIG, FALSE)
setDefault("tuneGrid", MVPA_CONFIG, NULL)
setDefault("method_params", MVPA_CONFIG, list())

MVPA_CONFIG$output <- makeOutputDir(MVPA_CONFIG$output)

MVPA_CONFIG$logFile <- file(paste0(MVPA_CONFIG$output, "/searchlight.log"), "w")


MVPA_CONFIG$maskVolume <- loadMask(MVPA_CONFIG)
#MVPA_CONFIG$maskVolume[,,c(1:11, 13:26)] <- 0

MVPA_CONFIG$full_design <- read.table(MVPA_CONFIG$table, header=TRUE, comment.char=";")


MVPA_CONFIG$train_subset <- loadSubset(MVPA_CONFIG$full_design, MVPA_CONFIG)
MVPA_CONFIG$train_design <- MVPA_CONFIG$full_design[MVPA_CONFIG$train_subset,]
MVPA_CONFIG$labels <- loadLabels(MVPA_CONFIG$train_design, MVPA_CONFIG)
MVPA_CONFIG$blockVar <- loadBlockColumn(MVPA_CONFIG, MVPA_CONFIG$train_design)
MVPA_CONFIG$train_datavec <- loadBrainData(MVPA_CONFIG, indices=which(MVPA_CONFIG$train_subset))
print(paste("subset contains", nrow(MVPA_CONFIG$train_design), "of", nrow(MVPA_CONFIG$full_design), "rows."))

MVPA_CONFIG$model <- loadModel(MVPA_CONFIG$method)
library(MVPA_CONFIG$model$library, character.only=TRUE)

dataset <- MVPADataset(MVPA_CONFIG$train_datavec, MVPA_CONFIG$labels, MVPA_CONFIG$maskVolume, MVPA_CONFIG$blockVar)
searchres <- searchlight(dataset$trainVec, dataset$Y, dataset$mask,dataset$blockVar, MVPA_CONFIG$radius, MVPA_CONFIG$method, ncores=MVPA_CONFIG$ncores)


lapply(1:length(searchres), function(i) {
  out <- paste0(MVPA_CONFIG$output, "/", names(searchres)[i], ".nii")
  writeVolume(searchres[[i]], out)  
})


