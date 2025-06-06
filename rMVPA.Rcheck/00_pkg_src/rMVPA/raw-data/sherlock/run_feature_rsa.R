#! /usr/bin/env Rscript

#library(rMVPA)
library(dplyr)
library(neuroim2)

#init_future(8, type="callr")

options(future.globals.maxSize= 8000*1024^2)


pdir <- "."

cropinfo <- read.table("movie_scan1_crop.tsv", header=TRUE)

#bp <- bidser::bids_project("/project/b/brad/dsets/sherlock",fmriprep=TRUE, prep_dir="derivatives")

scan <- paste0("sherlock_part1_small.nii.gz")

feats <- read.table("all_features_pca.csv", header=TRUE, sep=",")
feats <- feats %>% filter(part == "part1") %>% select(starts_with("low"))

# Generate new time indices
new_time <- seq(1, nrow(feats), length.out = 951)


featsd <- apply(feats, 2, function(vals) {
    approx(1:length(vals), vals, xout = new_time)$y
})

svec <- read_vec(scan)
cbeg <- cropinfo$Crop_from_beginning[4]
svec <- sub_vector(svec, (cbeg+1):(cbeg+951))

mask <- read_vol("group_bold_mask_small.nii")

rdes <- rMVPA::feature_rsa_design(F=featsd, labels=paste0("lab_", 1:951), k=99)
dset <- rMVPA::mvpa_dataset(svec, mask=mask)

block <- rep(1:10, length.out=951)
cval <- rMVPA::blocked_cross_validation(block)
rmod <- rMVPA::feature_rsa_model(dset, rdes, method="pca", crossval=cval, cache_pca=FALSE, nperm=100)


slight <- run_searchlight.feature_rsa_model(rmod, radius=12, method="standard")
