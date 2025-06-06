#! /usr/bin/env Rscript

.suppress <- suppressPackageStartupMessages
.suppress(library(neuroim))
.suppress(library(rMVPA))
.suppress(library(optparse))
.suppress(library(futile.logger))
.suppress(library(io))

option_list <- list(
  make_option(c("-t", "--test_design"), type="character", help="the file name of the evaluation design table (optional)"),
  make_option(c("-s", "--test_subset"), type="character", help="a quoted R expression used to subset the test_design data frame, e.g. 'id == 5001' -- Note: This will not subset the *image* data"),
  make_option(c("-f", "--folder"), type="character", help="the location of the output folder of the previously trained model"),
  make_option(c("-n", "--roinum"), type="character", help="the index of the region of interest (ROI) for which to make predictions (default is to make predicitons for all regions)"),
  make_option(c("-d", "--newdata"), type="character", help="the *new* image data for used for model predictions (should be a NIFTI formatted 3D or 4D image file)"),
  make_option(c("-o", "--output"), type="character", default="model_predictions.txt", help="the name of the output file containing model predictions. If 'test_design' is specified, predictions will be added as new columns design file.")
)

oparser <- OptionParser(usage = "MVPA_Predict.R [options]", option_list=option_list)
opt <- parse_args(oparser, positional_arguments=TRUE)
args <- opt$options

flog.info("command line args are ", args, capture=TRUE)

predictorName <- paste0(args$folder, "/predictor.RDS")
predictor <- readRDS(predictorName)
config <- qread(paste0(args$folder, "/config.yaml"))
testvec <- read_vec(args$newdata)


if (config$normalize) {
  flog.info("normalizing test vector")
  mask <- read_vol(config$mask)
  testvec <- normalizeSamples(testvec, mask)
}

outfile <- args$output
preds <- evaluateModel(predictor, testvec)

if (is.null(args$test_design)) {
  if (length(preds) == 1) {
    flog.info("Writing predictions to %s", args$output)
    p1 <- preds[[1]]
    oframe <- cbind(p1$class, p1$prob, rep(names(preds)[[1]], nrow(p1$prob)))
    names(oframe) <- c("pred_class", paste0("prob_", names(p1$prob)), "ROI")
    write.table(format(oframe, digits=2, scientific=FALSE, drop0trailing=TRUE), args$output, row.names=FALSE, quote=FALSE)
  } else {
    flog.info("Writing predictions to %s", args$output)
    oframe <- do.call(rbind, lapply(1:length(preds), function(i) {
      p1 <- preds[[i]]
      oframe <- cbind(p1$class, p1$prob, rep(names(preds)[[i]], nrow(p1$prob)))
      names(oframe) <- c("pred_class", paste0("prob_", names(p1$prob)), "ROI")
      oframe
    }))
    
    write.table(format(oframe,  digits=2, scientific=FALSE, drop0trailing=TRUE), args$output, row.names=FALSE,quote=FALSE)
  }
} else {
  
  flog.info("test data will be evaluated against labels stored in: %s", args$test_design)
  
  
  tdes <- read.table(args$test_design, header=TRUE)
  oframe <- do.call(rbind, lapply(1:length(preds), function(i) {
    cbind(tdes, constructPredictionFrame(preds[i]))
  }))
  
  oframe$pred_correct <- ifelse(oframe$pred_class == oframe[[config$label_column]], 1, 0)
  
  plabs <- sapply(oframe[[config$label_column]], function(lab) paste0("prob_", lab))
  oframe$prob_item <- sapply(1:length(plabs), function(i) oframe[i, plabs[[i]]])
  write.table(format(oframe,  digits=2, scientific=FALSE, drop0trailing=TRUE), args$output, row.names=FALSE,quote=FALSE)
  
  
  
}


