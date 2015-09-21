#! /usr/bin/env Rscript

.suppress <- suppressPackageStartupMessages
.suppress(library(neuroim))
.suppress(library(rMVPA))
.suppress(library(optparse))
.suppress(library(futile.logger))
.suppress(library(parallel))

option_list <- list(
  make_option(c("-d", "--design"), type="character", help="the name of the design file"),
  make_option(c("-p", "--polort"), type="numeric", help="the order of the polynomial for regressing out low-frequency trends (default is 0)", default=0),
  make_option(c("-z", "--zscore"), type="logical", action="store_true", help="z-score time series within each run (default is false)", default=FALSE),
  make_option(c("-b", "--block_column"), type="character", help="the name of the column in design file indicating the block number"),
  make_option(c("-o", "--output"), type="character", help="the name of the output image file, if missing a 'p' is prefixed to the input file"),
  make_option(c("-i", "--input"), type="character", help="the name of the input image file, which can be a concatenated set of image blocks"),
  make_option(c("-m", "--mask"), type="character", help="the name of the image mask used when computing global mean")
)


oparser <- OptionParser(usage = "prep_betas.R [options]", option_list=option_list)
opt <- parse_args(oparser, positional_arguments=TRUE)
args <- opt$options

flog.info("command line args are ", args, capture=TRUE)

inputFile <- args$input
outputFile <- if (is.null(args$output)) {
  paste0("p", args$input)
} else {
  args$ouput
}

if (!is.null(args$design)) {
  ## load design table
  design <- read.table(args$design, header=TRUE, comment.char=";")
} 

mask <- loadVolume(args$mask)
mask.idx <- which(mask>0)
mat <- as.matrix(loadVector(inputFile))
mat <- mat[mask.idx,]

global <- colMeans(mat)
run <- factor(design[[args$block_column]])
time <- unlist(lapply(split(run, run), function(r) seq_along(r)))

polort <- args$polort
zscore <- args$zscore

flog.info("number of blocks %s", length(levels(run)))
flog.info("time point per block", table(run), capture=TRUE)
flog.info("image matrix dim", dim(mat), capture=TRUE)
flog.info("zscore %s", zscore)
flog.info("polynomial order %s", polort)


res <- do.call(cbind, mclapply(1:nrow(mat), function(i) {
  vec <- mat[i,]
  lm.1 <- if (polort > 0) {
    lm(vec ~ global + run*poly(time,polort))
  } else {
    lm(vec ~ global + run)
  }
  
  if (zscore) {
    scale(resid(lm.1))
  } else {
    resid(lm.1)
  }
}))

svec <- SparseBrainVector(res, space(mask), as.logical(mask))
writeVector(svec, outputFile)
