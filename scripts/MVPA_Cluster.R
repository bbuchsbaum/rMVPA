#! /usr/bin/env Rscript

.suppress <- suppressPackageStartupMessages
.suppress(library(neuroim))
.suppress(library(rMVPA))
.suppress(library(optparse))
.suppress(library(futile.logger))

option_list <- list(
  make_option(c("-i", "--input"), type="character", help="the file name of the 4D image vecgtor to cluster"),
  make_option(c("-k", "--clusters"), type="numeric", help="the number of clusters to compute"),
  make_option(c("-s", "--shrink"), type="numeric", help="number of shrinkage iterations (default = 0)", default=0),
  make_option(c("-m", "--mask"), type="character", help="the file name of the volumetric image mask"),
  make_option(c("-e", "--exp"), type="numeric", help="the decay parameter for spatial penalty. Larger value increases the spatial 'pull' (default = .2)", default=.2),
  make_option(c("-o", "--output"), type="character", help="the name of the output file containing the clustered image")
)

oparser <- OptionParser(usage = "MVPA_Cluster.R [options]", option_list=option_list)
opt <- parse_args(oparser, positional_arguments=TRUE)
args <- opt$options

flog.info("command line args are ", args, capture=TRUE)

if (!file.exists(args$input)) {
  flog.error("cannot find input file %s", args$input)
  stop()
}

if (is.null(args$output)) {
  args$output <- paste0("clusout_K", args$k, ".nii")
}

if (!file.exists(args$mask)) {
  flog.error("cannot find mask file %s", args$mask)
  stop()
}


bvec <- loadVector(args$input)
mask <- loadVolume(args$mask)

res <- rMVPA:::slic_cluster(mask, bvec, args$clusters, decay=-(args$exp), nn=8, shrink=args$shrink)

writeVolume(res$clusvol, args$output)


