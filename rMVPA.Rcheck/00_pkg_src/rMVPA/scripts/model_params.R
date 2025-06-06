#! /usr/bin/env Rscript

.suppress <- suppressPackageStartupMessages

.suppress(library(rMVPA))
.suppress(library(futile.logger))
.suppress(library(caret))
.suppress(library(optparse))

option_list <- list(
  make_option(c("-m", "--model"), type="character", help="the name of the clasisifcation model")
)


oparser <- OptionParser(usage = "model_params.R [options]", option_list=option_list)
opt <- parse_args(oparser, positional_arguments=TRUE)
args <- opt$options

#flog.info("command line args are ", args, capture=TRUE)

print(rMVPA::loadModel(args$model)$parameters)
print("\\n")