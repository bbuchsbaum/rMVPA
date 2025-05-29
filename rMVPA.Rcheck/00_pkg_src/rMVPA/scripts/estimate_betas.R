#! /usr/bin/env Rscript



library(optparse)

suppressPackageStartupMessages(library("optparse"))

option_list <- list(make_option(c("-s", "--scans"), type="character", help="a file containing an ordered list of NIFTI 4D fMRI scans"),
                    make_option(c("-d", "--design"), type="character", help="design table with at least two columns named: 'onsets' and 'run' (and 'duration' is optional)"),
                    make_option(c("-l", "--duration"), type="numeric", default=1, help="the duration of each stimulus event (constant for all events; note: overrides 'duration' column in design file. default=1)"),  
                    make_option(c("-o", "--out"), type="character", default="betas", help="name of output file stem (default: betas)"),
                    make_option(c("-c", "--concatenate"), type="logical", default="TRUE", help="concatenate all beta estimates into one output file (default is TRUE)"),
                    make_option(c("-n", "--nuisance"), type="character", help="list of nuisance covariate files, one per run each on a separate line"),
                    make_option(c("-t", "--onsets"), type="character", help="name of onset column (default is 'onsets'", default="onsets"),
                    make_option(c("-r", "--runs"), type="character", help="optional name of file containing runs to subselect from design file"),
                    make_option(c("-p", "--polort"), type="numeric", default=3, help="number of polynomial regressors -- passed to as 'polort' 3dDeconvolve (default = 3)"),
                    make_option(c("-m", "--mask"), type="character", help="name of binary image mask"))
                   
oparser <- OptionParser(usage = "estimate_betas.R [options]", option_list=option_list)
opt <- parse_args(oparser, positional_arguments=TRUE)
args <- opt$options

if (!is.null(args$nuisance)) {
  nuisanceFiles <- scan(args$nuisance, "")
  nuisanceCovars <- lapply(nuisanceFiles, read.table, header=TRUE)
} else {
  nuisanceCovars = NULL
}

if (is.null(args$scans)) {
  stop("must supply a file containing list of scans with -s or --scans option")
}

scans <- scan(args$scans, "")
fexists <- sapply(scans, file.exists)

if (!all(fexists)) {
  stop(paste("could not find input file(s):", paste(scans[fexists], collapse=" ")))
}
       
  
if (is.null(args$design)) {
  stop("must supply a design table with -d or --design option")
}

design <- read.table(args$design, header=TRUE)

if (!is.null(args$runs)) {
  runset <- scan(args$runs)
  design <- subset(design, run %in% runset)
}


onsetColumn <- args$onsets

if (is.null(design[[onsetColumn]])) {
  stop(paste("design file must have column named:", onsetColumn))
}

if (is.null(design$run)) {
  stop("design file must have column named 'run'")
}


runs <- design$run
runLevs <- sort(unique(design$run))

message("Scanning runs are: ", runLevs)

if (!is.null(design$duration)) {
  message("found duration column in design file: will use for event duration vector")
  duration <- design$duration
} else {
  duration <- rep(as.numeric(args$duration), nrow(design))
  message("using constant event duration of: ", duration)
}

outfiles <- paste0(args$out, "_", runLevs, ".nii")


genStimFileArgs <- function(fnames) {
  
  paste(lapply(1:length(fnames), function(k) {
    paste("-stim_base", (k+1), "-stim_file", (k+1), fnames[k], "-stim_label", (k+1), paste0("nuisance_", k))
  }), collapse = " ")
}
  

for (run in 1:length(scans)) {

  nuisnames <- if (!is.null(nuisanceCovars)) {
    nuis <- nuisanceCovars[[run]]
    sapply(1:ncol(nuis), function(j) {
      oname <- paste0("nuisance_run_", run, "_reg#", j, ".1D")
      write(nuis[,j], file=oname, ncolumns=1)
      oname
    })
  } else {
    NULL
  }

  nstim <- if (!is.null(nuisanceCovars)) {
    length(nuisnames) + 1
  } else {
    1
  }

  onsets <- design[[onsetColumn]][design$run == runLevs[run]]
  rundurs <- duration[design$run == runLevs[run]]
  onsetMat <- cbind(onsets, rundurs)
  onsetStr <- sapply(1:nrow(onsetMat), function(i) paste0(onsetMat[i,1], ":", onsetMat[i,2]))
  

  onsname <- paste0("onsets_run", runLevs[run], ".1D")

  message(onsetStr)
  write(onsetStr, file=onsname, ncolumns=length(onsets))

  dmat <- paste0("R.xmat.", run, ".1D")

  

  cmd1 <- paste("/thalia/bblab/pkgs/afni_src/3dDeconvolve -input",
                scans[run],
                "-polort", args$polort, "-x1D", dmat,              
                "-x1D_stop -num_stimts", nstim,
                "-stim_times_IM 1", onsname, paste0("dmBLOCK"),
                "-stim_label 1 event",
                if (!is.null(nuisnames)) genStimFileArgs(nuisnames) else "")
                
   

  cmd2 <- paste("3dLSS -verb", "-prefix", outfiles[run], "-input", scans[run], "-mask", args$mask, "-matrix", dmat)
  print(paste("RUNNING:", cmd1))
  system(cmd1)
  print(paste("RUNNING", cmd2))
  system(cmd2)
}


conc_out <- paste0(args$out, "_all.nii.gz")

cmd <- paste("3dTcat -prefix", conc_out, paste(outfiles, collapse=" "))
system(cmd)
