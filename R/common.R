

#' @export
#' @import stringr
initializeConfiguration <- function(args) {
  
  if (!is.null(args$config)) {
    if (! file.exists(args$config)) {
      flog.error("cannot find configuration file: %s", args$config)
      stop()
    } else if (str_detect(args$config, "\\.yaml$")) {
      confyaml <- qread(args$config)
      config <- as.environment(confyaml)
    } else if (str_detect(args$config, "\\.[rR]")) {
      config <- new.env()
      source(args$config, config)
    }
  }
  
  config

  
}

#' @export
initializeStandardParameters <- function(config, args, analysisType) {
  #setArg("radius", config, args, 8)
  setArg("train_design", config, args, "mvpa_design.txt")
  setArg("test_design", config, args, NULL)
  setArg("train_data", config, args, "mvpa_design.txt")
  setArg("test_data", config, args, NULL)
  #setArg("type", config, args, "randomized")
  setArg("model", config, args, "corsim")
  setArg("pthreads", config, args, 1)
  setArg("label_column", config, args, "labels")
  setArg("output", config, args, paste0(analysisType, "_", config$labelColumn))
  setArg("block_column", config, args, "block")
  setArg("normalize", config, args, FALSE)
  #setDefault("autobalance", config, FALSE)
  setArg("tune_grid", config, args, NULL)
  #setDefault("method_params", config, list())
  #setArg("niter", config, args, 4)
  setArg("mask", config, args, NULL)
  config
}

#' @export
initializeData <- function(config) {
  
  if (!is.null(config$train_subset)) {
    config$train_datavec <- loadBrainData(config, "train_data", indices=which(config$train_subset))    
  } else {
    config$train_datavec <- loadBrainData(config, "train_data")  
  }
  
  
  if (!is.null(config$test_data)) {
    flog.info("loading test data: %s", config$test_data)
    if (!is.null(config$test_subset)) {
      config$test_datavec <- loadBrainData(config, "test_data", indices=which(config$test_subset))
    } else {
      config$test_datavec <- loadBrainData(config, "test_data")
    }
  }
  
  if (config$normalize) {
    flog.info("Normalizing: centering and scaling each volume of training data")
    norm_datavec <- do.call(cbind, eachVolume(config$train_datavec, function(x) scale(x), mask=config$maskVolume))
    config$train_datavec <- SparseBrainVector(norm_datavec, space(config$train_datavec), mask=config$maskVolume)
    
    if (!is.null(config$test_data)) {
      flog.info("Normalizing: centering and scaling each volume of test data")
      norm_datavec <- do.call(cbind, eachVolume(config$test_datavec, function(x) scale(x), mask=config$maskVolume))
      config$test_datavec <- SparseBrainVector(norm_datavec, space(config$test_datavec), mask=config$maskVolume)
    }
  }
  
  config

  
}

#' @export
initializeDesign <- function(config) {
  if (is.character(config$train_subset)) {
    config$train_subset <- eval(parse(text=config$train_subset))
  }
  
  if (is.character(config$test_subset)) {
    config$test_subset <- eval(parse(text=config$test_subset))
  }
  
  config$full_train_design <- read.table(config$train_design, header=TRUE, comment.char=";")
  config$train_subset <- loadSubset(config$full_train_design, config$train_subset)
  config$train_design <- config$full_train_design[config$train_subset,]
  config$labels <- loadLabels(config$train_design, config)
  config$block <- loadBlockColumn(config, config$train_design)
  
  flog.info(paste("training subset contains", nrow(config$train_design), "of", nrow(config$full_design), "rows."))
  
  if (!is.null(config$test_design)) {
    config$full_test_design <- read.table(config$test_design, header=TRUE, comment.char=";")
    config$test_subset <- loadSubset(config$full_test_design, config$test_subset)
    config$test_design <- config$full_test_design[config$test_subset,]
    config$testLabels <- loadLabels(config$test_design, config)     
    flog.info(paste("test subset contains", nrow(config$test_design), "of", nrow(config$full_test_design), "rows."))
    
  }
  
  config
  
}

#' @export
initializeTuneGrid <- function(args, config) {
  if (!is.null(args$tune_grid) && !args$tune_grid == "NULL") {
    params <- try(expand.grid(eval(parse(text=args$tune_grid))))
    if (inherits(params, "try-error")) {
      stop("could not parse tune_grid expresson: ", args$tune_grid)
    }
    flog.info("tuning grid is", params, capture=TRUE)
    config$tune_grid <- params
  } else if (!is.data.frame(config$tune_grid)) {
    params <- try(lapply(config$tune_grid, function(x) eval(parse(text=x))))
    if (inherits(params, "try-error")) {
      stop("could not parse tune_grid expresson: ", config$tune_grid)
    }
    
    config$tune_grid <- expand.grid(params)
    flog.info("tuning grid is", params, capture=TRUE)
  }
  
  config
}


#' @export
setDefault <- function(name, config, default) {
  if (is.null(config[[name]])) {
    config[[name]]<- default
  }
}

#' @export
setArg <- function(name, config, args, default) {
  if (is.null(config[[name]]) && is.null(args[[name]])) {
    config[[name]] <- default
  } else if (!is.null(args[[name]])) {
    config[[name]] <- args[[name]]
  } else if (is.null(config[[name]])) {
    config[[name]] <- default
  }    
}

#' @export
makeOutputDir <- function(dirname) {
  if (!file.exists(dirname)) {
    system(paste("mkdir", dirname))
    dirname
  } else {
    dirname <- paste(dirname, "+", sep="")
    Recall(dirname)
  }
}

#' @export
abort <- function(config, msg) {
  stop(msg)
}

#' @export
logit <- function(config, msg) {
  #writeLines(msg, config$logFile)
}

#' @export
loadModel <- function(name) {
  registry <- get("MVPAModels", .GlobalEnv)
  if (!is.null(registry[[name]])) {
    registry[[name]]       
  } else if (length(caret::getModelInfo(name)) > 0) {
    caret::getModelInfo(name)[[name]]    
  } else {
    abort(paste("unrecognized model: ", name))
  }
}

#' @export
loadMask <- function(config) {
  if (file.exists(config$mask)) {
    mask <- loadVolume(config$mask)
  } else {
    stop(paste("cannot find mask file named: ", config$mask))
  }
  
  mask
}

#' @export
loadDesign <- function(config, name) {
  if (!file.exists(config[[name]])) {
    abort("cannot find table named", config$table)
  } else {
    read.table(config[[name]], header=TRUE, comment.char=";")
  }
}

#' @export
loadLabels <- function(full_design, config) {
  if (is.null(full_design[[config$label_column]])) {
    abort(paste("Error: labelColumn", config$label_column, "not found"))
  } else {
    labels <- full_design[[config$label_column]]
  }
  labels
}

#' @export
loadSubset <- function(full_design, subset) {
  if (is.character(subset)) {
    subset <- eval(parse(text=subset))
  }
  
  keep <- if(is.null(subset)) rep(TRUE, nrow(full_design)) else {
    subexpr <- subset[[2]]
    keep <- eval(subexpr, full_design)
    if (sum(keep) == nrow(full_design)) {
      warning(paste("subset has same number of rows as full table"))
    }
    
    keep
  }
  
  keep
  
}

#' @export
loadBlockColumn <- function(config, design) {
  if (is.null(design[[config$block_column]])) {
    abort(paste("blockColumn variable named", config$blockColumn, "not found."))
  } else {  
    config$nfolds <- length(unique(design[[config$block_column]]))
    design[[config$block_column]]
  }
   
}

#' @export
#' 
loadBrainDataSequence <- function(fnames, config) {
  if (!all(file.exists(fnames))) {
    offenders <- fnames[!file.exists(fnames)]
    abort(config, paste("training data", offenders, "not found."))
  }
  
  vecmat <- do.call(rbind, lapply(fnames, function(fname) {
    flog.info("loading data file %s", fname)
    mat <- as.matrix(loadVector(fname, mask=config$maskVolume))
    flog.info("data file %s has %s voxels and %s samples", fname, ncol(mat), nrow(mat))
    mat
  }))
  
  SparseBrainVector(vecmat, space(config$maskVolume), mask=config$maskVolume)
}

#' @export
loadBrainData <- function(config, name, indices=NULL) {
  fname <- config[[name]]
  if (length(fname) > 1) {
    loadBrainDataSequence(fname, config)
  } else if (!file.exists(fname)) {
    abort(config, paste("training data", fname, "not found."))
  } else {
    flog.info("loading data file %s", fname)
    if (!is.null(indices)) {
      loadVector(fname, indices=indices, mask=config$maskVolume)
    } else {
      loadVector(fname, mask=config$maskVolume)
    }
    
  }
}


  