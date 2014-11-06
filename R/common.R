
#' @export
setDefault <- function(name, config, default) {
  if (is.null(config[[name]])) {
    config[[name]]<- default
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
  logFile <- config$logFile
  writeLines(msg, logFile)
  close(logFile)
  print(msg)
  stop()
}

#' @export
logit <- function(config, msg) {
  writeLines(msg, config$logFile)
}

#' @export
loadModel <- function(name) {
  registry <- get("MVPAModels", .GlobalEnv)
  if (!is.null(registry[[name]])) {
    registry[[name]]       
  } else if (length(caret::getModelInfo(name)) > 0) {
    caret::getModelInfo(name)[[name]]    
  } else {
    abort("unrecognized model: ", name)
  }
}

#' @export
loadMask <- function(config) {
  if (file.exists(config$mask)) {
    mask <- loadVolume(config$mask)
  } else {
    abort(paste("cannot find mask file named: ", config$mask))
  }
  
  mask
}

#' @export
loadDesign <- function(config) {
  if (!file.exists(config$table)) {
    abort("cannot find table named", config$table)
  } else {
    read.table(config$table, header=TRUE, comment.char=";")
  }
}

#' @export
loadLabels <- function(full_design, config) {
  if (is.null(full_design[[config$labelColumn]])) {
    abort(paste("Error: labelColumn", config$labelColumn, "not found"))
  } else {
    labels <- full_design[[config$labelColumn]]
  }
  labels
}

#' @export
loadSubset <- function(full_design, config) {
  keep <- if(is.null(config$subset)) rep(TRUE, nrow(full_design)) else {
    subexpr <- config$subset[[2]]
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
  if (is.null(design[[config$blockColumn]])) {
    abort(paste("blockColumn variable named", config$blockColumn, "not found."))
  } else {  
    config$nfolds <- length(unique(design[[config$blockColumn]]))
    design[[config$blockColumn]]
  }
   
}

#' @export
loadBrainData <- function(config, indices=NULL) {
  if (!file.exists(config$train_data)) {
    abort(config, paste("training data", config$train_data, "not found."))
  } else {
    #logit(paste("loading data file", config$train_data))
    if (!is.null(indices)) {
      loadVector(config$train_data, indices=indices, mask=config$maskVolume)
    } else {
      loadVector(config$train_data, mask=config$maskVolume)
    }
    
  }
}


  