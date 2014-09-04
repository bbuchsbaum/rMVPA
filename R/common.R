

setDefault <- function(name, config, default) {
  if (is.null(config[[name]])) {
    config[[name]]<- default
  }
}

makeOutputDir <- function(dirname) {
  if (!file.exists(dirname)) {
    system(paste("mkdir", dirname))
    dirname
  } else {
    dirname <- paste(dirname, "+", sep="")
    Recall(dirname)
  }
}

abort <- function(msg) {
  writeLines(msg, get("logFile", .GlobalEnv))
  close(logFile)
  print(msg)
  stop()
}


log <- function(msg) {
  writeLines(msg, get("logFile", .GlobalEnv))
}

loadMask <- function(config) {
  if (file.exists(config$mask)) {
    mask <- loadVolume(config$mask)
  } else {
    abort(paste("cannot find mask file named: ", config$mask))
  }
}

loadDesign <- function(config) {
  if (!file.exists(config$table)) {
    abort("cannot find table named", config$table)
  } else {
    read.table(config$table, header=TRUE, comment.char=";")
  }
}

loadLabels <- function(design, config) {
  if (is.null(design[[config$labelColumn]])) {
    abort(paste("Error: labelColumn", config$labelColumn, "not found"))
  } else {
    config$labels <- design[[config$labelColumn]]
  }
  config$labels
}

loadSubset <- function(design, config) {
  if(is.null(config$subset)) rep(TRUE, nrow(full_design)) else {
    subexpr <- config$subset[[2]]
    keep <- eval(subexpr, full_design)
    if (sum(keep) == nrow(full_design)) {
      warning(paste("subset has same number of rows as full table"))
    }
    
    keep
  }
}


  
}