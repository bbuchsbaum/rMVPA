


# initializeROIGrouping <- function(config) {
#   if (!is.null(config$roi_grouping)) {
#     roilist <- lapply(1:length(config$roi_grouping), function(i) {
#       grp <- config$roi_grouping[[i]]
#       idx <- which(config$ROIVolume %in% grp)
#       vol <- makeVolume(refvol=config$ROIVolume)
#       vol[idx] <- i
#       vol
#     })
#     
#     if (is.null(names(config$roi_grouping))) {
#       names(roilist) <- paste0("roi_group_", seq_along(config$roi_grouping))
#     } else {
#       names(roilist) <- names(config$roi_grouping)
#     }
#     config$ROIVolume <- roilist
#   } else {
#     config$ROIVolume <- list(config$ROIVolume)
#   }
# }


# initializeROISubset <- function(config) {
#   if (!is.null(config$roi_subset)) {
#     form <- try(eval(parse(text=config$roi_subset)))
#     if (inherits(form, "try-error")) {
#       flog.error("could not parse roi_subset parameter: %s", config$roi_subset)
#       stop()
#     }
#     
#     if (class(form) != "formula") {
#       flog.error("roi_subset argument must be an expression that starts with a ~ character")
#       stop()
#     }
#     
#     res <- as.logical(eval(form[[2]], list(x=config$ROIVolume)))
#     
#     config$ROIVolume[!res] <- 0
#     flog.info("roi_subset contains %s voxels", sum(config$ROIVolume > 0))
#   }
# }
# 
# 
# initMVPARegional
#  
# @param configFile an 'MVPA_Regional' configuration file
# @param args list of addiitonal override arguments
# @export
# initMVPARegional <- function(configFile, args=list(), verbose=FALSE) {
#   if (!verbose) {
#     flog.threshold(ERROR)
#   } else {
#     flog.threshold(DEBUG)
#   }
#   
#   config <- initializeConfiguration(list(config=configFile))
#   config <- initializeStandardParameters(config, args, "mvpa_regional")
# 
#   set_arg("savePredictors", config, args, FALSE)
#   
#   config <- initializeTuneGrid(args, config)
#   configParams <- as.list(config)
#   config <- initializeDesign(config)
#   
#   rowIndices <- which(config$train_subset)
#   config$ROIVolume <- loadMask(config)
#   
#   initializeROISubset(config)
#   initializeROIGrouping(config)
#   
#   parcellationVolume <- if (!is.null(config$parcellation)) {
#     loadVolume(config$parcellation)
#   }
#   
#   config$maskVolume <- as(Reduce("+", lapply(config$ROIVolume, function(roivol) as(roivol, "LogicalBrainVolume"))), "LogicalBrainVolume")
#   config <- initializeData(config)
#   
#   flog.info("number of training trials: %s", length(rowIndices))
#   
#   flog.info("max trial index: %s", max(rowIndices))
#   
#   flog.info("loading training data: %s", config$train_data)
#   
#   flog.info("mask contains %s voxels", sum(config$maskVolume))
#   
#   for (i in seq_along(config$ROIVolume)) {
#     rvol <- config$ROIVolume[[i]]
#     
#     flog.info("Region mask contains: %s ROIs", length(unique(rvol[rvol > 0])))
#     flog.info(paste("ROIs are for group ", i, "are:"), rvol, capture=TRUE)
#     
#   }
#   
#   flog.info("Running regional MVPA with parameters:", configParams, capture=TRUE)
#   
#   flog.info("With %s roi groups", length(config$ROIVolume))
#   
#   if (length(config$labels) != dim(config$train_datavec)[4]) {
#     flog.error("Number of volumes: %s must equal number of labels: %s", dim(config$train_datavec)[4], length(config$labels))
#     stop()
#   }
#   
#   featureSelector <- if (!is.null(config$feature_selector)) {
#     FeatureSelector(config$feature_selector$method, config$feature_selector$cutoff_type, as.numeric(config$feature_selector$cutoff_value))
#   }
#   
#   flog.info("feature selector: ", featureSelector, capture=TRUE)
#   
#   flog.info("bootstrap replications: ", config$bootstrap_replications, capture=TRUE)
#   
#   dataset <- MVPADataset(config$train_datavec, config$labels, config$maskVolume, config$block, config$test_datavec, 
#                          config$testLabels, modelName=config$model, tuneGrid=config$tune_grid,
#                          tuneLength=config$tune_length, testSplitVar=config$testSplitVar, testSplits=config$testSplits, 
#                          trainDesign=config$train_design,
#                          testDesign=config$test_design)
#   
#   for (varname in c("test_subset", "train_subset", "roi_subset", "split_by")) {
#     if (!is.null(configParams[[varname]]) && is(configParams[[varname]], "formula")) {
#       configParams[[varname]] <- Reduce(paste, deparse(configParams[[varname]]))
#     }
#   }
#   
#   for (lib in dataset$model$library) {
#     library(lib, character.only = TRUE)
#   }
#   
#   list(dataset=dataset, config=config)
#   
# }
# 
# 
# # @export
# initMVPASearchlight <- function(configFile, args=list(), verbose=FALSE) {
#   if (!verbose) {
#     flog.threshold(ERROR)
#   } else {
#     flog.threshold(DEBUG)
#   }
#   
#   
#   config <- initializeConfiguration(list(config=configFile))
#   config <- initializeStandardParameters(config, args, "mvpa_searchlight")
#   
#   set_arg("niter", config, args, 16)
#   set_arg("radius", config, args, 8)
#   set_arg("type", config, args, "randomized")
#   
#   config <- initializeTuneGrid(args, config)
#   configParams <- as.list(config)
#   config <- initializeDesign(config)
#   
#   config$maskVolume <- as(loadMask(config), "LogicalBrainVolume")
#   
#   
#   rowIndices <- which(config$train_subset)
#   config$ROIVolume <- loadMask(config)
#   
#   rowIndices <- which(config$train_subset)
#   flog.info("number of trials: %s", length(rowIndices))
#   flog.info("max trial index: %s", max(rowIndices))
#   flog.info("loading training data: %s", config$train_data)
#   flog.info("mask contains %s voxels", sum(config$maskVolume))
#   
#   config <- initializeData(config)
#   
#   flog.info("Running searchlight with parameters:", configParams, capture=TRUE)
#   
#   
#   dataset <- MVPADataset(config$train_datavec, 
#                          config$labels, 
#                          config$maskVolume, 
#                          config$block, 
#                          config$test_datavec, 
#                          config$testLabels, 
#                          modelName=config$model, 
#                          tuneGrid=config$tune_grid,
#                          tuneLength=config$tune_length, 
#                          testSplitVar=config$testSplitVar, 
#                          testSplits=config$testSplits,
#                          trainDesign=config$train_design,
#                          testDesign=config$test_design)
#   
#   for (lib in dataset$model$library) {
#     library(lib, character.only = TRUE)
#   }
#   
#   list(dataset=dataset, config=config)
# }




#' @import stringr
initialize_configuration <- function(args) {
  
  if (!is.null(args$config)) {
    if (!file.exists(args$config)) {
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


#' @noRd
initialize_standard_parameters <- function(config, args, analysisType) {
  set_arg("train_design", config, args, "mvpa_design.txt")
  set_arg("test_design", config, args, NULL)
  set_arg("train_data", config, args, "mvpa_design.txt")
  set_arg("test_data", config, args, NULL)
  set_arg("model", config, args, "corsim")
  set_arg("feature_selector", config, args, NULL)
  set_arg("pthreads", config, args, 1)
  set_arg("label_column", config, args, "labels")
  set_arg("skip_if_folder_exists", config, args, FALSE)
  set_arg("output", config, args, paste0(analysisType, "_", config$labelColumn))
  set_arg("block_column", config, args, "block")
  set_arg("normalize_samples", config, args, FALSE)
  set_arg("tune_grid", config, args, NULL)
  set_arg("mask", config, args, NULL)
  set_arg("class_metrics", config, args, TRUE)
  set_arg("split_by", config, args, NULL)
  set_arg("custom_performance", config, args, NULL)
  set_arg("test_label_column", config, args, NULL)
  set_arg("data_mode", config, args, "image")
  
  config
}

#' @noRd
normalize_image_samples <- function(bvec, mask) {
  norm_datavec <- do.call(cbind, eachVolume(bvec, function(x) scale(x)[,1], mask=mask))
  SparseBrainVector(norm_datavec, space(bvec), mask=mask)  
}


#' @noRd
normalize_surface_samples <- function(bvec, mask) {
  mat <- scale(bvec@data[indices(bvec), ,drop=FALSE])
  
  m2 <- matrix(0, length(nodes(bvec)), ncol(bvec@data))
  m2[indices(bvec),] <- mat
  
  BrainSurfaceVector(geometry(bvec), indices=indices(bvec), m2)
}

initialize_surface_data <- function(config) {
  if (!is.null(config$train_subset)) {
    indices <- which(config$train_subset)
    flog.info("length of training subset %s", length(indices))
  }
  
  train_surfaces <- load_surface_data(config, "train_data", colind=indices) 
  
  if (!is.null(config$test_data)) {
    flog.info("loading test surface data: %s", config$test_data)
    indices <- which(config$test_subset)
    flog.info("length of test subset %s", length(indices))
    
    test_surfaces <- if (!is.null(config$test_subset)) {
      load_surface_data(config, "test_data", colind=indices)
    } else {
      load_surface_data(config, "test_data")
    }
    
  } else {
    test_surfaces <- NULL
  }
  
  if (config$normalize_samples) {
    flog.info("Normalizing: centering and scaling each surface of training data")
    ret <- lapply(train_surfaces, normalize_surface_samples)
    names(ret) <- names(train_surfaces)
    train_surfaces <- ret
    
    if (!is.null(test_surfaces)) {
      flog.info("Normalizing: centering and scaling each surface of test data")
      ret <- lapply(test_surfaces, normalize_surface_samples)
      names(ret) <- names(test_surfaces)
      test_surfaces <- ret
    }
  }
  
  if (!is.null(test_surfaces) && !(length(train_surfaces) == length(test_surfaces))) {
    flog.info("number of training surfaces: %s", length(train_surfaces))
    flog.info("number of test surfaces: %s", length(test_surfaces))
    flog.error("number of training surface entries must equal number of test surface entries")
    stop()
  }
  
  ret <- if (!is.null(config$mask)) {
    flog.info("loading mask: %s ", config$mask)
    masksurf <- load_surface_mask(config$mask, train_surfaces)
    
    lapply(seq_along(train_surfaces), function(i) {
      mvpa_surface_dataset(train_surfaces[[i]], test_surfaces[[i]], name=names(train_surfaces)[i], mask=masksurf[[i]])
    })
    
  } else {
    lapply(seq_along(train_surfaces), function(i) {
      mvpa_surface_dataset(train_surfaces[[i]], test_surfaces[[i]], name=names(train_surfaces)[i])
    })
  }
  
  names(ret) <- names(train_surfaces)
  ret

 
}

initialize_feature_selection <- function(config) {
  if (!is.null(config$feature_selector)) {
    feature_selector(config$feature_selector$method, config$feature_selector$cutoff_type, as.numeric(config$feature_selector$cutoff_value))
  } else {
    NULL
  }
}


initialize_image_data <- function(config, mask) {
  if (!is.null(config$train_subset)) {
    indices <- which(config$train_subset)
    flog.info("length of training subset %s", length(indices))
  }
  
  mask_volume <- as(mask, "LogicalBrainVolume")
  
  train_datavec <- load_image_data(config, "train_data", mask_volume=mask_volume, indices=indices)    

  if (!is.null(config$test_data)) {
    flog.info("loading test data: %s", config$test_data)
    indices=which(config$test_subset)
    flog.info("length of test subset %s", length(indices))
    
    if (!is.null(config$test_subset)) {
      test_datavec <- load_image_data(config, "test_data", mask_volume=mask_volume, indices=indices)
    } else {
      test_datavec <- load_image_data(config, "test_data", mask_volume=mask_volume)
    }
  } else {
    test_datavec <- NULL
  }
  
  if (config$normalize_samples) {
    flog.info("Normalizing: centering and scaling each volume of training data")
    train_datavec <- normalize_image_samples(train_datavec, mask_volume)
    
    if (!is.null(test_datavec)) {
      flog.info("Normalizing: centering and scaling each volume of test data")
      test_datavec <- normalize_image_samples(test_datavec, mask_volume)
    }
  }
  
  mvpa_dataset(train_datavec, test_datavec, mask=mask)

}


initialize_design <- function(config) {
  if (is.character(config$train_subset)) {
    config$train_subset <- eval(parse(text=config$train_subset))
  }
  
  if (is.character(config$test_subset)) {
    config$test_subset <- eval(parse(text=config$test_subset))
  }
  

  ## full design
  config$full_train_design <- read.table(config$train_design, header=TRUE, comment.char=";")
  
  ## subset of training samples
  config$train_subset <- load_subset(config$full_train_design, config$train_subset)
  
  ## training design
  config$train_design <- config$full_train_design[config$train_subset,]
  
 
  flog.info(paste("training subset contains", nrow(config$train_design), "of", nrow(config$full_train_design), "rows."))
  
  if (!is.null(config$test_design) && is.null(config$test_data)) {
    flog.error("test_design %s is supplied with no test_data")
    stop()
  }
  
  if (is.null(config$test_design) && !is.null(config$test_data)) {
    flog.error("test_data %s is supplied with no test_design")
    stop()
  }
  
  
  #if (!is.null(config$test_subset) && is.null(config$test_design) && is.null(config$test_data)) {
  #  flog.info("test subset is taken from training design table")
  #  config$test_subset <- load_subset(config$full_train_design, config$test_subset)
  #  
  #  config$test_design <- config$full_train_design[config$test_subset,]
  #  config$full_test_design <- config$test_design
  #  config$testLabels <- loadLabels(config$test_design, config)   
  #}
  
  if (!is.null(config$test_design)) {
    has_test <- TRUE
    flog.info("test design %s is specified", config$test_design)
    config$full_test_design <- read.table(config$test_design, header=TRUE, comment.char=";")
    flog.info(paste("test design contains", nrow(config$full_test_design), "rows."))
    
    config$test_subset <- load_subset(config$full_test_design, config$test_subset)
    config$test_design <- config$full_test_design[config$test_subset,]
    
    flog.info(paste("test subset contains", nrow(config$test_design), "of", nrow(config$full_test_design), "rows."))
    
    #config$testLabels <- loadTestLabels(config$test_design, config)     
    #flog.info(paste("test subset contains", nrow(config$test_design), "of", nrow(config$full_test_design), "rows.")) 
    #flog.info(paste("first 10 test labels: ", head(config$testLabels, 10), capture=TRUE))
    
  } else {
    has_test <- FALSE
    flog.info("testing is via internal cross-validation")
    #config$testLabels <- config$labels
  }
  
  if (!is.null(config$split_by)) {
    flog.info("splitting performance metrics by %s: ", deparse(config$split_by))
  }
  
  
  if (has_test) {
    ## todo what if we don't have/need "block_column"
    mvpa_design(train_design=config$train_design, 
              y_train=config$label_column, 
              test_design=config$test_design, 
              y_test=config$test_label_column, 
              block_var=config$block_column, 
              split_by=config$split_by)
  } else {
    ## todo what if we don't have/need "block_column"
    mvpa_design(train_design=config$train_design, 
                y_train=config$label_column, 
                block_var=config$block_column, 
                split_by=config$split_by)
  }
   
  
  
}


initialize_tune_grid <- function(args, config) {
  if (!is.null(args$tune_grid) && !args$tune_grid == "NULL") {
    params <- try(expand.grid(eval(parse(text=args$tune_grid))))
    
    if (inherits(params, "try-error")) {
      stop("could not parse tune_grid expresson: ", args$tune_grid)
    }
    
    flog.info("tuning grid is", params, capture=TRUE)
    params
    
  } else if (!is.null(config$tune_grid) && !is.data.frame(config$tune_grid)) {
    params <- try(lapply(config$tune_grid, function(x) eval(parse(text=x))))
    if (inherits(params, "try-error")) {
      stop("could not parse tune_grid expresson: ", config$tune_grid)
    }
    
    flog.info("tuning grid is", params, capture=TRUE)
    expand.grid(params)
    
  } else if (is.data.frame(config$tune_grid)) {
    config$tune_grid
  } else {
    NULL
  }
  
  
}



set_default <- function(name, config, default) {
  if (is.null(config[[name]])) {
    config[[name]]<- default
  }
}


set_arg <- function(name, config, args, default) {
  if (is.null(config[[name]]) && is.null(args[[name]])) {
    config[[name]] <- default
  } else if (!is.null(args[[name]])) {
    config[[name]] <- args[[name]]
  } else if (is.null(config[[name]])) {
    config[[name]] <- default
  }    
}


make_output_dir <- function(dirname) {
  if (!file.exists(dirname)) {
    system(paste("mkdir", dirname))
    dirname
  } else {
    dirname <- paste(dirname, "+", sep="")
    Recall(dirname)
  }
}

initialize_crossval <- function(config, des=NULL) {
  cval <- if (is.null(config$cross_validation) && !is.null(des$block_var)) {
    flog.info("cross-validation type: cross validation using predefined blocking variable")
    blocked_cross_validation(des$block_var)
  } else if (!is.null(config$cross_validation)) {
    assertthat::assert_that(!is.null(des))
    cval <- config$cross_validation
    
    if (cval$name == "twofold" || cval$name == "two_fold" || cval$name == "two_fold_blocked_cross_validation") {
      flog.info("cross-validation type: twofold cross-validation.")
      if (is.null(cval$nreps)) {
        cval$nreps <- 10
      }
      flog.info("cross-validation reps: %s ", cval$nreps)
      twofold_blocked_cross_validation(block_var=des$block_var, nreps=cval$nreps)
    } else if (cval$name == "blocked" || cval$name == "blocked_cross_validation") {
      blocked_cross_validation(des$block_var)
    } else if (cval$name == "custom" || cval$name == "custom_cross_validation") {
      custom_cross_validation(cval$sample_set)
    } else {
      flog.error("unrecognized cross_validation type: %s", cval$name)
      stop()
    }
    
  } else {
    flog.info("cross-validation type: 5 fold cross-validation using random splits")
    kfold_cross_validation(nobs(des))
  }
  
  cval
}



#' load_model
#' @param name the name of the model
#' @examples load_model("sda")
#' @export
load_model <- function(name) {
  registry <- MVPAModels
  
  ret <- if (!is.null(registry[[name]])) {
    registry[[name]]   
  } else if (length(caret::getModelInfo(name)) > 0) {
    caret::getModelInfo(name)[[name]] 
  } else {
    stop(paste("unrecognized model: ", name))
  }
  
  ret$label <- name
  
  ret
  
}


load_mask <- function(config) {
  if (config$data_mode == "image") {
    loadVolume(config$mask)
  } else if (config$data_mode == "surface") {
    NULL
  }
}


load_design <- function(config, name) {
  if (!file.exists(config[[name]])) {
    futile.logger::flog.error(paste("cannot find table named: ", name))
    stop()
  } else {
    read.table(config[[name]], header=TRUE, comment.char=";")
  }
}

load_mvpa_model <- function(config, dataset, design, crossval, feature_selector) {
  mod <- load_model(config$model)
  mvp_mod <- mvpa_model(mod,dataset, design=design, 
                        model_type=config$model_type,
                        crossval=crossval,
                        feature_selector=feature_selector, 
                        tune_grid=config$tune_grid,
                        performance=config$performance,
                        class_metrics=config$class_metrics)
  
}


load_subset <- function(full_design, subset) {
  if (is.character(subset)) {
    if (substr(subset, 1,1) != "~") {
      subset <- paste0("~", subset)
    }   
    subset <- eval(parse(text=subset))
  } 
  
  keep <- if(is.null(subset)) rep(TRUE, nrow(full_design)) else {  
    subexpr <- subset[[2]]   
    keep <- eval(subexpr, full_design)
    if (sum(keep) == nrow(full_design)) {
      warning(paste("subset has same number of rows as full table"))
    }
    if (sum(keep) <= 1) {
      futile.error("train_subset %s results in design with only 1 row.", as.character(subexpr))
      stop()
    }
    keep
  }
  
  keep
  
}


load_image_data_series <- function(fnames, config, indices, mask_volume) {
  if (!all(file.exists(fnames))) {
    offenders <- fnames[!file.exists(fnames)]
    message(paste("data files", offenders, "not found."))
    stop()
  }
  
  
  ### TODO make more efficient. This loads in all data then subsets.
  vecmat <- do.call(rbind, lapply(seq_along(fnames), function(i) {
    fname <- fnames[i]
    flog.info("loading data file %s", fname)
    mat <- neuroim::as.matrix(loadVector(fname, mask=mask_volume))
    flog.info("data file %s has %s voxels and %s samples", fname, nrow(mat), ncol(mat))
    mat
  }))
  
  SparseBrainVector(vecmat[indices,], space(mask_volume), mask=mask_volume)
}


load_image_data <- function(config, name, mask_volume, indices=NULL) {
  fname <- config[[name]]
  if (length(fname) > 1) {
    load_image_data_series(fname, config, indices, mask_volume=mask_volume)
  } else if (!file.exists(fname)) {
    flog.error("datafile %s not found.", fname)
    stop()
  } else {
    flog.info("loading data file %s", fname)
    if (!is.null(indices)) {
      loadVector(fname, indices=indices, mask=mask_volume)
    } else {
      loadVector(fname, mask=mask_volume)
    }
    
  }
}

load_surface_mask <- function(masklist, train_surf) {
  sections <- names(train_surf)
  assert_that(all(names(sections) == names(masklist)))
  
  masksurfaces <- lapply(sections, function(section) {
    msurf <- neurosurf::loadSurface(train_surf[[section]]@geometry, masklist[[section]])
    flog.info("mask for %s has %s regions", section, length(unique(msurf@data)))
    msurf
  })
  
  names(masksurfaces) <- sections
  masksurfaces
}

load_surface_data <- function(config, name, nodeind=NULL, colind=NULL) {
  tdat <- config[[name]]
  sections <- names(tdat)
  
  flog.info("surface sections: ", sections, capture=TRUE)
  
  surfaces <- lapply(sections, function(section) {
    neurosurf::loadSurface(tdat[[section]]$geometry, tdat[[section]]$data, nodeind=nodeind, colind=colind)
  })
  
  names(surfaces) <- sections
  surfaces
    
}





  