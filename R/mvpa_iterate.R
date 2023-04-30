
#' @keywords internal
try_warning  <- function(expr) {
  warn <- err <- NULL
  value <- withCallingHandlers(
    tryCatch(expr, error=function(e) {
      err <<- e
      NULL
    }), warning=function(w) {
      warn <<- paste0(warn, str_trim(as.character(w)), sep=";")
      invokeRestart("muffleWarning")
    })
  list(value=value, warning=warn, error=err)
}


#' @keywords internal
wrap_result <- function(result_table, design, fit=NULL) {
 
  observed <- y_test(design)
  
  ## It could happen that not all design rows are actually tested, which is why we find the unqiue set of test indices
  testind <- unique(sort(unlist(result_table$test_ind)))
  
  if (is.factor(observed)) {
    prob <- matrix(0, length(testind), length(levels(observed)))
    colnames(prob) <- levels(observed)
  
    for (i in seq_along(result_table$probs)) {
      p <- as.matrix(result_table$probs[[i]])
      tind <- match(result_table$test_ind[[i]], testind)
      prob[tind,] <- prob[tind,] + p
    }
    
    ## probs must sum to one, can divide by sum.
    prob <- t(apply(prob, 1, function(vals) vals / sum(vals)))
    maxid <- max.col(prob)
    pclass <- levels(observed)[maxid]
  
    ## storing observed, testind, test_design 
    classification_result(observed[testind], pclass, prob, testind=testind, design$test_design, fit)
  } else {
    
      testind <- unique(sort(unlist(result_table$test_ind)))
      preds <- numeric(length(testind))
      
      for (i in seq_along(result_table$preds)) {
        #tind <- result_table$test_ind[[i]]
        tind <- match(result_table$test_ind[[i]], testind)
        preds[tind] <- preds[tind] + result_table$preds[[i]]
      }
      
      ## TODO check me
      counts <- table(sort(unlist(result_table$test_ind)))
      preds <- preds/counts
      regression_result(observed, preds, testind=testind, test_design=design$test_design, fit)
  }
}


#' external_crossval
#' @keywords internal
#' @importFrom stats predict
external_crossval <- function(roi, mspec, id, compute_performance=TRUE, return_fit=FALSE, permute=FALSE) {
  xtrain <- tibble::as_tibble(neuroim2::values(roi$train_roi), .name_repair=.name_repair)
 
  #dset <- mspec$dataset
  
  ytrain <- if (permute) {
    sample(y_train(mspec))
  } else {
    y_train(mspec)
  }
  
  ytest <- y_test(mspec)
  
  ind <- neuroim2::indices(roi$train_roi)
  
  #result <- try_warning(train_model(mspec, xtrain, ytrain, indices=ind, param=mspec$tune_grid, tune_reps=mspec$tune_reps))
  result <- try(train_model(mspec, xtrain, ytrain, indices=ind, 
                            param=mspec$tune_grid, 
                            tune_reps=mspec$tune_reps))
  
  if (inherits(result, "try-error")) {
    flog.warn("error fitting model %s : %s", id, attr(result, "condition")$message)
    ## error encountered, store error messages
    emessage <- if (is.null(attr(result, "condition")$message)) "" else attr(result, "condition")$message
    tibble::tibble(class=list(NULL), probs=list(NULL), y_true=list(ytest), 
                   fit=list(NULL), error=TRUE, error_message=emessage)
  } else {
  # if (!is.null(result$error)) {
  #   
  #   emessage <- if (!is.null(attr(result$error, "condition")$message)) {
  #     attr(result, "condition")$message
  #   } else {
  #     result$error
  #   }
  #   tibble::tibble(result=list(NULL), indices=list(ind), performance=list(NULL), id=id, 
  #                  error=TRUE, error_message=emessage, 
  #                  warning=!is.null(result$warning), 
  #                  warning_message=if (is.null(result$warning)) "~" else result$warning)
  # } else {
    pred <- predict(result, tibble::as_tibble(neuroim2::values(roi$test_roi), .name_repair=.name_repair), NULL)
    plist <- lapply(pred, list)
    plist$y_true <- list(ytest)
    plist$test_ind=list(as.integer(seq_along(ytest)))
  
    ret <- tibble::as_tibble(plist, .name_repair = .name_repair) 
  
    cres <- if (return_fit) {
      wrap_result(ret, mspec$design, result$fit)
    } else {
      wrap_result(ret, mspec$design)
    }
  
    if (compute_performance) {
      tibble::tibble(result=list(cres), indices=list(ind), 
                     performance=list(compute_performance(mspec, cres)), id=id, 
                 error=FALSE, error_message="~", 
                 warning=!is.null(result$warning), 
                 warning_message=if (is.null(result$warning)) "~" else result$warning)
    } else {
      tibble::tibble(result=list(cres), indices=list(ind), performance=list(NULL), id=id, 
                     error=FALSE, error_message="~", 
                     warning=!is.null(result$warning), 
                     warning_message=if (is.null(result$warning)) "~" else result$warning)
    }
  
  }
}
 

#' @keywords internal
#' @importFrom dplyr rowwise do bind_rows
#' @importFrom tibble as_tibble
internal_crossval <- function(roi, mspec, id, compute_performance=TRUE, return_fit=FALSE, permute=FALSE) {
  ## generate cross-validation samples
  ## this could be done outside the function??
  ## crossval_samples should really cache the indices rather than regenerate every iteration
  ## for large class sizes, this can create big objects for every voxel.

  
  samples <- if (!permute) {
    crossval_samples(mspec$crossval, tibble::as_tibble(neuroim2::values(roi$train_roi), .name_repair=.name_repair), y_train(mspec))
  } else {
    crossval_samples(mspec$crossval, tibble::as_tibble(neuroim2::values(roi$train_roi), .name_repair=.name_repair), sample(y_train(mspec)))
  }
  
  ## get ROI indices
  ind <- neuroim2::indices(roi$train_roi)
  
  ret <- samples %>% pmap(function(ytrain, ytest, train, test, .id) {
    ## fit model for cross-validation sample
    ##if (ncol(train) < 2) {
    ##  return(NULL)
    ##}  
    
    result <- try(train_model(mspec, tibble::as_tibble(train, .name_repair=.name_repair), ytrain, 
                                indices=ind, param=mspec$tune_grid, 
                                tune_reps=mspec$tune_reps))
      
      if (inherits(result, "try-error")) {
        flog.warn("error fitting model %s : %s", id, attr(result, "condition")$message)
        ## error encountered, store error messages
        emessage <- if (is.null(attr(result, "condition")$message)) "" else attr(result, "condition")$message
        tibble::tibble(class=list(NULL), probs=list(NULL), y_true=list(ytest), 
                       fit=list(NULL), error=TRUE, error_message=emessage)
      } else {
        pred <- predict(result, tibble::as_tibble(test, .name_repair=.name_repair), NULL)
        plist <- lapply(pred, list)
        plist$y_true <- list(ytest)
        plist$test_ind <- list(as.integer(test))
        plist$fit <- if (return_fit) list(result) else list(NULL)
        plist$error <- FALSE
        plist$error_message <- "~"
        tibble::as_tibble(plist, .name_repair=.name_repair) 
      }
  }) %>% purrr::discard(is.null) %>% dplyr::bind_rows()
  
  if (any(ret$error)) {
    emessage <- ret$error_message[which(ret$error)[1]]
    tibble::tibble(result=list(NULL), indices=list(ind), performance=list(NULL), 
                   error=TRUE, error_message=emessage)
  } else {
    cres <- if (return_fit) {
      predictor <- weighted_model(ret$fit)
      wrap_result(ret, mspec$design, predictor)
    } else {
      wrap_result(ret, mspec$design)
    }
    
    if (compute_performance) {
      tibble::tibble(result=list(cres), indices=list(ind), 
                   performance=list(compute_performance(mspec, cres)), 
                   id=id, error=FALSE, error_message="~")
    } else {
      tibble::tibble(result=list(cres), indices=list(ind), 
                     performance=list(NULL), 
                     id=id, error=FALSE, error_message="~")
    }
  }
  
}

#' @keywords internal
extract_roi <- function(sample, data) {
  r <- as_roi(sample,data)
  v <- neuroim2::values(r$train_roi)
  r <- try(filter_roi(r))
  if (inherits(r, "try-error") || ncol(v) < 2) {
    NULL
  } else {
    r
  }
}
  
  
#' MVPA Iteration for Voxel Sets with Parallelization
#'
#' This function fits a classification or regression model for each voxel set in a list using parallelization.
#' It can compute and store performance measures, return row-wise predictions, and return the model fit for each voxel set.
#'
#' @param mod_spec An object of class \code{mvpa_model} specifying the model.
#' @param vox_list A \code{list} of voxel indices or coordinates.
#' @param ids A \code{vector} of IDs for each voxel set (defaults to 1:length(vox_list)).
#' @param compute_performance A \code{logical} indicating whether to compute and store performance measures for each voxel set (defaults to TRUE).
#' @param return_predictions A \code{logical} indicating whether to return row-wise predictions for each voxel set (defaults to TRUE).
#' @param return_fits A \code{logical} indicating whether to return the model fit for each voxel set (defaults to FALSE).
#' @param batch_size An \code{integer} specifying the number of voxel sets to process in each batch (defaults to 10% of the total voxel sets).
#' @param permute A \code{logical} indicating whether to permute the labels (defaults to FALSE).
#' @param verbose A \code{logical} indicating whether to print progress messages (defaults to TRUE).
#'
#' @return A \code{data.frame} containing the results for each voxel set, including performance measures, predictions, and model fits, as specified by the input parameters.
#'
#' @details
#' This function utilizes parallel processing to speed up the process of fitting the specified model for each voxel set in a list.
#' The parallelization is achieved using the `furrr` package, which provides a parallel backend for the `purrr` package.
#' By default, it divides the voxel sets into batches and processes them in parallel according to the specified batch size.
#' The function provides options to control the return of performance measures, predictions, and model fits for each voxel set.
#'
#' @importFrom dplyr bind_rows
#' @importFrom furrr future_pmap
#' @importFrom purrr map
#' @export
mvpa_iterate <- function(mod_spec, vox_list, ids=1:length(vox_list), 
                         compute_performance=TRUE, 
                         return_predictions=TRUE,
                         return_fits=FALSE, 
                         batch_size=as.integer(.1*length(ids)),
                         permute=FALSE, verbose=TRUE) {

  assert_that(length(ids) == length(vox_list), 
              msg=paste("length(ids) = ", length(ids), "::", "length(vox_list) =", length(vox_list)))
  
  
  batch_size <- max(1, batch_size)
  nbatches <- as.integer(length(ids)/batch_size)
  batch_group <- sort(rep(1:nbatches, length.out=length(ids)))
  batch_ids <- split(1:length(ids), batch_group)
  rnums <- split(ids, batch_group)
  
  dset <- mod_spec$dataset
  tot <- length(ids)
  
  do_fun <- if (has_test_set(dset)) external_crossval else internal_crossval
  
  result <- purrr::map(1:length(batch_ids), function(i) {
    futile.logger::flog.info("mvpa_iterate: compute analysis for batch %s ...", i)
    sf <- get_samples(dset, vox_list[batch_ids[[i]]]) %>% mutate(.id=batch_ids[[i]], rnum=rnums[[i]])
    sf <- sf %>% rowwise() %>% mutate(roi=list(extract_roi(sample,dset))) %>% select(-sample)
    fut_mvpa(mod_spec, sf, tot, do_fun, verbose, permute, compute_performance, return_fits, return_predictions)
  }) %>% bind_rows()
  
  
  result
  
}

#' @keywords internal
#' @noRd
fut_mvpa <- function(mod_spec, sf, tot, do_fun, verbose, permute, compute_performance, return_fits, return_predictions) {
  mod_spec$dataset <- NULL
  gc()
  sf %>% furrr::future_pmap(function(.id, rnum, roi) {
    
    if (verbose && (as.numeric(.id) %% 100 == 0)) {
      perc <- as.integer(as.numeric(.id)/tot * 100)
      futile.logger::flog.info("mvpa_iterate: %s percent", perc)
    }
    
    result <- do_fun(roi, mod_spec, rnum, 
                     compute_performance=compute_performance,
                     return_fit=return_fits, permute=permute)
    
    if (!return_predictions) {
      result <- result %>% mutate(result = list(NULL))
    }
    
    result
  }) %>% purrr::discard(is.null) %>% dplyr::bind_rows()

}

#' @keywords internal
#' @importFrom Rfit rfit
run_rfit <- function(dvec, obj) {
  form <- paste("dvec", "~", paste(names(obj$design$model_mat), collapse = " + "))
  obj$design$model_mat$dvec <- dvec
  res <- Rfit::rfit(form, data=obj$design$model_mat)
  coef(res)[-1]
}


#' @keywords internal
#' @importFrom stats coef cor dist rnorm terms lm sd
run_lm <- function(dvec, obj) {
  form <- paste("dvec", "~", paste(names(obj$design$model_mat), collapse = " + "))
  vnames <- names(obj$design$model_mat)
  obj$design$model_mat$dvec <- dvec
 
  #browser()
  res <- lm(form, data=obj$design$model_mat)
  res <- coef(summary(res))[-1,3]
  names(res) <- vnames
  res
}

#' @keywords internal
run_cor <- function(dvec, obj, method) {
  res <- sapply(obj$design$model_mat, function(x) cor(dvec, x, method=method))
  names(res) <- names(obj$design$model_mat)
  res
}


#' Train an RSA Model
#'
#' This function trains an RSA (representational similarity analysis) model using the specified method and distance calculation.
#'
#' @param obj An object of class \code{rsa_model}.
#' @param train_dat The training data.
#' @param indices The indices of the training data.
#' @param wts Optional, the weights for the model training.
#' @param method The method used for model training. One of "lm" (linear regression), "rfit" (robust regression), "pearson" (Pearson correlation), or "spearman" (Spearman correlation). Default is "lm".
#' @param distmethod The method used for distance calculation. One of "pearson" (Pearson correlation) or "spearman" (Spearman correlation). Default is "pearson".
#' @param ... Additional arguments passed to the training method.
#' @return The trained model.
train_model.rsa_model <- function(obj, train_dat, indices, wts=NULL, method=c("lm", "rfit", "pearson", "spearman"), 
                                  distmethod=c("pearson", "spearman"), ...) {
  method <- match.arg(method)
  distmethod <- match.arg(distmethod)
  
  dtrain <- 1 - cor(t(train_dat), method=distmethod)
  dvec <- dtrain[lower.tri(dtrain)]
  
  if (!is.null(obj$design$include)) {
    dvec <- dvec[obj$design$include]
  }
  
  switch(method,
         rfit=run_rfit(dvec, obj),
         lm=run_lm(dvec,obj),
         pearson=run_cor(dvec,obj,"pearson"),
         spearman=run_cor(dvec,obj,"spearman"))
  
}


#' @importFrom neuroim2 indices values
do_rsa <- function(roi, mod_spec, rnum, method, distmethod) {
  xtrain <- tibble::as_tibble(neuroim2::values(roi$train_roi), .name_repair=.name_repair)
  ind <- indices(roi$train_roi)
  ret <- try(train_model(mod_spec, xtrain, ind, method=method, distmethod=distmethod))
  if (inherits(ret, "try-error")) {
    tibble::tibble(result=list(NULL), indices=list(ind), performance=list(ret), id=rnum, error=TRUE, error_message=attr(ret, "condition")$message)
  } else {
    tibble::tibble(result=list(NULL), indices=list(ind), performance=list(ret), id=rnum, error=FALSE, error_message="~")
  }
}


#' rsa_iterate
#'
#' Runs representational similarity analysis (RSA) for each voxel set in a list.
#'
#' @param mod_spec An object of class \code{rsa_model} specifying the RSA model.
#' @param vox_list A \code{list} of voxel indices or coordinates for each voxel set.
#' @param ids A \code{vector} of IDs for each voxel set (defaults to 1:length(vox_list)).
#' @param permute Logical, whether to permute the labels (defaults to FALSE).
#' @param regtype A character string specifying the analysis method. One of: \code{"pearson"}, \code{"spearman"}, \code{"lm"}, or \code{"rfit"} (defaults to "pearson").
#' @param distmethod A character string specifying the method used to compute distances between observations. One of: \code{"pearson"} or \code{"spearman"} (defaults to "spearman").
#' @importFrom dplyr do rowwise
#' @export
#' @inheritParams mvpa_iterate
rsa_iterate <- function(mod_spec, vox_list, ids=1:length(vox_list),  permute=FALSE, regtype=c("pearson", "spearman", "lm", "rfit"), 
                        distmethod=c("spearman", "pearson")) {
 
  distmethod <- match.arg(distmethod)
  regtype <- match.arg(regtype)
  
  message("regtype is:", regtype)
  
  assert_that(length(ids) == length(vox_list), msg=paste("length(ids) = ", length(ids), "::", "length(vox_list) =", length(vox_list)))
  sframe <- get_samples(mod_spec$dataset, vox_list)
  
  ## iterate over searchlights using parallel futures
  sf <- sframe %>% dplyr::mutate(rnum=ids) 
  fut_rsa(mod_spec,sf)
}


#' @keywords internal
fut_rsa <- function(mod_spec, sf) {
  mod_spec$dataset <- NULL
  gc()
  sf %>% dplyr::mutate(rnum=ids) %>% furrr::future_pmap(function(sample, rnum, .id) {
    ## extract_roi?
    do_rsa(as_roi(sample, mod_spec$dataset), mod_spec, rnum, method=regtype, distmethod=distmethod)
  }) %>% dplyr::bind_rows()
  
  
}


#' Train a MANOVA Model
#'
#' This function trains a multivariate analysis of variance (MANOVA) model using the specified design.
#'
#' @param obj An object of class \code{manova_model}.
#' @param train_dat The training data.
#' @param indices The indices of the training data.
#' @param ... Additional arguments passed to the training method.
#' @return A named numeric vector of -log(p-values) for each predictor in the MANOVA model.
#' @importFrom stats as.formula
train_model.manova_model <- function(obj, train_dat, indices, ...) {
  dframe <- obj$design$data
  dframe$response <- as.matrix(train_dat)
  form <- stats::as.formula(paste("response", paste(as.character(obj$design$formula), collapse='')))
 
  fres=ffmanova(form, data=dframe)
  pvals=fres$pValues
  names(pvals) <- sanitize(names(pvals))   
  lpvals <- -log(pvals)
  lpvals
}



#' @importFrom neuroim2 indices values
#' @keywords internal
do_manova <- function(roi, mod_spec, rnum) {
  xtrain <- tibble::as_tibble(neuroim2::values(roi$train_roi), .name_repair=.name_repair)
  ind <- indices(roi$train_roi)
  ret <- try(train_model(mod_spec, xtrain, ind))
  if (inherits(ret, "try-error")) {
    flog.warn("error fitting model %s : %s", rnum, attr(ret, "condition")$message)
    ## error encountered, store error messages
    emessage <- if (is.null(attr(ret, "condition")$message)) "" else attr(ret, "condition")$message
    tibble::tibble(result=list(NULL), indices=list(ind), performance=list(NULL), 
                   id=rnum, error=TRUE, error_message=emessage)
  } else {
      tibble::tibble(result=list(NULL), indices=list(ind), performance=list(ret), 
                     id=rnum, error=FALSE, error_message="~")
  }
}

#' MANOVA Iteration for Voxel Sets
#'
#' This function runs a MANOVA analysis for each of a list of voxel sets.
#'
#' @param mod_spec A \code{mvpa_model} object representing the model specification.
#' @param vox_list A \code{list} of voxel indices or coordinates.
#' @param ids A \code{vector} of IDs for each voxel set.
#' @param batch_size An \code{integer} specifying the batch size for processing (default is 10% of the total number of IDs).
#' @param permute A \code{logical} indicating whether to permute the voxel sets (default is FALSE).
#'
#' @return A \code{data.frame} containing the MANOVA results for each voxel set.
#'
#' @importFrom dplyr do rowwise
#' @export
manova_iterate <- function(mod_spec, vox_list, ids=1:length(vox_list),   batch_size=as.integer(.1*length(ids)), permute=FALSE) {
  assert_that(length(ids) == length(vox_list), msg=paste("length(ids) = ", length(ids), "::", "length(vox_list) =", length(vox_list)))
  futile.logger::flog.info("manova_iterate: extracting voxel ids")

  batch_size <- max(1, batch_size)
  nbatches <- as.integer(length(ids)/batch_size)
  batch_group <- sort(rep(1:nbatches, length.out=length(ids)))
  batch_ids <- split(1:length(ids), batch_group)
  rnums <- split(ids, batch_group)

  dset <- mod_spec$dataset
  ##mod_spec$dataset <- NULL

  tot <- length(ids)
  result <- purrr::map(1:length(batch_ids), function(i) {
    futile.logger::flog.info("manova_iterate: compute manovas ...")
    sf <- get_samples(mod_spec$dataset, vox_list[batch_ids[[i]]]) %>% mutate(.id=batch_ids[[i]], rnum=rnums[[i]])
    sf <- sf %>% rowwise() %>% mutate(roi=list(extract_roi(sample,dset))) %>% select(-sample)
    fut_manova(mod_spec, sf)
  }) %>% bind_rows()
  
  result

}

#'@keywords internal
fut_manova <- function(mod_spec, sf) {
  mod_spec$dataset <- NULL
  gc()
  sf %>% furrr::future_pmap(function(.id, rnum, roi) {
    do_manova(roi, mod_spec, rnum)
  }) %>% purrr::discard(is.null) %>% dplyr::bind_rows()
}


