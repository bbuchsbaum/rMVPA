
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
external_crossval <- function(roi, mspec, id, compute_performance=TRUE, return_fit=FALSE) {
  xtrain <- tibble::as_tibble(neuroim::values(roi$train_roi))
 
  dset <- mspec$dataset
  
  ytrain <- y_train(mspec)
  ytest <- y_test(mspec)
  
  ind <- neuroim::indices(roi$train_roi)
  
  result <- try_warning(train_model(mspec, xtrain, ytrain, indices=ind, param=mspec$tune_grid, tune_reps=mspec$tune_reps))
  
  if (!is.null(result$error)) {
    emessage <- attr(result, "condition")$message
    tibble::tibble(result=list(NULL), indices=list(ind), performance=list(NULL), id=id, 
                   error=TRUE, error_message=emessage, 
                   warning=!is.null(result$warning), warning_message=if (is.null(result$warning)) "~" else result$warning)
  } else {
   
    pred <- predict(result$value, tibble::as_tibble(neuroim::values(roi$test_roi)), NULL)
    plist <- lapply(pred, list)
    plist$y_true <- list(ytest)
    plist$test_ind=list(as.integer(seq_along(ytest)))
  
    ret <- tibble::as_tibble(plist) 
  
    cres <- if (return_fit) {
      wrap_result(ret, mspec$design, result$fit)
    } else {
      wrap_result(ret, mspec$design)
    }
  
    if (compute_performance) {
      tibble::tibble(result=list(cres), indices=list(ind), performance=list(compute_performance(mspec, cres)), id=id, 
                 error=FALSE, error_message="~", 
                 warning=!is.null(result$warning), warning_message=if (is.null(result$warning)) "~" else result$warning)
    } else {
      tibble::tibble(result=list(cres), indices=list(ind), performance=list(NULL), id=id, 
                     error=FALSE, error_message="~", 
                     warning=!is.null(result$warning), warning_message=if (is.null(result$warning)) "~" else result$warning)
    }
  }
}
 

#' @keywords internal
#' @importFrom dplyr rowwise do bind_rows
#' @importFrom tibble as_tibble
internal_crossval <- function(roi, mspec, id, compute_performance=TRUE, return_fit=FALSE) {
  
  ## generate cross-validation samples
  samples <- crossval_samples(mspec$crossval, tibble::as_tibble(values(roi$train_roi)), y_train(mspec))
  
  ## get ROI indices
  ind <- neuroim::indices(roi$train_roi)
  
  ret <- samples %>% pmap(function(ytrain, ytest, train, test, .id) {
      ## fit model for cross-validation sample
      result <- try(train_model(mspec, tibble::as_tibble(train), ytrain, 
                                indices=ind, param=mspec$tune_grid, 
                                tune_reps=mspec$tune_reps))
      
      if (inherits(result, "try-error")) {
        flog.warn("error fitting model %s : %s", id, attr(result, "condition")$message)
        ## error encountered, store error messages
        emessage <- attr(result, "condition")$message
        tibble::tibble(class=list(NULL), probs=list(NULL), y_true=list(ytest), 
                       fit=list(NULL), error=TRUE, error_message=emessage)
      } else {
        pred <- predict(result, tibble::as_tibble(test), NULL)
        plist <- lapply(pred, list)
        plist$y_true <- list(ytest)
        plist$test_ind <- list(as.integer(test))
        plist$fit <- if (return_fit) list(result) else list(NULL)
        plist$error <- FALSE
        plist$error_message <- "~"
        tibble::as_tibble(plist) 
      }
  }) %>% dplyr::bind_rows()
  

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
  
  
#' mvpa_iterate
#' 
#' Fit a classification/regression model for each voxel set in a list
#' 
#' @param mod_spec a class of type \code{mvpa_model}
#' @param vox_list a \code{list} of voxel indices/coordinates
#' @param ids a \code{vector} of ids for each voxel set
#' @param compute_performance compute and store performance measures for each voxel set
#' @param return_fits return the model fit for each voxel set?
#' @importFrom dplyr bind_rows
#' @export
mvpa_iterate <- function(mod_spec, vox_list, ids=1:length(vox_iter), compute_performance=TRUE, return_fits=FALSE) {
  assert_that(length(ids) == length(vox_list), 
              msg=paste("length(ids) = ", length(ids), "::", "length(vox_list) =", length(vox_list)))
  
  sframe <- get_samples(mod_spec$dataset, vox_list)

  do_fun <- if (has_test_set(mod_spec$dataset)) external_crossval else internal_crossval
  
  ### iterate over rows using parallel map with futures
  ret <- sframe %>% dplyr::mutate(rnum=ids) %>% furrr::future_pmap(function(sample, rnum, .id) {
    do_fun(as_roi(sample), mod_spec, .id, 
           compute_performance=compute_performance,
           return_fit=return_fits)
  }) %>% dplyr::bind_rows()
  
  ret
  
}

#' @importFrom Rfit rfit
run_rfit <- function(dvec, obj) {
  form <- paste("dvec", "~", paste(names(obj$design$model_mat), collapse = " + "))
  obj$design$model_mat$dvec <- dvec
  res <- Rfit::rfit(form, data=obj$design$model_mat)
  coef(res)[-1]
}


#' @keywords internal
run_lm <- function(dvec, obj) {
  form <- paste("dvec", "~", paste(names(obj$design$model_mat), collapse = " + "))
  vnames <- names(obj$design$model_mat)
  obj$design$model_mat$dvec <- dvec
 
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


#' @export
train_model.rsa_model <- function(obj, train_dat, indices, wts=NULL, method=c("lm", "rfit", "pearson", "spearman"), 
                                  distmethod=c("pearson", "spearman"), ...) {
  method <- match.arg(method)
  distmethod <- match.arg(distmethod)
  
  dtrain <- 1 - cor(t(train_dat), method=distmethod)
  dvec <- dtrain[lower.tri(dtrain)]
  
  if (! is.null(obj$design$include)) {
    dvec <- dvec[obj$design$include]
  }
  
  switch(method,
         rfit=run_rfit(dvec, obj),
         lm=run_lm(dvec,obj),
         pearson=run_cor(dvec,obj,"pearson"),
         spearman=run_cor(dvec,obj,"spearman"))
  
}


#' @importFrom neuroim indices values
do_rsa <- function(roi, mod_spec, rnum, method=method, distmethod=distmethod) {
  xtrain <- tibble::as_tibble(values(roi$train_roi))
  ind <- indices(roi$train_roi)
  ret <- train_model(mod_spec, xtrain, ind, method=method, distmethod=distmethod)
  tibble::tibble(result=list(NULL), indices=list(ind), performance=list(ret), id=rnum, error=FALSE, error_message="~")
}


#' rsa_iterate
#' 
#' Run RSA analysis for each of a list of voxels sets
#' 
#' @param mod_spec a class of type \code{mvpa_model}
#' @param vox_list a \code{list} of voxel indices/coordinates
#' @param ids a \code{vector} of ids for each voxel set
#' @param regtype the analysis method, one of: \code{lm}, \code{rfit}, \code{pearson}, \code{spearman}
#' @param distmethod the method used to computer distances between oservations, one of: \code{pearson}, \code{spearman}
#' @importFrom dplyr do rowwise
#' @export
rsa_iterate <- function(mod_spec, vox_list, ids=1:length(vox_iter), regtype=c("rfit", "lm", "pearson", "spearman"), 
                        distmethod=c("spearman", "pearson")) {
 
  distmethod <- match.arg(distmethod)
  regtype <- match.arg(regtype)
  
  print(paste("regtype", regtype))
  
  assert_that(length(ids) == length(vox_list), msg=paste("length(ids) = ", length(ids), "::", "length(vox_list) =", length(vox_list)))
  sframe <- get_samples(mod_spec$dataset, vox_list)
  
  ## iterate over searchlights using parallel futures
  ret <- sframe %>% dplyr::mutate(rnum=ids) %>% furrr::future_pmap(function(sample, rnum, .id) {
    do_rsa(as_roi(sample), mod_spec, rnum, regtype, distmethod)
  }) %>% dplyr::bind_rows()
}


