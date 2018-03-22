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

wrap_result <- function(result_table, design, fit=NULL) {

  observed <- y_test(design)
  
  if (is.factor(observed)) {
    prob <- matrix(0, length(observed), length(levels(observed)))
    colnames(prob) <- levels(observed)
  
    for (i in seq_along(result_table$probs)) {
      p <- as.matrix(result_table$probs[[i]])
      tind <- result_table$test_ind[[i]]
      prob[tind,] <- prob[tind,] + p
    }
  
    prob <- t(apply(prob, 1, function(vals) vals / sum(vals)))
    maxid <- max.col(prob)
    pclass <- levels(observed)[maxid]
  
    classification_result(observed, pclass, prob, design$test_design, fit)
  } else {
      preds <- numeric(length(observed))
      for (i in seq_along(result_table$preds)) {
        tind <- result_table$test_ind[[i]]
        preds[tind] <- result_table$preds[[i]]
      }
      
      counts <- table(unlist(result_table$test_ind))
      preds <- preds/counts
      regression_result(observed, preds, design$test_design, fit)
  }
}


#' external_crossval
external_crossval <- function(roi, mspec, id, return_fit=FALSE) {
  xtrain <- tibble::as_tibble(values(roi$train_roi))
 
  dset <- mspec$dataset
  
  ytrain <- y_train(mspec)
  ytest <- y_test(mspec)
  
  ind <- indices(roi$train_roi)
  
  result <- try_warning(train_model(mspec, xtrain, ytrain, indices=ind, param=mspec$tune_grid, tune_reps=mspec$tune_reps))
  
  if (!is.null(result$error)) {
    emessage <- attr(result, "condition")$message
    tibble::tibble(result=list(NULL), indices=list(ind), performance=list(NULL), id=id, 
                   error=TRUE, error_message=emessage, 
                   warning=!is.null(result$warning), warning_message=if (is.null(result$warning)) "~" else result$warning)
  } else {
   
    pred <- predict(result$value, tibble::as_tibble(values(roi$test_roi)), NULL)
    plist <- lapply(pred, list)
    plist$y_true <- list(ytest)
    plist$test_ind=list(as.integer(seq_along(ytest)))
  
    ret <- tibble::as_tibble(plist) 
  
    cres <- if (return_fit) {
      wrap_result(ret, mspec$design, result$fit)
    } else {
      wrap_result(ret, mspec$design)
    }
  
    tibble::tibble(result=list(cres), indices=list(ind), performance=list(compute_performance(mspec, cres)), id=id, 
                 error=FALSE, error_message="~", 
                 warning=!is.null(result$warning), warning_message=if (is.null(result$warning)) "~" else result$warning)
  }
}
 


#' @importFrom dplyr rowwise do
#' @importFrom tibble as_tibble
internal_crossval <- function(roi, mspec, id, return_fit=FALSE) {
  
  
  ## generate cross-validation samples
  samples <- crossval_samples(mspec$crossval, tibble::as_tibble(values(roi$train_roi)), y_train(mspec))
  
  ## get ROI indices
  ind <- indices(roi$train_roi)
  
  
  ret <- samples %>% dplyr::rowwise() %>% dplyr::do( {
    
      ## fit model for cross-validation sample
      result <- try(train_model(mspec, tibble::as_tibble(.$train), .$ytrain, indices=ind, param=mspec$tune_grid, tune_reps=mspec$tune_reps))
      
      if (inherits(result, "try-error")) {
        flog.warn("error fitting model %s : %s", id, attr(result, "condition")$message)
        ## error encountered, store error messages
        emessage <- attr(result, "condition")$message
        tibble::tibble(class=list(NULL), probs=list(NULL), y_true=list(.$ytest), fit=list(NULL), error=TRUE, error_message=emessage)
      } else {
        pred <- predict(result, tibble::as_tibble(.$test), NULL)
        plist <- lapply(pred, list)
        plist$y_true <- list(.$ytest)
        plist$test_ind <- list(as.integer(.$test))
        plist$fit <- if (return_fit) list(result) else list(NULL)
        plist$error <- FALSE
        plist$error_message <- "~"
        tibble::as_tibble(plist) 
      }
  })
  

  if (any(ret$error)) {
    emessage <- ret$error_message[which(ret$error)[1]]
    tibble::tibble(result=list(NULL), indices=list(ind), performance=list(NULL), error=TRUE, error_message=emessage)
  } else {
    cres <- if (return_fit) {
      predictor <- weighted_model(ret$fit)
      wrap_result(ret, mspec$design, predictor)
    } else {
      wrap_result(ret, mspec$design)
    }
    
    tibble::tibble(result=list(cres), indices=list(ind), performance=list(compute_performance(mspec, cres)), id=id, error=FALSE, error_message="~")
  }
}
  
  
#' mvpa_iterate
#' 
#' Fit a classification/regression model for each of a list of voxels sets
#' 
#' @param mod_spec a class of type \code{mvpa_model}
#' @param vox_list a \code{list} of voxel indices/coordinates
#' @param ids a \code{vector} of ids for each voxel set
#' @param return_fits return the model fit for each voxel set?
#' @importFrom dplyr do rowwise
#' @export
mvpa_iterate <- function(mod_spec, vox_list, ids=1:length(vox_iter), return_fits=FALSE) {
  assert_that(length(ids) == length(vox_list), msg=paste("length(ids) = ", length(ids), "::", "length(vox_list) =", length(vox_list)))
  
  
  sframe <- get_samples(mod_spec$dataset, vox_list)

  do_fun <- if (has_test_set(mod_spec$dataset)) external_crossval else internal_crossval

  ret <- sframe %>% dplyr::mutate(rnum=ids) %>% 
    dplyr::rowwise() %>% 
    #do(function(x) { flog.info("mvpa_iterate: %s ", .$rnum); x }) %>%
    dplyr::do(do_fun(as_roi(.$sample), mod_spec, .$rnum, return_fits))
  
  ret
  
}

#' @importFrom Rfit rfit
run_rfit <- function(dvec, obj) {
  form <- paste("dvec", "~", paste(names(obj$design$model_mat), collapse = " + "))
  obj$design$model_mat$dvec <- dvec
  res <- Rfit::rfit(form, data=obj$design$model_mat)
  coef(res)[-1]
}

run_lm <- function(dvec, obj) {
  form <- paste("dvec", "~", paste(names(obj$design$model_mat), collapse = " + "))
  vnames <- names(obj$design$model_mat)
  obj$design$model_mat$dvec <- dvec
 
  res <- lm(form, data=obj$design$model_mat)
  res <- coef(summary(res))[-1,3]
  names(res) <- vnames
  res
}

run_cor <- function(dvec, obj, method) {
  res <- sapply(obj$design$model_mat, function(x) cor(dvec, x, method=method))
  names(res) <- names(obj$design$model_mat)
  res
}


#' @export
train_model.rsa_model <- function(obj, train_dat, indices, wts=NULL, method=c("lm", "rfit", "pearson", "spearman"), 
                                  distmethod=c("pearson", "spearman")) {
  method <- match.arg(method)
  distmethod <- match.arg(distmethod)
  print(paste("train_model rsa method", method))
  
  
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



do_rsa <- function(roi, mod_spec, rnum, method=method, distmethod=distmethod) {
  print(paste("rsa method", method))
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
  
  ret <- sframe %>% dplyr::mutate(rnum=ids) %>% 
    dplyr::rowwise() %>% 
    dplyr::do(do_rsa(as_roi(.$sample), mod_spec, .$rnum, regtype, distmethod))
}


