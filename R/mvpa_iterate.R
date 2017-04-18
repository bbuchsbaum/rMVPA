
wrap_result <- function(result_table, design, fit=NULL) {

  observed <- y_test(design)
  
  if (is.factor(observed)) {
    prob <- matrix(0, length(observed), length(levels(observed)), dimnames=list(list(), levels(observed)))
    #colnames(prob) <- levels(observed)
  
    for (i in seq_along(result_table$probs)) {
      p <- as.matrix(result_table$probs[[i]])
      tind <- result_table$test_ind[[i]]
      prob[tind,] <- prob[tind,] + p
    }
  
    prob <- t(apply(prob, 1, function(vals) vals / sum(vals)))
    maxid <- apply(prob, 1, which.max)
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
  
  if (ncol(xtrain) < 2) {
    return(tibble::tibble())
  }
  
  dset <- mspec$dataset
  
  ytrain <- y_train(mspec)
  ytest <- y_test(mspec)
  ind <- indices(roi$train_roi)
  
  result <- try(train_model(mspec, xtrain, ytrain, indices=ind, param=mspec$tune_grid))
  
  if (inherits(result, "try-error")) {
    tibble::tibble(result=list(), indices=list(ind), performance=list(), id=id, 
                   error=TRUE, error_message=attr(result, "condition")$message)
  } else {
  
    pred <- predict(result, tibble::as_tibble(values(roi$test_roi)), NULL)
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
                 error=FALSE, error_message="~")
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
      result <- try(train_model(mspec, tibble::as_tibble(.$train), .$ytrain, indices=ind, param=mspec$tune_grid))
      
      if (inherits(result, "try-error")) {
        flog.debug("error: %s", attr(result, "condition")$message)
        ## error encountered, store error messages
        tibble(class=list(NULL), probs=list(NULL), y_true=list(.$ytest), fit=list(NULL), error=TRUE, error_message=attr(result, "condition")$message)
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
  assert_that(length(ids) == length(vox_list))
  
  sframe <- get_samples(mod_spec$dataset, vox_list)

  do_fun <- if (has_test_set(mod_spec$dataset)) external_crossval else internal_crossval

  ret <- sframe %>% dplyr::mutate(rnum=ids) %>% 
    dplyr::rowwise() %>% 
    dplyr::do(do_fun(as_roi(.$sample), mod_spec, .$rnum, return_fits))
  
  ret
  
}