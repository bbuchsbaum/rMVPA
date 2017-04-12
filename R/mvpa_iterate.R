
wrap_result <- function(result_table, dataset, fit=NULL) {
  observed <- y_test(dataset)
  
  if (is.factor(observed)) {
    prob <- matrix(0, length(observed), length(levels(observed)))
    colnames(prob) <- levels(observed)
  
    for (i in seq_along(result_table$probs)) {
      p <- as.matrix(result_table$probs[[i]])
      tind <- result_table$test_ind[[i]]
      prob[tind,] <- prob[tind,] + p
    }
  
    prob <- t(apply(prob, 1, function(vals) vals / sum(vals)))
    maxid <- apply(prob, 1, which.max)
    pclass <- levels(observed)[maxid]
    classification_result(observed, pclass, prob, dataset$design$test_design, fit)
  } else {
      preds <- numeric(length(observed))
      for (i in seq_along(result_table$preds)) {
        tind <- result_table$test_ind[[i]]
        preds[tind] <- result_table$preds[[i]]
      }
      
      counts <- table(unlist(result_table$test_ind))
      preds <- preds/counts
      regression_result(observed, preds, dataset$design$test_design, fit)
  }
}


#' extern_crossval
#' 
external_crossval <- function(roi, mspec, id, return_fit=FALSE) {
  xtrain <- tibble::as_tibble(values(roi$train_roi))
  
  if (ncol(xtrain) < 2) {
    return(tibble::tibble())
  }
  
  dset <- mspec$dataset
  
  ytrain <- y_train(dset)
  ytest <- y_test(dset)
  ind <- indices(roi$train_roi)
  
  result <- try(train_model(mspec, xtrain, ytrain, indices=ind, param=mspec$tune_grid))
  
  if (inherits(result, "try-error")) {
    warning(result)
    return(data.frame())
  }
  
  pred <- predict(result, tibble::as_tibble(values(roi$test_roi)), NULL)
  plist <- lapply(pred, list)
  plist$y_true <- list(ytest)
  plist$test_ind=list(as.integer(seq_along(ytest)))
  
  ret <- tibble::as_tibble(plist) 
  
  cres <- if (return_fit) {
    wrap_result(ret, dset, result$fit)
  } else {
    wrap_result(ret,dset)
  }
  
  tibble::tibble(result=list(cres), indices=list(ind), performance=list(compute_performance(mspec, cres)), id=id)
}
 


#' @importFrom dplyr rowwise do
#' @importFrom tibble as_tibble
internal_crossval <- function(roi, mspec, id, return_fit=FALSE) {
  
  samples <- crossval_samples(mspec$crossval, tibble::as_tibble(values(roi$train_roi)), y_train(mspec$dataset))
  ind <- indices(roi$train_roi)
  
  ret <- samples %>% dplyr::rowwise() %>% dplyr::do( {
    if (ncol(.$train) < 2) {
      data.frame()
    } else {
      result <- try(train_model(mspec, tibble::as_tibble(.$train), .$ytrain, indices=ind, param=mspec$tune_grid))
      
      if (inherits(result, "try-error")) {
        warning(result)
        data.frame()
      } else {
        pred <- predict(result, tibble::as_tibble(.$test), NULL)
        plist <- lapply(pred, list)
        plist$y_true <- list(.$ytest)
        plist$test_ind=list(as.integer(.$test))
        if (return_fit) {
          plist$fit <- list(result)
        }
      
        tibble::as_tibble(plist) 
      }
    }
    
  })
  
  if (nrow(ret) != nrow(samples)) {
    tibble::tibble()
  } else {
    cres <- if (return_fit) {
      predictor <- weighted_model(ret$fit)
      wrap_result(ret,mspec$dataset, predictor)
    } else {
      wrap_result(ret,mspec$dataset)
    }
    #tibble::tibble(result=list(cres), indices=list(ind), performance=list(as_tibble(as.list(performance(cres)))))
    tibble::tibble(result=list(cres), indices=list(ind), performance=list(compute_performance(mspec, cres)), id=id)
  }
}
  
  
#' mvpa_iterate
#' 
#' Fit a model for each of a list of voxels sets
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
    mutate(len=length(attr(sample[["vox"]], "indices"))) %>% 
    filter(len >= 2)
    dplyr::do(do_fun(as_roi(.$sample), mod_spec, .$rnum, return_fits))
  ret
  
}