
wrap_result <- function(result_table, dataset, predictor=NULL) {
  print(nrow(result_table))
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
    classification_result(observed, pclass, prob, dataset$design$test_design, predictor)
  } else {
      preds <- numeric(length(observed))
      for (i in seq_along(result_table$prediction)) {
        tind <- result_table$test_ind[[i]]
        preds[tind] <- result$prediction[[i]]$preds
      }
      
      counts <- table(unlist(result_table$test_ind))
      preds <- preds/counts
      regression_result(observed, preds, dataset$design$test_design, predictor)
  }
}


#' @importFrom dplyr rowwise, do
#' @importFrom tibble as_tibble
#' 
internal_crossval <- function(dset, roi, mspec, id, keep_predictor) {
  print(id)
  
  samples <- crossval_samples(mspec$crossval, tibble::as_tibble(values(roi$train_roi)), y_train(dset))
  ind <- indices(roi$train_roi)
  
  ret <- samples %>% dplyr::rowwise() %>% dplyr::do( {
    if (ncol(.$train) < 2) {
      data.frame()
    } else {
      result <- train_model(mspec, tibble::as_tibble(.$train), .$ytrain, indices=ind, param=mspec$tune_grid)
      pred <- predict(result, tibble::as_tibble(.$test), NULL)
      plist <- lapply(pred, list)
      plist$y_true <- list(.$ytest)
      plist$test_ind=list(as.integer(.$test))
      tibble::as_tibble(plist) 
    }
    
  })
  
  if (nrow(ret) == 0) {
    tibble::tibble()
  } else {
    cres <- wrap_result(ret,dset)
    #tibble::tibble(result=list(cres), indices=list(ind), performance=list(as_tibble(as.list(performance(cres)))))
    tibble::tibble(result=list(cres), indices=list(ind), performance=list(performance(cres)), id=id)
  }
}
  
  
#' @export
#' @importFrom dplyr do, rowwise
mvpa_iterate <- function(dset, vox_iter, mod_spec, ids=1:length(vox_iter)) {
  sframe <- get_samples(dset, vox_iter)

  ret <- sframe %>% dplyr::mutate(rnum=ids) %>% 
    dplyr::rowwise() %>% dplyr::do(internal_crossval(dset, as_roi(.$sample), mod_spec, .$rnum, FALSE))
  ## table(unlist(ret$indices))
  ret
  
}