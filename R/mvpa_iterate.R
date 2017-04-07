
wrap_result <- function(result_table, dataset, predictor=NULL) {
  observed <- if (is.null(dataset$test_data)) {
    dataset$design$y_train
  } else {
    dataset$design$y_test
  }
  
  if (is.factor(observed)) {
    prob <- matrix(0, length(observed), length(levels(observed)))
    colnames(prob) <- levels(observed)
  
    for (i in seq_along(result_table$probs)) {
      p <- as.matrix(result_table$probs[[i]])
      tind <- result_table$test_ind[[i]]
      prob[tind,] <- prob[tind,] + p
    }
  
    ## if there are multiple predictions on the same sample, this average over predictions
    
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
      
      counts <- table(unlist(result_table$test_ind[[i]]))
      preds <- preds/counts
      regression_result(observed, preds, testDesign, predictor)
  }
}


#' @importFrom dplyr rowwise
internal_crossval <- function(dset, roi, mspec, keep_predictor) {
  param <- mspec$tune_grid[1,]
  
  samples <- crossval_samples(mspec$crossval, as.data.frame(values(roi$train_roi)), y_train(dset))
  ind <- indices(roi$train_roi)
  
  ret <- samples %>% rowwise() %>% do( {
    result <- train_model(mspec, as_data_frame(.$train), .$ytrain, indices=ind, param=param)
    pred <- predict(result, as_data_frame(.$test), NULL)
    plist <- lapply(pred, list)
    plist$y_true <- list(.$ytest)
    plist$test_ind=list(as.integer(.$test))
    tibble::as_tibble(plist) 
    
  })
  
  cres <- wrap_result(ret,dset)
  tibble::tibble(result=list(cres), indices=list(ind))
}
  
  
#' @export
mvpa_iterate <- function(dset, vox_iter, mod_spec) {
  sframe <- get_samples(dset, vox_iter)
  ret <- sframe %>% rowwise() %>% do(internal_crossval(dset, as_roi(.$sample), mod_spec))
  ## table(unlist(ret$indices))
  ret
  
}