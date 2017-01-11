
internal_crossval <- function(dset, roi, mspec) {
  param <- mspec$tune_grid[1,]
  
  samples <- crossval_samples(mspec$crossval, as.data.frame(values(roi$train_roi)), y_train(dset))
  ind <- indices(roi$train_roi)
  
  ret <- samples %>% rowwise() %>% do( {
    result <- train_model(mspec, as_data_frame(.$train), .$ytrain, indices=ind, param=param)
    print(class(result))
    #browser()
    ## keep predictor?
    pred <- predict(result, as_data_frame(.$test), NULL)
    plist <- lapply(pred, list)
    plist$y_true <- list(.$ytest)
    plist$test_ind=list(as.integer(.$test))
    tibble::as_tibble(plist) 
    
  })
}
  
  
#' @export
mvpa_iterate <- function(dset, vox_iter, mod_spec) {
  sframe <- get_samples(dset, vox_iter)
  sframe %>% rowwise() %>% do(mvpa_crossval(as_roi(.$sample), mod_spec))
  
}