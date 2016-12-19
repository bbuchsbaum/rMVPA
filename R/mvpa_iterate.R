
mvpa_crossval <- function(dset, roi, mspec) {
  param <- mspec$tune_grid[1,]
  
  samples <- crossval_samples(mspec$crossval, as.data.frame(values(roi$train_roi)), y_train(dset))
  ind <- indices(roi$train_roi)
  
  ret <- samples %>% rowwise() %>% do( {
    browser()
    result <- train_model(mspec, as_data_frame(.$train), .$ytrain, indices=ind, param=param)
    print(class(result))
   
    tibble::as_tibble(predict(result, as_data_frame(.$test), NULL))
  })
}
  
  
#' @export
mvpa_iterate <- function(dset, vox_iter, mod_spec) {
  sframe <- get_samples(dset, vox_iter)
  sframe %>% rowwise() %>% do(mvpa_crossval(as_roi(.$sample), mod_spec))
  
}