
mvpa_crossval <- function(roi, mspec) {
  samples <- crossval_samples(mspec)
  samples %>% 
    train_model.model_spec <- function(obj, roi_train, Ytrain, param=NULL, wts=NULL)
}
  
  
#' @export
mvpa_iterate <- function(dset, vox_iter, mod_spec) {
  sframe <- get_samples(dset, vox_iter)
  sframe %>% rowwise() %>% do(mvpa_crossval(as_roi(.$sample), mod_spec))
  
}