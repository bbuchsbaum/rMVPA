
#' searchlight
#' 
#' @export
searchlight <- function(dset, radius=8, model=MVPA.lda_thomaz, ncores=1) {
  if (radius < 1 || radius > 100) {
    stop(paste("radius", radius, "outside allowable range (1-100)"))
  }
  
  if (length(folds) != length(labels)) {
    stop(paste("length of 'labels' must equal length of 'folds'", length(labels), "!=", length(folds)))
  }
  
  modelInstance <- createModel(dset, model)
  
  searchIter <- ihasNext(Searchlight(mask, radius))
  
  while (hasNext(searchIter)) {
    vox <- nextElem(searchIter) 
    fitModel(model, dset, vox, ncores)
    
    
    
    
  
  }
  
  
  
  
}