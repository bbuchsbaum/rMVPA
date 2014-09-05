
#' searchlight
#' @param dataset
#' @param radius the searchlight radus in mm
#' @param modelName the name of the classifcation model to be used
#' @param ncores the number of cores for parallel processign (default is 1)
#' @importFrom itertools ihasNext
#' @export
searchlight <- function(dset, radius=8, modelName="svmLinear", ncores=1) {
  if (radius < 1 || radius > 100) {
    stop(paste("radius", radius, "outside allowable range (1-100)"))
  }
  
  if (length(dset$blockVar) != length(dset$Y)) {
    stop(paste("length of 'labels' must equal length of 'cross validation blocks'", length(dset$Y), "!=", length(dset$blockVar)))
  }
  
  model <- loadModel(modelName)
  print(model) 
  
  searchIter <- itertools::ihasNext(Searchlight(dset$mask, radius))
  count <- 1
  while (itertools::hasNext(searchIter)) {
    print(count)
    vox <- nextElem(searchIter) 
    result <- fitModel(model, dset, vox, ncores)
    res <- performance(result)
    print(res)
    count <- count + 1
  }
   
}