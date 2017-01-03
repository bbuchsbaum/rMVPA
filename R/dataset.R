
roi_volume_matrix <- function(mat, refspace, indices, coords) {
  structure(mat,
            refspace=refspace,
            indices=indices,
            coords=coords,
            class=c("roi_volume_matrix", "matrix"))
            
}

roi_surface_matrix <- function(mat, refspace, indices, coords) {
  structure(mat,
            refspace=refspace,
            indices=indices,
            coords=coords,
            class=c("roi_surface_matrix", "matrix"))
  
}

has_test_set.mvpa_dataset <- function(obj) {
  !is.null(obj$design$y_test) 
}

y_train.mvpa_dataset <- function(obj) y_train(obj$design)

y_test.mvpa_dataset <- function(obj) y_test(obj$design)


#' @param train_data
#' @param test_data
#' @param mask
#' @param design
#' @importFrom assertthat assert_that
mvpa_dataset <- function(train_data,test_data=NULL, mask, design) {
  assert_that(inherits(design, "mvpa_design"))
  
  ret <- list(
    train_data=train_data,
    test_data=test_data,
    mask=mask,
    design=design
  )
  
  class(ret) <- c("mvpa_dataset", "list")
  ret
    
}

# 
# #' @importFrom assertthat assert_that
# mvpa_external_dataset <- function(train_data, test_data, mask, mvpades) {
#   assert_that(is(mvpades, "mvpa_external_design"))
#   ret <- list(
#     trainVec=trainVec,
#     testVec=testVec,
#     mask=mask,
#     design=mvpades
#   )
#   
#   class(ret) <- c("mvpa_internal_dataset", "mvpa_dataset", "list")
#   ret
#   
# }




# #' @export
# #' @import R6
# #' @import neuroim
# MVPADataset <- R6::R6Class("MVPADataset",
#                            public = list(
#                              trainVec=NULL,Y=NULL,mask=NULL, 
#                              blockVar=NULL,testVec=NULL,testY=NULL,
#                              parcellation=NULL,testSplitVar=NULL,testSplits=NULL, 
#                              trainDesign=NULL, testDesign=NULL,
#                              initialize=function(trainVec,
#                                                  Y,
#                                                  mask=NULL, 
#                                                  blockVar=NULL,
#                                                  testVec=NULL,
#                                                  testY=NULL,
#                                                  parcellation=NULL, 
#                                                  testSplitVar=NULL,
#                                                  testSplits=NULL, 
#                                                  trainDesign=NULL, 
#                                                  testDesign=NULL) {
#                                
#                                self$trainVec=trainVec
#                                self$Y=Y
#                                self$mask=mask
#                                self$blockVar=blockVar
#                                self$testVec=testVec
#                                self$testY=testY
#                                self$parcellation=parcellation
#                                self$testSplitVar=testSplitVar
#                                self$testSplits=testSplits
#                                self$trainDesign=testDesign
#                                self$testDesign=testDesign
#                                
#                                
#                              },
#                              searchlight = function(radius, method=c("randomized", "standard")) {
#                                if (method == "randomized") {
#                                  neuroim::RandomSearchlight(self$mask, radius)
#                                } else if (method == "standard") {
#                                  neuroim::Searchlight(self$mask, radius)
#                                } else {
#                                  stop(paste("unrecognized method: ", method))
#                                }
#                              },
#                              
#                              trainChunk = function(vox) {
#                                mat <- series(self$trainVec, vox)
#                                cds <- if (is.vector(vox)) {
#                                  cds <- indexToGrid(space(self$mask), vox)
#                                } else {
#                                  vox
#                                }
#                                neuroim::ROIVolume(space(self$mask), cds, data=mat)
#                                
#                              },
#                              
#                              testChunk = function(vox) {
#                                assertthat::assert_that(!is.null(self$testVec))
#                                cds <- if (is.vector(vox)) {
#                                  cds <- indexToGrid(space(self$mask), vox)
#                                } else {
#                                  vox
#                                }
#                                mat <- series(self$testVec, vox)
#                                neuroim::ROIVolume(space(self$mask), cds, data=mat)
#                              },
#                              
#                              convertScores = function(ids, resultTable) {
#                                ret <- lapply(1:ncol(resultTable), function(i) {
#                                  vol <- array(NA, dim(self$mask))   
#                                  vol[ids] <- resultTable[,i]
#                                  BrainVolume(vol, space(self$mask))
#                                })
#                                
#                                names(ret) <- colnames(resultTable)
#                                ret
#                              },
#                              
#                              averageOverIterations = function(vols) {
#                                nmetrics <- length(vols[[1]])
#                                res <- lapply(1:nmetrics, function(i) {
#                                  X <- do.call(cbind, lapply(vols, function(v) v[[i]]))
#                                  xmean <- rowMeans(X, na.rm=TRUE)
#                                  xmean[is.na(xmean)] <- 0
#                                  BrainVolume(xmean, space(self$mask))
#                                })
#                                
#                                names(res) <- names(vols[[1]])
#                                res
#                                
#                              }
#                              
#                              
#                              
#                              
#                              
#                              
#                              
#                            ))
