


#' @export
#' @import R6
#' @import neuroim
MVPADataset <- R6::R6Class("MVPADataset",
                           public = list(
                             trainVec=NULL,Y=NULL,mask=NULL, 
                             blockVar=NULL,testVec=NULL,testY=NULL,
                             parcellation=NULL,testSplitVar=NULL,testSplits=NULL, 
                             trainDesign=NULL, testDesign=NULL,
                             initialize=function(trainVec,
                                                 Y,
                                                 mask=NULL, 
                                                 blockVar=NULL,
                                                 testVec=NULL,
                                                 testY=NULL,
                                                 parcellation=NULL, 
                                                 testSplitVar=NULL,
                                                 testSplits=NULL, 
                                                 trainDesign=NULL, 
                                                 testDesign=NULL) {
                               
                               self$trainVec=trainVec
                               self$Y=Y
                               self$mask=mask
                               self$blockVar=blockVar
                               self$testVec=testVec
                               self$testY=testY
                               self$parcellation=parcellation
                               self$testSplitVar=testSplitVar
                               self$testSplits=testSplits
                               self$trainDesign=testDesign
                               self$testDesign=testDesign
                               
                               
                             },
                             searchlight = function(radius, method=c("randomized", "standard")) {
                               if (method == "randomized") {
                                 neuroim::RandomSearchlight(self$mask, radius)
                               } else if (method == "standard") {
                                 neuroim::Searchlight(self$mask, radius)
                               } else {
                                 stop(paste("unrecognized method: ", method))
                               }
                             },
                             
                             trainChunk = function(vox) {
                               mat <- series(self$trainVec, vox)
                               cds <- if (is.vector(vox)) {
                                 cds <- indexToGrid(space(self$mask), vox)
                               } else {
                                 vox
                               }
                               neuroim::ROIVolume(space(self$mask), cds, data=mat)
                               
                             },
                             
                             testChunk = function(vox) {
                               assertthat::assert_that(!is.null(self$testVec))
                               cds <- if (is.vector(vox)) {
                                 cds <- indexToGrid(space(self$mask), vox)
                               } else {
                                 vox
                               }
                               mat <- series(self$testVec, vox)
                               neuroim::ROIVolume(space(self$mask), cds, data=mat)
                             },
                             
                             convertScores = function(ids, resultTable) {
                               ret <- lapply(1:ncol(resultTable), function(i) {
                                 vol <- array(NA, dim(self$mask))   
                                 vol[ids] <- resultTable[,i]
                                 BrainVolume(vol, space(self$mask))
                               })
                               
                               names(ret) <- colnames(resultTable)
                               ret
                             },
                             
                             averageOverIterations = function(vols) {
                               nmetrics <- length(vols[[1]])
                               res <- lapply(1:nmetrics, function(i) {
                                 X <- do.call(cbind, lapply(vols, function(v) v[[i]]))
                                 xmean <- rowMeans(X, na.rm=TRUE)
                                 xmean[is.na(xmean)] <- 0
                                 BrainVolume(xmean, space(self$mask))
                               })
                               
                               names(res) <- names(vols[[1]])
                               res
                               
                             }
                             
                             
                             
                             
                             
                             
                             
                           ))

#' @export
#' @import neuroim
#' @import R6
MVPASurfaceDataset <- R6::R6Class("MVPASurfaceDataset",
                                  public = list(
                                    trainVec=NULL, Y=NULL,
                                    mask=NULL, blockVar=NULL,testVec=NULL,
                                    testY=NULL,parcellation=NULL, 
                                    testSplitVar=NULL,testSplits=NULL, 
                                    trainDesign=NULL, testDesign=NULL,
                                    initialize=function(trainVec,
                                                        Y,
                                                        mask=NULL, 
                                                        blockVar=NULL,
                                                        testVec=NULL,
                                                        testY=NULL,
                                                        parcellation=NULL, 
                                                        testSplitVar=NULL,
                                                        testSplits=NULL, 
                                                        trainDesign=NULL, 
                                                        testDesign=NULL) {
                                      
                                      self$trainVec=trainVec
                                      self$Y=Y
                                      self$mask=mask
                                      self$blockVar=blockVar
                                      self$testVec=testVec
                                      self$testY=testY
                                      self$parcellation=parcellation
                                      self$testSplitVar=testSplitVar
                                      self$testSplits=testSplits
                                      self$trainDesign=testDesign
                                      self$testDesign=testDesign
                                      
                                    },
                                    
                                    searchlight = function(radius, method=c("randomized", "standard")) {
                                      if (method == "randomized") {
                                        neuroim::RandomSurfaceSearchlight(self$trainVec@geometry, radius, self$mask)
                                      } else if (method == "standard") {
                                        neuroim::SurfaceSearchlight(self$trainVec@geometry, radius, self$mask)
                                      } else {
                                        stop(paste("unrecognized method: ", method))
                                      }
                                    },
                                    
                                    trainChunk = function(ids) {
                                      mat <- series(self$trainVec, ids)
                                      neuroim::ROISurface(geometry=self$trainVec@geometry,indices=ids,data=as.matrix(mat))
                                    },
                                    
                                    testChunk = function(ids) {
                                      assertthat::assert_that(!is.null(self$testVec))
                                      mat <- series(self$testVec, vox)
                                      neuroim::ROISurface(geometry=self$testVec@geometry,indices=ids,data=as.matrix(mat))
                                    },
                                    
                                    convertScores = function(ids, resultTable) {
                                      assert_that(length(ids) == nrow(resultTable))
                                      ord <- order(ids)
                                      df1 <- as.data.frame(cbind(ids[ord], resultTable[ord,]))
                                      names(df1) <- c("node", colnames(resultTable))
                                      df1
                                    },
                                    averageOverIterations = function(tbls) {
                                      #allnodes <- sort(unique(unlist(lapply(tbls, function(x) x[,1]))))
                                      
                                      
                                      nodes <- tbls[[1]][,1]
                                      nmetrics <- ncol(tbls[[1]]) - 1
                                      scores <- lapply(tbls, function(x) x[,-1])
                                      
                                      out <- do.call(cbind, lapply(1:nmetrics, function(i) {
                                        X <- do.call(cbind, lapply(scores, function(sc) sc[,i]))
                                        xmean <- rowMeans(X, na.rm=TRUE)
                                        xmean[is.na(xmean)] <- 0
                                        xmean
                                      }))
                                      
                                      out <- as.data.frame(cbind(nodes, out))
                                      colnames(out) <- c("node", colnames(scores[[1]]))
                                      out
                                      
                                    }
                                    
                                    
                                  ))


#' @export
#' @import R6
BaseModel <- R6::R6Class(
  "BaseModel",
  public = list(
    model_name = NA,
    
    initialize = function(name) {
      self$model_name <- name
    },
    
    run = function(dataset, vox, crossVal, featureSelector = NULL) {
      stop("unimplemented")
    },
    
    combineResults = function(resultList) {
      stop("unimplemented ")
    },
    
    saveResults = function(results, folder) {
      stop("unimplemented")
    }
    
  )
)

#' @export
#' @import R6
CaretModelWrapper <- R6::R6Class(
  "CaretModelWrapper",
  inherit = BaseModel,
  public = list(
    model = NA,
    model_name = NULL,
    tuneGrid = NULL,
    customPerformance = NULL,
    
    initialize = function(model, tuneGrid, customPerformance=NULL) { 
      self$model_name <- model$label
      self$model <- model
      
      if (!missing(tuneGrid))
        self$tuneGrid = tuneGrid
      
      if (!is.null(customPerformance)) {
        assert_that(is.function(customPerformance)) 
        self$customPerformance = customPerformance
      }
      
      for (lib in self$model$library) {
        library(lib, character.only = TRUE)
      }
    },
    
    fit = function(x, y, ...) {
      self$model$fit(x,y, ...)
    },
    
    grid = function(x, y, len) {
      self$model$grid(xy,y,len)
    },
    
    predict = function(modelFit,newdata) {
      self$model$predict(modelFit, newdata)
    },
    
    prob = function(modelFit, newdata) {
      self$model$prob(modelFit,newdata)
    },
    
    run = function(dataset, roi, crossVal, featureSelector = NULL, subIndices=NULL) {
      result <- mvpa_crossval(dataset, roi, crossVal, self$model, self$tuneGrid, featureSelector, subIndices)
      
      observed <- if (is.null(dataset$testVec)) {
        if (is.null(subIndices)) dataset$Y else dataset$Y[subIndices] 
      } else {
        dataset$testY
      }
      
      testDesign <- if (is.null(dataset$testVec)) {
        if (is.null(subIndices)) dataset$trainDesign else dataset$trainDesign[subIndices,] 
      } else {
        dataset$testDesign
      }
      
      predictor <-
        if (length(result$predictor) > 1) {
          WeightedPredictor(result$predictor)
        } else {
          result$predictor
        }
      
      if (is.factor(observed)) {
        
        prob <- matrix(0, length(observed), length(levels(dataset$Y)))
        colnames(prob) <- levels(dataset$Y)
        
        for (i in seq_along(result$prediction)) {
          p <- as.matrix(result$prediction[[i]]$probs)
          tind <- result$testIndices[[i]]
          prob[tind,] <- prob[tind,] + p
        }
        
        prob <- t(apply(prob, 1, function(vals) vals / sum(vals)))
        maxid <- apply(prob, 1, which.max)
        pclass <- levels(dataset$Y)[maxid]
        classification_result(observed, pclass, prob, testDesign, predictor)
        
      } else {
        
        preds <- numeric(length(observed))
        for (i in seq_along(result$prediction)) {
          testInd <- result$testIndices[[i]]
          preds[testInd] <- result$prediction[[i]]$preds
        }
        
        regression_result(observed, preds, testDesign, predictor)
        
      }
    },
    
    performance = function(result, vox, splitList=NULL, classMetrics=NULL) {
      
    },
    
    combineResults = function(resultList) {
      rois <- sapply(resultList, function(res) attr(res, "ROINUM"))
      
      predictorList <- lapply(resultList, function(res) {
        MVPAVoxelPredictor(res$predictor, attr(res, "vox"))
      })
      
      predFrame <- if (is.factor(resultList[[1]]$observed)) {
        as.data.frame(do.call(rbind, lapply(resultList, function(res) {
          data.frame(ROI=rep(attr(res, "ROINUM"), length(res$observed)), observed=res$observed, pred=res$predicted, correct=as.character(res$observed) == as.character(res$predicted), prob=res$prob)
        })))
      } else {
        as.data.frame(do.call(rbind, lapply(resultList, function(res) {
          data.frame(ROI=rep(attr(res, "ROINUM"), length(res$observed)), observed=res$observed, pred=res$predicted)
        })))
      }
      
      ret <- list(predictor=ListPredictor(predictorList, rois), predictions=predFrame)
      class(ret) <- c("ClassificationResultList", "list")
      ret
    },
    
    saveResults = function(results, folder) {
      write.table(format(results$predictions,  digits=2, scientific=FALSE, drop0trailing=TRUE), 
                  paste0(paste0(folder, "/prediction_table.txt")), row.names=FALSE, quote=FALSE)  
      if (!is.null(results$predictor)) {
        saveRDS(results$predictor, paste0(folder, "/predictor.RDS"))
      }
    }
  )
)


