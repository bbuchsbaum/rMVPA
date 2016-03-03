
#' @export
#' @import R6
MVPADataset <- R6::R6Class("MVPADataset",
                           public = list(
                             trainVec=NULL,
                             Y=NULL,
                             mask=NULL, 
                             blockVar=NULL,
                             testVec=NULL,
                             testY=NULL,
                             parcellation=NULL, 
                             testSplitVar=NULL,
                             testSplits=NULL, 
                             trainDesign=NULL, 
                             testDesign=NULL,
                      
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
                             
                    
                           }))





#' @export
#' @import R6
BaseModel <- R6::R6Class("BaseModel",
                         public = list(
                           model_name=NA,
                           
                           initialize = function(name) {
                             self$model_name <- name
                           },
                           
                           #fit = function(x, y, ...) { stop("unimplemented") },
                           
                           #predict = function(modelFit,newdata) { stop("unimplemented") },
                           
                           #prob = function(modelFit, newdata) { stop("unimplemnted") },
                           
                           run = function(dataset, vox, crossVal, featureSelector = NULL) { stop("unimplemented") },
                           
                           combineResults = function(resultList) { stop("unimplemented ") },
                           
                           saveResults = function(results, folder) { stop("unimplemented") }
                           
                         ))

#' @export
#' @import R6
CaretModelWrapper <- R6::R6Class(
  "CaretModelWrapper",
  inherit = BaseModel,
  public = list(
    model = NA,
    tuneGrid = NULL,
    initialize = function(model, tuneGrid) {
      self$model_name <- model$label
      self$model <- model
      if (!missing(tuneGrid))
        self$tuneGrid = tuneGrid
      
      for (lib in self$model$library) {
        library(lib, character.only = TRUE)
      }
    },
    
    fit = function(x, y, ...) {
      self$model$fit(x@data,y, ...)
    },
    
    grid = function(x, y, len) {
      self$model$grid(xy,y,len)
    },
    
    predict = function(modelFit,newdata) {
      self$model$predict(modelFit, newdata@data)
    },
    
    prob = function(modelFit, newdata) {
      self$model$prob(modelFit,newdata@data)
    },
    
    run = function(dataset, vox, crossVal, featureSelector = NULL, subIndices=NULL) {
      result <- mvpa_crossval(dataset, vox, crossVal, self$model, self$tuneGrid, featureSelector, subIndices)
      
      observed <- if (is.null(dataset$testVec)) {
        if (is.null(subIndices)) dataset$Y else dataset$Y[subIndices] 
      } else {
        dataset$testY
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
          testInd <- result$testIndices[[i]]
          prob[testInd,] <- prob[testInd,] + p
        }
       
      
        prob <- t(apply(prob, 1, function(vals) vals / sum(vals)))
        maxid <- apply(prob, 1, which.max)
        pclass <- levels(dataset$Y)[maxid]
        classificationResult(observed, pclass, prob, predictor)
      } else {
        
        preds <- numeric(length(observed))
        for (i in seq_along(result$prediction)) {
          testInd <- result$testIndices[[i]]
          preds[testInd] <- result$prediction[[i]]$preds
        }
        
        classificationResult(observed, preds, NULL, predictor)
        
      
    }
  },
    
    performance = function(result) {
      
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


# #' @export
# RSAModel <- R6::R6Class(
#   "RSAModel",
#   inherit = BaseModel,
#   
#   public = list(
#     model_name = "RSA",
#     simFun = NA,
#     contrastMatrix=NULL,
#     
#     initialize = function(simFun = cor, contrastMatrix=NULL) {
#       self$simFun = simFun
#     },
#     
#     run = function(dataset, vox, crossVal, featureSelector = NULL) {
#       patternSimilarity(dataset, vox, self$simFun, self$contrastMatrix)
#     },
#     
#     
#     performance = function(simResult) {
#       list(avgContrast=simResult$avgContrast, effContrast=simResult$avgContrast/simResult$sdContrast)
#     },
#     
#     
#     combineResults = function(resultList) {
#       rois <- sapply(resultList, function(res) attr(res, "ROINUM"))
#       
#       corMatList <- lapply(resultList, function(res) {
#            res$corMat   
#       })
#       
#       lapply(resultList, function(res) {
#         data.frame(ROI=rep(attr(res, "ROINUM"), length(res$observed)), observed=res$observed, pred=res$predicted, correct=as.character(res$observed) == as.character(res$predicted), prob=res$prob)
#         
#       
#       
#       
#       
#          
# 
#       #   names(simMatList) <- rois
#       #   names(simWithinList) <- rois
#       #   
#       #   ret <- list(simMatList=simMatList, simWithinList=simWithinList)
#       #   
#       #   class(ret) <- c("SimilarityResultList", "list")
#       #   ret
#       # }
#       # 
#     }
#   )
# )



