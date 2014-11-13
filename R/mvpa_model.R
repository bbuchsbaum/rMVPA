

#' @export
#' @import MASS
corsimFit <- function(x, y, method, robust) {
  estimator <- if (robust) {
    function(vec)  {
      h <- try(huber(vec))
      if (inherits(h, "try-error")) {
        median(vec)
      } else {
        h$mu
      }
    }      
  } else {
    mean
  }
  
  
  list(conditionMeans = splitReduce(as.matrix(x), y, estimator), levs=levels(y), method=method, robust=robust)
}

#' @export
predict.corsimFit <- function(modelFit, newData) {
  res <- sapply(1:nrow(newData), function(i) {
    pattern <- newData[i,]
    which.max(cor(pattern, t(modelFit$conditionMeans), method=modelFit$method))
  })
  
  modelFit$levs[res]
}

#' @export
prob.corsimFit <- function(modelFit, newData) {
  scores <- cor(t(newData), t(modelFit$conditionMeans), method=modelFit$method)
  
  mc <- scores[cbind(1:nrow(scores), max.col(scores, ties.method = "first"))]
  probs <- exp(scores - mc)
  probs <- zapsmall(probs/rowSums(probs))
  colnames(probs) <- modelFit$levs
  probs
}

#' @export
MVPAModels <- list()

MVPAModels$pca_lda <- list(type = "Classification", 
                               library = "MASS", 
                               loop = NULL, 
                               parameters=data.frame(parameters="ncomp", class="numeric", labels="ncomp"),
                               grid=function(x, y, len = 5) {
                                 data.frame(ncomp=1:len)
                               },
                               fit=function(x, y, wts, param, lev, last, weights, classProbs, ...) {
                                 pres <- prcomp(x, scale=TRUE)
                                 lda(pres$x[,1:param$ncomp])
                               },
                                 
                               predict=function(modelFit, newdata, preProc = NULL, submodels = NULL) {
                                 modelFit$lev[predict(modelFit, as.matrix(newdata))$class]
                               },
                               prob=function(modelFit, newdata, preProc = NULL, submodels = NULL) {
                                 predict(modelFit, as.matrix(newdata))$posterior
                               })

MVPAModels$nearestMean <- list(type = "Classification", 
                          library = "klaR", 
                          loop = NULL, 
                          parameters=data.frame(parameters="gamma", class="numeric", labels="gamma"),
                          grid=function(x, y, len = NULL) {
                            data.frame(gamma=.1)
                          },
                          fit=function(x, y, wts, param, lev, last, weights, classProbs, ...) klaR::nm(x,y, param$gamma),
                          predict=function(modelFit, newdata, preProc = NULL, submodels = NULL) {
                            modelFit$lev[predict(modelFit, as.matrix(newdata))$class]
                          },
                          prob=function(modelFit, newdata, preProc = NULL, submodels = NULL) {
                            predict(modelFit, as.matrix(newdata))$posterior
                          })
    
MVPAModels$corsim <- list(type = "Classification", 
                    library = "rMVPA", 
                    loop = NULL, 
                    parameters=data.frame(parameters=c("method", "robust"), class=c("character", "logical"), label=c("correlation type: pearson, spearman, or kendall", "mean or huber")),
                    grid=function(x, y, len = NULL) if (len == 1) { data.frame(method="pearson", robust=FALSE) } else { expand.grid(method=c("pearson", "spearman", "kendall"), robust=c(TRUE, FALSE)) },
                    fit=function(x, y, wts, param, lev, last, weights, classProbs, ...) corsimFit(x,y, as.character(param$method), param$robust),
                    predict=function(modelFit, newdata, preProc = NULL, submodels = NULL) predict.corsimFit(modelFit, as.matrix(newdata)),
                    prob=function(modelFit, newdata, preProc = NULL, submodels = NULL) {
                      prob.corsimFit(modelFit, as.matrix(newdata))
                    })


MVPAModels$sda_notune <- list(type = "Classification", 
                 library = "sda", 
                 loop = NULL, 
                 parameters=data.frame(parameters="parameter", class="character", label="parameter"),
                 grid=function(x, y, len = NULL) data.frame(parameter="none"),
                 fit=function(x, y, wts, param, lev, last, weights, classProbs, ...) sda::sda(Xtrain=as.matrix(x), L=y, verbose=FALSE, ...),
                 predict=function(modelFit, newdata, preProc = NULL, submodels = NULL) predict(modelFit, as.matrix(newdata), verbose=FALSE)$class,
                 prob=function(modelFit, newdata, preProc = NULL, submodels = NULL) {
                    predict(modelFit, as.matrix(newdata),verbose=FALSE)$posterior
                 })

MVPAModels$sda_ranking <- list(type = "Classification", 
                              library = "sda", 
                              loop = NULL, 
                              parameters=data.frame(parameters="parameter", class="character", label="parameter"),
                              grid=function(x, y, len = NULL) data.frame(parameter="none"),
                              fit=function(x, y, wts, param, lev, last, weights, classProbs, ...) {
              
                                x <- as.matrix(x)                             
                                ind <- if (ncol(x) > 20) {
                                  rank <- sda::sda.ranking(Xtrain=x, L=y, fdr=TRUE, verbose=FALSE, ...)
                                  hcind <- which.max(rank[,"HC"])
                                  
                                  keep.ind <- if (length(hcind) < 2) {
                                    1:2
                                  } else {
                                    hcind                               
                                  }                                                                
                                  rank[keep.ind,"idx"]
                                } else if (ncol(x) <= 3) {
                                  1:ncol(x)
                      
                                } else {
                                  rank <- sda::sda.ranking(Xtrain=x, L=y, fdr=FALSE, verbose=FALSE, ...)
                                  rank[1:(ncol(x)/2), "idx"]
                                }
                                
                                if (length(ind) < 2) {
                                  browser()
                                }
                                fit <- sda::sda(Xtrain=x[,ind,drop=FALSE], L=y, verbose=FALSE)
                                attr(fit, "keep.ind") <- ind
                                fit
                              },
                              predict=function(modelFit, newdata, preProc = NULL, submodels = NULL) {                        
                                predict(modelFit, as.matrix(newdata[,attr(modelFit, "keep.ind"), drop=FALSE]), verbose=FALSE)$class
                              },
                              prob=function(modelFit, newdata, preProc = NULL, submodels = NULL) {
                                predict(modelFit, as.matrix(newdata[,attr(modelFit, "keep.ind"), drop=FALSE]),verbose=FALSE)$posterior
                              })


              
MVPAModels$lda_thomaz <- list(type = "Classification", 
                        library = "sparsediscrim", 
                        loop = NULL, 
                        parameters=data.frame(parameters="parameter", class="character", label="parameter"),
                        grid=function(x, y, len = NULL) data.frame(parameter="none"),
                        fit=function(x, y, wts, param, lev, last, weights, classProbs, ...) { lda_thomaz(x,y, ...) },
                        predict=function(modelFit, newdata, preProc = NULL, submodels = NULL) { predict(modelFit, as.matrix(newdata))$class },
                        prob=function(modelFit, newdata, preProc = NULL, submodels = NULL) { 
                          scores <- t(predict(modelFit, newdata)$scores)
                          mc <- scores[cbind(1:nrow(scores), max.col(scores, ties.method = "first"))]
                          probs <- exp(scores - mc)
                          zapsmall(probs/rowSums(probs))
                        })


#' create an \code{MVPA} instance
#' @param trainVec
#' @param Y
#' @param mask
#' @param blockVar
#' @export
MVPADataset <- function(trainVec, Y, mask, blockVar) {
  ret <- list(
    trainVec=trainVec,
    Y=Y,
    mask=mask,
    blockVar=blockVar)
  
  class(ret) <- c("MVPADataset", "list")
  ret
}

#' create an \code{TwoWayClassification} instance
#' @param vox
#' @param model
#' @param observed
#' @param predicted
#' @param probs
#' @export
TwoWayClassificationResult <- function(vox, model, observed, predicted, probs, modelFits=NULL) {
  ret <- list(vox=vox,
              model=model,
              observed=observed,
              predicted=predicted,
              probs=as.matrix(probs),
              modelFits=modelFits)
  
  class(ret) <- c("TwoWayClassificationResult", "list")
  ret
  
}

#' create an \code{TwoWayClassification} instance
#' @param vox
#' @param model
#' @param observed
#' @param predicted
#' @param probs
#' @export
MultiWayClassificationResult <- function(vox, model, observed, predicted, probs, modelFits=NULL) {
  ret <- list(vox=vox,
              model=model,
              observed=observed,
              predicted=predicted,
              probs=as.matrix(probs),
              modelFits=modelFits)
  
  class(ret) <- c("MultiWayClassificationResult", "list")
  ret
  
}


#' @export
performance <- function(x,...) {
  UseMethod("performance")
}

#' @export
performance.TwoWayClassificationResult <- function(x) {
  c(Accuracy=sum(x$observed == x$predicted)/length(x$observed), AUC=Metrics::auc(x$observed == levels(x$observed)[2], x$probs[,2]))
}

#' @export
performance.MultiWayClassificationResult <- function(x) {
  obs <- as.character(x$observed)
  
  aucres <- sapply(1:ncol(x$prob), function(i) {
    lev <- colnames(x$prob)[i]
    pos <- obs == lev
    pclass <- x$prob[,i]
    pother <- rowMeans(x$prob[,-i])
    Metrics::auc(as.numeric(pos), pclass - pother)
  })
  
  names(aucres) <- paste0("AUC_", colnames(x$prob))

  c(Accuracy=sum(obs == as.character(x$predicted))/length(obs), Combined_AUC=mean(aucres), aucres)
}


#' create a classification model
createModel <- function(method) {
  ret <- list(method=method)
  class(ret) <- c(method, "list")
  ret
}


trainModel <- function(x, ...) {
  UseMethod("trainModel")
}

#' @export
crossval <- function(X, Y, foldSplit, method, ncores=2, tuneGrid=NULL, tuneLength=1, fast=TRUE, finalFit=FALSE) {
  
  if (is.null(tuneGrid)) {
    tuneGrid <- method$grid(X, Y, tuneLength)
  }
  
 
  res <- #parallel::mclapply(foldSplit, function(fidx) {
    
  
  lapply(seq_along(foldSplit), function(find) {
    fidx <- foldSplit[[find]]
    Xtrain <- X[-fidx,]
    Ytrain <- Y[-fidx]
    Xtest <- X[fidx,]   
    
    
    ret <- if (nrow(tuneGrid) == 1 && fast) {
      ## fast fit
      fit <- method$fit(Xtrain, Ytrain, NULL, tuneGrid, classProbs=TRUE)
      #fit <- caret::train(Xtrain, Ytrain, method=method, trControl=ctrl, tuneGrid=tuneGrid)
      probs <- method$prob(fit, newdata=Xtest)
      #probs <- predict(fit, Xtest, type="prob")
      cpred <- apply(probs,1, which.max)
      cpred <- levels(Ytrain)[cpred]
      list(class=cpred, probs=probs, fit=fit)
    } else {
      
      foldSplitM <- foldSplit[-find]
      index <- relist(1:length(unlist(foldSplitM)), foldSplitM)
      allInd <- 1:nrow(Xtrain)
      index <- lapply(index, function(i) allInd[-i])
       
      ctrl <- if (nrow(tuneGrid) == 1) {
        caret::trainControl("none", verboseIter=TRUE, classProbs=TRUE,allowParallel=FALSE)
      } else {
        caret::trainControl("cv", verboseIter=TRUE, classProbs=TRUE,allowParallel=FALSE, index=index)
      }
       
      fit <- caret::train(as.data.frame(Xtrain), Ytrain, method=method, trControl=ctrl, tuneGrid=tuneGrid)
      list(class=predict(fit, newdata=Xtest), probs=predict(fit, newdata=Xtest, type="prob"), fit=fit)
    }
       
  })#, mc.cores=ncores)
  
  final <- if (finalFit) {
    ## fit final model to whole data set
    ctrl <- if (nrow(tuneGrid) == 1) {
      caret::trainControl("none", verboseIter=TRUE, classProbs=TRUE,allowParallel=FALSE)
    } else {
      index <- relist(1:length(unlist(foldSplit)), foldSplit)
      allInd <- 1:nrow(X)
      index <- lapply(index, function(i) allInd[-i])
      caret::trainControl("cv", verboseIter=TRUE, classProbs=TRUE,allowParallel=FALSE, index=index)
    }
    
    caret::train(X, Y, method=method, trControl=ctrl, tuneGrid=tuneGrid)
  } else {
    NULL
  }
  
  probMat <- do.call(rbind, lapply(res, "[[", "probs"))
  predClass <- unlist(lapply(res, "[[", "class"))
  modelFits <- lapply(res, "[[", "fit")
 
  if (is.factor(Y)) {
    list(class=predClass, prob=probMat, modelFits=modelFits, finalFit=final)
  } else {
    stop("regression not supported yet")
  }
}

#' @import neuroim
#' @export
fitMVPAModel <- function(model, bvec, Y, blockVar, voxelGrid, ncores=2, tuneGrid=NULL, tuneLength=1, fast=TRUE, finalFit=FALSE) {
  
  M <- series(bvec, voxelGrid) 
  
  if (ncol(M) < 2) {
    stop("feature matrix must have at least two columns")
  }
  
  hasVariance <- which(apply(M, 2, sd) > 0)
  M <- M[, hasVariance]
  
  if (ncol(M) < 2) {
    stop("feature matrix must have at least two columns with nonzero variance")
  }
  
  voxelGrid <- voxelGrid[hasVariance, ]
  
  result <- crossval(as.matrix(M), Y, split(1:length(blockVar), blockVar), model, ncores, tuneGrid,tuneLength,fast=fast,finalFit=finalFit)  
  
  if (is.factor(Y) && length(levels(Y)) == 2) {
    TwoWayClassificationResult(voxelGrid, model, Y, result$class, result$prob, result$modelFits)
  } else if (is.factor(Y) && length(levels(Y)) >= 3) {
    MultiWayClassificationResult(voxelGrid, model, Y, result$class, result$prob, result$modelFits)
  } else {
    stop("regression not supported yet.")
  }
}

trainModel.default <- function(model, dset, voxelGrid, ncores=2) {
  M <- series(dset$trainVec, voxelGrid)
  fitlist <- crossval(M, dset$Y, dset$blockVar, model)
}



trainModel.lda_strimmer <- function(dset, voxelGrid, ncores=2) {
  M <- series(dset$trainVec, voxelGrid)
  ctrl <- caret::trainControl("none")
  fit <- caret::train(x$x, x$y, method=MVPA.sda, trControl=ctrl, tuneLength=2)
}

trainModel.lda_thomaz <- function(x) {
  ctrl <- caret::trainControl("none")
  fit <- caret::train(x$x, x$y, method=MVPA.lda_thomaz, trControl=ctrl, tuneLength=0)
}
