

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
                            
                                 pres <- prcomp(as.matrix(x), scale=TRUE)

                                 lda.fit <- lda(pres$x[,1:param$ncomp, drop=FALSE], y)
                                 attr(lda.fit, "ncomp") <- param$ncomp
                                 attr(lda.fit, "pcfit") <- pres
                                 lda.fit
                               },
                                 
                               predict=function(modelFit, newdata, preProc = NULL, submodels = NULL) {
                                 compind <- seq(1, attr(modelFit, "ncomp"))
                                 
                                 pcfit <- attr(modelFit, "pcfit")
                                 colnames(newdata) <- rownames(pcfit$rotation)
                                 pcx <- predict(pcfit, newdata)[,compind,drop=FALSE]
                                 predict(modelFit, pcx)$class
                               },
                               prob=function(modelFit, newdata, preProc = NULL, submodels = NULL) {
                                 compind <- seq(1, attr(modelFit, "ncomp"))
                                 pcfit <- attr(modelFit, "pcfit")
                                 colnames(newdata) <- rownames(pcfit$rotation)
                                 pcx <- predict(pcfit, newdata)[,compind,drop=FALSE]
                                 predict(modelFit, pcx)$posterior                              
                               })


MVPAModels$gpca_lda <- list(type = "Classification", 
                               library = c("sGPCA", "MASS"), 
                               loop = NULL, 
                               parameters=data.frame(parameters=c("ncomp", "theta"), class=c("numeric", "numeric"), labels=c("number of PCs", "smoothing")),
                               grid=function(x, y, len = NULL) {
                                 data.frame(ncomp=5, theta=5)
                               },
                               fit=function(x, y, wts, param, lev, last, weights, classProbs, ...) {
                                 args <- list(...)
                                 vox <- args$vox
                                 xs <- scale(x, scale=FALSE)
                                 xc <- attr(xs, "scaled:center")
                                 
                                 R <- Exp.cov(vox,theta=param$theta)
                                 er <- eigen(R,only.values=TRUE)
                                 R <- R/max(er$values)
                                 fit <- gpca(xs,diag(nrow(x)),R,K=param$ncomp)
                          
                                 lda.fit <-lda(fit$U, y)
                                 attr(lda.fit, "centroid") <- xc
                                 attr(lda.fit, "pcfit") <- fit
                                 attr(lda.fit, "ncomp") <- param$ncomp
                                 lda.fit
                               },
                               predict=function(modelFit, newdata, preProc = NULL, submodels = NULL) {
                                 compind <- seq(1, attr(modelFit, "ncomp"))
                                 pcx <- sweep(newdata, 2, attr(modelFit, "centroid")) %*% attr(modelFit, "pcfit")$V
                                 predict(modelFit, pcx)$class 
                               },
                               prob=function(modelFit, newdata, preProc = NULL, submodels = NULL) {
                                 compind <- seq(1, attr(modelFit, "ncomp"))
                                 pcx <- sweep(newdata, 2, attr(modelFit, "centroid")) %*% attr(modelFit, "pcfit")$V
                                 predict(modelFit, pcx)$posterior
                               })

MVPAModels$liblinear <- list(type = "Classification", 
                               library = "LiblineaR", 
                               loop = NULL, 
                               parameters=data.frame(parameters=c("type", "cost"), class=c("numeric", "numeric"), labels=c("model type", "cost of constraints violation")),
                               grid=function(x, y, len = NULL) {
                                 data.frame(type=0, cost=heuristicC(x))
                               },
                               fit=function(x, y, wts, param, lev, last, weights, classProbs, ...) LiblineaR::LiblineaR(x,y,param$type, param$cost),
                               predict=function(modelFit, newdata, preProc = NULL, submodels = NULL) {
                                 predict(modelFit, as.matrix(newdata))$predictions
                               },
                               prob=function(modelFit, newdata, preProc = NULL, submodels = NULL) {
                                 predict(modelFit, as.matrix(newdata), prob=TRUE)
                               })


MVPAModels$nearestMean <- list(type = "Classification", 
                          library = "klaR", 
                          loop = NULL, 
                          parameters=data.frame(parameters="gamma", class="numeric", labels="gamma"),
                          grid=function(x, y, len = NULL) {
                            data.frame(gamma=.3)
                          },
                          fit=function(x, y, wts, param, lev, last, weights, classProbs, ...) klaR::nm(x,y, param$gamma),
                          predict=function(modelFit, newdata, preProc = NULL, submodels = NULL) {
                            predict(modelFit, as.matrix(newdata))$class
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
                          scores <- -t(predict(modelFit, newdata)$scores)
                          mc <- scores[cbind(1:nrow(scores), max.col(scores, ties.method = "first"))]
                          probs <- exp(scores - mc)
                          zapsmall(probs/rowSums(probs))
                        })


#' create an \code{MVPA} instance
#' @param trainVec
#' @param Y
#' @param mask
#' @param blockVar
#' @param testVec
#' @param testY
#' @param modelName
#' @param tuneGrid
#' @export
MVPADataset <- function(trainVec, Y, mask, blockVar, testVec, testY, modelName="corsim", tuneGrid=NULL) {
  
  model <- loadModel(modelName)
  
  testSets <- split(1:length(blockVar), blockVar)
  trainSets <- invertFolds(testSets, length(Y))
  
  #index <- invertFolds(testSets[-blockIndex], nrow(Xtrain))   
  
  ret <- list(
    trainVec=trainVec,
    Y=Y,
    mask=mask,
    blockVar=blockVar,
    testVec=testVec,
    testY=testY,
    model=model,
    tuneGrid=tuneGrid,
    testSets=testSets,
    trainSets=trainSets)
   
    
  
  class(ret) <- c("MVPADataset", "list")
  ret
}

NullResult <- function(vox, model, observed) {
  ret <- list(vox=vox,
              model=model,
              observed=observed)
  class(ret) <- c("NullResult", "list")
  ret
  
}

#' create an \code{TwoWayClassification} instance
#' @param vox
#' @param model
#' @param observed
#' @param predicted
#' @param probs
#' @param finalFit
#' @export
TwoWayClassificationResult <- function(vox, model, observed, predicted, probs, finalFit=NULL) {
  ret <- list(vox=vox,
              model=model,
              observed=observed,
              predicted=predicted,
              probs=as.matrix(probs),
              finalFit=finalFit)
  
  class(ret) <- c("TwoWayClassificationResult", "list")
  ret
  
}

#' create an \code{TwoWayClassification} instance
#' @param vox
#' @param model
#' @param observed
#' @param predicted
#' @param probs
#' @param finalFit
#' @export
MultiWayClassificationResult <- function(vox, model, observed, predicted, probs, finalFit=NULL) {
  ret <- list(vox=vox,
              model=model,
              observed=observed,
              predicted=predicted,
              probs=as.matrix(probs),
              finalFit=finalFit)
  
  class(ret) <- c("MultiWayClassificationResult", "list")
  ret
  
}


#' @export
RawModel <- function(model, Xtrain, Ytrain, Xtest, Ytest, tuneGrid) {
  fit <- model$fit(Xtrain, Ytrain, NULL, tuneGrid, classProbs=TRUE)  
  ret <- list(
              model=model,
              Xtrain=Xtrain,
              Ytrain=Ytrain,
              Xtest=Xtest,
              Ytest=Ytest,
              tuneGrid=tuneGrid,
              modelFit=fit)
  
  class(ret) <- c("RawModel", "list")
  ret
}

#' @export
RawPredictor <- function(fit, model) {
  ret <- list(fit=fit,
              model=model)
  
  class(ret) <- c("RawPredictor", "list")
  ret
}

#' @export
CaretPredictor <- function(fit) {
  ret <- list(fit=fit)
  
  class(ret) <- c("CaretPredictor", "list")
  ret
}

#' @export
MVPAPredictor <- function(predictor, voxelGrid) {
  ret <- list(predictor=predictor, voxelGrid=voxelGrid)
  
  class(ret) <- c("MVPAPredictor", "list")
  ret
}


#' @export
CaretModel <- function(model, Xtrain, Ytrain, Xtest, Ytest, tuneGrid, tuneControl,...) {
  fit <- caret::train(Xtrain, Ytrain, method=model, trControl=tuneControl, tuneGrid=tuneGrid, ...)
  ret <- list(
              model=model,
              Xtrain=Xtrain,
              Ytrain=Ytrain,
              Xtest=Xtest,
              Ytest=Ytest,
              tuneGrid=tuneGrid,
              tuneControl=tuneControl,
              modelFit=fit)
  
  class(ret) <- c("CaretModel", "list")
  ret
}

#' @export
evaluateModel <- function(x, newdata) {
  UseMethod("evaluateModel")
}

#' @export
asPredictor <- function(x, voxelGrid) {
  UseMethod("asPredictor")
}

#' @export
asPredictor.RawModel <- function(x, voxelGrid) {
  MVPAPredictor(RawPredictor(x$modelFit, x$model), voxelGrid)
}

#' @export
asPredictor.CaretModel <- function(x, voxelGrid, brainSpace) {
  MVPAPredictor(CaretPredictor(x$modelFit), voxelGrid)
}


#' @export
evaluateModel.RawModel <- function(x, newdata=NULL) {
  if (is.null(newdata)) {
    newdata=x$Xtest
  }
  
  probs <- x$model$prob(x$modelFit, newdata=x$Xtest)     
  cpred <- apply(probs,1, which.max)
  cpred <- levels(x$Ytrain)[cpred]
  list(class=cpred, probs=probs)
}

#' @export
evaluateModel.MVPAPredictor <- function(x, newdata) {
  if (is.null(newdata)) {
    stop("newdata cannot be null")
  }
  
  X <- series(newdata,x$voxelGrid)
  evaluateModel(x$predictor, X)
}


#' @export
evaluateModel.RawPredictor <- function(x, newdata=NULL) {
  if (is.null(newdata)) {
    stop("newdata cannot be null")
  }
  
  probs <- x$model$prob(x$fit, newdata=newdata)     
  cpred <- apply(probs,1, which.max)
  cpred <- levels(x$Ytrain)[cpred]
  list(class=cpred, probs=probs)
}


#' @export
evaluateModel.CaretModel <- function(x, newdata=NULL) {
  if (is.null(newdata)) {
    newdata=x$Xtest
  }
  
  list(class=predict(x$modelFit, newdata=x$Xtest), probs=predict(x$modelFit, newdata=x$Xtest, type="prob"))
}

#' @export
evaluateModel.CaretPredictor <- function(x, newdata=NULL) {
  if (is.null(newdata)) {
    stop("newdata cannot be null")
  }
  
  list(class=predict(x$fit, newdata=newdata), probs=predict(x$fit, newdata=newdata, type="prob"))
}

#' @export
evaluateModelList <- function(modelList, newdata=NULL) {
  results <- lapply(modelList, evaluateModel, newdata)
  probMat <- do.call(rbind, lapply(results, "[[", "probs"))
  predClass <- unlist(lapply(results, "[[", "class"))
  list(class=predClass, probs=probMat)
}


#' @export
performance <- function(x,...) {
  UseMethod("performance")
}



#' @export
performance.TwoWayClassificationResult <- function(x,...) {
  
  ncorrect <- sum(x$observed == x$predicted)
  ntotal <- length(x$observed)
  maxClass <- max(table(x$observed))
  
  out <- binom.test(ncorrect,
                    ntotal,
                    p = maxClass/ntotal,
                    alternative = "greater")
  
  c(ZAccuracy=-qnorm(out$p.value), Accuracy=ncorrect/ntotal, AUC=Metrics::auc(x$observed == levels(x$observed)[2], x$probs[,2])-.5)
}

#' @export
performance.MultiWayClassificationResult <- function(x,...) {
  obs <- as.character(x$observed)
  
  ncorrect <- sum(obs == x$predicted)
  ntotal <- length(obs)
  maxClass <- max(table(obs))
  
  out <- binom.test(ncorrect,
                    ntotal,
                    p = maxClass/ntotal,
                    alternative = "greater")
  
  aucres <- sapply(1:ncol(x$prob), function(i) {
    lev <- colnames(x$prob)[i]
    pos <- obs == lev
    pclass <- x$prob[,i]
    pother <- rowMeans(x$prob[,-i])
    Metrics::auc(as.numeric(pos), pclass - pother)-.5
  })
  
  names(aucres) <- paste0("AUC_", colnames(x$prob))

  c(ZAccuracy=-qnorm(out$p.value), Accuracy=sum(obs == as.character(x$predicted))/length(obs), Combined_AUC=mean(aucres), aucres)
}



.noneControl <- caret::trainControl("none", verboseIter=TRUE, classProbs=TRUE)
.cvControl <- caret::trainControl("cv", verboseIter=TRUE, classProbs=TRUE)  


#' @export
invertFolds <- function(foldSplit, len) {
  index <- relist(1:length(unlist(foldSplit)), foldSplit)
  allInd <- 1:len
  index <- lapply(index, function(i) allInd[-i])
}

#' @export
fitFinalModel <- function(Xtrain, Ytrain,  method, Xtest=NULL, Ytest=NULL, tuneGrid=NULL, tuneControl=NULL, ...) {  
 
  if (is.null(tuneGrid)) {
    tuneGrid <- method$grid(Xtrain, Ytrain, 1)
  }
  
  if (nrow(tuneGrid) == 1) {
    RawModel(method, Xtrain, Ytrain, Xtest, Ytest, tuneGrid)   
  } else {
    if (is.null(tuneControl)) {
      tuneControl <- .cvControl
    }
    CaretModel(method, Xtrain, Ytrain, Xtest, Ytest, tuneGrid, tuneControl, ...)
  }
}
  
#' @export
crossValidate <- function(X, Y, trainSets, testSets, model, tuneGrid, fast=TRUE, finalFit=FALSE, ncores=1, ...) {
 
  if (!is.factor(Y)) {
    stop("regression not supported yet") 
  }
   
  .lapply <- if (ncores > 1) {
      function(sets, FUN, ...) {
        parallel::mclapply(sets, FUN, ..., mc.cores=ncores)
      }
    } else {
      lapply      
    }
  
  fittedModels <- 
    .lapply(seq_along(testSets), function(blockIndex) {
      testIndices <- testSets[[blockIndex]]
      Xtrain <- X[-testIndices,]
      
      ret <- if (nrow(tuneGrid) == 1 && fast) {
        RawModel(model, Xtrain, Y[-testIndices], X[testIndices,], Y[testIndices], tuneGrid)        
      } else { 
          
        if (nrow(tuneGrid) == 1) {
          CaretModel(model, Xtrain, Y[-testIndices], X[testIndices,], Y[testIndices], tuneGrid, .noneControl)
        } else {
  
          index <- invertFolds(testSets[-blockIndex], nrow(Xtrain)) 
          ctrl <- caret::trainControl("cv", verboseIter=TRUE, classProbs=TRUE, index=index)
          CaretModel(model, Xtrain, Y[-testIndices], X[testIndices,], Y[testIndices], tuneGrid, ctrl)
        }
        
      }
    })
  
  
  final <- if (finalFit) {
    ## fit final model to whole data set
    ctrl <- if (nrow(tuneGrid) == 1) {
      .noneControl
    } else {
      caret::trainControl("cv", verboseIter=TRUE, classProbs=TRUE, index=trainSets)
    }    
    
    CaretModel(model, X, Y, NULL, NULL, tuneGrid, ctrl)
  } 
  
  result <- evaluateModelList(fittedModels)
  list(class=result$class, prob=result$prob, observed=Y, modelFits=fittedModels, finalFit=final)
   
}

#' @import neuroim
#' @export
fitMVPAModel <- function(dataset, voxelGrid, tuneLength=1, fast=TRUE, finalFit=FALSE, ncores=1) {
  
  M <- series(dataset$trainVec, voxelGrid) 
  
  if (ncol(M) < 2) {
    stop("feature matrix must have at least two columns: returning NullResult")
  }
  

  hasVariance <- which(apply(M, 2, sd, na.rm=TRUE) > 0)
  M <- M[, hasVariance, drop=FALSE]
  
  hasNA <- apply(M,1, function(x) any(is.na(x)))
  numNAs <- sum(hasNA)
  
  if (numNAs > 0) {
    stop("training data has NA values, aborting")
  } 
  
  tuneGrid <- if (is.null(dataset$tuneGrid)) {
    tuneGrid <- dataset$model$grid(M, dataset$Y, tuneLength)
  } else {
    dataset$tuneGrid
  }
  
  if (ncol(M) < 2) {
    stop("feature matrix must have at least two columns with nonzero variance")
  }
  
  voxelGrid <- voxelGrid[hasVariance, ]
  
  result <- if (is.null(dataset$testVec)) {
    crossValidate(M, dataset$Y, dataset$trainSets, dataset$testSets, dataset$model, tuneGrid, tuneLength, testVec=dataset$testVec, testY=dataset$testY, fast=fast,finalFit=finalFit, ncores=ncores)  
  } else {
    Xtest <- series(dataset$testVec, voxelGrid)
    modelFit <- if (nrow(tuneGrid) == 1) {
      fitFinalModel(M, dataset$Y,  dataset$model, Xtest, dataset$testY, tuneGrid)
    } else {      
      ctrl <- caret::trainControl("cv", verboseIter=TRUE, classProbs=TRUE, index=dataset$trainSets) 
      fitFinalModel(M, dataset$Y,  dataset$model, dataset$Xtest, dataset$testY, tuneGrid, tuneControl=ctrl)     
    }
    
    preds <- evaluateModel(modelFit)
    list(class=preds$class, prob=preds$prob, observed=dataset$testY, modelFits=NULL, finalFit=modelFit)
  }
  
  if (is.factor(dataset$Y) && length(levels(dataset$Y)) == 2) {
    TwoWayClassificationResult(voxelGrid, dataset$model, result$observed, result$class, result$prob, result$finalFit)
  } else if (is.factor(dataset$Y) && length(levels(dataset$Y)) >= 3) {
    MultiWayClassificationResult(voxelGrid, dataset$model, result$observed, result$class, result$prob, result$finalFit)
  } else {
    stop("regression not supported yet.")
  }
}


