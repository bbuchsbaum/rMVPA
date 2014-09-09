


corsimFit <- function(x, y) {
  list(conditionMeans = splitReduce(x, y, mean), levs=levels(y))
}

predict.corsimFit <- function(modelFit, newData) {
  res <- sapply(1:nrow(newData), function(i) {
    pattern <- newData[i,]
    which.max(cor(pattern, t(modelFit$conditionMeans)))
  })
  
  modelFit$levs[res]
}

prob.corsimFit <- function(modelFit, newData) {
  scores <- cor(t(newData), t(modelFit$conditionMeans))
  
  mc <- scores[cbind(1:nrow(scores), max.col(scores, ties.method = "first"))]
  probs <- exp(scores - mc)
  probs <- zapsmall(probs/rowSums(probs))
  colnames(probs) <- modelFit$levs
  probs
}

MVPAModels <- list()

MVPAModels$pca_lda <- list(type = "Classification", 
                               library = "MASS", 
                               loop = NULL, 
                               parameters=data.frame(parameters="ncomp", class="numeric", labels="ncomp"),
                               grid=function(x, y, len = 5) {
                                 data.frame(ncomp=1:len)
                               },
                               fit=function(x, y, wts, param, lev, last, weights, classProbs, ...) {
                                 browser()
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
                          fit=function(x, y, wts, param, lev, last, weights, classProbs, ...) nm(x,y, param$gamma),
                          predict=function(modelFit, newdata, preProc = NULL, submodels = NULL) {
                            modelFit$lev[predict(modelFit, as.matrix(newdata))$class]
                          },
                          prob=function(modelFit, newdata, preProc = NULL, submodels = NULL) {
                            predict(modelFit, as.matrix(newdata))$posterior
                          })
    
MVPAModels$corsim <- list(type = "Classification", 
                    library = "rMVPA", 
                    loop = NULL, 
                    parameters=data.frame(parameters="none", class="numeric", labels="none"),
                    grid=function(x, y, len = NULL) data.frame(none=1),
                    fit=function(x, y, wts, param, lev, last, weights, classProbs, ...) corsimFit(x,y),
                    predict=function(modelFit, newdata, preProc = NULL, submodels = NULL) predict.corsimFit(modelFit, as.matrix(newdata)),
                    prob=function(modelFit, newdata, preProc = NULL, submodels = NULL) {
                      prob.corsimFit(modelFit, as.matrix(newdata))
                    })


MVPAModels$sda_notune <- list(type = "Classification", 
                 library = "sda", 
                 loop = NULL, 
                 parameters=data.frame(parameters="lambda", class="numeric", labels="lambda"),
                 grid=function(x, y, len = NULL) data.frame(lambda=0),
                 fit=function(x, y, wts, param, lev, last, weights, classProbs, ...) sda(Xtrain=as.matrix(x), L=y, verbose=FALSE, ...),
                 predict=function(modelFit, newdata, preProc = NULL, submodels = NULL) predict(modelFit, as.matrix(newdata), verbose=FALSE),
                 prob=function(modelFit, newdata, preProc = NULL, submodels = NULL) {
                    predict(modelFit, as.matrix(newdata),verbose=FALSE)$posterior
                 })

MVPAModels$sda_ranking <- list(type = "Classification", 
                              library = "sda", 
                              loop = NULL, 
                              parameters=data.frame(parameters="lambda", class="numeric", labels="lambda"),
                              grid=function(x, y, len = NULL) data.frame(lambda=0),
                              fit=function(x, y, wts, param, lev, last, weights, classProbs, ...) {
                                x <- as.matrix(x)
                                
                                ind <- if (ncol(x) > 20) {
                                  rank <- sda.ranking(Xtrain=x, L=y, fdr=TRUE, verbose=FALSE, ...)
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
                                  rank <- sda.ranking(Xtrain=x, L=y, fdr=FALSE, verbose=FALSE, ...)
                                  rank[1:(ncol(x)/2), "idx"]
                                }
                                
                                if (length(ind) < 2) {
                                  browser()
                                }
                                fit <- sda(Xtrain=x[,ind,drop=FALSE], L=y, verbose=FALSE)
                                attr(fit, "keep.ind") <- ind
                                fit
                              },
                              predict=function(modelFit, newdata, preProc = NULL, submodels = NULL) {                        
                                predict(modelFit, as.matrix(newdata[,attr(modelFit, "keep.ind"), drop=FALSE]), verbose=FALSE)
                              },
                              prob=function(modelFit, newdata, preProc = NULL, submodels = NULL) {
                                predict(modelFit, as.matrix(newdata[,attr(modelFit, "keep.ind"), drop=FALSE]),verbose=FALSE)$posterior
                              })


                


MVPAModels$lda_thomaz <- list(type = "Classification", 
                        library = "sparsediscrim", 
                        loop = NULL, 
                        parameters=data.frame(parameters="prior", class="numeric", labels="prior"),
                        grid=function(x, y, len = NULL) { data.frame(prior=.5) },
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
MVPADataset <- function(trainVec, Y, mask, blockVar) {
  ret <- list(
    trainVec=trainVec,
    Y=Y,
    mask=mask,
    blockVar=blockVar)
  
  class(ret) <- c("MVPADataset", "list")
  ret
}


TwoWayClassificationResult <- function(observed, predicted, probs) {
  ret <- list(observed=observed,
              predicted=predicted,
              probs=probs)
  
  class(ret) <- c("TwoWayClassificationResult", "list")
  ret
  
}

MultiWayClassificationResult <- function(observed, predicted, probs) {
  ret <- list(observed=observed,
              predicted=predicted,
              probs=as.matrix(probs))
  
  
  class(ret) <- c("MultiWayClassificationResult", "list")
  ret
  
}

performance <- function(x,...) {
  UseMethod("performance")
}



performance.TwoWayClassificationResult <- function(x) {
  c(acc=sum(x$observed == x$predicted)/length(x$observed), acc=Metrics::auc(x$observed == levels(x$observed)[2], x$probs[,2]))
}

performance.MultiWayClassificationResult <- function(x) {
  obs <- as.character(x$observed)
  
  aucres <- sapply(1:ncol(x$prob), function(i) {
    lev <- colnames(x$prob)[i]
    pos <- obs == lev
    pclass <- x$prob[,i]
    pother <- rowMeans(x$prob[,-i])
    Metrics::auc(as.numeric(pos), pclass - pother)
  })
  
  names(aucres) <- paste0("auc_", colnames(x$prob))

  c(accuracy=sum(obs == as.character(x$predicted))/length(obs), mauc=mean(aucres), aucres)
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


crossval <- function(X, Y, foldSplit, method, ncores=2, tuneGrid=NULL, tuneLength=1) {
  if (is.null(tuneGrid)) {
    tuneGrid <- method$grid(tuneLength)
  }
  
  if (nrow(tuneGrid) == 1) {
    ctrl <- caret::trainControl("none", verboseIter=TRUE, classProbs=TRUE)
  } else {
    ctrl <- caret::trainControl("cv", verboseIter=TRUE, classProbs=TRUE)
  }
  
  
  
  res <- #parallel::mclapply(foldSplit, function(fidx) {
  lapply(foldSplit, function(fidx) {
    Xtrain <- X[-fidx,]
    Ytrain <- Y[-fidx]
    Xtest <- X[fidx,]   
    
    if (nrow(tuneGrid) == 1) {
      ## fast fit
      fit <- method$fit(Xtrain, Ytrain, NULL, tuneGrid, classProbs=TRUE)
      probs <- method$prob(fit, newdata=Xtest)
      cpred <- apply(probs,1, which.max)
      cpred <- levels(Ytrain)[cpred]
      data.frame(class=cpred, probs)
    } else {
      fit <- caret::train(Xtrain, Ytrain, method=method, trControl=ctrl, tuneGrid=tuneGrid)
      data.frame(class=predict(fit, newdata=Xtest), predict(fit, newdata=Xtest, type="prob"))
    }
  })#, mc.cores=ncores)
  

  
  ret <- do.call(rbind, res)
  if (is.factor(Y) && length(levels(Y)) == 2) {
    TwoWayClassificationResult(Y, ret[,1], ret[,2:3])
  } else if (is.factor(Y) && length(levels(Y)) > 2) {
    MultiWayClassificationResult(Y, ret[,1], ret[,2:ncol(ret)])
  }
}

fitModel <- function(model, dset, voxelGrid, ncores=2, tuneGrid=NULL, tuneLength=NULL) {
  M <- series(dset$trainVec, voxelGrid)
  
  result <- crossval(as.matrix(M), dset$Y, split(1:length(dset$blockVar), dset$blockVar), model, ncores, tuneGrid,tuneLength)
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
