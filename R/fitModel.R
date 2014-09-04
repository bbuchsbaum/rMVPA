

MVPA.lda_strimmer <- list(type = "Classification", 
                 library = "sda", 
                 loop = NULL, 
                 parameters=data.frame(parameters="lambda", class="numeric", labels="lambda"),
                 grid=function(x, y, len = NULL) data.frame(lambda=seq(.2, .9, length.out=len)),
                 fit=function(x, y, wts, param, lev, last, weights, classProbs, ...) sda(Xtrain=as.matrix(x), L=y, ...),
                 predict=function(modelFit, newdata, preProc = NULL, submodels = NULL) predict(modelFit, as.matrix(newdata))$class,
                 prob=function(modelFit, newdata, preProc = NULL, submodels = NULL) {
                    predict(modelFit, as.matrix(newdata))$posterior
                })


MVPA.lda_thomaz <- list(type = "Classification", 
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




#' create a classification model
createModel <- function(method) {
  ret <- list(method=method)
  class(ret) <- c(method, "list")
  ret
}


trainModel <- function(x, ...) {
  UseMethod("trainModel")
}


crossval <- function(X, Y, foldSplit, method, ncores=2) {
  ctrl <- caret::trainControl("none")
  res <- lapply(foldSplit, function(fidx) {
    Xtrain <- X[-fidx,]
    Ytrain <- Y[-fidx]
    Xtest <- X[fidx,]   
    fit <- caret::train(Xtrain, Ytrain, method=method, trControl=ctrl, tuneLength=1)
    cbind(class=predict(fit, newdata=Xtest), predict(fit, newdata=Xtest, type="prob"))
  })
  
  ret <- do.call(rbind, res)
}

fitModel <- function(model, dset, voxelGrid, ncores=2) {
  M <- series(dset$trainVec, voxelGrid)
  fitlist <- crossval(M, dset$Y, split(1:length(dset$blockVar), dset$blockVar), model)
}

trainModel.default <- function(model, dset, voxelGrid, ncores=2) {
  M <- series(dset$trainVec, voxelGrid)
  fitlist <- crossval(M, dset$Y, dset$blockVar, model)
}



trainModel.lda_strimmer <- function(dset, voxelGrid, ncores=2) {
  M <- series(dset$trainVec, voxelGrid)
  ctrl <- caret::trainControl("none")
  fit <- caret::train(x$x, x$y, method=MVPA.sda, trControl=ctrl, tuneLength=1)
}

trainModel.lda_thomaz <- function(x) {
  ctrl <- caret::trainControl("none")
  fit <- caret::train(x$x, x$y, method=MVPA.lda_thomaz, trControl=ctrl, tuneLength=0)
}
