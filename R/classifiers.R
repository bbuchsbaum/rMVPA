
#' @export
MVPAModels <- list()


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
MVPAModels$wsrf <- list(type = "Classification", 
                             library = "wsrf", 
                             label="wsrf",
                             loop = NULL, 
                             parameters=data.frame(parameters=c("mtry"), class=c("numeric"), labels=c("number of variables at each split")),
                             grid=function(x, y, len = NULL) {
                               data.frame(mtry=8)
                             },
                             fit=function(x, y, wts, param, lev, last, weights, classProbs, ...) {
                               form <- as.formula(paste(y, "~ ."))
                               model.wsrf.1 <- wsrf(form, data=as.data.frame(x))
                             },
                             predict=function(modelFit, newdata, preProc = NULL, submodels = NULL) {
                               predict(modelFit, as.matrix(newdata), type="response")
                             },
                             prob=function(modelFit, newdata, preProc = NULL, submodels = NULL) {
                               predict(modelFit, as.matrix(newdata), type="prob")
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
                               predict(modelFit, as.matrix(newdata), prob=TRUE)$probabilities
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

MVPAModels$xgboost <- list(type = "Classification", 
                              library = "xgboost", 
                              loop = NULL, 
                              parameters=data.frame(parameters=c("nrounds", "max.depth", "eta"), class=c("numeric", "numeric", "numeric"), 
                                                    label=c("max number of iterations", "maximum depth of tree", "step size")),
                              grid=function(x, y, len = NULL) expand.grid(nrounds = seq(1,len),
                                                                          max.depth = seq(1,len),
                                                                          eta=1),
                              fit=function(x, y, wts, param, lev, last, weights, classProbs, ...) {
                                if (length(levels(y)) > 2) {
                                  stop("xgboost does not support multiclass classification")
                                }
                                y0 <- ifelse(y == levels(y)[1], 1, 0)
                                param$objective <- "binary:logistic"
                                xgfit <- xgboost(x,y0,params=param, nrounds=param$nrounds)
                                mfit <- list(xgfit=xgfit, levs=levels(y))
                                
                              },
                              predict=function(modelFit, newdata, preProc = NULL, submodels = NULL) {
                                p <- predict(modelFit$xgfit, as.matrix(newdata))
                                ifelse(p > .5, modelFit$levs[1], modelFit$levs[2])
                              },
                              prob=function(modelFit, newdata, preProc = NULL, submodels = NULL) {
                                p <- predict(modelFit$xgfit, as.matrix(newdata))
                                ret <- cbind(p, 1-p)
                                colnames(ret) <- modelFit$levs
                                ret
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
                              fit=function(x, y, wts, param, lev, last, weights, classProbs, ...) { sparsediscrim::lda_thomaz(x,y, ...) },
                              predict=function(modelFit, newdata, preProc = NULL, submodels = NULL) { predict(modelFit, as.matrix(newdata))$class },
                              prob=function(modelFit, newdata, preProc = NULL, submodels = NULL) { 
                                scores <- -t(predict(modelFit, newdata)$scores)
                                mc <- scores[cbind(1:nrow(scores), max.col(scores, ties.method = "first"))]
                                probs <- exp(scores - mc)
                                zapsmall(probs/rowSums(probs))
                              })


MVPAModels$pls_rf <- list(label = "PLS-RF",
                                    library = c("pls", "randomForest"),
                                    type = "Regression",
                                    ## Tune over both parameters at the same time
                                    parameters = data.frame(parameter = c('ncomp', 'mtry'),
                                                            class = c("numeric", 'numeric'),
                                                            label = c('#Components',
                                                                      '#Randomly Selected Predictors')),
                                    grid = function(x, y, len = NULL) {
                                      grid <- expand.grid(ncomp = seq(1, min(ncol(x) - 1, len), by = 1),
                                                          mtry = 1:len)
                                      ## We can't have mtry > ncomp
                                      grid <- subset(grid, mtry <= ncomp)
                                    },
                                    loop = NULL,
                                    fit = function(x, y, wts, param, lev, last, classProbs, ...) {
                                      ## First fit the pls model, generate the training set scores,
                                      ## then attach what is needed to the random forest object to 
                                      ## be used later
                                      
                                      ## plsr only has a formula interface so create one data frame
                                      dat <- if(!is.data.frame(x)) x <- as.data.frame(x)
                                      dat$y <- y
                                      pre <- plsr(y~ ., data = dat, ncomp = param$ncomp)
                                      scores <- predict(pre, x, type = "scores")
                                      colnames(scores) <- paste("score", 1:param$ncomp, sep = "")
                                      mod <- randomForest(scores, y, mtry = param$mtry, ...)
                                      mod$projection <- pre$projection
                                      mod
                                    },
                                    predict = function(modelFit, newdata, submodels = NULL) {
                                      ## Now apply the same scaling to the new samples
                                      scores <- as.matrix(newdata)  %*% modelFit$projection
                                      colnames(scores) <- paste("score", 1:ncol(scores), sep = "")
                                      scores <- as.data.frame(scores)
                                      ## Predict the random forest model
                                      predict(modelFit, scores)
                                    },
                                    prob = NULL,
                                    varImp = NULL,
                                    predictors = function(x, ...) rownames(x$projection),
                                    levels = function(x) x$obsLevels,
                                    sort = function(x) x[order(x[,1]),])




