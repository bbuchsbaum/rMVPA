
## use R6 to wrap models
## 




#' @export
MVPAModels <- new.env()


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
    "mean"
  }
  
  if (estimator == "mean") {
    list(conditionMeans=groupMeans(X, 1, Y), levs=levels(y), method=method, robust=robust)
  } else {
    list(conditionMeans = splitReduce(as.matrix(x), y, estimator), levs=levels(y), method=method, robust=robust)
  }
}


predict_corsimFit <- function(modelFit, newData) {
  cres <- cor(t(newData), t(modelFit$conditionMeans), method=modelFit$method)
  res <- max.col(cres)
  modelFit$levs[res]
}


prob_corsimFit <- function(modelFit, newData) {
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
                           label="pca_lda",
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

MVPAModels$clusterSVM <- list(type = "Classification", 
                           library = "SwarmSVM", 
                           loop = NULL, 
                           label="clusterSVM",
                           parameters=data.frame(parameters=c("K", "lambda", "cost"), 
                                                 class=c("numeric", "numeric", "numeric"), 
                                                 labels=c("number of clusters", "global regulzariaztion", "svm cost")),
                           grid=function(x, y, len = 5) {
                             expand.grid(K=1:len, lambda=seq(0,2, length.out=len), cost=1)
                           },
                           fit=function(x, y, wts, param, lev, last, weights, classProbs, ...) {
                             SwarmSVM::clusterSVM(x=as.matrix(x), y=y, centers=param$K, lambda=param$lambda, cost=param$cost, type=0)
                           },
                           
                           predict=function(modelFit, newdata, preProc = NULL, submodels = NULL) {
                             pred = predict(modelFit, as.matrix(newdata))$predictions
                           },
                           prob=function(modelFit, newdata, preProc = NULL, submodels = NULL) {
                             pred = predict(modelFit, as.matrix(newdata), prob=TRUE)$probabilities
                                                        
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
                               predict(modelFit, as.matrix(newdata), prob=TRUE)$probabilities
                             })


MVPAModels$nearestMean <- list(type = "Classification", 
                               library = "klaR", 
                               label="nearestMean",
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

MVPAModels$corclass <- list(type = "Classification", 
                          library = "rMVPA", 
                          label="corclass",
                          loop = NULL, 
                          parameters=data.frame(parameters=c("method", "robust"), class=c("character", "logical"), label=c("correlation type: pearson, spearman, or kendall", "mean or huber")),
                          grid=function(x, y, len = NULL) if (len == 1) { data.frame(method="pearson", robust=FALSE) } else { expand.grid(method=c("pearson", "spearman", "kendall"), robust=c(TRUE, FALSE)) },
                          fit=function(x, y, wts, param, lev, last, weights, classProbs, ...) corsimFit(x,y, as.character(param$method), param$robust),
                          predict=function(modelFit, newdata, preProc = NULL, submodels = NULL) predict_corsimFit(modelFit, as.matrix(newdata)),
                          prob=function(modelFit, newdata, preProc = NULL, submodels = NULL) {
                            prob_corsimFit(modelFit, as.matrix(newdata))
                          })


MVPAModels$sda_notune <- list(type = "Classification", 
                              library = "sda", 
                              label="sda_notune",
                              loop = NULL, 
                              parameters=data.frame(parameters="parameter", class="character", label="parameter"),
                              grid=function(x, y, len = NULL) data.frame(parameter="none"),
                              fit=function(x, y, wts, param, lev, last, weights, classProbs, ...) {                       
                                sda::sda(Xtrain=as.matrix(x), L=y, verbose=FALSE, ...)
                              },
                              predict=function(modelFit, newdata, preProc = NULL, submodels = NULL) predict(modelFit, as.matrix(newdata), verbose=FALSE)$class,
                              prob=function(modelFit, newdata, preProc = NULL, submodels = NULL) {
                                predict(modelFit, as.matrix(newdata),verbose=FALSE)$posterior
                              })


MVPAModels$sda_boot <- list(type = "Classification", 
                              library = "sda", 
                              label="sda_boot",
                              loop = NULL, 
                              parameters=data.frame(parameters=c("reps"), class=c("numeric"), label=c("number of bootstap resamples")),
                              grid=function(x, y, len = NULL) data.frame(reps=10, frac=1),
                              fit=function(x, y, wts, param, lev, last, weights, classProbs, ...) {    
                                x <- as.matrix(x)
                                mfits <- lapply(seq(1, param$reps), function(i) {
                                  message("fitting sda model ", i)
                                  row.idx <- sample(1:nrow(x), nrow(x), replace=TRUE)
                                  fit <- sda::sda(Xtrain=x[row.idx,], L=y[row.idx], verbose=FALSE, ...)
                                
                                })
                                
                              },
                              predict=function(modelFit, newdata, preProc = NULL, submodels = NULL) {
                                preds <- lapply(modelFit, function(fit) {
                                  predict(fit, as.matrix(newdata), verbose=FALSE)$posterior
                                })
                                
                                prob <- preds[!sapply(preds, function(x) is.null(x))]
                                pfinal <- Reduce("+", prob)/length(prob)
                                
                                colnames(pfinal)[apply(pfinal, 1, which.max)]
                              },
                                
                                
                              prob=function(modelFit, newdata, preProc = NULL, submodels = NULL) {
                                preds <- lapply(modelFit, function(fit) {
                                  predict(fit, as.matrix(newdata), verbose=FALSE)$posterior
                                })
                                
                                prob <- preds[!sapply(preds, function(x) is.null(x))]
                                pfinal <- Reduce("+", prob)/length(prob)
                                
                              })

#' @export
#' @import memoise
#' @import sda
memo_rank <- memoise(function(X, L,fdr) {
  sda::sda.ranking(X,L,fdr=fdr,verbose=FALSE)
})



MVPAModels$glmnet_opt <- list(type = "Classification", 
                              library = c("c060", "glmnet"), 
                              label="glmnet_opt",
                              loop = NULL, 
                              parameters=data.frame(parameters="parameter", class="character", label="parameter"),
                              grid=function(x, y, len = NULL) data.frame(parameter="none"),
                              fit=function(x, y, wts, param, lev, last, weights, classProbs, ...) {
                                bounds <- t(data.frame(alpha=c(0, 1)))
                                fam <- if (length(levels(y) > 2)) {
                                  "multinomial"
                                } else {
                                  "binomial"
                                }
                                colnames(bounds)<-c("lower","upper")
                                x <- as.matrix(x)                          
                                fit <- epsgo(Q.func="tune.glmnet.interval",
                                             bounds=bounds,
                                             parms.coding="none",
                                             seed = 1234,
                                             show="none",
                                             fminlower = -100,
                                             x = x, y = y, family = fam,
                                             type.min = "lambda.1se",
                                             foldid=createFolds(y,k=5,list=FALSE),
                                             type.measure = "mse")
                                sfit <- summary(fit, verbose=FALSE)
                                fit <- glmnet(x,y, family=fam,alpha=sfit$opt.alpha, nlambda=20)
                                fit$opt_lambda <- sfit$opt.lambda
                                fit$opt_alpha <- sfit$opt.alpha
                                fit
                              },
                              predict=function(modelFit, newdata, preProc = NULL, submodels = NULL) {                        
                                out <- predict(modelFit, as.matrix(newdata), s=modelFit$opt_lambda, type="class")
                                if (is.matrix(out)) {
                                  out[,1]
                                } else {
                                  out
                                }
                              },
                              prob=function(modelFit, newdata, preProc = NULL, submodels = NULL) {
                                obsLevels <- if("classnames" %in% names(modelFit)) modelFit$classnames else NULL
                                if(length(obsLevels) == 2) {
                                  probs <- as.vector(probs)
                                  probs <- as.data.frame(cbind(1-probs, probs))
                                  colnames(probs) <- modelFit$obsLevels
                                } else {
                                  probs <- as.data.frame(probs[,,1,drop = FALSE])
                                  names(probs) <- modelFit$obsLevels
                                }
                                
                                probs
                              })



MVPAModels$sparse_sda <- list(type = "Classification", 
                               library = "sda", 
                               label="sparse_sda",
                               loop = NULL, 
                               parameters=data.frame(parameters=c("frac", "lambda"), class=c("numeric", "numeric"), label=c("fraction of features to keep (frac > 0 a frac <= 1)", "lambda")),
                               grid=function(x, y, len = NULL) expand.grid(frac=seq(.1,1,length.out=len), lambda=seq(.01,.99,length.out=len)),
                               fit=function(x, y, wts, param, lev, last, weights, classProbs, ...) {
                                 
                                 x <- as.matrix(x)                          
                                 nkeep <- max(param$frac * ncol(x),1)
                                 print(nkeep)
                                 rank <- memo_rank(x, L=y, fdr=FALSE)
                                 ind <- rank[,"idx"][1:nkeep]
                                 
                                 fit <- sda::sda(Xtrain=x[,ind,drop=FALSE], L=y, lambda=param$lambda, verbose=FALSE)
                                 attr(fit, "keep.ind") <- ind
                                 fit
                               },
                               predict=function(modelFit, newdata, preProc = NULL, submodels = NULL) {                        
                                 predict(modelFit, as.matrix(newdata[,attr(modelFit, "keep.ind"), drop=FALSE]), verbose=FALSE)$class
                               },
                               prob=function(modelFit, newdata, preProc = NULL, submodels = NULL) {
                                 predict(modelFit, as.matrix(newdata[,attr(modelFit, "keep.ind"), drop=FALSE]),verbose=FALSE)$posterior
                               })


MVPAModels$sda_ranking <- list(type = "Classification", 
                               library = "sda", 
                               label="sda_ranking",
                               loop = NULL, 
                               parameters=data.frame(parameters="parameter", class="character", label="parameter"),
                               grid=function(x, y, len = NULL) data.frame(parameter="none"),
                               fit=function(x, y, wts, param, lev, last, weights, classProbs, ...) {
                                 
                                 x <- as.matrix(x)                             
                                 ind <- if (ncol(x) > 2) {
                                   ## rank features
                                   rank <- sda::sda.ranking(Xtrain=x, L=y, fdr=TRUE, verbose=FALSE, ...)
                                   
                                   #thresh based on higher criticism
                                   hcind <- which.max(rank[,"HC"])
                                   
                                   
                                   keep.ind <- if (hcind < 2) {
                                     seq(1, min(ncol(x), 2))
                                   } else {
                                     1:hcind                               
                                   }   
                                   rank[keep.ind,"idx"]
                                } else {
                                  ## keep everything
                                  1:ncol(x)
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

MVPAModels$mgsda <- list(type = "Classification", 
                              library = "MGSDA", 
                              loop = NULL, 
                              parameters=data.frame(parameters="lambda", class="numeric", label="sparsity penalty"),
                              grid=function(x, y, len = NULL) data.frame(lambda=seq(.001, .99, by=len)),
                              fit=function(x, y, wts, param, lev, last, weights, classProbs, ...) { 
                                ycodes <- as.integer(y)
                                V <- MGSDA::dLDA(x,ycodes, lambda= param$lambda, ...) 
                                modelFit$ycodes <- ycodes
                                modelFit$V <- V
                                modelFit$xtrain <- xtrain
                                modelFit$ytrain <- y
                                modelFit
                              },
                              predict=function(modelFit, newdata, preProc = NULL, submodels = NULL) {  
                                preds <- MSGDA::classifiyV(modelFit$xtrain, modelFit$ycodes, newdata, modelFit$V)
                                levels(modelFit$ytrain)[preds]
                              },
                              prob=function(modelFit, newdata, preProc = NULL, submodels = NULL) { 
                               
                              })

MVPAModels$lda_thomaz <- list(type = "Classification", 
                              library = "sparsediscrim", 
                              loop = NULL, 
                              parameters=data.frame(parameters="parameter", class="character", label="parameter"),
                              grid=function(x, y, len = NULL) data.frame(parameter="none"),
                              fit=function(x, y, wts, param, lev, last, weights, classProbs, ...) { sparsediscrim:::lda_thomaz(x,y, ...) },
                              predict=function(modelFit, newdata, preProc = NULL, submodels = NULL) { sparsediscrim:::predict.lda_thomaz(modelFit, as.matrix(newdata))$class },
                              prob=function(modelFit, newdata, preProc = NULL, submodels = NULL) { 
                                scores <- -t(sparsediscrim:::predict.lda_thomaz(modelFit, newdata)$scores)
                                mc <- scores[cbind(1:nrow(scores), max.col(scores, ties.method = "first"))]
                                probs <- exp(scores - mc)
                                zapsmall(probs/rowSums(probs))
                              })

MVPAModels$hdrda <- list(type = "Classification", 
                              library = "sparsediscrim", 
                              label="hdrda",
                              loop = NULL, 
                              parameters=data.frame(parameters=c("lambda", "gamma"), class=c("numeric", "numeric"), label=c("HRDA pooling parameter", "shrinkage parameter")),
                              grid=function(x, y, len = NULL) expand.grid(lambda=seq(.99, .001, length.out=len), gamma=seq(.001, .99, length.out=len)),
                              fit=function(x, y, wts, param, lev, last, weights, classProbs, ...) { sparsediscrim:::hdrda(x,y, lambda=param$lambda, gamma=param$gamma,...) },
                              predict=function(modelFit, newdata, preProc = NULL, submodels = NULL) { sparsediscrim:::predict.hdrda(modelFit, as.matrix(newdata))$class },
                              prob=function(modelFit, newdata, preProc = NULL, submodels = NULL) { 
                                posterior <- sparsediscrim:::predict.hdrda(modelFit, newdata)$posterior
                                t(apply(posterior,1, function(x) x/sum(x)))
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



# MVPAModels$gmd <- list(type = "Classification", 
#                        library = "GMD", 
#                        label="histogram distance",
#                        loop=NULL,
#                        parameters=data.frame(parameters="k", class="integer", label="number of nearest neighbors"),
#                        grid=function(x, y, len = NULL) data.frame(k=5),
#                        fit=function(x, y, wts, param, lev, last, weights, classProbs, ...) {
#                          
#                          #breaks <- gbreaks(as.vector(x),25)
#                          #xh <- as.mhist(lapply(1:nrow(x), function(i) ghist(x[i,], breaks=breaks)))
#                          #xh <- apply(x,1,function(vals) ghist(vals, n=20), simplify=FALSE)
#                          #ind <- expand.grid(i=1:length(xh), j=1:length(xh))
#                          #Dmat <- mapply(function(i,j) {
#                          #  .Call("gmd0", xh[[i]], xh[[j]], 0, PACKAGE = "GMD")
#                          #}, ind[,1], ind[,2])
#                          list(x=x, y=y, levs=levels(y), k=param$k, breaks = gbreaks(as.vector(x),25))
#                        },
#                        predict=function(modelFit, newdata, preProc = NULL, submodels = NULL) {       
#                          ret <- unlist(lapply(1:nrow(newdata), function(i) {
#                            vals <- newdata[i,]
#                            hvals <- ghist(vals, breaks=modelFit$breaks)
#                            D <- apply(modelFit$x, 1, function(tvals) {
#                              #.Call("gmd0", hvals, ghist(tvals, breaks=modelFit$breaks), 0, PACKAGE = "GMD")
#                              gmdp(hvals, ghist(tvals, breaks=modelFit$breaks), sliding=FALSE)
#                            })
#                            
#                            #ind <- order(D)[1:modelFit$k]
#                            minD <- min(D)
#                            agg <- aggregate(D ~ modelFit$y, FUN=median)
#                            ord <- order(agg[,2])
#                            modelFit$levs[ord[1]]
#                            #names(which.max(table(modelFit$y[ind])))[1]
#                          }))
#                          
#                        },
#                        prob=function(modelFit, newdata, preProc = NULL, submodels = NULL) {
#                          
#                          #predict(modelFit, as.matrix(newdata[,attr(modelFit, "keep.ind"), drop=FALSE]),verbose=FALSE)$posterior
#                        })
