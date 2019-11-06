
## use R6 to wrap models
## 

colHuber <- function(x,k=1.5, tol=1e-04) {
  mu <- matrixStats::colMedians(x)
  s <- matrixStats::colMads(x)
  n <- nrow(x)
  sapply(1:length(mu), function(i) {
    repeat {
      yy <- pmin(pmax(mu[i] - k * s[i], x[,i]), mu[i] + k * s[i])
      mu1 <- sum(yy)/n
      if (abs(mu[i] - mu1) < tol * s[i]) 
        break
      mu <- mu1
    }
    mu
  })
}




#' @noRd
MVPAModels <- new.env()



#' @importFrom MASS huber
#' @importFrom stats median
#' @noRd
corsimFit <- function(x, y, method, robust) {
  estimator <- if (robust) {
    function(vec)  {
      h <- try(MASS::huber(vec))
      if (inherits(h, "try-error")) {
        median(vec)
      } else {
        h$mu
      }
    }      
  } else {
    "mean"
  }
  
  if (identical("mean", estimator)) {
    list(conditionMeans=group_means(x, 1, y), levs=levels(y), method=method, robust=robust)
  } else {
    list(conditionMeans = neuroim2::split_reduce(as.matrix(x), y, estimator), levs=levels(y), method=method, robust=robust)
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
  #probs <- exp(scores - mc)
  probs <- exp(sweep(scores, 1, mc, "-"))
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
                             scmat <- scale(as.matrix(x))
                             pres <- svd::propack.svd(scmat, neig=param$ncomp)
                             #pres <- prcomp(as.matrix(x), scale=TRUE)
                             #lda.fit <- lda(pres$x[,1:param$ncomp, drop=FALSE], y)
                             lda.fit <- lda(pres$u[, 1:param$ncomp, drop=FALSE], y)
                             attr(lda.fit, "ncomp") <- param$ncomp
                             attr(lda.fit, "pcfit") <- pres
                             attr(lda.fit, "center") <- attr(scmat, "scaled:center")
                             attr(lda.fit, "scale") <- attr(scmat, "scaled:scale")
                             lda.fit
                           },
                           
                           predict=function(modelFit, newdata, preProc = NULL, submodels = NULL) {
                             compind <- seq(1, attr(modelFit, "ncomp"))
                             
                             pcfit <- attr(modelFit, "pcfit")
                             #colnames(newdata) <- rownames(pcfit$rotation)
                             #pcx <- predict(pcfit, newdata)[,compind,drop=FALSE]
                             pcx <- scale(newdata, attr(pcfit, "center"), attr(pcfit, "scale")) %*% pcfit$v
                             predict(modelFit, pcx)$class
                           },
                           prob=function(modelFit, newdata, preProc = NULL, submodels = NULL) {
                             compind <- seq(1, attr(modelFit, "ncomp"))
                             pcfit <- attr(modelFit, "pcfit")
                             #colnames(newdata) <- rownames(pcfit$rotation)
                             #pcx <- predict(pcfit, newdata)[,compind,drop=FALSE]
                             pcx <- scale(newdata, attr(modelFit, "center"), attr(modelFit, "scale")) %*% pcfit$v
                             predict(modelFit, pcx)$posterior                              
                           })




MVPAModels$corclass <- list(type = "Classification", 
                          library = "rMVPA", 
                          label="corclass",
                          loop = NULL, 
                          parameters=data.frame(parameters=c("method", "robust"), class=c("character", "logical"), label=c("correlation type: pearson, spearman, or kendall", "mean or huber")),
                          grid=function(x, y, len = NULL) if (is.null(len) || len == 1) { data.frame(method="pearson", robust=FALSE) } else { expand.grid(method=c("pearson", "spearman", "kendall"), 
                                                                                                                                          robust=c(TRUE, FALSE)) },
                          fit=function(x, y, wts, param, lev, last, weights, classProbs, ...) corsimFit(x,y, as.character(param$method), param$robust),
                          
                          predict=function(modelFit, newdata, preProc = NULL, submodels = NULL) predict_corsimFit(modelFit, as.matrix(newdata)),
                          prob=function(modelFit, newdata, preProc = NULL, submodels = NULL) {
                            prob_corsimFit(modelFit, as.matrix(newdata))
                          })

MVPAModels$corsim <- MVPAModels$corclass

MVPAModels$sda_notune <- list(type = "Classification", 
                              library = "sda", 
                              label="sda_notune",
                              loop = NULL, 
                              parameters=data.frame(parameters="parameter", class="character", label="parameter"),
                              grid=function(x, y, len = NULL) data.frame(parameter="none"),
                              fit=function(x, y, wts, param, lev, last, weights, classProbs, ...) {                       
                                sda::sda(Xtrain=as.matrix(x), L=y, verbose=FALSE, ...)
                              },
                              predict=function(modelFit, newdata, preProc = NULL, submodels = NULL) {
                                
                                predict(modelFit, as.matrix(newdata), verbose=FALSE)$class
                              },
                              prob=function(modelFit, newdata, preProc = NULL, submodels = NULL) {
            
                                predict(modelFit, as.matrix(newdata),verbose=FALSE)$posterior
                              })


MVPAModels$sda_boot <- list(type = "Classification", 
                              library = "sda", 
                              label="sda_boot",
                              loop = NULL, 
                              parameters=data.frame(parameters=c("reps", "frac", "lambda_min", "lambda_max"), 
                                                    class=c("numeric", "numeric", "numeric", "numeric"), 
                                                    label=c("number of bootstap resamples", "fraction of features to select", "min lambda", "max lambda")),
                              grid=function(x, y, len = NULL) data.frame(reps=10, frac=1, lambda_min=.01, lambda_max=.8),
                              fit=function(x, y, wts, param, lev, last, weights, classProbs, ...) {   
                                
                                x <- as.matrix(x)
                                mfits <- list()
                                count <- 1
                                failures <- 0
                                
                                split_y <- split(seq_along(y), y)
                                lambda <- seq(param$lambda_min, param$lambda_max, length.out=param$reps)
                                
                                if (!all(sapply(split_y, function(yy) length(yy > 0)))) {
                                  stop("every factor level in 'y' must have at least 1 instance, cannot run sda_boot")
                                }
                                
                                assertthat::assert_that(param$frac > 0, msg="sda_boot: 'frac' parameter must be greater than 0")
                                assertthat::assert_that(param$reps > 0, msg="sda_boot: 'reps' parameter must be greater than 0")
                                
                                while (count <= param$reps) {
                                  #message("fitting sda model ", count)
                                  
                                  ysam <- lapply(split_y, function(idx) if (length(idx) == 1) idx else sample(idx, length(idx), replace=TRUE))
                                  row.idx <- sort(unlist(ysam))
                                  
                                  ret <- if (param$frac > 0 && param$frac < 1) {
                                    nkeep <- max(param$frac * ncol(x),1)
                                    ind <- sample(1:ncol(x), nkeep)
                                    fit <- sda::sda(Xtrain=x[row.idx,ind,drop=FALSE], L=y[row.idx], lambda=lambda[count], verbose=FALSE,...)
                                    attr(fit, "keep.ind") <- ind
                                    fit
                                  } else {
                                    fit <- try(sda::sda(Xtrain=x[row.idx,], L=y[row.idx], lambda=lambda[count], verbose=FALSE, ...))
                                    attr(fit, "keep.ind") <- 1:ncol(x)
                                    fit
                                  }
                                  
                                  if (!inherits(ret, "try-error")) {
                                    mfits[[count]] <- ret
                                    count <- count + 1
                                  } else {
                                    message("sda model fit error, retry attempt: ", failures + 1)
                                    failures <- failures + 1
                                    if (failures > 10) {
                                      return(ret)
                                    }
                                  }
                                }
                              
                          
                                ret <- list(fits=mfits)
                                class(ret) <- "sda_boot"
                                ret
                                
                              },
                              predict=function(modelFit, newdata, preProc = NULL, submodels = NULL) {
                                
                                preds <- lapply(modelFit$fits, function(fit) {
                                  ind <- attr(fit, "keep.ind")
                                  predict(fit, as.matrix(newdata)[,ind], verbose=FALSE)$posterior
                                })
                                
                                prob <- preds[!sapply(preds, function(x) is.null(x))]
                                pfinal <- Reduce("+", prob)/length(prob)
                                
                                colnames(pfinal)[apply(pfinal, 1, which.max)]
                            
                              },
                                
                                
                              prob=function(modelFit, newdata, preProc = NULL, submodels = NULL) {
                
                                preds <- lapply(modelFit$fits, function(fit) {
                                  ind <- attr(fit, "keep.ind")
                                  predict(fit, as.matrix(newdata)[,ind], verbose=FALSE)$posterior
                                })
                                
                                prob <- preds[!sapply(preds, function(x) is.null(x))]
                                Reduce("+", prob)/length(prob)
                                
                              })
MVPAModels$lda_thomaz_boot <- list(type = "Classification", 
                            library = "sparsediscrim", 
                            label="lda_thomaz_boot",
                            loop = NULL, 
                            parameters=data.frame(parameters=c("reps", "frac"), 
                                                  class=c("numeric", "numeric"), 
                                                  label=c("number of bootstap resamples", "fraction of features to select")),
                            grid=function(x, y, len = NULL) data.frame(reps=10, frac=1),
                            fit=function(x, y, wts, param, lev, last, weights, classProbs, ...) {   
                              
                              x <- as.matrix(x)
                              mfits <- list()
                              count <- 1
                              failures <- 0
                              
                              split_y <- split(seq_along(y), y)
                              
                              if (!all(sapply(split_y, function(yy) length(yy > 0)))) {
                                stop("every factor level in 'y' must have at least 1 instance, cannot run lda_thomaz_boot")
                              }
                              
                              assertthat::assert_that(param$frac > 0, msg="lda_thomaz_boot: 'frac' parameter must be greater than 0")
                              
                              while (count <= param$reps) {
                                #message("fitting sda model ", count)
                                
                                ysam <- lapply(split_y, function(idx) if (length(idx) == 1) idx else sample(idx, length(idx), replace=TRUE))
                                row.idx <- sort(unlist(ysam))
                                
                                ret <- if (param$frac > 0 && param$frac < 1) {
                                  nkeep <- max(param$frac * ncol(x),1)
                                  ind <- sample(1:ncol(x), nkeep)
                                  fit <- sparsediscrim::lda_thomaz(x=x[row.idx,ind,drop=FALSE], y[row.idx])
                                  attr(fit, "keep.ind") <- ind
                                  fit
                                } else {
                                  fit <- try(sparsediscrim::lda_thomaz(x[row.idx,,drop=FALSE], y[row.idx]))
                                  attr(fit, "keep.ind") <- 1:ncol(x)
                                  fit
                                }
                                
                                if (!inherits(ret, "try-error")) {
                                  mfits[[count]] <- ret
                                  count <- count + 1
                                } else {
                                  message("lda_thomaz model fit error, retry attempt: ", failures + 1)
                                  failures <- failures + 1
                                  if (failures > 10) {
                                    return(ret)
                                  }
                                }
                              }
                              
                              
                              ret <- list(fits=mfits)
                              class(ret) <- "lda_thomaz_boot"
                              ret
                              
                            },
                            predict=function(modelFit, newdata, preProc = NULL, submodels = NULL) {
                              
                              
                              preds <- lapply(modelFit$fits, function(fit) {
                                ind <- attr(fit, "keep.ind")
                                
                                scores <- -t(sparsediscrim:::predict.lda_thomaz(fit, newdata[,ind])$scores)
                                mc <- scores[cbind(1:nrow(scores), max.col(scores, ties.method = "first"))]
                                probs <- exp(scores - mc)
                                zapsmall(probs/rowSums(probs))
                              })
                              
                              prob <- preds[!sapply(preds, function(x) is.null(x))]
                              pfinal <- Reduce("+", prob)/length(prob)
                              
                              colnames(pfinal)[apply(pfinal, 1, which.max)]
                              
                            },
                            
                            
                            prob=function(modelFit, newdata, preProc = NULL, submodels = NULL) {
                              preds <- lapply(modelFit$fits, function(fit) {
                                ind <- attr(fit, "keep.ind")
                                
                                scores <- -t(sparsediscrim:::predict.lda_thomaz(fit, newdata[,ind])$scores)
                                mc <- scores[cbind(1:nrow(scores), max.col(scores, ties.method = "first"))]
                                probs <- exp(scores - mc)
                                zapsmall(probs/rowSums(probs))
                              })
                              
                              prob <- preds[!sapply(preds, function(x) is.null(x))]
                              pfinal <- Reduce("+", prob)/length(prob)
                            
                              
                            })


#' @importFrom memoise memoise
#' @importFrom sda sda.ranking
#' @noRd
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
                               parameters=data.frame(parameters=c("frac", "lambda"), class=c("numeric", "numeric"), 
                                                     label=c("fraction of features to keep (frac > 0 a frac <= 1)", "lambda")),
                               grid=function(x, y, len = NULL) expand.grid(frac=seq(.1,1,length.out=len), lambda=seq(.01,.99,length.out=len)),
                               fit=function(x, y, wts, param, lev, last, weights, classProbs, ...) {
                                 print(param)
                                 x <- as.matrix(x)                          
                                 nkeep <- max(param$frac * ncol(x),1)
                                 
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
