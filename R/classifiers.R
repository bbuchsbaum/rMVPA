#' @keywords internal
#' @noRd
colHuber <- function(x, k=1.5, tol=1e-04) {
  mu <- matrixStats::colMedians(x)
  s <- matrixStats::colMads(x)
  n <- nrow(x)
  sapply(seq_along(mu), function(i) {
    repeat {
      yy <- pmin(pmax(mu[i] - k * s[i], x[, i]), mu[i] + k * s[i])
      mu1 <- sum(yy)/n
      if (abs(mu[i] - mu1) < tol * s[i])
        break
      mu[i] <- mu1
    }
    mu[i]
  })
}


#' Pre-defined MVPA Classification Models
#'
#' An environment containing custom classification models for MVPA analysis.
#'
#' @format An environment with the following models:
#' \describe{
#'   \item{corclass}{Correlation-based classifier using template matching with options (pearson, spearman, kendall)}
#'   \item{corsim}{Alias for corclass}
#'   \item{sda_notune}{Shrinkage Discriminant Analysis (SDA) without parameter tuning}
#'   \item{sda_boot}{SDA with bootstrap resampling and feature selection}
#'   \item{glmnet_opt}{Elastic net classifier (glmnet) with optimized alpha/lambda via EPSGO}
#'   \item{sparse_sda}{SDA with sparsity constraints and feature selection}
#'   \item{sda_ranking}{SDA with feature ranking and selection via higher criticism}
#'   \item{mgsda}{Multi-Group Sparse Discriminant Analysis}
#'   \item{lda_thomaz}{Modified LDA classifier for high-dimensional data}
#'   \item{hdrda}{High-Dimensional Regularized Discriminant Analysis}
#' }
#'
#' @details
#' Models are accessed via \code{load_model(name)}. Each implements caret's \code{fit}, \code{predict}, and \code{prob} methods.
#'
#' @examples
#' # Load simple SDA classifier
#' model <- load_model("sda_notune")
#'
#' # Load correlation classifier
#' model <- load_model("corclass")
#'
#' @export
MVPAModels <- new.env()

#' @importFrom MASS huber
#' @importFrom stats median
#' @noRd
corsimFit <- function(x, y, method, robust) {
  estimator <- if (robust) {
    function(vec) {
      h <- try(MASS::huber(vec), silent=TRUE)
      if (inherits(h, "try-error")) median(vec) else h$mu
    }
  } else {
    "mean"
  }

  lev <- levels(y)
  if (identical("mean", estimator)) {
    list(conditionMeans=group_means(x, 1, y), levs=lev, method=method, robust=robust)
  } else {
    list(conditionMeans = neuroim2::split_reduce(as.matrix(x), y, estimator), levs=lev, method=method, robust=robust)
  }
}

#' @keywords internal
#' @noRd
predict_corsimFit <- function(modelFit, newData) {
  cres <- cor(t(newData), t(modelFit$conditionMeans), method=modelFit$method)
  res <- max.col(cres)
  modelFit$levs[res]
}

#' @keywords internal
#' @noRd
prob_corsimFit <- function(modelFit, newData) {
  scores <- cor(t(newData), t(modelFit$conditionMeans), method=modelFit$method)
  mc <- scores[cbind(seq_len(nrow(scores)), max.col(scores, ties.method = "first"))]
  probs <- exp(sweep(scores, 1, mc, "-"))
  probs <- zapsmall(probs / rowSums(probs))
  colnames(probs) <- modelFit$levs
  probs
}

# PCA + LDA model
# Store lev after training so predictions can reference classes.
#' @keywords internal
#' @noRd
MVPAModels$pca_lda <- list(
  type = "Classification",
  library = "MASS",
  loop = NULL,
  label="pca_lda",
  parameters=data.frame(parameters="ncomp", class="numeric", labels="ncomp"),
  grid=function(x, y, len = 5) data.frame(ncomp=1:len),

  fit=function(x, y, wts, param, lev, last, weights, classProbs, ...) {
    scmat <- scale(as.matrix(x))
    pres <- svd::propack.svd(scmat, neig=param$ncomp)
    lda.fit <- lda(pres$u[, 1:param$ncomp, drop=FALSE], y)
    attr(lda.fit, "ncomp") <- param$ncomp
    attr(lda.fit, "pcfit") <- pres
    attr(lda.fit, "center") <- attr(scmat, "scaled:center")
    attr(lda.fit, "scale") <- attr(scmat, "scaled:scale")
    attr(lda.fit, "obsLevels") <- lev
    lda.fit
  },

  predict=function(modelFit, newdata, preProc = NULL, submodels = NULL) {
    compind <- seq_len(attr(modelFit, "ncomp"))
    pcfit <- attr(modelFit, "pcfit")
    pcx <- scale(newdata, attr(modelFit, "center"), attr(modelFit, "scale")) %*% pcfit$v
    predict(modelFit, pcx[, compind, drop=FALSE])$class
  },

  prob=function(modelFit, newdata, preProc = NULL, submodels = NULL) {
    compind <- seq_len(attr(modelFit, "ncomp"))
    pcfit <- attr(modelFit, "pcfit")
    pcx <- scale(newdata, attr(modelFit, "center"), attr(modelFit, "scale")) %*% pcfit$v
    predict(modelFit, pcx[, compind, drop=FALSE])$posterior
  }
)

# corclass and corsim
# Ensure lev is known from fit if needed. corsimFit stores lev internally.
#' @keywords internal
#' @noRd
MVPAModels$corclass <- list(
  type = "Classification",
  library = "rMVPA",
  label="corclass",
  loop = NULL,
  parameters=data.frame(parameters=c("method", "robust"), class=c("character", "logical"),
                        label=c("correlation type: pearson, spearman, or kendall", "mean or huber")),
  grid=function(x, y, len = NULL) {
    if (is.null(len) || len == 1) {
      data.frame(method="pearson", robust=FALSE)
    } else {
      expand.grid(method=c("pearson","spearman","kendall"), robust=c(TRUE,FALSE))
    }
  },

  fit=function(x, y, wts, param, lev, last, weights, classProbs, ...) {
    corsimFit(x, y, as.character(param$method), param$robust)
  },

  predict=function(modelFit, newdata, preProc = NULL, submodels = NULL) {
    predict_corsimFit(modelFit, as.matrix(newdata))
  },

  prob=function(modelFit, newdata, preProc = NULL, submodels = NULL) {
    prob_corsimFit(modelFit, as.matrix(newdata))
  }
)

#' @keywords internal
#' @noRd
MVPAModels$corsim <- MVPAModels$corclass

# sda_notune
# Add levels storage
#' @keywords internal
#' @noRd
MVPAModels$sda_notune <- list(
  type = "Classification",
  library = "sda",
  label="sda_notune",
  loop = NULL,
  parameters=data.frame(parameters="parameter", class="character", label="parameter"),
  grid=function(x, y, len = NULL) data.frame(parameter="none"),

  fit=function(x, y, wts, param, lev, last, weights, classProbs, ...) {
    m <- sda::sda(Xtrain=as.matrix(x), L=y, verbose=FALSE, ...)
    m$obsLevels <- lev
    m
  },

  predict=function(modelFit, newdata, preProc = NULL, submodels = NULL) {
    predict(modelFit, as.matrix(newdata), verbose=FALSE)$class
  },

  prob=function(modelFit, newdata, preProc = NULL, submodels = NULL) {
    predict(modelFit, as.matrix(newdata), verbose=FALSE)$posterior
  }
)

# sda_boot
# Store lev and ensure predictions label correctly
#' @keywords internal
#' @noRd
MVPAModels$sda_boot <- list(
  type = "Classification",
  library = "sda",
  label="sda_boot",
  loop = NULL,
  parameters=data.frame(parameters=c("reps","frac","lambda_min","lambda_max"),
                        class=c("numeric","numeric","numeric","numeric"),
                        label=c("boot reps","frac features","min lambda","max lambda")),
  grid=function(x, y, len = NULL) data.frame(reps=10, frac=1, lambda_min=.01, lambda_max=.8),

  fit=function(x, y, wts, param, lev, last, weights, classProbs, ...) {
    x <- as.matrix(x)
    mfits <- list()
    count <- 1
    failures <- 0

    split_y <- split(seq_along(y), y)
    lambda <- seq(param$lambda_min, param$lambda_max, length.out=param$reps)

    if (!all(sapply(split_y, function(yy) length(yy) > 0))) {
      stop("Every level in 'y' must have at least one instance for sda_boot.")
    }

    assertthat::assert_that(param$frac > 0, msg="sda_boot: 'frac' must be > 0")
    assertthat::assert_that(param$reps > 0, msg="sda_boot: 'reps' must be > 0")

    while (count <= param$reps) {
      ysam <- lapply(split_y, function(idx) if (length(idx)==1) idx else sample(idx, length(idx), replace=TRUE))
      row.idx <- sort(unlist(ysam))

      ret <- if (param$frac > 0 && param$frac < 1) {
        nkeep <- max(param$frac * ncol(x),1)
        ind <- sample(seq_len(ncol(x)), nkeep)
        fit <- try(sda::sda(Xtrain=x[row.idx,ind,drop=FALSE], L=y[row.idx], lambda=lambda[count], verbose=FALSE, ...),
                   silent=TRUE)
        attr(fit, "keep.ind") <- ind
        fit
      } else {
        fit <- try(sda::sda(Xtrain=x[row.idx,], L=y[row.idx], lambda=lambda[count], verbose=FALSE, ...),
                   silent=TRUE)
        attr(fit, "keep.ind") <- seq_len(ncol(x))
        fit
      }

      if (!inherits(ret, "try-error")) {
        mfits[[count]] <- ret
        count <- count + 1
      } else {
        message("sda model fit error, retry: ", failures + 1)
        failures <- failures + 1
        if (failures > 10) return(ret)
      }
    }

    out <- list(fits=mfits, lev=lev)
    class(out) <- "sda_boot"
    out
  },

  predict=function(modelFit, newdata, preProc = NULL, submodels = NULL) {
    preds <- lapply(modelFit$fits, function(fit) {
      ind <- attr(fit, "keep.ind")
      predict(fit, as.matrix(newdata)[,ind,drop=FALSE], verbose=FALSE)$posterior
    })

    prob <- preds[!sapply(preds, is.null)]
    pfinal <- Reduce("+", prob)/length(prob)
    modelFit$lev[apply(pfinal, 1, which.max)]
  },

  prob=function(modelFit, newdata, preProc = NULL, submodels = NULL) {
    preds <- lapply(modelFit$fits, function(fit) {
      ind <- attr(fit, "keep.ind")
      predict(fit, as.matrix(newdata)[,ind,drop=FALSE], verbose=FALSE)$posterior
    })
    prob <- preds[!sapply(preds, is.null)]
    pfinal <- Reduce("+", prob)/length(prob)
    colnames(pfinal) <- modelFit$lev
    pfinal
  }
)


# lda_thomaz_boot
# Store lev similarly, ensure correct predictions and probs
#' @keywords internal
#' @noRd
MVPAModels$lda_thomaz_boot <- list(
  type = "Classification",
  library = "sparsediscrim",
  label="lda_thomaz_boot",
  loop = NULL,
  parameters=data.frame(parameters=c("reps","frac"), class=c("numeric","numeric"),
                        label=c("bootstrap reps","fraction of features")),
  grid=function(x, y, len = NULL) data.frame(reps=10, frac=1),

  fit=function(x, y, wts, param, lev, last, weights, classProbs, ...) {
    x <- as.matrix(x)
    mfits <- list()
    count <- 1
    failures <- 0

    split_y <- split(seq_along(y), y)
    if (!all(sapply(split_y, function(yy) length(yy) > 0))) {
      stop("All levels in 'y' must have ≥1 instance for lda_thomaz_boot.")
    }
    assertthat::assert_that(param$frac > 0, msg="lda_thomaz_boot: 'frac' must be > 0")

    while (count <= param$reps) {
      ysam <- lapply(split_y, function(idx) if (length(idx)==1) idx else sample(idx, length(idx), replace=TRUE))
      row.idx <- sort(unlist(ysam))

      ret <- if (param$frac > 0 && param$frac < 1) {
        nkeep <- max(param$frac * ncol(x),1)
        ind <- sample(seq_len(ncol(x)), nkeep)
        fit <- try(sparsediscrim::lda_thomaz(x=x[row.idx,ind,drop=FALSE], y[row.idx]), silent=TRUE)
        attr(fit, "keep.ind") <- ind
        fit
      } else {
        fit <- try(sparsediscrim::lda_thomaz(x[row.idx,,drop=FALSE], y[row.idx]), silent=TRUE)
        attr(fit, "keep.ind") <- seq_len(ncol(x))
        fit
      }

      if (!inherits(ret, "try-error")) {
        mfits[[count]] <- ret
        count <- count + 1
      } else {
        message("lda_thomaz model fit error, retry: ", failures + 1)
        failures <- failures + 1
        if (failures > 10) return(ret)
      }
    }

    out <- list(fits=mfits, lev=lev)
    class(out) <- "lda_thomaz_boot"
    out
  },

  predict=function(modelFit, newdata, preProc = NULL, submodels = NULL) {
    preds <- lapply(modelFit$fits, function(fit) {
      ind <- attr(fit, "keep.ind")
      scores <- -t(predict(fit, newdata[,ind])$scores)
      mc <- scores[cbind(seq_len(nrow(scores)), max.col(scores, ties.method = "first"))]
      probs <- exp(scores - mc)
      zapsmall(probs/rowSums(probs))
    })

    prob <- preds[!sapply(preds, is.null)]
    pfinal <- Reduce("+", prob)/length(prob)
    modelFit$lev[apply(pfinal, 1, which.max)]
  },

  prob=function(modelFit, newdata, preProc = NULL, submodels = NULL) {
    preds <- lapply(modelFit$fits, function(fit) {
      ind <- attr(fit, "keep.ind")
      scores <- -t(predict(fit, newdata[,ind])$scores)
      mc <- scores[cbind(seq_len(nrow(scores)), max.col(scores))]
      probs <- exp(scores - mc)
      zapsmall(probs/rowSums(probs))
    })

    prob <- preds[!sapply(preds, is.null)]
    pfinal <- Reduce("+", prob)/length(prob)
    # pfinal has no colnames here since lda_thomaz doesn't store them directly.
    # Normally, need levels. If needed: colnames(pfinal) <- modelFit$lev
    pfinal
  }
)

# memoized ranking for sda methods
#' @importFrom memoise memoise
#' @importFrom sda sda.ranking
#' @noRd
memo_rank <- memoise::memoise(function(X, L, fdr) {
  sda::sda.ranking(X,L,fdr=fdr,verbose=FALSE)
})

# glmnet_opt
# Store lev, fix prob calculation
#' @keywords internal
#' @noRd
MVPAModels$glmnet_opt <- list(
  type = "Classification",
  library = c("c060","glmnet"),
  label="glmnet_opt",
  loop = NULL,
  parameters=data.frame(parameters="parameter", class="character", label="parameter"),

  grid=function(x, y, len = NULL) data.frame(parameter="none"),

  fit=function(x, y, wts, param, lev, last, weights, classProbs, ...) {
    x <- as.matrix(x)
    bounds <- t(data.frame(alpha=c(0,1)))
    colnames(bounds)<-c("lower","upper")
    fam <- if (length(lev) > 2) "multinomial" else "binomial"

    # Generate fold assignments for cv.glmnet
    foldid <- create_mvpa_folds(y, k = 5, list = FALSE, seed = 1234)
    fit <- epsgo(Q.func="tune.glmnet.interval",
                 bounds=bounds,
                 parms.coding="none",
                 seed=1234,
                 show="none",
                 fminlower=-100,
                 x=x, y=y, family=fam,
                 type.min="lambda.1se",
                 foldid=foldid,
                 type.measure="mse")
    sfit <- summary(fit, verbose=FALSE)
    modelFit <- glmnet(x,y,family=fam,alpha=sfit$opt.alpha, nlambda=20)
    modelFit$opt_lambda <- sfit$opt.lambda
    modelFit$opt_alpha <- sfit$opt.alpha
    modelFit$obsLevels <- lev
    modelFit
  },

  predict=function(modelFit, newdata, preProc = NULL, submodels = NULL) {
    out <- predict(modelFit, as.matrix(newdata), s=modelFit$opt_lambda, type="class")
    if (is.matrix(out)) out[,1] else out
  },

  prob=function(modelFit, newdata, preProc = NULL, submodels = NULL) {
    lev <- modelFit$obsLevels
    fam <- if (length(lev) > 2) "multinomial" else "binomial"
    probs <- predict(modelFit, as.matrix(newdata), s=modelFit$opt_lambda, type="response")
    if (fam=="binomial") {
      # Binary
      probs <- as.vector(probs)
      probs <- data.frame(probs_1 = 1 - probs, probs_2 = probs)
      colnames(probs) <- lev
    } else {
      # Multinomial returns a 3D array: (N, classes, 1)
      # Convert to data.frame
      probs <- as.data.frame(probs[,,1,drop=FALSE])
      colnames(probs) <- lev
    }
    probs
  }
)

# sparse_sda
# Store lev if needed, no explicit predict labeling needed except from posterior
#' @keywords internal
#' @noRd
MVPAModels$sparse_sda <- list(
  type = "Classification",
  library = "sda",
  label="sparse_sda",
  loop = NULL,
  parameters=data.frame(parameters=c("frac","lambda"), class=c("numeric","numeric"),
                        label=c("frac features","lambda")),
  grid=function(x, y, len = NULL) expand.grid(frac=seq(.1,1,length.out=len), lambda=seq(.01,.99,length.out=len)),

  fit=function(x, y, wts, param, lev, last, weights, classProbs, ...) {
    x <- as.matrix(x)
    nkeep <- max(param$frac * ncol(x),1)
    rank <- memo_rank(x, L=y, fdr=FALSE)
    ind <- rank[,"idx"][1:nkeep]

    fit <- sda::sda(Xtrain=x[,ind,drop=FALSE], L=y, lambda=param$lambda, verbose=FALSE)
    attr(fit, "keep.ind") <- ind
    fit$obsLevels <- lev
    fit
  },

  predict=function(modelFit, newdata, preProc = NULL, submodels = NULL) {
    predict(modelFit, as.matrix(newdata[,attr(modelFit, "keep.ind"), drop=FALSE]), verbose=FALSE)$class
  },

  prob=function(modelFit, newdata, preProc = NULL, submodels = NULL) {
    predict(modelFit, as.matrix(newdata[,attr(modelFit, "keep.ind"), drop=FALSE]),verbose=FALSE)$posterior
  }
)


# sda_ranking
# Store lev
#' @keywords internal
#' @noRd
MVPAModels$sda_ranking <- list(
  type = "Classification",
  library = "sda",
  label="sda_ranking",
  loop = NULL,
  parameters=data.frame(parameters="parameter", class="character", label="parameter"),
  grid=function(x, y, len = NULL) data.frame(parameter="none"),

  fit=function(x, y, wts, param, lev, last, weights, classProbs, ...) {
    x <- as.matrix(x)
    if (ncol(x) > 2) {
      rank <- sda::sda.ranking(Xtrain=x, L=y, fdr=TRUE, verbose=FALSE, ...)
      hcind <- which.max(rank[,"HC"])
      keep.ind <- if (hcind < 2) seq(1, min(ncol(x), 2)) else 1:hcind
      ind <- rank[keep.ind,"idx"]
    } else {
      ind <- seq_len(ncol(x))
    }

    fit <- sda::sda(Xtrain=x[,ind,drop=FALSE], L=y, verbose=FALSE)
    attr(fit, "keep.ind") <- ind
    fit$obsLevels <- lev
    fit
  },

  predict=function(modelFit, newdata, preProc = NULL, submodels = NULL) {
    predict(modelFit, as.matrix(newdata[,attr(modelFit, "keep.ind"), drop=FALSE]), verbose=FALSE)$class
  },

  prob=function(modelFit, newdata, preProc = NULL, submodels = NULL) {
    predict(modelFit, as.matrix(newdata[,attr(modelFit, "keep.ind"), drop=FALSE]),verbose=FALSE)$posterior
  }
)


# mgsda
# Fix MSGDA::classifiyV to MGSDA::classifyV and store lev
# Define modelFit as a list after training
#' @keywords internal
#' @noRd
MVPAModels$mgsda <- list(
  type = "Classification",
  library = "MGSDA",
  loop = NULL,
  parameters=data.frame(parameters="lambda", class="numeric", label="sparsity penalty"),

  grid=function(x, y, len = NULL) data.frame(lambda=seq(.001, .99, length.out=len)),

  fit=function(x, y, wts, param, lev, last, weights, classProbs, ...) {
    ycodes <- as.integer(y)
    V <- MGSDA::dLDA(as.matrix(x), ycodes, lambda=param$lambda, ...)
    modelFit <- list(
      V = V,
      ycodes = ycodes,
      xtrain = as.matrix(x),
      ytrain = y,
      obsLevels = lev
    )
    modelFit
  },

  predict=function(modelFit, newdata, preProc = NULL, submodels = NULL) {
    preds <- MGSDA::classifyV(modelFit$xtrain, modelFit$ycodes, as.matrix(newdata), modelFit$V)
    modelFit$obsLevels[preds]
  },

  prob=function(modelFit, newdata, preProc = NULL, submodels = NULL) {
    # Not implemented in original code. If probabilities are not supported, return NULL or implement if possible.
    # Here we return NULL to avoid errors.
    NULL
  }
)

# lda_thomaz
# Store lev
#' @keywords internal
#' @noRd
MVPAModels$lda_thomaz <- list(
  type = "Classification",
  library = "sparsediscrim",
  loop = NULL,
  parameters=data.frame(parameters="parameter", class="character", label="parameter"),

  grid=function(x, y, len = NULL) data.frame(parameter="none"),

  fit=function(x, y, wts, param, lev, last, weights, classProbs, ...) {
    fit <- sparsediscrim::lda_thomaz(as.matrix(x), y, ...)
    fit$obsLevels <- lev
    fit
  },

  predict=function(modelFit, newdata, preProc = NULL, submodels = NULL) {
    # Returns a list with element $class
    predict(modelFit, as.matrix(newdata))$class
  },

  prob=function(modelFit, newdata, preProc = NULL, submodels = NULL) {
    p <- predict(modelFit, as.matrix(newdata), type="prob")
    # p is posterior probabilities with columns corresponding to classes in modelFit$obsLevels
    if (!is.null(modelFit$obsLevels)) colnames(p) <- modelFit$obsLevels
    p
  }
)

# hdrda
# Store lev
#' @keywords internal
#' @noRd
MVPAModels$hdrda <- list(
  type = "Classification",
  library = "sparsediscrim",
  label="hdrda",
  loop = NULL,
  parameters=data.frame(parameters=c("lambda","gamma"), class=c("numeric","numeric"),
                        label=c("HRDA pooling","shrinkage parameter")),

  grid=function(x, y, len = NULL) expand.grid(lambda=seq(.99, .001, length.out=len), gamma=seq(.001, .99, length.out=len)),

  fit=function(x, y, wts, param, lev, last, weights, classProbs, ...) {
    fit <- sparsediscrim::hdrda(as.matrix(x), y, lambda=param$lambda, gamma=param$gamma, ...)
    fit$obsLevels <- lev
    fit
  },

  predict=function(modelFit, newdata, preProc = NULL, submodels = NULL) {
    predict(modelFit, as.matrix(newdata))$class
  },

  prob=function(modelFit, newdata, preProc = NULL, submodels = NULL) {
    posterior <- predict(modelFit, as.matrix(newdata))$posterior
    posterior <- t(apply(posterior,1,function(x) x/sum(x)))
    if (!is.null(modelFit$obsLevels)) colnames(posterior) <- modelFit$obsLevels
    posterior
  }
)


#' @importFrom stats var dnorm
#' @keywords internal
#' @noRd
MVPAModels$naive_bayes <- list(
  type = "Classification",
  library = NULL,
  label = "naive_bayes",
  loop = NULL,
  parameters = data.frame(
    parameter = "parameter",
    class = "character",
    label = "parameter"
  ),
  grid = function(x, y, len = NULL) {
    data.frame(parameter = "none")
  },

  fit = function(x, y, wts, param, lev, last, weights, classProbs, ...) {
    x <- as.matrix(x)
    y <- factor(y)
    classes <- levels(y)
    n_class <- length(classes)

    log_priors <- log(table(y) / length(y))
    obs_by_class <- split(seq_len(nrow(x)), y)
    mus <- vars <- matrix(
      NA,
      nrow = n_class,
      ncol = ncol(x),
      dimnames = list(classes, NULL)
    )

    for (k in seq_along(classes)) {
      idx <- obs_by_class[[classes[k]]]
      samples <- x[idx, , drop = FALSE]
      mus[k, ] <- colMeans(samples)

      nk <- length(idx)
      vars_k <- apply(samples, 2, var) * (nk - 1) / nk
      zero_var <- vars_k <= .Machine$double.eps
      if (any(zero_var)) {
        nz <- vars_k[vars_k > .Machine$double.eps]
        eps <- if (length(nz) > 0) min(nz) * 1e-6 else 1e-10
        vars_k[zero_var] <- max(eps, 1e-10)
        warning(
          sprintf(
            "Naive Bayes: %d features had zero variance for class '%s'. Replaced with %g.",
            sum(zero_var), classes[k], max(eps, 1e-10)
          )
        )
      }
      vars[k, ] <- vars_k
    }

    modelFit <- list(
      mus = mus,
      vars = vars,
      log_class_priors = as.numeric(log_priors),
      classes = classes,
      obsLevels = classes
    )
    attr(modelFit, "obsLevels") <- classes
    attr(modelFit, "problemType") <- "Classification"
    modelFit
  },

  predict = function(modelFit, newdata, preProc = NULL, submodels = NULL) {
    newdata <- as.matrix(newdata)
    log_post <- calculate_log_posteriors(modelFit, newdata)
    factor(modelFit$classes[max.col(log_post)], levels = modelFit$classes)
  },

  prob = function(modelFit, newdata, preProc = NULL, submodels = NULL) {
    newdata <- as.matrix(newdata)
    log_post <- calculate_log_posteriors(modelFit, newdata)
    probs <- t(apply(log_post, 1, function(row) {
      max_log <- max(row)
      exp(row - max_log) / sum(exp(row - max_log))
    }))
    colnames(probs) <- modelFit$classes
    probs
  }
)

# Helper function to calculate log posterior probabilities (used by predict and prob)
# This avoids code duplication
calculate_log_posteriors <- function(modelFit, newdata) {
  mus <- modelFit$mus
  vars <- modelFit$vars
  log_priors <- modelFit$log_class_priors
  classes <- modelFit$classes
  n_class <- length(classes)

  if (ncol(newdata) != ncol(mus)) {
    stop(
      sprintf(
        "Feature dimension mismatch: model trained with %d features, newdata has %d",
        ncol(mus), ncol(newdata)
      )
    )
  }

  log_posts <- matrix(NA, nrow = nrow(newdata), ncol = n_class,
                       dimnames = list(NULL, classes))

  for (k in seq_along(classes)) {
    mu <- mus[k, ]
    var <- vars[k, ]
    ll <- sapply(seq_along(mu), function(j) {
      dnorm(newdata[, j], mean = mu[j], sd = sqrt(var[j]), log = TRUE)
    })
    ll[!is.finite(ll)] <- -1e100
    log_posts[, k] <- rowSums(ll) + log_priors[k]
  }

  all_inf <- apply(log_posts, 1, function(row) all(!is.finite(row)))
  if (any(all_inf)) {
    log_posts[all_inf, ] <- log(1 / n_class)
    warning(sprintf(
      "%d samples had non-finite log posteriors for all classes. Assigned equal probability.",
      sum(all_inf)
    ))
  }

  log_posts
}

