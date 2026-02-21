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
#'   \item{spacenet_tvl1}{Spatially-regularized sparse linear model with TV-L1 penalty for global whole-brain analysis}
#' }
#'
#' @details
#' Models are accessed via \code{load_model(name)}. Each model specification includes \code{fit}, \code{predict}, and \code{prob} methods.
#'
#' The \code{spacenet_tvl1} model follows the SpaceNet formulation used in Nilearn.
#'
#' @references
#' Gramfort, A., Thirion, B., & Varoquaux, G. (2013).
#' \emph{Identifying predictive regions from fMRI with TV-L1 prior}.
#' Pattern Recognition in Neuroimaging (PRNI), IEEE.
#' https://inria.hal.science/hal-00839984
#'
#' @return An environment containing registered MVPA model specifications.
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

#' @keywords internal
#' @noRd
quiet_sda <- function(...) {
  res <- NULL
  invisible(capture.output(res <- sda::sda(...)))
  res
}

#' @keywords internal
#' @noRd
quiet_sda_ranking <- function(...) {
  res <- NULL
  invisible(capture.output(res <- sda::sda.ranking(...)))
  res
}

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
  X <- as.matrix(newData)
  M <- as.matrix(modelFit$conditionMeans)
  if (identical(modelFit$method, "pearson")) {
    # Fast row-wise Pearson correlation via standardization + tcrossprod
    p <- ncol(X)
    if (is.null(p) || p < 2L) {
      cres <- cor(t(X), t(M), method = "pearson")
    } else {
      Xc <- sweep(X, 1, rowMeans(X), "-")
      Mc <- sweep(M, 1, rowMeans(M), "-")
      Xs <- sweep(Xc, 1, pmax(matrixStats::rowSds(X), .Machine$double.eps), "/")
      Ms <- sweep(Mc, 1, pmax(matrixStats::rowSds(M), .Machine$double.eps), "/")
      cres <- tcrossprod(Xs, Ms) / (p - 1)
    }
  } else {
    cres <- cor(t(X), t(M), method = modelFit$method)
  }
  res <- max.col(cres)
  factor(modelFit$levs[res], levels = modelFit$levs)
}

#' @keywords internal
#' @noRd
prob_corsimFit <- function(modelFit, newData) {
  X <- as.matrix(newData)
  M <- as.matrix(modelFit$conditionMeans)
  if (identical(modelFit$method, "pearson")) {
    p <- ncol(X)
    if (is.null(p) || p < 2L) {
      scores <- cor(t(X), t(M), method = "pearson")
    } else {
      Xc <- sweep(X, 1, rowMeans(X), "-")
      Mc <- sweep(M, 1, rowMeans(M), "-")
      Xs <- sweep(Xc, 1, pmax(matrixStats::rowSds(X), .Machine$double.eps), "/")
      Ms <- sweep(Mc, 1, pmax(matrixStats::rowSds(M), .Machine$double.eps), "/")
      scores <- tcrossprod(Xs, Ms) / (p - 1)
    }
  } else {
    scores <- cor(t(X), t(M), method = modelFit$method)
  }
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

#' @keywords internal
#' @noRd
.dual_lda_row_softmax <- function(scores) {
  row_max <- apply(scores, 1, max)
  exp_scores <- exp(sweep(scores, 1, row_max, "-"))
  exp_scores / rowSums(exp_scores)
}

#' @keywords internal
#' @noRd
.dual_lda_fit_core <- function(x, y, gamma) {
  x <- as.matrix(x)
  y <- as.factor(y)

  if (length(levels(y)) < 2L) {
    stop("dual_lda requires at least two classes.")
  }
  if (!is.finite(gamma) || gamma <= 0) {
    stop("dual_lda: `gamma` must be a positive finite scalar.")
  }

  classes <- levels(y)
  n <- nrow(x)
  p <- ncol(x)
  k <- length(classes)

  class_counts <- as.numeric(table(y)[classes])
  priors <- class_counts / sum(class_counts)

  means <- matrix(0, nrow = k, ncol = p, dimnames = list(classes, NULL))
  for (i in seq_len(k)) {
    idx <- which(y == classes[i])
    means[i, ] <- colMeans(x[idx, , drop = FALSE])
  }

  class_id <- match(y, classes)
  centered <- x - means[class_id, , drop = FALSE]
  sigma <- crossprod(centered) + diag(gamma, p)
  chol_sigma <- chol(sigma)

  means_t <- t(means)
  inv_sigma_means <- backsolve(chol_sigma, forwardsolve(t(chol_sigma), means_t))
  lin_const <- -0.5 * colSums(means_t * inv_sigma_means) + log(pmax(priors, .Machine$double.eps))

  list(
    classes = classes,
    means = means,
    inv_sigma_means = inv_sigma_means,
    lin_const = lin_const,
    priors = priors,
    gamma = gamma,
    obsLevels = classes
  )
}

#' @keywords internal
#' @noRd
.dual_lda_prob_core <- function(modelFit, newdata) {
  x <- as.matrix(newdata)
  scores <- x %*% modelFit$inv_sigma_means
  scores <- sweep(scores, 2, modelFit$lin_const, "+")
  probs <- .dual_lda_row_softmax(scores)
  colnames(probs) <- modelFit$classes
  probs
}

#' @keywords internal
#' @noRd
.dual_lda_predict_core <- function(modelFit, newdata) {
  probs <- .dual_lda_prob_core(modelFit, newdata)
  factor(modelFit$classes[max.col(probs)], levels = modelFit$classes)
}

#' @keywords internal
#' @noRd
MVPAModels$dual_lda <- list(
  type = "Classification",
  library = "rMVPA",
  loop = NULL,
  label = "dual_lda",
  parameters = data.frame(parameters = "gamma", class = "numeric", labels = "gamma"),

  grid = function(x, y, len = NULL) {
    if (is.null(len) || len <= 1L) {
      data.frame(gamma = 1e-2)
    } else {
      data.frame(gamma = 10 ^ seq(-4, 1, length.out = len))
    }
  },

  fit = function(x, y, wts, param, lev, last, weights, classProbs, ...) {
    fit <- .dual_lda_fit_core(x, y, gamma = as.numeric(param$gamma))
    fit$obsLevels <- if (!is.null(lev)) lev else fit$classes
    fit
  },

  predict = function(modelFit, newdata, preProc = NULL, submodels = NULL) {
    .dual_lda_predict_core(modelFit, newdata)
  },

  prob = function(modelFit, newdata, preProc = NULL, submodels = NULL) {
    .dual_lda_prob_core(modelFit, newdata)
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
    m <- quiet_sda(Xtrain=as.matrix(x), L=y, verbose=FALSE, ...)
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
        fit <- try(quiet_sda(Xtrain=x[row.idx,ind,drop=FALSE], L=y[row.idx], lambda=lambda[count], verbose=FALSE, ...),
                   silent=TRUE)
        attr(fit, "keep.ind") <- ind
        fit
      } else {
        fit <- try(quiet_sda(Xtrain=x[row.idx,], L=y[row.idx], lambda=lambda[count], verbose=FALSE, ...),
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
      stop("All levels in 'y' must have >= 1 instance for lda_thomaz_boot.")
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
  quiet_sda_ranking(X,L,fdr=fdr,verbose=FALSE)
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

    fit <- quiet_sda(Xtrain=x[,ind,drop=FALSE], L=y, lambda=param$lambda, verbose=FALSE)
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
      rank <- quiet_sda_ranking(Xtrain=x, L=y, fdr=TRUE, verbose=FALSE, ...)
      hcind <- which.max(rank[,"HC"])
      keep.ind <- if (hcind < 2) seq(1, min(ncol(x), 2)) else 1:hcind
      ind <- rank[keep.ind,"idx"]
    } else {
      ind <- seq_len(ncol(x))
    }

    fit <- quiet_sda(Xtrain=x[,ind,drop=FALSE], L=y, verbose=FALSE)
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


# -----------------------------------------------------------------------------
# SpaceNet (TV-L1 via single-loop PDHG)
# -----------------------------------------------------------------------------

#' @keywords internal
#' @noRd
.spacenet_soft_threshold <- function(x, thresh) {
  sign(x) * pmax(abs(x) - thresh, 0)
}

#' @keywords internal
#' @noRd
.spacenet_sigmoid <- function(z) {
  out <- numeric(length(z))
  pos <- z >= 0
  out[pos] <- 1 / (1 + exp(-z[pos]))
  ez <- exp(z[!pos])
  out[!pos] <- ez / (1 + ez)
  out
}

#' @keywords internal
#' @noRd
.spacenet_softmax <- function(eta) {
  eta <- as.matrix(eta)
  row_max <- apply(eta, 1L, max)
  ex <- exp(sweep(eta, 1L, row_max, "-"))
  rs <- rowSums(ex)
  rs[!is.finite(rs) | rs <= 0] <- 1
  ex / rs
}

#' @keywords internal
#' @noRd
.spacenet_y01 <- function(y, positive = NULL) {
  if (is.factor(y)) {
    lev <- levels(y)
    if (length(lev) < 2L) {
      return(rep(0, length(y)))
    }
    pos <- if (is.null(positive)) lev[2L] else as.character(positive)
    return(as.numeric(as.character(y) == pos))
  }

  y_num <- as.numeric(y)
  vals <- sort(unique(y_num[is.finite(y_num)]))
  if (length(vals) == 0L) return(rep(0, length(y_num)))
  if (length(vals) == 1L) return(as.numeric(y_num == vals[1L]))
  if (length(vals) == 2L) return(as.numeric(y_num == vals[2L]))
  as.numeric(y_num > mean(vals))
}

#' @keywords internal
#' @noRd
.spacenet_spectral_norm_squared <- function(X, n_iter = 20L) {
  X <- as.matrix(X)
  n <- ncol(X)
  if (n == 0L) return(0)
  v <- rnorm(n)
  v <- v / max(sqrt(sum(v * v)), .Machine$double.eps)
  for (i in seq_len(n_iter)) {
    v <- as.vector(crossprod(X, X %*% v))
    nv <- sqrt(sum(v * v))
    if (!is.finite(nv) || nv <= .Machine$double.eps) break
    v <- v / nv
  }
  w <- as.vector(X %*% v)
  sum(w * w)
}

#' @keywords internal
#' @noRd
.spacenet_screen_support <- function(X, y, screening_percentile = 100,
                                     min_features = 20L,
                                     feature_ids = NULL) {
  p <- ncol(X)
  if (p == 0L || is.null(p)) return(integer(0))
  if (screening_percentile >= 100 || p <= 100) return(seq_len(p))

  keep_k <- max(min_features, ceiling((screening_percentile / 100) * p))
  keep_k <- min(keep_k, p)

  y_num <- if (is.factor(y)) .spacenet_y01(y) else as.numeric(y)

  yc <- y_num - mean(y_num)
  xc <- sweep(X, 2L, colMeans(X), "-")
  scores <- abs(as.vector(crossprod(xc, yc)))
  scores[!is.finite(scores)] <- 0

  if (!is.null(feature_ids) && anyDuplicated(feature_ids) > 0L) {
    feature_ids <- as.integer(feature_ids)
    fid <- as.character(feature_ids)
    agg <- tapply(scores, fid, function(v) {
      vv <- v[is.finite(v)]
      if (length(vv) == 0L) 0 else max(vv)
    })
    agg[!is.finite(agg)] <- 0

    group_sizes <- table(fid)
    basis_count <- max(as.integer(group_sizes), 1L)
    min_groups <- max(1L, ceiling(min_features / basis_count))

    keep_groups <- max(min_groups, ceiling((screening_percentile / 100) * length(agg)))
    keep_groups <- min(keep_groups, length(agg))

    keep_ids <- names(sort(agg, decreasing = TRUE))[seq_len(keep_groups)]
    return(which(fid %in% keep_ids))
  }

  ord <- order(scores, decreasing = TRUE)
  ord[seq_len(keep_k)]
}

#' @keywords internal
#' @noRd
.spacenet_make_alpha_grid <- function(X, y, l1_ratio,
                                      n_alphas = 12L,
                                      alpha_min_ratio = 1e-3,
                                      is_classif = FALSE) {
  X <- as.matrix(X)
  n <- nrow(X)
  if (n <= 1L || ncol(X) == 0L) {
    return(rep(1e-2, n_alphas))
  }

  Xc <- sweep(X, 2L, colMeans(X), "-")
  if (is_classif) {
    y01 <- .spacenet_y01(y)
    g0 <- as.vector(crossprod(Xc, 0.5 - y01)) / n
  } else {
    yc <- as.numeric(y) - mean(as.numeric(y))
    g0 <- -as.vector(crossprod(Xc, yc)) / n
  }

  denom <- max(as.numeric(l1_ratio), 1e-3)
  alpha_max <- max(abs(g0)) / denom
  if (!is.finite(alpha_max) || alpha_max <= 0) {
    alpha_max <- 1e-2
  }
  alpha_min <- max(alpha_max * alpha_min_ratio, 1e-8)
  exp(seq(log(alpha_max), log(alpha_min), length.out = n_alphas))
}

#' @keywords internal
#' @noRd
.spacenet_edges_from_feature_ids <- function(feature_ids, spatial_mask = NULL) {
  feature_ids <- as.integer(feature_ids)
  p <- length(feature_ids)
  empty <- list(edges = matrix(integer(0), nrow = 0L, ncol = 2L), d_norm2 = 0)

  if (p <= 1L || is.null(spatial_mask) || !inherits(spatial_mask, "NeuroVol")) {
    return(empty)
  }
  dims <- dim(spatial_mask)[1:3]
  if (length(dims) < 3L || any(!is.finite(dims))) {
    return(empty)
  }
  full_len <- prod(dims)
  valid <- feature_ids >= 1L & feature_ids <= full_len
  if (!all(valid)) {
    feature_ids <- feature_ids[valid]
    p <- length(feature_ids)
    if (p <= 1L) return(empty)
  }

  uniq_ids <- sort(unique(feature_ids))
  n_unique <- length(uniq_ids)
  if (n_unique <= 1L) return(empty)

  idx_map <- integer(full_len)
  idx_map[uniq_ids] <- seq_len(n_unique)

  coords <- arrayInd(uniq_ids, .dim = dims)
  offsets <- rbind(
    c(1L, 0L, 0L),
    c(0L, 1L, 0L),
    c(0L, 0L, 1L)
  )

  all_edges <- vector("list", nrow(offsets))
  edge_count <- 0L

  for (k in seq_len(nrow(offsets))) {
    nbr <- sweep(coords, 2L, offsets[k, ], "+")
    in_bounds <-
      nbr[, 1L] >= 1L & nbr[, 1L] <= dims[1L] &
      nbr[, 2L] >= 1L & nbr[, 2L] <= dims[2L] &
      nbr[, 3L] >= 1L & nbr[, 3L] <= dims[3L]

    if (!any(in_bounds)) next

    src_u <- which(in_bounds)
    nbr_lin <- (nbr[in_bounds, 1L] - 1L) +
      dims[1L] * (nbr[in_bounds, 2L] - 1L) +
      dims[1L] * dims[2L] * (nbr[in_bounds, 3L] - 1L) + 1L

    dst_u <- idx_map[nbr_lin]
    ok <- dst_u > 0L
    if (!any(ok)) next

    e <- cbind(src_u[ok], dst_u[ok])
    keep <- e[, 1L] < e[, 2L]
    if (!any(keep)) next

    edge_count <- edge_count + 1L
    all_edges[[edge_count]] <- e[keep, , drop = FALSE]
  }

  if (edge_count == 0L) {
    return(empty)
  }

  base_edges <- do.call(rbind, all_edges[seq_len(edge_count)])
  if (is.null(base_edges) || nrow(base_edges) == 0L) {
    return(empty)
  }

  id_groups <- split(seq_len(p), factor(feature_ids, levels = uniq_ids))
  max_copies <- max(lengths(id_groups), 1L)

  if (max_copies <= 1L) {
    col_map <- vapply(id_groups, function(ix) ix[1L], integer(1))
    edges <- cbind(col_map[base_edges[, 1L]], col_map[base_edges[, 2L]])
  } else {
    idx_mat <- matrix(NA_integer_, nrow = n_unique, ncol = max_copies)
    for (i in seq_len(n_unique)) {
      ix <- sort(id_groups[[i]])
      idx_mat[i, seq_along(ix)] <- ix
    }

    edge_list <- vector("list", max_copies)
    edge_count <- 0L
    for (b in seq_len(max_copies)) {
      src <- idx_mat[base_edges[, 1L], b]
      dst <- idx_mat[base_edges[, 2L], b]
      ok <- !is.na(src) & !is.na(dst)
      if (!any(ok)) next
      e <- cbind(src[ok], dst[ok])
      keep <- e[, 1L] < e[, 2L]
      if (!any(keep)) next
      edge_count <- edge_count + 1L
      edge_list[[edge_count]] <- e[keep, , drop = FALSE]
    }

    if (edge_count == 0L) return(empty)
    edges <- do.call(rbind, edge_list[seq_len(edge_count)])
  }

  if (is.null(edges) || nrow(edges) == 0L) return(empty)
  edges <- unique(edges)

  deg <- tabulate(c(edges[, 1L], edges[, 2L]), nbins = p)
  d_norm2 <- 2 * max(deg, 1L)

  list(edges = edges, d_norm2 = d_norm2)
}

#' @keywords internal
#' @noRd
.spacenet_pdhg_solver <- function(X, y, alpha, l1_ratio,
                                  edges,
                                  d_norm2,
                                  loss = c("mse", "logistic"),
                                  max_iter = 300L,
                                  tol = 1e-4,
                                  data_lipschitz = NULL,
                                  init = NULL,
                                  sigma = 1.0,
                                  theta = 1.0) {
  loss <- match.arg(loss)
  X <- as.matrix(X)
  n <- nrow(X)
  p <- ncol(X)

  lambda1 <- alpha * l1_ratio
  lambda_tv <- alpha * (1 - l1_ratio)

  if (is.null(data_lipschitz)) {
    spec2 <- .spacenet_spectral_norm_squared(X)
    data_lipschitz <- if (loss == "logistic") {
      0.25 * (spec2 / max(n, 1L)) + 0.25
    } else {
      spec2 / max(n, 1L)
    }
  }
  data_lipschitz <- as.numeric(data_lipschitz)[1L]
  if (!is.finite(data_lipschitz) || data_lipschitz < 0) {
    data_lipschitz <- 0
  }

  tau <- 0.99 / (data_lipschitz + sigma * d_norm2 + 1e-12)

  w <- if (!is.null(init$w) && length(init$w) == p) {
    as.numeric(init$w)
  } else {
    numeric(p)
  }
  wbar <- if (!is.null(init$wbar) && length(init$wbar) == p) {
    as.numeric(init$wbar)
  } else {
    w
  }
  b <- if (loss == "logistic") {
    if (!is.null(init$b)) as.numeric(init$b) else 0
  } else {
    0
  }
  bbar <- if (loss == "logistic") {
    if (!is.null(init$bbar)) as.numeric(init$bbar) else b
  } else {
    0
  }
  p_dual <- if (!is.null(edges) && nrow(edges) > 0L) {
    if (!is.null(init$p_dual) && length(init$p_dual) == nrow(edges)) {
      as.numeric(init$p_dual)
    } else {
      numeric(nrow(edges))
    }
  } else {
    numeric(0)
  }

  rel_hist <- numeric(max_iter)

  for (it in seq_len(max_iter)) {
    w_old <- w
    b_old <- b

    btp <- numeric(p)
    if (length(p_dual) > 0L && lambda_tv > 0) {
      diff_bar <- wbar[edges[, 2L]] - wbar[edges[, 1L]]
      p_dual <- p_dual + sigma * diff_bar
      p_dual <- pmax(-lambda_tv, pmin(lambda_tv, p_dual))

      btp[edges[, 1L]] <- btp[edges[, 1L]] - p_dual
      btp[edges[, 2L]] <- btp[edges[, 2L]] + p_dual
    } else if (length(p_dual) > 0L) {
      p_dual[] <- 0
    }

    if (loss == "logistic") {
      eta <- as.vector(X %*% w + b)
      prob <- .spacenet_sigmoid(eta)
      diff <- prob - y
      grad_w <- as.vector(crossprod(X, diff)) / max(n, 1L)
      grad_b <- mean(diff)
      w <- w - tau * (grad_w + btp)
      if (lambda1 > 0) {
        w <- .spacenet_soft_threshold(w, tau * lambda1)
      }
      b <- b - tau * grad_b
    } else {
      resid <- as.vector(X %*% w - y)
      grad_w <- as.vector(crossprod(X, resid)) / max(n, 1L)
      w <- w - tau * (grad_w + btp)
      if (lambda1 > 0) {
        w <- .spacenet_soft_threshold(w, tau * lambda1)
      }
      b <- 0
    }

    wbar <- w + theta * (w - w_old)
    if (loss == "logistic") {
      bbar <- b + theta * (b - b_old)
    }

    delta <- c(w - w_old, b - b_old)
    denom <- max(1, sqrt(sum(c(w_old, b_old)^2)))
    rel <- sqrt(sum(delta * delta)) / denom
    rel_hist[it] <- rel
    if (is.finite(rel) && rel < tol) {
      rel_hist <- rel_hist[seq_len(it)]
      break
    }
  }

  list(
    w = w,
    b = b,
    init = list(
      w = w,
      wbar = wbar,
      b = b,
      bbar = bbar,
      p_dual = p_dual,
      tau = tau,
      sigma = sigma,
      theta = theta,
      lipschitz = data_lipschitz
    ),
    rel_history = rel_hist
  )
}

#' @keywords internal
#' @noRd
.spacenet_cv_path <- function(X,
                              y,
                              feature_ids,
                              spatial_mask,
                              l1_ratio,
                              n_alphas,
                              alpha_min_ratio,
                              screening_percentile,
                              max_iter,
                              tol,
                              inner_folds,
                              is_classif) {
  n <- nrow(X)
  p <- ncol(X)

  alpha_grid <- .spacenet_make_alpha_grid(
    X, y,
    l1_ratio = l1_ratio,
    n_alphas = n_alphas,
    alpha_min_ratio = alpha_min_ratio,
    is_classif = is_classif
  )

  if (inner_folds < 2L || n < 6L) {
    return(list(alpha_grid = alpha_grid, cv_loss = rep(NA_real_, length(alpha_grid)), best_alpha = alpha_grid[1L]))
  }

  k <- min(inner_folds, n)
  fold_y <- if (is_classif) {
    factor(ifelse(.spacenet_y01(y) > 0, "pos", "neg"))
  } else {
    y
  }
  fold_tests <- create_mvpa_folds(fold_y, k = k, list = TRUE, seed = 1234)

  fold_loss <- matrix(NA_real_, nrow = length(alpha_grid), ncol = length(fold_tests))

  for (fi in seq_along(fold_tests)) {
    test_idx <- as.integer(fold_tests[[fi]])
    train_idx <- setdiff(seq_len(n), test_idx)
    if (length(train_idx) < 2L || length(test_idx) < 1L) next

    X_train <- X[train_idx, , drop = FALSE]
    X_test <- X[test_idx, , drop = FALSE]
    y_train <- y[train_idx]
    y_test <- y[test_idx]

    support <- .spacenet_screen_support(
      X_train, y_train,
      screening_percentile = screening_percentile,
      feature_ids = feature_ids
    )
    if (length(support) == 0L) support <- seq_len(p)

    X_train <- X_train[, support, drop = FALSE]
    X_test <- X_test[, support, drop = FALSE]
    fids_fold <- feature_ids[support]

    x_mean <- colMeans(X_train)
    X_train <- sweep(X_train, 2L, x_mean, "-")
    X_test <- sweep(X_test, 2L, x_mean, "-")

    if (is_classif) {
      y_train_proc <- .spacenet_y01(y_train)
      y_test_proc <- .spacenet_y01(y_test)
      loss_name <- "logistic"
    } else {
      y_mean <- mean(as.numeric(y_train))
      y_train_proc <- as.numeric(y_train) - y_mean
      y_test_proc <- as.numeric(y_test)
      loss_name <- "mse"
    }

    edge_info <- .spacenet_edges_from_feature_ids(fids_fold, spatial_mask)
    spec2 <- .spacenet_spectral_norm_squared(X_train)
    data_L <- if (is_classif) {
      0.25 * (spec2 / max(nrow(X_train), 1L)) + 0.25
    } else {
      spec2 / max(nrow(X_train), 1L)
    }

    init <- NULL
    for (ai in seq_along(alpha_grid)) {
      alpha <- alpha_grid[ai]
      fit <- .spacenet_pdhg_solver(
        X = X_train,
        y = y_train_proc,
        alpha = alpha,
        l1_ratio = l1_ratio,
        edges = edge_info$edges,
        d_norm2 = edge_info$d_norm2,
        loss = loss_name,
        max_iter = max(40L, as.integer(max_iter)),
        tol = 2 * tol,
        data_lipschitz = data_L,
        init = init
      )
      init <- fit$init

      if (is_classif) {
        eta <- as.vector(X_test %*% fit$w + fit$b)
        pr <- .spacenet_sigmoid(eta)
        pr <- pmin(pmax(pr, 1e-8), 1 - 1e-8)
        fold_loss[ai, fi] <- -mean(y_test_proc * log(pr) + (1 - y_test_proc) * log(1 - pr))
      } else {
        pred <- as.vector(X_test %*% fit$w + y_mean)
        fold_loss[ai, fi] <- mean((pred - y_test_proc)^2)
      }
    }
  }

  mean_loss <- rowMeans(fold_loss, na.rm = TRUE)
  if (all(!is.finite(mean_loss))) {
    best_idx <- 1L
  } else {
    mean_loss[!is.finite(mean_loss)] <- Inf
    best_idx <- which.min(mean_loss)
  }

  list(alpha_grid = alpha_grid, cv_loss = mean_loss, best_alpha = alpha_grid[best_idx])
}

#' @keywords internal
#' @noRd
.spacenet_fit_binary <- function(X, y01, feature_ids, spatial_mask,
                                 l1_ratio, n_alphas, alpha_min_ratio,
                                 screening_percentile, max_iter, tol,
                                 inner_folds) {
  X <- as.matrix(X)
  p <- ncol(X)
  y01 <- as.numeric(y01)

  cv_path <- .spacenet_cv_path(
    X = X,
    y = y01,
    feature_ids = feature_ids,
    spatial_mask = spatial_mask,
    l1_ratio = l1_ratio,
    n_alphas = n_alphas,
    alpha_min_ratio = alpha_min_ratio,
    screening_percentile = screening_percentile,
    max_iter = max_iter,
    tol = tol,
    inner_folds = inner_folds,
    is_classif = TRUE
  )

  alpha_opt <- as.numeric(cv_path$best_alpha)
  support <- .spacenet_screen_support(
    X, y01,
    screening_percentile = screening_percentile,
    feature_ids = feature_ids
  )
  if (length(support) == 0L) support <- seq_len(p)

  X_fit <- X[, support, drop = FALSE]
  fids_fit <- feature_ids[support]

  x_mean <- colMeans(X_fit)
  X_centered <- sweep(X_fit, 2L, x_mean, "-")

  edge_info <- .spacenet_edges_from_feature_ids(fids_fit, spatial_mask)
  spec2 <- .spacenet_spectral_norm_squared(X_centered)
  data_L <- 0.25 * (spec2 / max(nrow(X_centered), 1L)) + 0.25

  final_fit <- .spacenet_pdhg_solver(
    X = X_centered,
    y = y01,
    alpha = alpha_opt,
    l1_ratio = l1_ratio,
    edges = edge_info$edges,
    d_norm2 = edge_info$d_norm2,
    loss = "logistic",
    max_iter = max_iter,
    tol = tol,
    data_lipschitz = data_L,
    init = NULL
  )

  beta <- numeric(p)
  beta[support] <- final_fit$w
  intercept <- as.numeric(final_fit$b - sum(x_mean * final_fit$w))

  list(
    beta = beta,
    intercept = intercept,
    support = support,
    alpha_opt = alpha_opt,
    alpha_grid = cv_path$alpha_grid,
    cv_loss = cv_path$cv_loss
  )
}

#' @keywords internal
#' @noRd
MVPAModels$spacenet_tvl1 <- list(
  type = c("Regression", "Classification"),
  library = "rMVPA",
  label = "spacenet_tvl1",
  loop = NULL,
  parameters = data.frame(
    parameters = c("l1_ratio", "n_alphas", "alpha_min_ratio", "screening_percentile", "max_iter", "tol", "inner_folds"),
    class = c("numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric"),
    label = c("L1 Ratio", "# Alphas", "Alpha Min Ratio", "Screening Percentile", "Max Iter", "Tolerance", "Inner CV Folds")
  ),
  grid = function(x, y, len = NULL) {
    data.frame(
      l1_ratio = 0.5,
      n_alphas = if (is.null(len) || len <= 1) 10 else max(6, as.integer(len)),
      alpha_min_ratio = 1e-3,
      screening_percentile = 20,
      max_iter = 150,
      tol = 1e-4,
      inner_folds = 3
    )
  },
  fit = function(x, y, wts, param, lev, last, weights, classProbs,
                 feature_ids = NULL, spatial_mask = NULL, ...) {
    X <- as.matrix(x)
    n <- nrow(X)
    p <- ncol(X)

    if (p < 2L || n < 4L) {
      stop("spacenet_tvl1: requires at least 4 samples and 2 features.")
    }

    is_classif <- is.factor(y)

    feature_ids <- if (is.null(feature_ids)) seq_len(p) else as.integer(feature_ids)
    if (length(feature_ids) != p) {
      feature_ids <- seq_len(p)
    }

    l1_ratio <- as.numeric(param$l1_ratio)
    l1_ratio <- max(0, min(1, ifelse(length(l1_ratio) == 0L || is.na(l1_ratio), 0.5, l1_ratio[1L])))
    n_alphas <- max(4L, as.integer(param$n_alphas[1L]))
    alpha_min_ratio <- as.numeric(param$alpha_min_ratio[1L])
    if (!is.finite(alpha_min_ratio) || alpha_min_ratio <= 0 || alpha_min_ratio >= 1) {
      alpha_min_ratio <- 1e-3
    }
    screening_percentile <- as.numeric(param$screening_percentile[1L])
    if (!is.finite(screening_percentile)) screening_percentile <- 100
    screening_percentile <- max(1, min(100, screening_percentile))
    max_iter <- max(50L, as.integer(param$max_iter[1L]))
    tol <- as.numeric(param$tol[1L])
    if (!is.finite(tol) || tol <= 0) tol <- 1e-4
    inner_folds <- max(2L, as.integer(param$inner_folds[1L]))

    if (is_classif) {
      n_classes <- length(lev)
      loss_name <- "logistic"

      if (n_classes <= 1L) {
        stop("spacenet_tvl1: classification requires at least 2 classes.")
      }

      if (n_classes == 2L) {
        bin_fit <- .spacenet_fit_binary(
          X = X,
          y01 = .spacenet_y01(y, positive = lev[2L]),
          feature_ids = feature_ids,
          spatial_mask = spatial_mask,
          l1_ratio = l1_ratio,
          n_alphas = n_alphas,
          alpha_min_ratio = alpha_min_ratio,
          screening_percentile = screening_percentile,
          max_iter = max_iter,
          tol = tol,
          inner_folds = inner_folds
        )

        beta <- bin_fit$beta
        intercept <- bin_fit$intercept
        support <- bin_fit$support
        alpha_opt <- bin_fit$alpha_opt
        alpha_grid <- bin_fit$alpha_grid
        cv_loss <- bin_fit$cv_loss
      } else {
        beta <- matrix(0, nrow = p, ncol = n_classes, dimnames = list(NULL, lev))
        intercept <- numeric(n_classes)
        names(intercept) <- lev
        support <- vector("list", n_classes)
        alpha_opt <- numeric(n_classes)
        names(alpha_opt) <- lev
        alpha_grid <- vector("list", n_classes)
        names(alpha_grid) <- lev
        cv_loss <- vector("list", n_classes)
        names(cv_loss) <- lev

        for (ci in seq_len(n_classes)) {
          cls <- lev[ci]
          bin_fit <- .spacenet_fit_binary(
            X = X,
            y01 = .spacenet_y01(y, positive = cls),
            feature_ids = feature_ids,
            spatial_mask = spatial_mask,
            l1_ratio = l1_ratio,
            n_alphas = n_alphas,
            alpha_min_ratio = alpha_min_ratio,
            screening_percentile = screening_percentile,
            max_iter = max_iter,
            tol = tol,
            inner_folds = inner_folds
          )
          beta[, ci] <- bin_fit$beta
          intercept[ci] <- bin_fit$intercept
          support[[ci]] <- bin_fit$support
          alpha_opt[ci] <- bin_fit$alpha_opt
          alpha_grid[[ci]] <- bin_fit$alpha_grid
          cv_loss[[ci]] <- bin_fit$cv_loss
        }
      }
    } else {
      cv_path <- .spacenet_cv_path(
        X = X,
        y = y,
        feature_ids = feature_ids,
        spatial_mask = spatial_mask,
        l1_ratio = l1_ratio,
        n_alphas = n_alphas,
        alpha_min_ratio = alpha_min_ratio,
        screening_percentile = screening_percentile,
        max_iter = max_iter,
        tol = tol,
        inner_folds = inner_folds,
        is_classif = FALSE
      )

      alpha_opt <- as.numeric(cv_path$best_alpha)
      alpha_grid <- cv_path$alpha_grid
      cv_loss <- cv_path$cv_loss

      support <- .spacenet_screen_support(
        X, y,
        screening_percentile = screening_percentile,
        feature_ids = feature_ids
      )
      if (length(support) == 0L) support <- seq_len(p)

      X_fit <- X[, support, drop = FALSE]
      fids_fit <- feature_ids[support]

      x_mean <- colMeans(X_fit)
      X_centered <- sweep(X_fit, 2L, x_mean, "-")

      y_offset <- mean(as.numeric(y))
      y_proc <- as.numeric(y) - y_offset
      loss_name <- "mse"

      edge_info <- .spacenet_edges_from_feature_ids(fids_fit, spatial_mask)
      spec2 <- .spacenet_spectral_norm_squared(X_centered)
      data_L <- spec2 / max(nrow(X_centered), 1L)

      final_fit <- .spacenet_pdhg_solver(
        X = X_centered,
        y = y_proc,
        alpha = alpha_opt,
        l1_ratio = l1_ratio,
        edges = edge_info$edges,
        d_norm2 = edge_info$d_norm2,
        loss = loss_name,
        max_iter = max_iter,
        tol = tol,
        data_lipschitz = data_L,
        init = NULL
      )

      beta <- numeric(p)
      beta[support] <- final_fit$w
      intercept <- as.numeric(y_offset - sum(x_mean * final_fit$w))
    }

    out <- list(
      beta = beta,
      intercept = intercept,
      support = support,
      alpha_opt = alpha_opt,
      alpha_grid = alpha_grid,
      cv_loss = cv_loss,
      l1_ratio = l1_ratio,
      screening_percentile = screening_percentile,
      loss = loss_name,
      obsLevels = if (is_classif) lev else NULL
    )
    class(out) <- c("spacenet_fit", "list")
    out
  },
  predict = function(modelFit, newdata, preProc = NULL, submodels = NULL) {
    X <- as.matrix(newdata)
    beta <- modelFit$beta
    intercept <- modelFit$intercept

    if (!is.null(modelFit$obsLevels)) {
      if (is.matrix(beta) && ncol(beta) > 1L) {
        eta <- sweep(X %*% beta, 2L, intercept, "+")
        pred <- colnames(eta)[max.col(eta)]
        factor(pred, levels = modelFit$obsLevels)
      } else {
        eta <- as.vector(X %*% as.numeric(beta) + as.numeric(intercept))
        pred <- ifelse(eta >= 0, modelFit$obsLevels[2L], modelFit$obsLevels[1L])
        factor(pred, levels = modelFit$obsLevels)
      }
    } else {
      as.vector(X %*% as.numeric(beta) + as.numeric(intercept))
    }
  },
  prob = function(modelFit, newdata, preProc = NULL, submodels = NULL) {
    if (is.null(modelFit$obsLevels)) {
      stop("spacenet_tvl1 probabilities are only available for classification.")
    }
    X <- as.matrix(newdata)
    beta <- modelFit$beta
    intercept <- modelFit$intercept

    if (is.matrix(beta) && ncol(beta) > 1L) {
      eta <- sweep(X %*% beta, 2L, intercept, "+")
      probs <- .spacenet_softmax(eta)
      colnames(probs) <- modelFit$obsLevels
      probs
    } else {
      eta <- as.vector(X %*% as.numeric(beta) + as.numeric(intercept))
      p2 <- .spacenet_sigmoid(eta)
      probs <- cbind(1 - p2, p2)
      colnames(probs) <- modelFit$obsLevels
      probs
    }
  }
)
