#' Create a Feature-Based RSA Design
#'
#' Creates a design for feature-based Representational Similarity Analysis (RSA).
#' You can either supply a similarity matrix S (and optionally select dimensions)
#' or directly supply a feature matrix F.
#'
#' @param S A symmetric similarity matrix representing the feature space relationships. 
#'          If NULL, you must supply F.
#' @param F A feature space matrix (observations by features). If supplied, this overrides S and k.
#' @param labels Vector of labels corresponding to the rows/columns of S or observations of F.
#' @param k Integer specifying the number of feature dimensions to retain when using S. If 0 (default),
#'          automatically determines dimensions using eigenvalue threshold > 1.
#'
#' @return A \code{feature_rsa_design} object (S3 class) containing:
#'   \describe{
#'     \item{S}{The input similarity matrix (if used)}
#'     \item{F}{Feature space projection matrix}
#'     \item{labels}{Vector of observation labels}
#'   }
#'
#' @details
#' If F is supplied directly, it is used as is.
#' If only S is supplied, the design computes eigen decomposition of S to find a suitable 
#' feature space (F).
#'
#' @export
feature_rsa_design <- function(S=NULL, F=NULL, labels, k=0) {
  assertthat::assert_that(!is.null(labels))
  
  if (!is.null(F)) {
    # If F is provided, trust user supplied features
    assertthat::assert_that(is.matrix(F))
    assertthat::assert_that(nrow(F) == length(labels))
    ret <- list(S=S, F=F, labels=labels)
  } else {
    # Must have S
    assertthat::assert_that(!is.null(S))
    assertthat::assert_that(is.matrix(S))
    assertthat::assert_that(nrow(S) == length(labels))
    assertthat::assert_that(isSymmetric(S))
    
    S <- (S + t(S))/2
    
    if (k == 0) {
      eres <- eigen(S)
      k <- max(which(eres$values > 1))
      k <- max(k, 2)
      F <- eres$vectors[, 1:k, drop=FALSE]
    } else {
      assertthat::assert_that(k > 0 && k <= nrow(S))
      eres <- eigen(S)
      F <- eres$vectors[, 1:k, drop=FALSE]
    }
    ret <- list(S=S, F=F, labels=labels)
  }
  
  class(ret) <- "feature_rsa_design"
  ret
}


#' Create a Feature-Based RSA Model
#'
#' Creates a model for feature-based Representational Similarity Analysis (RSA) that relates neural patterns
#' (X) to a predefined feature space (F).
#'
#' @param dataset An \code{mvpa_dataset} object containing the neural data (X).
#' @param design A \code{feature_rsa_design} object specifying the feature space (F).
#' @param method Character string specifying the analysis method. One of:
#'   \describe{
#'     \item{scca}{Sparse Canonical Correlation Analysis}
#'     \item{pls}{Partial Least Squares}
#'     \item{pca}{Principal Component Analysis + linear regression on F}
#'   }
#' @param crossval Optional cross-validation specification.
#'
#' @return A \code{feature_rsa_model} object (S3 class).
#'
#' @details
#' Feature RSA models analyze the relationship between neural patterns X and a predefined feature space F.
#' Methods:
#'   - scca: Finds canonical correlations between X and F and reconstructs F from X via these canonical directions.
#'   - pls: Uses partial least squares regression to predict F from X.
#'   - pca: PCA on X, then regress F_sc ~ PC(X_sc) for prediction.
#'
#' @export
feature_rsa_model <- function(dataset,
                              design,
                              method=c("scca", "pls", "pca"),
                              crossval=NULL) {
  
  method <- match.arg(method)
  
  assertthat::assert_that(inherits(dataset, "mvpa_dataset"))
  assertthat::assert_that(inherits(design, "feature_rsa_design"))
  
  if (is.null(crossval) && !is.null(design$block_var)) {
    crossval <- blocked_cross_validation(design$block_var)
  }
  
  assertthat::assert_that(!is.null(crossval))
  
  ret <- list(method=method,
              dataset=dataset,
              design=design,
              crossval=crossval)
  
  class(ret) <- "feature_rsa_model"
  ret
}


# Helper to standardize data
.standardize <- function(X) {
  cm <- colMeans(X)
  csd <- apply(X,2,sd)
  csd[csd==0] <- 1
  X_sc <- scale(X, center=cm, scale=csd)
  list(X_sc=X_sc, mean=cm, sd=csd)
}

# Given scca result, predict F from X:
# Steps:
# 1) Standardize X with training stats
# 2) Compute CCAX = X_sc %*% WX
# 3) CCAY = CCAX * lambda (canonical correlations)
# 4) Y_sc = CCAY %*% t(WY)
# 5) Unscale Y_sc to get Y_pred
.predict_scca <- function(model, X_new) {
  tm <- model$trained_model
  Xsc <- sweep(sweep(X_new,2,model$scca_x_mean,"-"),2,model$scca_x_sd,"/")
  CCAX_new <- Xsc %*% tm$WX
  CCAY_new <- sweep(CCAX_new,2,tm$lambda,"*")
  Y_sc <- CCAY_new %*% t(tm$WY)
  Y_pred <- sweep(sweep(Y_sc,2,model$scca_f_sd,"*"),2,model$scca_f_mean,"+")
  Y_pred
}


# For PCA method
.predict_pca <- function(model, X_new) {
  Xsc <- sweep(sweep(X_new,2,model$pca_x_mean,"-"),2,model$pca_x_sd,"/")
  PC_new <- Xsc %*% model$pcarot
  F_sc_pred <- cbind(1, PC_new) %*% model$pca_coefs
  F_pred <- sweep(sweep(F_sc_pred,2,model$pca_f_sd,"*"),2,model$pca_f_mean,"+")
  F_pred
}


#' @export
train_model.feature_rsa_model <- function(obj, train_dat, ytrain, ...) {
  X <- as.matrix(train_dat)
  Fsub <- as.matrix(ytrain)  # ytrain is the already subsetted portion of F
  
  if (obj$method == "pls") {
    obj$trained_model <- pls::plsr(Fsub ~ X, scale=TRUE)
  } else if (obj$method == "scca") {
    sx <- .standardize(X)
    sf <- .standardize(Fsub)
    scca_res <- whitening::scca(sx$X_sc, sf$X_sc, scale=FALSE)
    obj$trained_model <- scca_res
    obj$scca_x_mean <- sx$mean
    obj$scca_x_sd <- sx$sd
    obj$scca_f_mean <- sf$mean
    obj$scca_f_sd <- sf$sd
  } else if (obj$method == "pca") {
    sx <- .standardize(X)
    sf <- .standardize(Fsub)
    pca_res <- prcomp(sx$X_sc, scale.=FALSE)
    PC_train <- pca_res$x
    PC_train_i <- cbind(1, PC_train)
    coefs <- solve(t(PC_train_i) %*% PC_train_i, t(PC_train_i) %*% sf$X_sc)
    obj$trained_model <- pca_res
    obj$pcarot <- pca_res$rotation
    obj$pca_x_mean <- sx$mean
    obj$pca_x_sd <- sx$sd
    obj$pca_f_mean <- sf$mean
    obj$pca_f_sd <- sf$sd
    obj$pca_coefs <- coefs
  }
  obj
}


#' Predict Method for Feature RSA Model
#'
#' @param object The trained feature RSA model
#' @param newdata New data
#' @param ... Additional args
#' @return Predicted values in the feature space
#' @export
predict_model.feature_rsa_model <- function(object, newdata, ...) {
  X <- as.matrix(newdata)
  method <- object$method
  if (method == "pls") {
    pred_arr <- predict(object$trained_model, newdata=data.frame(X))
    ncomp <- object$trained_model$ncomp
    pred <- pred_arr[,,ncomp, drop=FALSE]
    pred <- pred[,,1]
    return(pred)
  } else if (method == "scca") {
    return(.predict_scca(object, X))
  } else if (method == "pca") {
    return(.predict_pca(object, X))
  }
}


#' Evaluate Method for Feature RSA Model
#'
#' @param object The feature RSA model
#' @param predicted Predicted values
#' @param observed Observed values
#' @param ... Additional args
#' @return Performance metrics
#' @export
evaluate_model.feature_rsa_model <- function(object, predicted, observed, ...) {
  cors <- diag(cor(predicted, observed))
  mse <- mean((predicted - observed)^2)
  
  list(
    correlations = cors,
    mse = mse
  )
}


#' @export
y_train.feature_rsa_model <- function(object) {
  object$F
}

format_result.feature_rsa_model <- function(obj, result, error_message=NULL, context, ...) {
  if (!is.null(error_message)) {
    return(tibble::tibble(
      observed=list(NULL),
      predicted=list(NULL),
      error=TRUE,
      error_message=error_message
    ))
  } else {
    # Predict on test data
    testX <- tibble::as_tibble(context$test, .name_repair=.name_repair)
    pred <- predict_model(obj, testX)
    
    # observed is ytest
    observed <- as.matrix(context$ytest)
    
    # Evaluate
    perf <- evaluate_model(obj, pred, observed)
    
    # Return a tibble
    # Store predicted and observed for optional inspection
    tibble::tibble(
      observed=list(observed),
      predicted=list(pred),
      performance=list(perf),
      error=FALSE,
      error_message="~"
    )
  }
}

#' Print Method for Feature RSA Model
#'
#' @param x The feature RSA model
#' @param ... Additional args
#' @export
print.feature_rsa_model <- function(x, ...) {
  cat("Feature RSA Model\n")
  cat("Method:", x$method, "\n")
  cat("Number of features:", ncol(x$design$F), "\n")
  cat("Number of observations:", nrow(x$design$F), "\n")
}

#' Summary Method for Feature RSA Model
#'
#' @param object The feature RSA model
#' @param ... Additional args
#' @export
summary.feature_rsa_model <- function(object, ...) {
  print(object)
  if (!is.null(object$trained_model)) {
    cat("\nModel Performance:\n")
    print(object$performance)
  }
}


#' Run regional RSA analysis on a specified Feature RSA model
#'
#' This function runs a regional analysis using a feature RSA model and region mask.
#'
#' @param model_spec A \code{feature_rsa_model} object.
#' @param region_mask A mask representing different brain regions.
#' @param coalesce_design_vars If TRUE, merges design variables into prediction table.
#' @param processor A custom processor function for ROIs. If NULL, uses defaults.
#' @param verbose Print progress messages.
#' @param ... Additional arguments
#' 
#' @details
#' This integrates `feature_rsa_model` with the MVPA framework, similar to `run_regional.mvpa_model`.
#'
#' @export
run_regional.feature_rsa_model <- function(model_spec, region_mask, coalesce_design_vars=FALSE, processor=NULL, 
                                           verbose=FALSE, ...) {  
  prepped <- prep_regional(model_spec, region_mask)
  
  # uses mvpa_iterate
  results <- mvpa_iterate(model_spec, prepped$vox_iter, ids=prepped$region_set, processor=processor, verbose=verbose, ...)
  
  perf <- if (model_spec$dataset$compute_performance) comp_perf(results, region_mask) else list(vols=list(), perf_mat=tibble::tibble())
  
  prediction_table <- if (model_spec$dataset$return_predictions) {
    combine_regional_results(results) 
  } else {
    NULL
  }
  
  if (coalesce_design_vars && !is.null(prediction_table)) {
    prediction_table <- coalesce_join(prediction_table, test_design(model_spec$design), 
                                      by=".rownum")
  }
  
  fits <- if (model_spec$dataset$return_fits) {
    lapply(results$result, "[[", "predictor")
  } else {
    NULL
  }
  
  regional_mvpa_result(model_spec=model_spec, performance_table=perf$perf_mat, 
                       prediction_table=prediction_table, vol_results=perf$vols, fits=fits)
}
