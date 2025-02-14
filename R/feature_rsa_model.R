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
  
  create_model_spec("feature_rsa_model", dataset=dataset, design=design, method=method, crossval=crossval, 
                    compute_performance=TRUE, return_fits=FALSE)
}


#' @noRd
.standardize <- function(X) {
  cm <- colMeans(X)
  csd <- apply(X,2,sd)
  csd[csd==0] <- 1
  X_sc <- scale(X, center=cm, scale=csd)
  list(X_sc=X_sc, mean=cm, sd=csd)
}

#' @noRd
.predict_scca <- function(model, X_new) {
  tm <- model$trained_model
  Xsc <- sweep(sweep(X_new,2,model$scca_x_mean,"-"),2,model$scca_x_sd,"/")
  xcoef <- t(tm$WX)
  ycoef <- t(tm$WY)
  Lambda_q <- diag(tm$lambda)
  y_inv      <- corpcor::pseudoinverse(ycoef)
  Yhat <- Xsc %*% xcoef %*% Lambda_q %*% y_inv
  Yhat <- sweep(sweep(Yhat, 2, model$scca_f_sd, "*"), 2, model$scca_f_mean, "+")
}

.predict_scca2 <- function(model, X_new) {
  tm <- model$trained_model
  
  # Standardize X
  Xsc <- sweep(sweep(X_new,2,model$scca_x_mean,"-"),2,model$scca_x_sd,"/")
  
  # Grab the canonical weights
  WX <- tm$WX  # shape: (#canonical, ncol(X_sc)) or (2 x 10)
  WY <- tm$WY  # shape: (#canonical, ncol(F_sc)) or (2 x 2)
  
  # Canonical correlations
  Lambda_q <- diag(tm$lambda)  # (2 x 2)
  
  Yhat_sc <- Xsc %*% t(WX) %*% Lambda_q %*% t(WY)
  
  # Unscale to original F space
  Yhat <- sweep(sweep(Yhat_sc, 2, model$scca_f_sd, "*"), 2, model$scca_f_mean, "+")
  Yhat
}


#' @noRd
.predict_pca <- function(model, X_new) {
  Xsc <- sweep(sweep(X_new,2,model$pca_x_mean,"-"),2,model$pca_x_sd,"/")
  PC_new <- Xsc %*% model$pcarot
  F_sc_pred <- cbind(1, PC_new) %*% model$pca_coefs
  F_pred <- sweep(sweep(F_sc_pred,2,model$pca_f_sd,"*"),2,model$pca_f_mean,"+")
  F_pred
}


#' @export
train_model.feature_rsa_model <- function(obj, train_dat, ytrain, indices, ...) {
  X <- as.matrix(train_dat)
  Fsub <- as.matrix(ytrain)  # ytrain is the already subsetted portion of F
  result <- list()
  if (obj$method == "pls") {
    df_train <- as.data.frame(X)
    if(is.null(colnames(df_train))) {
      colnames(df_train) <- paste0("V", 1:ncol(X))
    }
    fit <- pls::plsr(I(Fsub) ~ ., data = df_train, scale = TRUE)
    result$trained_model <- fit
  } else if (obj$method == "scca") {
    sx <- .standardize(X)
    sf <- .standardize(Fsub)
    scca_res <- whitening::scca(sx$X_sc, sf$X_sc, scale=FALSE)
    result$trained_model <- scca_res
    result$scca_x_mean <- sx$mean
    result$scca_x_sd <- sx$sd
    result$scca_f_mean <- sf$mean
    result$scca_f_sd <- sf$sd
  } else if (obj$method == "pca") {
    sx <- .standardize(X)
    sf <- .standardize(Fsub)
    pca_res <- prcomp(sx$X_sc, scale.=FALSE)
    PC_train <- pca_res$x
    # Truncate the number of principal components to avoid singularity; use up to 10 PCs
    k <- min(ncol(PC_train), 10)
    PC_train_subset <- PC_train[, 1:k, drop=FALSE]
    # Fit a multivariate linear model using lm(), predicting standardized F from the truncated PCs
    df <- as.data.frame(PC_train_subset)
    fit <- lm(sf$X_sc ~ ., data = df)
    coefs <- coef(fit)  # The coefficients matrix (Intercept and slopes) for each response variable
    result$trained_model <- pca_res
    # Store only the first k principal rotation vectors for use in prediction
    result$pcarot <- pca_res$rotation[, 1:k, drop=FALSE]
    result$pca_x_mean <- sx$mean
    result$pca_x_sd <- sx$sd
    result$pca_f_mean <- sf$mean
    result$pca_f_sd <- sf$sd
    result$pca_coefs <- coefs
  }
  result
}


#' Predict Method for Feature RSA Model
#'
#' @param object The trained feature RSA model
#' @param newdata New data
#' @param ... Additional args
#' @return Predicted values in the feature space
#' @export
predict_model.feature_rsa_model <- function(object, fit, newdata, ...) {
  X <- as.matrix(newdata)
  method <- object$method
  if (method == "pls") {
    df_new <- as.data.frame(X)
    if(is.null(colnames(df_new))) {
      colnames(df_new) <- paste0("V", 1:ncol(X))
    }
    pred_arr <- predict(fit$trained_model, newdata = df_new, ncomp = fit$ncomp)
    pred <- pred_arr[,, fit$ncomp, drop = FALSE]
    pred <- drop(pred)
    return(pred)
  } else if (method == "scca") {
    return(.predict_scca(fit, X))
  } else if (method == "pca") {
    return(.predict_pca(fit, X))
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
y_train.feature_rsa_model <- function(obj) {
  obj$design$F
}

#' @export
y_train.feature_rsa_design <- function(obj) {
  obj$F
}

#' @export
format_result.feature_rsa_model <- function(obj, result, error_message=NULL, context, ...) {
  if (!is.null(error_message)) {
    return(tibble::tibble(observed=list(NULL), predicted=list(NULL), error=TRUE, error_message=error_message))
  } else {
    Xtest <- tibble::as_tibble(context$test)
    pred <- predict_model.feature_rsa_model(obj, result, Xtest)
    
    observed <- as.matrix(context$ytest)
    perf <- evaluate_model.feature_rsa_model(obj, pred, observed)
    
    tibble::tibble(
      observed=list(observed),
      predicted=list(pred),
      performance=list(perf),
      error=FALSE,
      error_message="~"
    )
  }
}



#' Merge Multiple Results for Feature RSA Model
#'
#' @param obj A \code{feature_rsa_model} object
#' @param result_set A data frame of results from cross-validation folds
#' @param indices The voxel indices used (may not be relevant for feature_rsa_model)
#' @param id An identifier for the merged result (e.g., ROI id)
#' @param ... Additional arguments
#' @return A tibble with merged results
#' @export
merge_results.feature_rsa_model <- function(obj, result_set, indices, id, ...) {
  if (any(result_set$error)) {
    emessage <- result_set$error_message[which(result_set$error)[1]]
    return(tibble::tibble(
      result=list(NULL), indices=list(indices),
      performance=list(NULL), id=id,
      error=TRUE, error_message=emessage
    ))
  }
  
  # Combine observed/predicted across folds
  observed_list <- result_set$observed
  predicted_list <- result_set$predicted
  combined_observed <- do.call(rbind, observed_list)
  combined_predicted <- do.call(rbind, predicted_list)
  
  # Evaluate performance on combined data
  perf <- evaluate_model.feature_rsa_model(obj, combined_predicted, combined_observed)
  perf <- list(mse=perf$mse, correlation=mean(perf$correlations))
  
  # Store combined result
  combined_result <- list(
    observed=combined_observed,
    predicted=combined_predicted
  )
  
  
  tibble::tibble(
    result=list(combined_result),
    indices=list(indices),
    performance=list(perf),
    id=id,
    error=FALSE,
    error_message="~"
  )
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
 

  perf <- if (model_spec$compute_performance) comp_perf(results, region_mask) else list(vols=list(), perf_mat=tibble::tibble())
  
  prediction_table <- if (model_spec$return_predictions) {
    combine_regional_results(results)
  } else {
    NULL
  }

  if (coalesce_design_vars && !is.null(prediction_table)) {
    prediction_table <- coalesce_join(prediction_table, test_design(model_spec$design),
                                      by=".rownum")
  }

  fits <- if (model_spec$return_fits) {
    lapply(results$result, "[[", "predictor")
  } else {
    NULL
  }

  regional_mvpa_result(model_spec=model_spec, performance_table=perf$perf_mat,
                       prediction_table=prediction_table, vol_results=perf$vols, fits=fits)
}

#' Print Method for Feature RSA Design
#'
#' @param x A feature_rsa_design object.
#' @param ... Additional arguments (ignored).
#' @export
print.feature_rsa_design <- function(x, ...) {
  # Create a border line for styling
  border <- crayon::bold(crayon::cyan(strrep("=", 50)))
  
  # Header
  cat(border, "\n")
  cat(crayon::bold(crayon::cyan("          Feature RSA Design          \n")))
  cat(border, "\n\n")
  
  # Extract key details
  n_obs <- nrow(x$F)
  n_feat <- ncol(x$F)
  
  # Print number of observations and feature dimensions
  cat(crayon::bold(crayon::green("Number of Observations: ")), n_obs, "\n")
  cat(crayon::bold(crayon::green("Feature Dimensions:     ")), n_feat, "\n")
  
  # Indicate whether a similarity matrix was provided
  if (!is.null(x$S)) {
    cat(crayon::bold(crayon::magenta("Similarity Matrix:      ")), "Provided\n")
  } else {
    cat(crayon::bold(crayon::magenta("Similarity Matrix:      ")), 
        "Not provided (using feature matrix F directly)\n")
  }
  
  # Print first few labels
  n_labels <- length(x$labels)
  n_to_print <- min(5, n_labels)
  label_str <- paste(x$labels[1:n_to_print], collapse = ", ")
  if (n_labels > n_to_print) {
    label_str <- paste0(label_str, ", ...")
  }
  cat(crayon::bold(crayon::yellow("Labels (first few):   ")), label_str, "\n")
  
  # Footer
  cat("\n", border, "\n")
}

#' Print Method for Feature RSA Model
#'
#' @param x A feature_rsa_model object.
#' @param ... Additional arguments (ignored).
#' @export
print.feature_rsa_model <- function(x, ...) {
  # Create a border line for styling
  border <- crayon::bold(crayon::cyan(strrep("=", 50)))
  
  # Header
  cat(border, "\n")
  cat(crayon::bold(crayon::cyan("          Feature RSA Model           \n")))
  cat(border, "\n\n")
  
  # Display the method used (e.g., scca, pls, or pca)
  cat(crayon::bold(crayon::green("Method: ")), x$method, "\n")
  
  # Check if the design component is present to extract dimensions
  if (!is.null(x$design)) {
    n_obs <- nrow(x$design$F)
    n_feat <- ncol(x$design$F)
  } else {
    n_obs <- "Unknown"
    n_feat <- "Unknown"
  }
  
  cat(crayon::bold(crayon::green("Number of Observations: ")), n_obs, "\n")
  cat(crayon::bold(crayon::green("Feature Dimensions:     ")), n_feat, "\n")
  
  # Indicate training status
  if (!is.null(x$trained_model)) {
    cat(crayon::bold(crayon::magenta("Status: ")), "Trained model available\n")
  } else {
    cat(crayon::bold(crayon::magenta("Status: ")), "Model not yet trained\n")
  }
  
  # Display cross-validation status
  if (!is.null(x$crossval)) {
    cat(crayon::bold(crayon::yellow("Cross-Validation: ")), "Configured\n")
  } else {
    cat(crayon::bold(crayon::yellow("Cross-Validation: ")), "Not configured\n")
  }
  
  # Footer
  cat("\n", border, "\n")
}

