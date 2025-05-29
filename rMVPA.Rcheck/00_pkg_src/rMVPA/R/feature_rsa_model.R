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
#'          automatically determines dimensions using eigenvalue threshold > 1 (minimum 2 dimensions kept).
#'          This parameter is ignored if F is supplied directly (k becomes ncol(F)).
#' @param max_comps Initial upper limit for the number of components to be derived from the
#'          feature space F by subsequent `feature_rsa_model` methods (PCA, PLS, SCCA).
#'          This value is automatically capped by the final feature dimensionality `k`. Default 10.
#' @param block_var Optional blocking variable for cross-validation. If provided and
#'          `crossval` is `NULL` in `feature_rsa_model`, a blocked cross-validation
#'          scheme will be generated using this vector.
#'
#' @return A \code{feature_rsa_design} object (S3 class) containing:
#'   \describe{
#'     \item{S}{The input similarity matrix (if used)}
#'     \item{F}{Feature space projection matrix (k dimensions)}
#'     \item{labels}{Vector of observation labels}
#'     \item{k}{The final number of feature dimensions used}
#'     \item{max_comps}{The upper limit on components (<= k)}
#'     \item{block_var}{Optional blocking variable for cross-validation}
#'   }
#'
#' @details
#' This function defines the feature space representation for the analysis.
#' If F is supplied directly, it is used as-is, and `k` becomes `ncol(F)`.
#' If only S is supplied, an eigen decomposition of S is performed.
#' `k` determines how many eigenvectors form the feature matrix F. If `k=0`,
#' dimensions with eigenvalues > 1 are kept (minimum 2).
#' `max_comps` sets an upper bound for the number of components that model-fitting
#' methods (like PCA, PLS, SCCA in `feature_rsa_model`) can use, and it cannot
#' exceed the final feature dimensionality `k`.
#'
#' @export
feature_rsa_design <- function(S=NULL, F=NULL, labels, k=0, max_comps=10, block_var=NULL) {
  assertthat::assert_that(!is.null(labels))
  
  if (!is.null(F)) {
    assertthat::assert_that(is.matrix(F))
    assertthat::assert_that(nrow(F) == length(labels))
    k_value <- ncol(F)
    max_comps <- min(max_comps, k_value)
    
    ret <- list(
      S=S,
      F=F,
      labels=labels,
      k=k_value,
      max_comps=max_comps,
      block_var=block_var
    )
  } else {
    assertthat::assert_that(!is.null(S))
    assertthat::assert_that(is.matrix(S))
    assertthat::assert_that(nrow(S) == length(labels))
    assertthat::assert_that(isSymmetric(S))
    
    S <- (S + t(S)) / 2
    eres <- eigen(S)
    vals <- eres$values
    vecs <- eres$vectors
    
    if (k == 0) {
      k <- sum(vals > 1)
      k <- max(k, 2)
    } else {
      assertthat::assert_that(k > 0 && k <= nrow(S))
    }
    F <- vecs[, seq_len(k), drop=FALSE]
    max_comps <- min(max_comps, k)
    
    ret <- list(
      S=S,
      F=F,
      labels=labels,
      k=k,
      max_comps=max_comps,
      block_var=block_var
    )
  }
  
  class(ret) <- "feature_rsa_design"
  ret
}


#' Create a Feature-Based RSA Model
#'
#' Creates a model for feature-based Representational Similarity Analysis (RSA) that relates neural patterns
#' (X) to a predefined feature space (F).
#'
#' @param dataset An \code{mvpa_dataset} object containing the neural data (\code{X}).
#' @param design A \code{feature_rsa_design} object specifying the feature space (\code{F})
#'   and including the component limit (`max_comps`).
#' @param method Character string specifying the analysis method. One of:
#'   \describe{
#'     \item{scca}{Sparse Canonical Correlation Analysis relating X and F.}
#'     \item{pls}{Partial Least Squares regression predicting X from F.}
#'     \item{pca}{Principal Component Analysis on F, followed by regression predicting X from the PCs.}
#'     \item{glmnet}{Elastic net regression predicting X from F using glmnet with multivariate Gaussian response.}
#'   }
#' @param crossval Optional cross-validation specification.
#' @param cache_pca Logical, if TRUE and method is "pca", cache the PCA decomposition of the
#'   feature matrix F across cross-validation folds involving the same training rows.
#'   Defaults to FALSE.
#' @param alpha Numeric value between 0 and 1, only used when method="glmnet". Controls the elastic net
#'   mixing parameter: 1 for lasso (default), 0 for ridge, values in between for a mixture.
#'   Defaults to 0.5 (equal mix of ridge and lasso).
#' @param cv_glmnet Logical, if TRUE and method="glmnet", use cv.glmnet to automatically select the
#'   optimal lambda value via cross-validation. Defaults to FALSE.
#' @param lambda Optional numeric value or sequence of values, only used when method="glmnet" and
#'   cv_glmnet=FALSE. Specifies the regularization parameter. If NULL (default), a sequence will be 
#'   automatically determined by glmnet.
#' @param nperm Integer, number of permutations to run for statistical testing of model performance
#'   metrics after merging cross-validation folds. Default 0 (no permutation testing).
#' @param permute_by DEPRECATED. Permutation is always done by shuffling rows of the predicted matrix.
#' @param save_distributions Logical, if TRUE and nperm > 0, save the full null distributions
#'   from the permutation test. Defaults to FALSE.
#' @param ... Additional arguments (currently unused).
#'
#' @return A \code{feature_rsa_model} object (S3 class).
#'
#' @details
#' Feature RSA models analyze how well a feature matrix \code{F} (defined in the `design`)
#' relates to neural data \code{X}. The `max_comps` parameter, inherited from the `design` object,
#' sets an upper limit on the number of components used:
#'   - \strong{pca}: Performs PCA on \code{F}. `max_comps` limits the number of principal components
#'     (selected by variance explained) used to predict \code{X}. Actual components used: `min(max_comps, available_PCs)`.
#'   - \strong{pls}: Performs PLS regression predicting \code{X} from \code{F}. `max_comps` sets the
#'     maximum number of PLS components to compute. Actual components used may be fewer based on the PLS algorithm.
#'   - \strong{scca}: Performs SCCA between \code{X} and \code{F}. `max_comps` limits the number of
#'     canonical components retained (selected by correlation strength). Actual components used: `min(max_comps, effective_components)`.
#'   - \strong{glmnet}: Performs elastic net regression predicting \code{X} from \code{F} using the glmnet package
#'     with multivariate Gaussian response family. The regularization (lambda) can be automatically selected via cross-validation
#'     if cv_glmnet=TRUE. The alpha parameter controls the balance between L1 (lasso) and L2 (ridge) regularization.
#'
#' **Performance Metrics** (computed by `evaluate_model` after cross-validation):
#'   - `mean_correlation`: Average correlation between predicted and observed patterns for corresponding trials/conditions (diagonal of the prediction-observation correlation matrix).
#'   - `cor_difference`: The `mean_correlation` minus the average off-diagonal correlation (`mean_correlation` - `off_diag_correlation`). Measures how much better the model predicts the correct trial/condition compared to incorrect ones.
#'   - `mean_rank_percentile`: Average percentile rank of the diagonal correlations. For each condition, ranks how well the model's prediction correlates with the correct observed pattern compared to incorrect patterns. Values range from 0 to 1, with 0.5 expected by chance and 1 indicating perfect discrimination.
#'   - `voxel_correlation`: Correlation between the vectorized predicted and observed data matrices across all trials and voxels.
#'   - `mse`: Mean Squared Error between predicted and observed values.
#'   - `r_squared`: Proportion of variance in the observed data explained by the predicted data.
#'   - `p_*`, `z_*`: If `nperm > 0`, permutation-based p-values and z-scores for the above metrics, assessing significance against a null distribution generated by shuffling predicted trial labels.
#'
#' The number of components actually used (`ncomp`) for the region/searchlight is also included in the performance output.
#'
#' @export
feature_rsa_model <- function(dataset,
                               design,
                               method = c("scca", "pls", "pca", "glmnet"),
                               crossval = NULL,
                               cache_pca = FALSE,
                               alpha = 0.5,
                               cv_glmnet = FALSE,
                               lambda = NULL,
                               nperm = 0,
                               permute_by = c("features", "observations"),
                               save_distributions = FALSE,
                               ...) {
  
  method <- match.arg(method)
  permute_by <- match.arg(permute_by)
  assertthat::assert_that(inherits(dataset, "mvpa_dataset"))
  assertthat::assert_that(inherits(design, "feature_rsa_design"))
  
  # Additional validation for dataset dimensions
  mask_dims <- dim(dataset$mask)[1:3]
  total_voxels <- prod(mask_dims)
  active_voxels <- sum(dataset$mask > 0)
  
  if (total_voxels <= 1) {
    stop("Invalid dataset for feature_rsa_model: Only 1 voxel detected (dimensions ",
         paste(mask_dims, collapse="×"),
         "). Feature RSA analysis requires multiple voxels.")
  }
  
  if (active_voxels <= 1) {
    stop("Invalid dataset for feature_rsa_model: Only ", active_voxels,
         " active voxel(s) in mask. Feature RSA analysis requires multiple active voxels.")
  }
  
  if (is.null(crossval)) {
    if (!is.null(design$block_var)) {
      crossval <- blocked_cross_validation(design$block_var)
    } else {
      stop("crossval must be provided or design must include block_var")
    }
  }
  assertthat::assert_that(!is.null(crossval))
  
  # GLMNet specific validation
  if (method == "glmnet") {
    # Check if glmnet is available
    if (!requireNamespace("glmnet", quietly = TRUE)) {
      stop("Package 'glmnet' is required for glmnet method. Please install it with: install.packages('glmnet')")
    }
    
    # Validate alpha parameter
    assertthat::assert_that(is.numeric(alpha) && alpha >= 0 && alpha <= 1,
                           msg = "alpha must be between 0 and 1 (inclusive)")
    
    # Validate lambda if provided
    if (!is.null(lambda)) {
      assertthat::assert_that(is.numeric(lambda) && all(lambda > 0),
                             msg = "lambda must be positive")
    }
  }
  
  max_comps <- design$max_comps
  
  model_spec <- create_model_spec(
    "feature_rsa_model", 
    dataset = dataset,
    design  = design, 
    method  = method, 
    crossval= crossval, 
    compute_performance = TRUE, 
    return_fits = FALSE
  )
  
  # Single "max_comps" in use
  model_spec$max_comps <- max_comps
  
  # new logic for caching
  model_spec$cache_pca <- cache_pca
  
  # We'll store PCA objects in a named list:
  # each entry is keyed by a "hash" of the training row indices
  model_spec$pca_cache <- list()
  
  # GLMNet parameters
  if (method == "glmnet") {
    model_spec$alpha <- alpha
    model_spec$cv_glmnet <- cv_glmnet
    model_spec$lambda <- lambda
  }
  
  model_spec$nperm <- nperm
  model_spec$permute_by <- permute_by
  model_spec$save_distributions <- save_distributions
  
  model_spec
}


#' @noRd
.hash_row_indices <- function(idx) {
  # Sort them so that permutations of the same set produce the same key
  idx_sorted <- sort(idx)
  # Convert to string
  paste0(idx_sorted, collapse="_")
  # Or use a real hash function, e.g., digest::digest(idx_sorted)
}

#' @noRd
.standardize <- function(X) {
  cm <- colMeans(X)
  csd <- apply(X, 2, sd)
  csd[csd == 0] <- 1
  X_sc <- scale(X, center=cm, scale=csd)
  list(X_sc = X_sc, mean = cm, sd = csd)
}



#' @noRd
.predict_scca <- function(model, F_new) {
  # SCCA prediction predicts X from F
  # Check if training failed to find components (ncomp=0)
  if (is.null(model$ncomp) || model$ncomp < 1) {
      futile.logger::flog.warn(".predict_scca: SCCA training found no components (ncomp=%s). Using mean brain pattern fallback for prediction.", 
                                ifelse(is.null(model$ncomp), "NULL", model$ncomp))
      
      # Fallback: Predict the mean brain pattern from training
      n_test <- nrow(F_new)
      mean_x_train <- model$scca_x_mean
      n_voxels <- length(mean_x_train)
      
      if (n_test <= 0 || n_voxels <= 0) {
         stop("Cannot create fallback prediction: invalid dimensions (n_test=%d, n_voxels=%d)", n_test, n_voxels)
      }
      
      # Create matrix repeating the mean training pattern for each test sample
      X_pred_fallback <- matrix(mean_x_train, nrow = n_test, ncol = n_voxels, byrow = TRUE)
      return(X_pred_fallback)
  }
  
  # --- Proceed with standard SCCA prediction if ncomp >= 1 ---
  
  # First check dimensions and ensure compatibility
  F_new <- as.matrix(F_new)
  n_features_new <- ncol(F_new)
  n_features_expected <- length(model$scca_f_mean)
  tm <- model$trained_model
  ncomp <- model$ncomp
  
  # Strictly enforce dimensional consistency
  if (n_features_new != n_features_expected) {
    stop(sprintf("Feature matrix dimension mismatch for SCCA: expected %d features but got %d. This indicates a data inconsistency between training and testing.",
                n_features_expected, n_features_new))
  }
  
  # Standard case: dimensions match
  # Standardize features using the stored means and standard deviations
  Fsc <- sweep(sweep(F_new, 2, model$scca_f_mean, "-"), 2, model$scca_f_sd, "/")
  
  # Use canonical directions from F to predict X
  fcoef <- t(tm$WY)[, 1:ncomp, drop=FALSE]
  xcoef <- t(tm$WX)[, 1:ncomp, drop=FALSE]
  x_inv <- corpcor::pseudoinverse(xcoef)
  
  # Project features to canonical space
  F_canonical <- Fsc %*% fcoef

  # -- Incorporate canonical correlations for prediction --
  # Get the canonical correlations corresponding to the selected components
  canonical_corrs <- tm$lambda[1:ncomp]
  if (is.null(canonical_corrs) || length(canonical_corrs) != ncomp) {
      stop(sprintf("Prediction error: Could not retrieve the expected %d canonical correlations.", ncomp))
  }
  
  # Scale F canonical variates by the canonical correlations to get predicted U variates
  # U_pred = V_test * Lambda
  U_pred <- sweep(F_canonical, 2, canonical_corrs, "*")
  
  # Map predicted U variates back to standardized X space using pseudoinverse of X weights
  # Xhat_sc = U_pred * Wx^+
  Xhat <- U_pred %*% x_inv
  # --------------------------------------------------------

  # Xhat <- F_canonical %*% x_inv # OLD direct mapping
  Xhat <- sweep(sweep(Xhat, 2, model$scca_x_sd, "*"), 2, model$scca_x_mean, "+")
  return(Xhat)
}


#' @noRd
.predict_pca <- function(model, F_new) {
  # F_new is test features (subset for that fold).
  F_new <- as.matrix(F_new)
  
  # Standardize with the same mean/sd from training
  f_means <- model$pca_f_mean
  f_sds   <- model$pca_f_sd
  
  if (length(f_means) != ncol(F_new)) {
    stop("Mismatch in the number of features for PCA prediction vs training.")
  }
  
  Fsc <- sweep(sweep(F_new, 2, f_means, "-"), 2, f_sds, "/")
  
  # Project onto the fold-trained PCA rotation
  PC_new <- Fsc %*% model$pcarot
  
  # Predict X in standardized space
  X_sc_pred <- cbind(1, PC_new) %*% model$pca_coefs  # dim => [nTest, nVoxelsInRegion (or columns in X)]
  
  # Undo standardization of X
  x_means <- model$pca_x_mean
  x_sds   <- model$pca_x_sd
  
  X_pred <- sweep(sweep(X_sc_pred, 2, x_sds, "*"), 2, x_means, "+")
  X_pred
}


#' @importFrom glmnet cv.glmnet glmnet 
#' @noRd
.predict_glmnet <- function(model, F_new) {
  # F_new is test features (subset for that fold)
  F_new <- as.matrix(F_new)
  
  # Check feature dimensions
  if (ncol(F_new) != length(model$glmnet_f_mean)) {
    stop(sprintf("Feature matrix dimension mismatch for GLMNet: expected %d features but got %d.",
                 length(model$glmnet_f_mean), ncol(F_new)))
  }
  
  
  # Standardize features using the stored means and standard deviations
  Fsc <- sweep(sweep(F_new, 2, model$glmnet_f_mean, "-"), 2, model$glmnet_f_sd, "/")
  
  # Make predictions with the trained glmnet model
  # For mgaussian family, predict returns a list with one matrix per response
  preds_mat <- if (model$cv_glmnet) {
    # Use the CV-selected lambda
    drop(predict(model$trained_model, newx = Fsc, s = "lambda.min"))
  } else {
    # Use the single lambda or the first lambda in sequence
    drop(predict(model$trained_model, newx = Fsc, s = model$lambda_used))
  }


  
  # For mgaussian, preds_scaled is a list where each element is a matrix
  # with number of rows = nrow(F_new) and columns = 1 (or nlambda if using multiple lambda values)
  # We need to convert this to a matrix of the same dimensions as X
  
  # Determine the number of response variables (voxels)
  nvox <- ncol(preds_mat)
  
  # Unstandardize: convert from z-scores back to original scale
  X_pred <- sweep(sweep(preds_mat, 2, model$glmnet_x_sd, "*"), 2, model$glmnet_x_mean, "+")
  
  return(X_pred)
}


#' @export
predict_model.feature_rsa_model <- function(object, fit, newdata, ...) {
  # Check if the 'fit' object contains an error from the training stage
  if (!is.null(fit$error)) {
     error_msg <- sprintf("predict_model: Cannot predict, training failed with error: %s", fit$error)
     futile.logger::flog.error(error_msg)
     stop(error_msg) # Stop prediction if training failed
  }
  
  # Check if trained_model is missing, even if no explicit error was set
  if (is.null(fit$trained_model) && object$method != "scca") { # SCCA might proceed with ncomp=0
      error_msg <- sprintf("predict_model (%s): 'trained_model' is missing in the fit object provided. Cannot predict.", object$method)
      futile.logger::flog.error(error_msg)
      stop(error_msg)
  }

 
  method <- object$method
  F_new  <- as.matrix(newdata)
  
  # Basic check for newdata dimensions
  if (nrow(F_new) < 1) {
      stop("predict_model: newdata (F_new) has 0 rows.")
  }

  # Wrap the entire prediction logic in tryCatch
  predictions <- tryCatch({
    if (method == "pls") {
      
      # Retrieve necessary components from 'fit'
      pls_model <- fit$trained_model
      f_mean <- fit$pls_f_mean
      f_sd <- fit$pls_f_sd
      x_mean <- fit$pls_x_mean
      x_sd <- fit$pls_x_sd
      ncomp_to_use <- fit$ncomp

      # Validate fit components
      if (is.null(pls_model) || is.null(f_mean) || is.null(f_sd) || is.null(x_mean) || is.null(x_sd) || is.null(ncomp_to_use)) {
        stop("predict_model (PLS): Missing essential components (model, means, sds, ncomp) in the fit object.")
      }
      if (ncomp_to_use < 1) {
         # This case should ideally be handled by training setting ncomp appropriately,
         # but predict.plsr might behave unexpectedly with ncomp=0.
         # Let's return mean prediction similar to SCCA? Or error out?
         # For now, let's error, assuming training should prevent ncomp=0 unless intended.
         stop(sprintf("predict_model (PLS): Number of components to use (ncomp=%d) is less than 1.", ncomp_to_use))
      }

      expected_cols <- length(f_mean)
      actual_cols <- ncol(F_new)
      if (actual_cols != expected_cols) {
        stop(sprintf("predict_model (PLS): Feature column mismatch. Model expects %d columns, newdata has %d.",
                     expected_cols, actual_cols))
      }
      
      # Standardize newdata using training parameters
      sf_test <- scale(F_new, center = f_mean, scale = f_sd)
      # Check for NaNs after scaling (can happen if sd was zero and scale didn't handle it perfectly)
      if (any(is.nan(sf_test))) {
         stop("predict_model (PLS): NaNs generated during standardization of newdata. Check feature variances.")
      }

      # Make predictions using pls::predict
      # preds_raw is expected to be the prediction using ncomp_to_use components,
      # potentially in a 3D array [obs, resp, 1]
      preds_raw <- predict(pls_model, newdata = sf_test, ncomp = ncomp_to_use)
      
      # Drop any singleton dimensions (e.g., the component dimension if size 1)
      preds_sc <- drop(preds_raw)

      # Ensure preds_sc is a matrix even if only one observation/response
      if (!is.matrix(preds_sc)) {
         preds_sc <- matrix(preds_sc, nrow=nrow(F_new), ncol=length(x_mean))
      }

      # Un-standardize: X = X_sc * sd + mean
      preds <- sweep(sweep(preds_sc, 2, x_sd, "*"), 2, x_mean, "+")
      return(preds)
      
    } else if (method == "scca") {
      # SCCA prediction logic is in .predict_scca
      # Pass the entire 'fit' object which contains model and standardization params
      return(.predict_scca(fit, F_new)) 
      
    } else if (method == "pca") {
      # PCA prediction logic is in .predict_pca
      # Pass the entire 'fit' object
      return(.predict_pca(fit, F_new))
      
    } else if (method == "glmnet") {
      # GLMNet prediction logic is in .predict_glmnet
      # Pass the entire 'fit' object
      return(.predict_glmnet(fit, F_new))
      
    } else {
      stop(paste("Unknown method in predict_model.feature_rsa_model:", method))
    }
  }, error = function(e) {
      # Catch any error from the specific prediction methods or checks
      error_msg <- sprintf("predict_model (%s): Prediction failed - %s", method, e$message)
      futile.logger::flog.error(error_msg)
      # Re-throw the error so format_result can catch it
      stop(error_msg)
  })
  
  # Final check on prediction output
  if (is.null(predictions) || !is.matrix(predictions)) {
     error_msg <- sprintf("predict_model (%s): Prediction result is NULL or not a matrix. Check internal prediction logic.", method)
     futile.logger::flog.error(error_msg)
     stop(error_msg)
  }
   if (nrow(predictions) != nrow(F_new)) {
     error_msg <- sprintf("predict_model (%s): Prediction result has %d rows, but expected %d (matching newdata).", 
                         method, nrow(predictions), nrow(F_new))
     futile.logger::flog.error(error_msg)
     stop(error_msg)
  }
  
  return(predictions)
}


#' @noRd
#' @keywords internal
#' Helper that performs permutation testing for Feature RSA
#' 
#' @param observed Matrix of observed data
#' @param predicted Matrix of predicted data
#' @param nperm Number of permutations
#' @param save_distributions Logical, whether to save all permutation distributions
#' @param mean_cor Scalar: the observed mean correlation
#' @param cor_difference Scalar: the observed correlation difference (mean_cor - off_diag_cor)
#' @param mean_rank_percentile Scalar: the observed mean percentile rank of diagonal correlations
#' @param voxel_cor Scalar: the observed voxel correlation
#' @param mse Scalar: the observed MSE
#' @param r_squared Scalar: the observed R^2
#' @param cor_temporal_means Scalar: the observed correlation of temporal means
#' @param mean_voxelwise_temporal_cor Scalar: the observed mean voxelwise temporal correlation
#' @param cors Vector of diagonal correlations
#' @return A list with p-values, z-scores, and optionally a list of permutations
.perm_test_feature_rsa <- function(observed,
                                   predicted,
                                   nperm,
                                   save_distributions,
                                   mean_cor,
                                   cor_difference,
                                   mean_rank_percentile,
                                   voxel_cor,
                                   mse,
                                   r_squared,
                                   cor_temporal_means,
                                   mean_voxelwise_temporal_cor,
                                   cors)
{
  message("Performing permutation tests with ", nperm, " permutations... (feature_rsa_model)")
  
  n_rows <- nrow(predicted)
  n_cols <- ncol(predicted)
  rss <- sum((observed - predicted)^2)
  tss <- sum((observed - mean(observed))^2)

  # Counters for how many permutations are "better" than the true model
  count_better_mean_corr  <- 0
  count_better_cor_diff <- 0
  count_better_rank_perc <- 0
  count_better_voxel_corr <- 0
  count_better_mse        <- 0
  count_better_r_squared  <- 0
  count_better_cor_temp_means <- 0 # New
  count_better_mean_vox_temp_cor <- 0 # New

  # For mean & SD (to calculate z-scores)
  sum_mean_corr    <- 0; sum_sq_mean_corr <- 0
  sum_cor_diff    <- 0; sum_sq_cor_diff <- 0
  sum_rank_perc    <- 0; sum_sq_rank_perc <- 0
  sum_voxel_corr    <- 0; sum_sq_voxel_corr <- 0
  sum_mse    <- 0; sum_sq_mse <- 0
  sum_r_squared    <- 0; sum_sq_r_squared <- 0
  sum_cor_temp_means <- 0; sum_sq_cor_temp_means <- 0 # New
  sum_mean_vox_temp_cor <- 0; sum_sq_mean_vox_temp_cor <- 0 # New

  # Optionally store entire distributions
  if (save_distributions) {
    perm_mean_corr     <- numeric(nperm)
    perm_cor_diff      <- numeric(nperm)
    perm_rank_perc     <- numeric(nperm)
    perm_voxel_corr    <- numeric(nperm)
    perm_mse           <- numeric(nperm)
    perm_r_squared     <- numeric(nperm)
    perm_correlations  <- matrix(NA, nrow=nperm, ncol=length(cors))
    perm_cor_temp_means <- numeric(nperm) # New
    perm_mean_vox_temp_cor <- numeric(nperm) # New
  }
  
  # Pre-calculate observed spatial means for efficiency in loop
  mean_obs_across_space <- rowMeans(observed, na.rm = TRUE)
  
  for (i in seq_len(nperm)) {
    # Permute row order of predicted
    perm_idx <- sample(n_rows)
    perm_pred <- predicted[perm_idx, , drop=FALSE]

    # --- Compute metrics for this permutation ---
    # Standard metrics
    perm_cormat <- tryCatch(cor(perm_pred, observed), error=function(e) matrix(NA, nrow=n_rows, ncol=n_rows))
    perm_cors <- diag(perm_cormat)
    pmc  <- mean(perm_cors, na.rm=TRUE)
    
    n <- nrow(perm_cormat)
    perm_off_diag <- NA_real_
    if (n > 1 && sum(!is.na(perm_cormat)) > sum(!is.na(perm_cors))) {
        perm_off_diag <- (sum(perm_cormat, na.rm=TRUE) - sum(perm_cors, na.rm=TRUE)) / (n*n - n)
    }
    pcd <- pmc - perm_off_diag  # Correlation difference
    
    perm_ranks <- rep(NA_real_, n)
    if (n > 1) {
        for (j in 1:n) {
            condition_cors <- perm_cormat[j, ]
            n_valid_cors <- sum(!is.na(condition_cors))
            if (n_valid_cors > 1) {
                perm_ranks[j] <- (sum(condition_cors <= condition_cors[j], na.rm = TRUE) - 1) / (n_valid_cors - 1)
            }
        }
    }
    prp <- mean(perm_ranks, na.rm=TRUE)  # Mean rank percentile
    
    pvc  <- tryCatch(cor(as.vector(perm_pred), as.vector(observed)), error=function(e) NA_real_)
    pmse <- mean((perm_pred - observed)^2, na.rm=TRUE)
    prsq <- NA_real_
    if (tss > 0) {
      rss_perm <- sum((observed - perm_pred)^2, na.rm=TRUE)
      prsq <- 1 - rss_perm / tss
    }

    # Correlation of Temporal Means for permutation
    pctm <- NA_real_
    if (nrow(perm_pred) > 1) {
       mean_perm_pred_across_space <- rowMeans(perm_pred, na.rm = TRUE)
       pctm <- tryCatch(cor(mean_perm_pred_across_space, mean_obs_across_space), error=function(e) NA_real_)
    }

    # Mean Voxelwise Temporal Correlation for permutation
    pmvtc <- NA_real_
    if (nrow(perm_pred) > 1 && n_cols > 0) {
       perm_voxel_cors <- numeric(n_cols)
       for (k in 1:n_cols) {
           perm_voxel_cors[k] <- tryCatch(cor(observed[, k], perm_pred[, k]), error=function(e) NA)
       }
       pmvtc <- mean(perm_voxel_cors, na.rm = TRUE)
    }

    # --- Update counters --- 
    if (!is.na(pmc) && !is.na(mean_cor) && pmc >= mean_cor) count_better_mean_corr <- count_better_mean_corr + 1
    if (!is.na(pcd) && !is.na(cor_difference) && pcd >= cor_difference) count_better_cor_diff <- count_better_cor_diff + 1
    if (!is.na(prp) && !is.na(mean_rank_percentile) && prp >= mean_rank_percentile) count_better_rank_perc <- count_better_rank_perc + 1
    if (!is.na(pvc) && !is.na(voxel_cor) && pvc >= voxel_cor) count_better_voxel_corr <- count_better_voxel_corr + 1
    if (!is.na(pmse) && !is.na(mse) && pmse <= mse) count_better_mse <- count_better_mse + 1 # Lower is better for MSE
    if (!is.na(prsq) && !is.na(r_squared) && prsq >= r_squared) count_better_r_squared <- count_better_r_squared + 1
    if (!is.na(pctm) && !is.na(cor_temporal_means) && pctm >= cor_temporal_means) count_better_cor_temp_means <- count_better_cor_temp_means + 1 # New
    if (!is.na(pmvtc) && !is.na(mean_voxelwise_temporal_cor) && pmvtc >= mean_voxelwise_temporal_cor) count_better_mean_vox_temp_cor <- count_better_mean_vox_temp_cor + 1 # New

    # --- Sums for z-scores --- 
    if (!is.na(pmc)) { sum_mean_corr    <- sum_mean_corr + pmc; sum_sq_mean_corr <- sum_sq_mean_corr + pmc^2 }
    if (!is.na(pcd)) { sum_cor_diff    <- sum_cor_diff + pcd; sum_sq_cor_diff <- sum_sq_cor_diff + pcd^2 }
    if (!is.na(prp)) { sum_rank_perc    <- sum_rank_perc + prp; sum_sq_rank_perc <- sum_sq_rank_perc + prp^2 }
    if (!is.na(pvc)) { sum_voxel_corr    <- sum_voxel_corr + pvc; sum_sq_voxel_corr <- sum_sq_voxel_corr + pvc^2 }
    if (!is.na(pmse)) { sum_mse    <- sum_mse + pmse; sum_sq_mse <- sum_sq_mse + pmse^2 }
    if (!is.na(prsq)) { sum_r_squared    <- sum_r_squared + prsq; sum_sq_r_squared <- sum_sq_r_squared + prsq^2 }
    if (!is.na(pctm)) { sum_cor_temp_means <- sum_cor_temp_means + pctm; sum_sq_cor_temp_means <- sum_sq_cor_temp_means + pctm^2 } # New
    if (!is.na(pmvtc)) { sum_mean_vox_temp_cor <- sum_mean_vox_temp_cor + pmvtc; sum_sq_mean_vox_temp_cor <- sum_sq_mean_vox_temp_cor + pmvtc^2 } # New

    # Possibly store the full permutation distribution
    if (save_distributions) {
      perm_mean_corr[i]    <- pmc
      perm_cor_diff[i]     <- pcd
      perm_rank_perc[i]    <- prp
      perm_voxel_corr[i]   <- pvc
      perm_mse[i]          <- pmse
      perm_r_squared[i]    <- prsq
      perm_correlations[i,] <- perm_cors
      perm_cor_temp_means[i] <- pctm # New
      perm_mean_vox_temp_cor[i] <- pmvtc # New
    }
  }

  # --- Compute p-values --- 
  # Calculate effective N for each metric (number of permutations where metric was computable)
  # Needs stored distributions if save_distributions=TRUE
  n_eff_list <- list()
  if (save_distributions) {
     n_eff_list$mean_corr <- sum(!is.na(perm_mean_corr))
     n_eff_list$cor_diff <- sum(!is.na(perm_cor_diff))
     n_eff_list$rank_perc <- sum(!is.na(perm_rank_perc))
     n_eff_list$voxel_corr <- sum(!is.na(perm_voxel_corr))
     n_eff_list$mse <- sum(!is.na(perm_mse))
     n_eff_list$r_squared <- sum(!is.na(perm_r_squared))
     n_eff_list$cor_temp_means <- sum(!is.na(perm_cor_temp_means))
     n_eff_list$mean_vox_temp_cor <- sum(!is.na(perm_mean_vox_temp_cor))
  } else {
     # If not saving distributions, estimate N_eff based on initial valid check (less precise)
     n_eff_list <- lapply(list(mean_cor, cor_difference, mean_rank_percentile, voxel_cor, mse, r_squared, cor_temporal_means, mean_voxelwise_temporal_cor), function(x) if(is.na(x)) 0 else nperm)
     names(n_eff_list) <- c("mean_corr", "cor_diff", "rank_perc", "voxel_corr", "mse", "r_squared", "cor_temp_means", "mean_vox_temp_cor")
  }
  
  # Helper for p-value calculation
  safe_p <- function(count_better, n_eff) {
       if (n_eff > 0) (count_better + 1) / (n_eff + 1) else NA_real_
  }

  p_mean_corr  <- safe_p(count_better_mean_corr, n_eff_list$mean_corr)
  p_cor_diff   <- safe_p(count_better_cor_diff, n_eff_list$cor_diff)
  p_rank_perc  <- safe_p(count_better_rank_perc, n_eff_list$rank_perc)
  p_voxel_corr <- safe_p(count_better_voxel_corr, n_eff_list$voxel_corr)
  p_mse        <- safe_p(count_better_mse, n_eff_list$mse)
  p_r_squared  <- safe_p(count_better_r_squared, n_eff_list$r_squared)
  p_cor_temp_means <- safe_p(count_better_cor_temp_means, n_eff_list$cor_temp_means)
  p_mean_vox_temp_cor <- safe_p(count_better_mean_vox_temp_cor, n_eff_list$mean_vox_temp_cor)

  # --- Compute means and SDs of permutation distributions --- 
  safe_mean_sd <- function(sum_val, sum_sq_val, n_eff) {
      if (n_eff > 0) {
          mean_perm = sum_val / n_eff
          var_perm = max(0, (sum_sq_val / n_eff) - mean_perm^2)
          sd_perm = sqrt(var_perm)
      } else {
          mean_perm = NA_real_
          sd_perm = NA_real_
      }
      list(mean = mean_perm, sd = sd_perm)
  }
  
  stats_mean_corr <- safe_mean_sd(sum_mean_corr, sum_sq_mean_corr, n_eff_list$mean_corr)
  stats_cor_diff <- safe_mean_sd(sum_cor_diff, sum_sq_cor_diff, n_eff_list$cor_diff)
  stats_rank_perc <- safe_mean_sd(sum_rank_perc, sum_sq_rank_perc, n_eff_list$rank_perc)
  stats_voxel_corr <- safe_mean_sd(sum_voxel_corr, sum_sq_voxel_corr, n_eff_list$voxel_corr)
  stats_mse <- safe_mean_sd(sum_mse, sum_sq_mse, n_eff_list$mse)
  stats_r_squared <- safe_mean_sd(sum_r_squared, sum_sq_r_squared, n_eff_list$r_squared)
  stats_cor_temp_means <- safe_mean_sd(sum_cor_temp_means, sum_sq_cor_temp_means, n_eff_list$cor_temp_means) # New
  stats_mean_vox_temp_cor <- safe_mean_sd(sum_mean_vox_temp_cor, sum_sq_mean_vox_temp_cor, n_eff_list$mean_vox_temp_cor) # New
  
  mean_perm_mean_corr <- stats_mean_corr$mean; sd_perm_mean_corr <- stats_mean_corr$sd
  mean_perm_cor_diff <- stats_cor_diff$mean; sd_perm_cor_diff <- stats_cor_diff$sd
  mean_perm_rank_perc <- stats_rank_perc$mean; sd_perm_rank_perc <- stats_rank_perc$sd
  mean_perm_voxel_corr <- stats_voxel_corr$mean; sd_perm_voxel_corr <- stats_voxel_corr$sd
  mean_perm_mse <- stats_mse$mean; sd_perm_mse <- stats_mse$sd
  mean_perm_r_squared <- stats_r_squared$mean; sd_perm_r_squared <- stats_r_squared$sd
  mean_perm_cor_temp_means <- stats_cor_temp_means$mean; sd_perm_cor_temp_means <- stats_cor_temp_means$sd # New
  mean_perm_mean_vox_temp_cor <- stats_mean_vox_temp_cor$mean; sd_perm_mean_vox_temp_cor <- stats_mean_vox_temp_cor$sd # New
  
  # --- z-scores --- 
  eps <- .Machine$double.eps
  safe_z <- function(observed_val, mean_perm, sd_perm, lower_is_better=FALSE) {
      if (is.na(observed_val) || is.na(mean_perm) || is.na(sd_perm)) return(NA_real_)
      # Ensure sd is not effectively zero
      sd_use <- max(sd_perm, eps)
      if (lower_is_better) {
         (mean_perm - observed_val) / sd_use # e.g., for MSE
      } else {
         (observed_val - mean_perm) / sd_use # For correlations, R^2 etc.
      }
  }
  
  z_mean_corr  <- safe_z(mean_cor, mean_perm_mean_corr, sd_perm_mean_corr)
  z_cor_diff   <- safe_z(cor_difference, mean_perm_cor_diff, sd_perm_cor_diff)
  z_rank_perc  <- safe_z(mean_rank_percentile, mean_perm_rank_perc, sd_perm_rank_perc)
  z_voxel_corr <- safe_z(voxel_cor, mean_perm_voxel_corr, sd_perm_voxel_corr)
  z_mse        <- safe_z(mse, mean_perm_mse, sd_perm_mse, lower_is_better=TRUE) 
  z_r_squared  <- safe_z(r_squared, mean_perm_r_squared, sd_perm_r_squared)
  z_cor_temp_means <- safe_z(cor_temporal_means, mean_perm_cor_temp_means, sd_perm_cor_temp_means) # New
  z_mean_vox_temp_cor <- safe_z(mean_voxelwise_temporal_cor, mean_perm_mean_vox_temp_cor, sd_perm_mean_vox_temp_cor) # New

  out <- list(
    p_values = c(mean_correlation = p_mean_corr,
                 cor_difference = p_cor_diff,
                 mean_rank_percentile = p_rank_perc,
                 voxel_correlation= p_voxel_corr,
                 mse = p_mse,
                 r_squared = p_r_squared,
                 cor_temporal_means = p_cor_temp_means,
                 mean_voxelwise_temporal_cor = p_mean_vox_temp_cor),
    z_scores = c(mean_correlation = z_mean_corr,
                 cor_difference = z_cor_diff,
                 mean_rank_percentile = z_rank_perc,
                 voxel_correlation= z_voxel_corr,
                 mse = z_mse,
                 r_squared = z_r_squared,
                 cor_temporal_means = z_cor_temp_means,
                 mean_voxelwise_temporal_cor = z_mean_vox_temp_cor)
  )

  if (save_distributions) {
    out$permutation_distributions <- list(
      mean_correlation  = perm_mean_corr,
      cor_difference    = perm_cor_diff,
      mean_rank_percentile = perm_rank_perc,
      voxel_correlation = perm_voxel_corr,
      mse               = perm_mse,
      r_squared         = perm_r_squared,
      correlations      = perm_correlations,
      cor_temporal_means = perm_cor_temp_means,
      mean_voxelwise_temporal_cor = perm_mean_vox_temp_cor
    )
  }
  
  out
}



#' Evaluate model performance for feature RSA
#'
#' Computes correlation-based metrics (diag correlation, mean correlation, voxel correlation),
#' MSE, R^2, and optionally performs permutation tests (via a helper function).
#'
#' @param object The feature RSA model
#' @param predicted Matrix of predicted values (from feature space F to voxel space X)
#' @param observed Matrix of observed values (actual voxel space X)
#' @param nperm Number of permutations for statistical testing (default: 0, no permutation)
#' @param save_distributions Logical indicating whether to save full permutation distributions
#' @param ... Additional arguments
#'
#' @return A list containing performance metrics and optional permutation results
#' @export
evaluate_model.feature_rsa_model <- function(object,
                                             predicted,
                                             observed,
                                             nperm = 0,
                                             save_distributions = FALSE,
                                             ...) 
{
  observed  <- as.matrix(observed)
  predicted <- as.matrix(predicted)

  
  
  # Check for constant predictions (zero variance) which cause issues
  if (any(apply(predicted, 2, stats::sd) == 0) || any(apply(observed, 2, stats::sd) == 0)) {
    warning("evaluate_model: Predictions or observed data have zero variance in some columns. Correlation metrics may be NA.")
  }
  
  if (ncol(observed) != ncol(predicted)) {
    stop(sprintf("Mismatch in columns: predicted has %d, observed has %d.", 
                 ncol(predicted), ncol(observed)))
  }
  
  # Base RSA metrics
  cormat     <- cor(predicted, observed)
  cors       <- diag(cormat)
  mean_cor   <- mean(cors, na.rm = TRUE) # Add na.rm = TRUE for robustness
  
  # Calculate mean of off-diagonal correlations
  n <- nrow(cormat)
  off_diag_cors <- (sum(cormat, na.rm = TRUE) - sum(cors, na.rm = TRUE)) / (n*n - n) # Add na.rm
  
  # New metric: mean diagonal correlation minus mean off-diagonal correlation
  # This measures how much better the model predicts the correct condition
  # compared to incorrect conditions.
  cor_difference <- mean_cor - off_diag_cors
  
  # Calculate rank percentile for each condition
  ranks <- numeric(n)
  for (i in 1:n) {
    # Get correlations for the ith predicted pattern with all observed patterns
    condition_cors <- cormat[i, ]
    # Compute percentile rank of diagonal correlation among all other correlations
    # Higher correlation = better rank; adjust to 0-1 scale excluding self-comparison
    ranks[i] <- (sum(condition_cors <= condition_cors[i], na.rm = TRUE) - 1) / (sum(!is.na(condition_cors)) - 1) # Handle NAs
  }
  mean_rank_percentile <- mean(ranks, na.rm = TRUE) # Add na.rm
  voxel_cor  <- cor(as.vector(predicted), as.vector(observed))
  mse        <- mean((predicted - observed)^2, na.rm=TRUE) # Add na.rm
  rss        <- sum((observed - predicted)^2, na.rm=TRUE)  # Add na.rm
  tss        <- sum((observed - mean(observed, na.rm=TRUE))^2, na.rm=TRUE) # Add na.rm
  r_squared  <- if (tss == 0) NA else 1 - (rss / tss) # Handle zero total sum of squares
  
  # Log warnings if key metrics are NA/NaN
  if (!is.finite(mean_cor)) futile.logger::flog.warn("Mean correlation is NA/NaN.")
  if (!is.finite(cor_difference)) futile.logger::flog.warn("Correlation difference is NA/NaN.")
  if (!is.finite(mean_rank_percentile)) futile.logger::flog.warn("Mean rank percentile is NA/NaN.")
  if (!is.finite(voxel_cor)) futile.logger::flog.warn("Voxel correlation is NA/NaN.")
  if (!is.finite(r_squared)) futile.logger::flog.warn("R-squared is NA/NaN.")
  
  # --- Calculate Correlation of Temporal Means (Spatial Averages) ---
  cor_temporal_means <- NA_real_
  if (nrow(observed) > 1 && nrow(predicted) > 1) { # Need >1 observation
      mean_obs_across_space <- tryCatch(rowMeans(observed, na.rm = TRUE), error=function(e) NULL)
      mean_pred_across_space <- tryCatch(rowMeans(predicted, na.rm = TRUE), error=function(e) NULL)
      if (!is.null(mean_obs_across_space) && !is.null(mean_pred_across_space)) {
          cor_temporal_means <- tryCatch(cor(mean_obs_across_space, mean_pred_across_space), error=function(e) NA_real_)
      } else {
         warning("evaluate_model: Could not compute rowMeans for cor_temporal_means.")
      }
  } else {
       warning("evaluate_model: Cannot calculate cor_temporal_means with < 2 observations.")
  }
  if (!is.finite(cor_temporal_means)) cor_temporal_means <- NA_real_ # Ensure NA if calculation failed
  # ------------------------------------------------------------------

  # --- Calculate Mean Voxelwise Temporal Correlation ---
  mean_voxelwise_temporal_cor <- NA_real_
  if (nrow(observed) > 1 && nrow(predicted) > 1 && ncol(observed) > 0) { # Need >1 observation and >0 voxels
      num_voxels <- ncol(observed)
      voxel_cors <- numeric(num_voxels)
      for (i in 1:num_voxels) {
          voxel_cors[i] <- tryCatch(cor(observed[, i], predicted[, i]), error = function(e) NA)
      }
      mean_voxelwise_temporal_cor <- mean(voxel_cors, na.rm = TRUE)
  } else {
       warning("evaluate_model: Cannot calculate mean_voxelwise_temporal_cor with < 2 observations or 0 voxels.")
  }
   if (!is.finite(mean_voxelwise_temporal_cor)) mean_voxelwise_temporal_cor <- NA_real_ # Ensure NA if calculation failed
  # ------------------------------------------------------


  perm_results <- NULL
  # Placeholder for incremental correlation - calculation requires comparing
  # results from models trained on feature subsets, which is not possible
  # within the evaluation of a single model run.
  # incremental_correlation <- NA_real_  # <<< REMOVE THIS COMPLETELY

  if (nperm > 0) {
    perm_results <- .perm_test_feature_rsa(
      observed = observed,
      predicted = predicted,
      nperm = nperm,
      save_distributions = save_distributions,
      mean_cor = mean_cor,
      cor_difference = cor_difference,
      mean_rank_percentile = mean_rank_percentile,
      voxel_cor = voxel_cor,
      mse = mse,
      r_squared = r_squared,
      # incremental_correlation = incremental_correlation, # Pass the calculated value # <<< REMOVE THIS
      cor_temporal_means = cor_temporal_means, # Pass new metric 1
      mean_voxelwise_temporal_cor = mean_voxelwise_temporal_cor, # Pass new metric 2
      cors = cors
    )
  }
  
  list(
    correlations        = cors,
    mean_correlation    = mean_cor,
    off_diag_correlation= off_diag_cors,
    cor_difference      = cor_difference,
    mean_rank_percentile = mean_rank_percentile,
    voxel_correlation   = voxel_cor,
    mse                 = mse,
    r_squared           = r_squared,
    cor_temporal_means = cor_temporal_means, # Add new metric 1
    mean_voxelwise_temporal_cor = mean_voxelwise_temporal_cor, # Add new metric 2
    permutation_results = perm_results
  )
}


#' @export
train_model.feature_rsa_model <- function(obj, X, y, indices, ...) {
  
  # X: brain data (samples x voxels)
  # y: should be the Feature Matrix F (samples x features)
  Fsub <- y
  
  result <- list(method=obj$method, design=obj$design)
  
  # Check for minimum data size
  if (nrow(X) < 2 || ncol(X) < 1 || nrow(Fsub) < 2 || ncol(Fsub) < 1) {
    error_msg <- sprintf("Insufficient data for training (X dims: %d x %d, F dims: %d x %d). Requires at least 2 samples and 1 voxel/feature.", 
                         nrow(X), ncol(X), nrow(Fsub), ncol(Fsub))
    futile.logger::flog.error(error_msg)
    result$error <- error_msg
    return(result)
  }
  
  # ---- PLS Train ----
  if (obj$method == "pls") {
    pls_res <- tryCatch({
      # Check for near-zero variance *before* standardization attempt
      var_X <- apply(X, 2, var, na.rm = TRUE)
      var_F <- apply(Fsub, 2, var, na.rm = TRUE)
      if (any(var_X < 1e-10) || any(var_F < 1e-10)) {
        stop("Near zero variance detected in X or F before standardization.")
      }
      
      sx <- .standardize(X)
      sf <- .standardize(Fsub)
      
      # Check variance *after* standardization (shouldn't happen if .standardize handles sd=0, but double-check)
      if (any(sx$sd < 1e-10) || any(sf$sd < 1e-10)) {
         stop("Near zero variance detected after standardization.")
      }
      
      max_k_possible <- min(nrow(sf$X_sc) - 1, ncol(sf$X_sc))
      k <- min(obj$max_comps, max_k_possible)
      
      if (k < 1) {
        stop(sprintf("Calculated number of PLS components (%d) is less than 1 (max_comps: %d, max_possible: %d).", 
                     k, obj$max_comps, max_k_possible))
      }
      
      # Fit PLS model
      pls::plsr(sx$X_sc ~ sf$X_sc, ncomp = k, scale = FALSE, validation = "none")
      
    }, error = function(e) {
      error_msg <- sprintf("train_model (PLS): Error during training - %s", e$message)
      futile.logger::flog.error(error_msg)
      # Return an error indicator object
      list(error = error_msg) 
    })

    # Check if tryCatch returned an error object
    if (!is.null(pls_res$error)) {
       result$error <- pls_res$error
       return(result)
    }

    # Store necessary results ONLY if successful
    result$trained_model <- pls_res
    result$pls_x_mean    <- sx$mean # Need sx/sf from the try block scope
    result$pls_x_sd      <- sx$sd
    result$pls_f_mean    <- sf$mean
    result$pls_f_sd      <- sf$sd
    result$ncomp         <- k       # Need k from the try block scope

  } else if (obj$method == "scca") {
    # ---- SCCA Train ----
    scca_res <- tryCatch({
       sx <- .standardize(X)
       sf <- .standardize(Fsub)
       if (any(sx$sd < 1e-8) || any(sf$sd < 1e-8)) {
          stop("Zero variance detected after standardization.")
       }
       # Store standardization details temporarily within this scope
       list(res=whitening::scca(sx$X_sc, sf$X_sc, scale=FALSE), sx=sx, sf=sf)
    }, error = function(e) {
      error_msg <- sprintf("train_model (SCCA): SCCA execution error - %s", e$message)
      futile.logger::flog.error(error_msg)
      list(error = error_msg)
    })

    # Check if tryCatch returned an error object
    if (!is.null(scca_res$error)) {
      result$error <- scca_res$error
      return(result)
    }
    
    # Extract results if successful
    scca_fit <- scca_res$res
    sx <- scca_res$sx
    sf <- scca_res$sf
    
    effective_ncomp <- if (!is.null(scca_fit$lambda)) sum(abs(scca_fit$lambda) > 1e-6) else 0
    ncomp <- min(effective_ncomp, obj$max_comps)
    
    # Store standardization info regardless (needed for potential fallback in predict)
    result$scca_x_mean <- sx$mean
    result$scca_x_sd   <- sx$sd
    result$scca_f_mean <- sf$mean
    result$scca_f_sd   <- sf$sd
    result$ncomp       <- ncomp 
        
    if (ncomp < 1) {
      futile.logger::flog.info("train_model (SCCA): No effective canonical components found (effective: %d, max_comps: %d). Prediction will use mean fallback.", 
                               effective_ncomp, obj$max_comps)
      # Still store the (potentially empty) scca_fit
      result$trained_model <- scca_fit
    } else {
      result$trained_model <- scca_fit
    }
    
  } else if (obj$method == "pca") {
    #browser()
    #
    # -- PCA with possible caching --
    #
    pca_result <- tryCatch({
        key <- .hash_row_indices(indices)
        cached_pca <- NULL
        if (isTRUE(obj$cache_pca) && !is.null(obj$pca_cache[[key]])) {
          cached_pca <- obj$pca_cache[[key]]
          sx <- .standardize(X) # Still standardize X for this ROI
          pca_res <- cached_pca$pca_res
          sf_mean <- cached_pca$f_mean
          sf_sd   <- cached_pca$f_sd
        } else {
          sx <- .standardize(X)
          sf <- .standardize(Fsub)
          if (any(sf$sd < 1e-10)) { # Check variance before prcomp
             stop("Zero variance detected in features (F) before PCA.")
          }
          pca_res <- prcomp(sf$X_sc, scale.=FALSE)
          sf_mean <- sf$mean
          sf_sd   <- sf$sd
          # Store if caching is enabled
          if (isTRUE(obj$cache_pca)) {
            obj$pca_cache[[key]] <- list(
              pca_res = pca_res, f_mean = sf_mean, f_sd = sf_sd
            )
          }
        }
        
        available_k <- ncol(pca_res$x)
        k <- min(obj$max_comps, available_k)
        if (k < 1) {
          stop(sprintf("No principal components available (k=%d) after PCA (max_comps: %d, available: %d).", 
                      k, obj$max_comps, available_k))
        }
        PC_train_subset <- pca_res$x[, seq_len(k), drop=FALSE]
        
        df_pcs <- as.data.frame(PC_train_subset)
        if (nrow(df_pcs) <= k) {
          stop(sprintf("Insufficient data for PCA regression (samples: %d, components: %d). Need more samples than components.", 
                       nrow(df_pcs), k))
        }
        
        fit <- lm(sx$X_sc ~ ., data=df_pcs)
        
        # Return all necessary components from successful run
        list(
          pcarot     = pca_res$rotation[, seq_len(k), drop=FALSE],
          pca_f_mean = sf_mean,
          pca_f_sd   = sf_sd,
          pca_coefs  = coef(fit),
          pca_x_mean = sx$mean,
          pca_x_sd   = sx$sd,
          ncomp      = k,
          trained_model = fit # Add the lm fit object here
        )
        
    }, error = function(e) {
        error_msg <- sprintf("train_model (PCA): Error during training - %s", e$message)
        futile.logger::flog.error(error_msg)
        list(error = error_msg)
    })

    # Check if tryCatch returned an error object
    if (!is.null(pca_result$error)) {
        result$error <- pca_result$error
        return(result)
    }
    
    # Assign results if successful
    result <- c(result, pca_result) # Merge the list of results
    
  } else if (obj$method == "glmnet") {
    #
    # ---- GLMNet Train ----
    #
    glm_result <- tryCatch({
        # Standardize X and F
        sx <- .standardize(X)
        sf <- .standardize(Fsub)
        
        if (any(sx$sd < 1e-10) || any(sf$sd < 1e-10)) { # Check variance
             stop("Zero variance detected in X or F after standardization.")
        }

        if (nrow(sx$X_sc) < 2 || nrow(sf$X_sc) < 2) {
          stop(sprintf("Insufficient observations for GLMNet (X: %d, F: %d). Requires >= 2.", nrow(X), nrow(Fsub)))
        }
        
        lambda_to_use <- obj$lambda
        cv_results <- NULL # Placeholder for CV output
        cv_error <- NULL # Placeholder for CV specific error
        
        # Determine if CV should run
        run_cv <- isTRUE(obj$cv_glmnet)
        
        if (run_cv) {
          n_obs <- nrow(sf$X_sc)
          if (n_obs < 3) { # cv.glmnet default nfolds=10 requires >=3
              futile.logger::flog.warn("train_model (GLMNet CV): Too few observations (%d) for reliable CV. Skipping CV.", n_obs)
              run_cv <- FALSE
          } else {
             foldid <- tryCatch({
                 # Use default k-fold (typically 10), let cv.glmnet handle if n_obs < nfolds
                 # Using internal cv.glmnet fold generation might be more robust
                 NULL 
             }, error = function(e) {
                 futile.logger::flog.warn("train_model (GLMNet CV): Error creating fold IDs - %s. Skipping CV.", e$message)
                 run_cv <<- FALSE # Modify run_cv in the outer scope
                 NULL
             })
             
             if (run_cv) { # Check again if foldid creation failed
                cv_fit <- tryCatch({
                    glmnet::cv.glmnet(
                      x = sf$X_sc, 
                      y = sx$X_sc,
                      family = "mgaussian",
                      alpha = obj$alpha,
                      lambda = obj$lambda, # Pass user lambda if specified
                      foldid = foldid,    # Pass NULL to let cv.glmnet create folds
                      standardize = FALSE,
                      intercept = TRUE
                    )
                }, error = function(e) {
                    cv_error <<- sprintf("cv.glmnet failed: %s", e$message) # Assign to outer scope
                    futile.logger::flog.warn("train_model (GLMNet CV): %s. Fitting with standard glmnet instead.", cv_error)
                    run_cv <<- FALSE # Modify run_cv in the outer scope
                    NULL # Return NULL to indicate CV failure
                })
                
                if (run_cv && !is.null(cv_fit)) { # If CV succeeded
                    lambda_to_use <- cv_fit$lambda.min
                    cv_results <- cv_fit # Store CV results
                }
             }
          }
        }
        
        # Fit standard glmnet (either as fallback or primary)
        final_fit <- glmnet::glmnet(
          x = sf$X_sc, 
          y = sx$X_sc,
          family = "mgaussian",
          alpha = obj$alpha,
          lambda = lambda_to_use, # Use CV lambda if available, otherwise obj$lambda
          standardize = FALSE,
          intercept = TRUE
        )
        
        # Determine lambda used for prediction
        lambda_used_for_pred <- if (run_cv && !is.null(cv_results)) {
           lambda_to_use # lambda.min from successful CV
        } else if (!is.null(final_fit$lambda)) {
           final_fit$lambda[1] # First lambda if multiple were fit (e.g., obj$lambda=NULL)
        } else {
           NA # Should not happen if fit succeeded
        }

        # Calculate ncomp proxy
        ncomp_proxy <- NA
        if (!is.null(final_fit) && !is.null(lambda_used_for_pred) && is.finite(lambda_used_for_pred)) {
           coefs <- tryCatch(glmnet::coef.glmnet(final_fit, s = lambda_used_for_pred), error=function(e) NULL)
           if (!is.null(coefs) && is.list(coefs)) { # mgaussian returns a list
              nonzero_count <- sapply(coefs, function(cm) sum(as.matrix(cm[-1,]) != 0)) # Exclude intercept
              ncomp_proxy <- round(mean(nonzero_count))
           } else {
              futile.logger::flog.warn("train_model (GLMNet): Could not extract coefficients to calculate ncomp proxy.")
           }
        } else {
           futile.logger::flog.warn("train_model (GLMNet): Could not determine lambda used or fit failed; cannot calculate ncomp proxy.")
        }

        # Return results
        list(
          trained_model = final_fit,
          glmnet_x_mean = sx$mean,
          glmnet_x_sd   = sx$sd,
          glmnet_f_mean = sf$mean,
          glmnet_f_sd   = sf$sd,
          cv_glmnet     = (run_cv && !is.null(cv_results)), # True only if CV ran *and* succeeded
          cv_results    = cv_results, # Store CV object if it succeeded
          cv_error      = cv_error,   # Store CV error message if it occurred
          lambda_used   = lambda_used_for_pred,
          ncomp         = ncomp_proxy
        )
        
    }, error = function(e) {
        # Catch errors from standardization or the final glmnet fit
        error_msg <- sprintf("train_model (GLMNet): Error during training - %s", e$message)
        futile.logger::flog.error(error_msg)
        list(error = error_msg)
    })
    
    # Check if tryCatch returned an error object
    if (!is.null(glm_result$error)) {
        result$error <- glm_result$error
        return(result)
    }
    
    # Log CV error if it occurred but didn't stop the process
    if (!is.null(glm_result$cv_error)) {
       # This was already logged as warning, but good to have in final result list too?
       # Maybe add it to the result list itself
       result$cv_warning <- glm_result$cv_error
    }
    
    # Assign results if successful
    result <- c(result, glm_result)

  } else {
    # This case should ideally not be reached if method is matched earlier
    error_msg <- paste("Unknown method in train_model.feature_rsa_model:", obj$method)
    futile.logger::flog.error(error_msg)
    result$error <- error_msg
    return(result)
  }

  
  # Check for NULL trained_model just in case
  if (is.null(result$trained_model)) {
     error_msg <- sprintf("train_model (%s): Training finished but 'trained_model' is NULL. This indicates an unexpected issue.", obj$method)
     futile.logger::flog.error(error_msg)
     result$error <- error_msg
     # Ensure ncomp is NA if model is NULL
     if (!"ncomp" %in% names(result)) result$ncomp <- NA 
  }
  
  # Ensure ncomp exists in the result list, set to NA if missing
  if (!"ncomp" %in% names(result)) {
      futile.logger::flog.warn("train_model (%s): 'ncomp' was not set during training. Setting to NA.", obj$method)
      result$ncomp <- NA_real_
  }
  
  return(result) # Return the final result list
}


#' @export
y_train.feature_rsa_model <- function(obj) {
  obj$design$F  # Features are used as predictors (y in training function)
}

#' @export
y_train.feature_rsa_design <- function(obj) {
  obj$F  # Features are used as predictors
}

#' @export
format_result.feature_rsa_model <- function(obj, result, error_message=NULL, context, ...) {
  
  if (!is.null(error_message)) {
    return(tibble::tibble(
      observed      = list(NULL), 
      predicted     = list(NULL),
      result        = list(NULL),
      performance   = list(NULL),
      error         = TRUE, 
      error_message = error_message
    ))
  }
  
  Xobs  <- as.data.frame(context$test)
  Ftest <- as.matrix(context$ytest)
  
  Xpred <- tryCatch({
    predict_model.feature_rsa_model(obj, result, Ftest)
  }, error=function(e) {
       # Log the specific error
       futile.logger::flog.warn("format_result: Prediction failed - %s", e$message)
       # Return a list indicating error and the message
       list(error = TRUE, message = paste("Prediction failed:", e$message))
    })
  
  # Check if the prediction step returned an error list
  if (is.list(Xpred) && !is.null(Xpred$error) && Xpred$error) {
    return(tibble::tibble(
      observed      = list(NULL), 
      predicted     = list(NULL),
      result        = list(NULL),
      performance   = list(NULL),
      error         = TRUE, 
      error_message = Xpred$message # Use the captured error message
    ))
  }
  
  # Check if Xpred is NULL for any other unexpected reason (shouldn't happen ideally)
  if (is.null(Xpred)) {
     return(tibble::tibble(
      observed      = list(NULL), 
      predicted     = list(NULL),
      result        = list(NULL),
      performance   = list(NULL),
      error         = TRUE, 
      error_message = "Prediction returned NULL unexpectedly."
    ))
  }
  
  # Evaluate WITHOUT permutations at the fold level
  perf <- evaluate_model.feature_rsa_model(
    object = obj,
    predicted = Xpred,
    observed  = Xobs,
    nperm = 0  # no permutation here
  )
  
  # Get ncomp from the first fold's result (assuming it's consistent)
  ncomp_used <- result$ncomp
  
  # Summarize
  perf_mat <- matrix(
    c(perf$mean_correlation,
      perf$cor_difference,
      perf$mean_rank_percentile,
      perf$voxel_correlation,
      perf$mse,
      perf$r_squared,
      perf$cor_temporal_means, # Add new metric 1
      perf$mean_voxelwise_temporal_cor, # Add new metric 2
      ncomp_used),
    nrow = 1,
    ncol = 9,
    dimnames = list(NULL, c("mean_correlation", "cor_difference", "mean_rank_percentile", "voxel_correlation", "mse", "r_squared", "cor_temporal_means", "mean_voxelwise_temporal_cor", "ncomp"))
  )
  
  tibble::tibble(
    observed    = list(Xobs),
    predicted   = list(Xpred),
    result      = list(result),
    performance = list(perf_mat),
    error       = FALSE,
    error_message = "~"
  )
}

#' @export
merge_results.feature_rsa_model <- function(obj, result_set, indices, id, ...) {
 
  if (any(result_set$error)) {
    emessage <- result_set$error_message[ which(result_set$error)[1] ]
    return(
      tibble::tibble(
        result       = list(NULL),
        indices      = list(indices),
        performance  = list(NULL),
        id           = id,
        error        = TRUE,
        error_message= emessage
      )
    )
  }
  
  observed_list  <- result_set$observed
  predicted_list <- result_set$predicted

  
  
  # Get the list of results from each fold (contains ncomp)
  fold_results_list <- result_set$result
  
  combined_observed  <- do.call(rbind, observed_list)
  combined_predicted <- do.call(rbind, predicted_list)
  
  # Extract ncomp from the first fold's result 
  # (result_set$result is a list, first element is result from fold 1)
  # Add robustness in case ncomp is NULL or NA
  ncomp_val <- fold_results_list[[1]]$ncomp
  ncomp_used <- if (!is.null(ncomp_val) && is.finite(ncomp_val)) {
    ncomp_val
  } else {
    NA # Use NA if ncomp wasn't properly recorded
  }
  
  # Now we do permutations (if nperm>0 in the model spec)
  perf <- evaluate_model.feature_rsa_model(
    object    = obj,
    predicted = combined_predicted,
    observed  = combined_observed,
    nperm     = obj$nperm,
    save_distributions = obj$save_distributions
  )
  
  # Collate results
  base_metrics <- c(
    perf$mean_correlation,
    perf$cor_difference,
    perf$mean_rank_percentile,
    perf$voxel_correlation,
    perf$mse,
    perf$r_squared,
    # perf$incremental_correlation, # REMOVED
    perf$cor_temporal_means, # Add new metric 1
    perf$mean_voxelwise_temporal_cor, # Add new metric 2
    ncomp_used
  )
  base_names <- c(
    "mean_correlation", "cor_difference", "mean_rank_percentile", 
    "voxel_correlation", "mse", "r_squared", 
    # "incremental_correlation", # REMOVED
    "cor_temporal_means", "mean_voxelwise_temporal_cor", # Add names here
    "ncomp"
  )

  if (is.null(perf$permutation_results)) {
      perf_values <- base_metrics
      perf_names <- base_names
  } else {
      perm_p_values <- perf$permutation_results$p_values
      perm_z_scores <- perf$permutation_results$z_scores
      
      # Ensure order matches p_values/z_scores structure in .perm_test_feature_rsa
      # Dynamically get names to be robust
      p_names <- paste0("p_", names(perm_p_values))
      z_names <- paste0("z_", names(perm_z_scores))

      perf_values <- c(base_metrics, perm_p_values, perm_z_scores)
      perf_names <- c(base_names, p_names, z_names)
  }
  
  perf_mat <- matrix(
      perf_values,
      nrow = 1,
      ncol = length(perf_values),
      dimnames = list(NULL, perf_names)
  )
  
  # Remove any potential columns that are all NA (handles case where incremental_corr placeholders might be NA)
  # Though they are explicitly set to NA_real_, this adds robustness
  perf_mat <- perf_mat[, colSums(is.na(perf_mat)) < nrow(perf_mat), drop = FALSE]

  tibble::tibble(
    result      = list(NULL),
    indices     = list(indices),
    performance = list(perf_mat),
    id          = id,
    error       = FALSE,
    error_message = "~"
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
  
  # Display the single max_comps limit stored in the design
  # This is the upper limit for components derived from the feature space (F)
  # for *any* method used in feature_rsa_model.
  cat(crayon::bold(crayon::blue("Max Components Limit:   ")),
      if(!is.null(x$max_comps)) x$max_comps else "Not explicitly set (using default)", "\n")

  # Indicate if a blocking variable was supplied
  if (!is.null(x$block_var)) {
    cat(crayon::bold(crayon::blue("Blocking Variable:      ")), "Provided (", length(unique(x$block_var)), " blocks)\n")
  } else {
    cat(crayon::bold(crayon::blue("Blocking Variable:      ")), "None\n")
  }
  
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
  
  # Display component limit
  comp_limit <- if (!is.null(x$max_comps)) {
    x$max_comps
  } else if (!is.null(x$design$max_comps)) {
    x$design$max_comps
  } else {
    "Default"
  }
  cat(crayon::bold(crayon::blue("Max components limit:   ")), comp_limit, "\n")

  if (x$method == "pca" && !is.null(x$precomputed_pca)) {
    cat(crayon::bold(crayon::green("PCA Optimization:       ")),
        "PCA precomputed for efficiency\n")
  } else if (x$method == "glmnet") {
    cat(crayon::bold(crayon::blue("Elastic Net alpha:      ")), x$alpha, "\n")
    cat(crayon::bold(crayon::blue("Cross-validate lambda:  ")), 
        ifelse(isTRUE(x$cv_glmnet), "Yes", "No"), "\n")
    
    if (!is.null(x$lambda)) {
      lambda_str <- if (length(x$lambda) > 3) {
        paste0(paste(x$lambda[1:3], collapse=", "), ", ...")
      } else {
        paste(x$lambda, collapse=", ")
      }
      cat(crayon::bold(crayon::blue("Lambda values:         ")), lambda_str, "\n")
    }
  }
  
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


