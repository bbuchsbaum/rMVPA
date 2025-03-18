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
#' @param max_comps Maximum number of components to use (for PCA/PLS/SCCA). Default 10.
#'
#' @return A \code{feature_rsa_design} object (S3 class) containing:
#'   \describe{
#'     \item{S}{The input similarity matrix (if used)}
#'     \item{F}{Feature space projection matrix}
#'     \item{labels}{Vector of observation labels}
#'     \item{k}{The number of feature dimensions}
#'     \item{max_comps}{Maximum number of components for any method}
#'   }
#'
#' @details
#' If F is supplied directly, it is used as-is.
#' If only S is supplied, the design computes an eigen decomposition of S 
#' (keeping factors with eigenvalue > 1, unless you specify \code{k}). 
#'
#' @export
feature_rsa_design <- function(S=NULL, F=NULL, labels, k=0, max_comps=10) {
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
      max_comps=max_comps
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
      max_comps=max_comps
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
#' @param design A \code{feature_rsa_design} object specifying the feature space (\code{F}).
#' @param method Character string specifying the analysis method. One of:
#'   \describe{
#'     \item{scca}{Sparse Canonical Correlation Analysis}
#'     \item{pls}{Partial Least Squares}
#'     \item{pca}{Principal Component Analysis on feature space}
#'   }
#' @param crossval Optional cross-validation specification.
#'
#' @return A \code{feature_rsa_model} object (S3 class).
#'
#' @details
#' Feature RSA models analyze how well a feature matrix \code{F} can predict
#' neural data \code{X}. For example:
#'   - \strong{pca}: PCA on \code{F}, then regress \code{X} on those principal components.
#'   - \strong{pls}: partial least squares regression of \code{X} on \code{F}.
#'   - \strong{scca}: sparse CCA of \code{X} and \code{F}.
#'
#' @export
feature_rsa_model <- function(dataset,
                               design,
                               method = c("scca", "pls", "pca"),
                               crossval = NULL,
                               cache_pca = FALSE,
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
  
  if (is.null(crossval) && !is.null(design$block_var)) {
    crossval <- blocked_cross_validation(design$block_var)
  }
  assertthat::assert_that(!is.null(crossval))
  
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
  # SCCA prediction now predicts X from F (reversed from previous implementation)
  # First check dimensions and ensure compatibility
  F_new <- as.matrix(F_new)
  n_features_new <- ncol(F_new)
  n_features_expected <- length(model$scca_f_mean)
  tm <- model$trained_model
  ncomp <- model$ncomp
  
  # Check if dimensions match between F_new and the stored means/SDs
  if (n_features_new != n_features_expected) {
    # Dimensions don't match, log warning
    warning(sprintf("Feature matrix dimension mismatch for SCCA: expected %d features but got %d. Using only the common features.",
                   n_features_expected, n_features_new))
    
    # Special case: If we have a single feature but model expects more
    if (n_features_new == 1 && n_features_expected > 1) {
      # For this case, we need to create a synthetic feature matrix
      # by duplicating the single feature to match expected dimensions
      F_expanded <- matrix(rep(F_new, n_features_expected),
                          nrow = nrow(F_new),
                          ncol = n_features_expected)
      
      # Standardize the expanded features
      Fsc <- sweep(sweep(F_expanded, 2, model$scca_f_mean, "-"), 2, model$scca_f_sd, "/")
      
      # Use canonical directions from F to predict X
      fcoef <- t(tm$WY)[, 1:ncomp, drop=FALSE]
      xcoef <- t(tm$WX)[, 1:ncomp, drop=FALSE]
      Lambda_q <- diag(tm$lambda[1:ncomp])
      x_inv <- corpcor::pseudoinverse(xcoef)
      
      # Predict X from F using canonical correlations
      Xhat <- Fsc %*% fcoef %*% Lambda_q %*% x_inv
      Xhat <- sweep(sweep(Xhat, 2, model$scca_x_sd, "*"), 2, model$scca_x_mean, "+")
      return(Xhat)
    }
    
    # General case: Use only common features
    common_features <- min(n_features_new, n_features_expected)
    
    # Adjust the dimensions of all relevant objects
    if (common_features < n_features_new) {
      # If F_new has more columns than expected, subset it
      F_new <- F_new[, 1:common_features, drop = FALSE]
    }
    
    # Ensure we only use the means and SDs for common features
    f_means_to_use <- model$scca_f_mean[1:common_features]
    f_sds_to_use <- model$scca_f_sd[1:common_features]
    
    # We need to create a projection matrix that maps from our reduced feature space
    # to the canonical components
    if (common_features < ncol(tm$WY)) {
      # Use only the weights corresponding to the common features
      WY_reduced <- tm$WY[1:common_features, , drop=FALSE]
      
      # Standardize features using the stored means and standard deviations
      Fsc <- sweep(sweep(F_new, 2, f_means_to_use, "-"), 2, f_sds_to_use, "/")
      
      # Create a projection directly to canonical space
      # This bypasses the full feature space
      canonical_proj <- Fsc %*% WY_reduced
      
      # Now use the canonical projections to predict X
      xcoef <- t(tm$WX)[, 1:ncomp, drop=FALSE]
      Lambda_q <- diag(tm$lambda[1:ncomp])
      x_inv <- corpcor::pseudoinverse(xcoef)
      
      # Predict X directly from canonical projections
      Xhat <- canonical_proj %*% Lambda_q %*% x_inv
      Xhat <- sweep(sweep(Xhat, 2, model$scca_x_sd, "*"), 2, model$scca_x_mean, "+")
      return(Xhat)
    }
  }
  
  # Standard case: dimensions match
  # Standardize features using the stored means and standard deviations
  Fsc <- sweep(sweep(F_new, 2, model$scca_f_mean, "-"), 2, model$scca_f_sd, "/")
  
  # Use canonical directions from F to predict X
  fcoef <- t(tm$WY)[, 1:ncomp, drop=FALSE]
  xcoef <- t(tm$WX)[, 1:ncomp, drop=FALSE]
  Lambda_q <- diag(tm$lambda[1:ncomp])
  x_inv <- corpcor::pseudoinverse(xcoef)
  
  # Predict X from F using canonical correlations
  Xhat <- Fsc %*% fcoef %*% Lambda_q %*% x_inv
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


#' Predict Method for Feature RSA Model
#'
#' @param object The trained feature RSA model (method, etc.)
#' @param fit The fitted model object from training_model()
#' @param newdata New feature data matrix (rows = test obs, cols = features)
#' @param ... Additional args
#' @return Predicted brain activity \code{X} values (matrix).
#' @export
#' @export
predict_model.feature_rsa_model <- function(object, fit, newdata, ...) {
  method <- object$method
  F_new  <- as.matrix(newdata)
  
  if (method == "pls") {
    # Wrap new feature matrix in a data frame with the SAME name 'Fsub'
    test_df <- data.frame(Fsub = I(newdata))
    # The 'I(...)' calls 'AsIs' wrapper so that 'newdata' remains matrix-like.
    # This matches how the training used data.frame(X, Fsub).
    
    # Now call predict() with the same formula model
    preds <- predict(
      fit$trained_model,
      newdata = test_df,
      ncomp   = fit$ncomp
    )
    
    
   
    # Typically a 3D array [nTest, nResponseCols, nComps].
    # If ncomp=1, it might be [nTest, nResponseCols].
    # So we 'drop' the extra dimension:
    preds <- drop(preds)
    
    return(preds)
    
  } else if (method == "scca") {
    #
    # ---- SCCA Predict ----
    #
    return(.predict_scca(fit, F_new))
    
  } else if (method == "pca") {
    #
    # ---- PCA Predict ----
    #
    return(.predict_pca(fit, F_new))
    
    } else {
    stop("Unknown method in feature_rsa_model.")
  }
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
#' @param spatial_specificity Scalar: the observed spatial specificity (mean_cor - off_diag_cor)
#' @param voxel_cor Scalar: the observed voxel correlation
#' @param mse Scalar: the observed MSE
#' @param r_squared Scalar: the observed R^2
#' @param cors Vector of diagonal correlations
#' @return A list with p-values, z-scores, and optionally a list of permutations
.perm_test_feature_rsa <- function(observed,
                                   predicted,
                                   nperm,
                                   save_distributions,
                                   mean_cor,
                                   spatial_specificity,
                                   voxel_cor,
                                   mse,
                                   r_squared,
                                   cors)
{
  message("Performing permutation tests with ", nperm, " permutations... (feature_rsa_model)")
  
  n_rows <- nrow(predicted)
  rss <- sum((observed - predicted)^2)
  tss <- sum((observed - mean(observed))^2)

  # Counters for how many permutations are "better" than the true model
  count_better_mean_corr  <- 0
  count_better_spatial_spec <- 0
  count_better_voxel_corr <- 0
  count_better_mse        <- 0
  count_better_r_squared  <- 0

  # For mean & SD (to calculate z-scores)
  sum_mean_corr    <- 0
  sum_sq_mean_corr <- 0
  sum_spatial_spec    <- 0
  sum_sq_spatial_spec <- 0
  sum_voxel_corr    <- 0
  sum_sq_voxel_corr <- 0
  sum_mse    <- 0
  sum_sq_mse <- 0
  sum_r_squared    <- 0
  sum_sq_r_squared <- 0

  # Optionally store entire distributions
  if (save_distributions) {
    perm_mean_corr     <- numeric(nperm)
    perm_spatial_spec  <- numeric(nperm)
    perm_voxel_corr    <- numeric(nperm)
    perm_mse           <- numeric(nperm)
    perm_r_squared     <- numeric(nperm)
    perm_correlations  <- matrix(NA, nrow=nperm, ncol=length(cors))
  }
  
  for (i in seq_len(nperm)) {
    # Permute row order of predicted
    perm_idx <- sample(n_rows)
    perm_pred <- predicted[perm_idx, , drop=FALSE]

    # Compute metrics for this permutation
    perm_cormat <- cor(perm_pred, observed)
    perm_cors <- diag(perm_cormat)
    pmc  <- mean(perm_cors)
    
    # Calculate off-diagonal correlations for permutation
    n <- nrow(perm_cormat)
    perm_off_diag <- (sum(perm_cormat) - sum(perm_cors)) / (n*n - n)
    pss <- pmc - perm_off_diag  # Spatial specificity
    
    pvc  <- cor(as.vector(perm_pred), as.vector(observed))
    pmse <- mean((perm_pred - observed)^2)
    prsq <- 1 - sum((observed - perm_pred)^2) / tss

    # Update counters
    if (pmc >= mean_cor) count_better_mean_corr       <- count_better_mean_corr + 1
    if (pss >= spatial_specificity) count_better_spatial_spec <- count_better_spatial_spec + 1
    if (pvc >= voxel_cor) count_better_voxel_corr     <- count_better_voxel_corr + 1
    if (pmse <= mse) count_better_mse                 <- count_better_mse + 1
    if (prsq >= r_squared) count_better_r_squared     <- count_better_r_squared + 1

    # Sums for z-scores
    sum_mean_corr    <- sum_mean_corr + pmc
    sum_sq_mean_corr <- sum_sq_mean_corr + pmc^2
    sum_spatial_spec    <- sum_spatial_spec + pss
    sum_sq_spatial_spec <- sum_sq_spatial_spec + pss^2
    sum_voxel_corr    <- sum_voxel_corr + pvc
    sum_sq_voxel_corr <- sum_sq_voxel_corr + pvc^2
    sum_mse    <- sum_mse + pmse
    sum_sq_mse <- sum_sq_mse + pmse^2
    sum_r_squared    <- sum_r_squared + prsq
    sum_sq_r_squared <- sum_sq_r_squared + prsq^2

    # Possibly store the full permutation distribution
    if (save_distributions) {
      perm_mean_corr[i]    <- pmc
      perm_spatial_spec[i] <- pss
      perm_voxel_corr[i]   <- pvc
      perm_mse[i]          <- pmse
      perm_r_squared[i]    <- prsq
      perm_correlations[i,] <- perm_cors
    }
  }

  # Compute p-values
  p_mean_corr  <- count_better_mean_corr / nperm
  p_spatial_spec <- count_better_spatial_spec / nperm
  p_voxel_corr <- count_better_voxel_corr / nperm
  p_mse        <- count_better_mse / nperm
  p_r_squared  <- count_better_r_squared / nperm

  # Compute means and SDs of permutation distributions
  mean_perm_mean_corr   <- sum_mean_corr / nperm
  sd_perm_mean_corr     <- sqrt((sum_sq_mean_corr / nperm) - mean_perm_mean_corr^2)
  mean_perm_spatial_spec <- sum_spatial_spec / nperm
  sd_perm_spatial_spec   <- sqrt((sum_sq_spatial_spec / nperm) - mean_perm_spatial_spec^2)
  mean_perm_voxel_corr  <- sum_voxel_corr / nperm
  sd_perm_voxel_corr    <- sqrt((sum_sq_voxel_corr / nperm) - mean_perm_voxel_corr^2)
  mean_perm_mse         <- sum_mse / nperm
  sd_perm_mse           <- sqrt((sum_sq_mse / nperm) - mean_perm_mse^2)
  mean_perm_r_squared   <- sum_r_squared / nperm
  sd_perm_r_squared     <- sqrt((sum_sq_r_squared / nperm) - mean_perm_r_squared^2)

  # z-scores
  eps <- .Machine$double.eps
  z_mean_corr  <- (mean_cor  - mean_perm_mean_corr)  / (sd_perm_mean_corr  + eps)
  z_spatial_spec <- (spatial_specificity - mean_perm_spatial_spec) / (sd_perm_spatial_spec + eps)
  z_voxel_corr <- (voxel_cor - mean_perm_voxel_corr) / (sd_perm_voxel_corr + eps)
  z_mse        <- (mean_perm_mse - mse) / (sd_perm_mse + eps)
  z_r_squared  <- (r_squared - mean_perm_r_squared) / (sd_perm_r_squared + eps)

  out <- list(
    p_values = c(mean_correlation = p_mean_corr,
                 spatial_specificity = p_spatial_spec,
                 voxel_correlation= p_voxel_corr,
                 mse = p_mse,
                 r_squared = p_r_squared),
    z_scores = c(mean_correlation = z_mean_corr,
                 spatial_specificity = z_spatial_spec,
                 voxel_correlation= z_voxel_corr,
                 mse = z_mse,
                 r_squared = z_r_squared)
  )

  if (save_distributions) {
    out$permutation_distributions <- list(
      mean_correlation  = perm_mean_corr,
      spatial_specificity = perm_spatial_spec,
      voxel_correlation = perm_voxel_corr,
      mse               = perm_mse,
      r_squared         = perm_r_squared,
      correlations      = perm_correlations
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
  
  if (ncol(observed) != ncol(predicted)) {
    stop(sprintf("Mismatch in columns: predicted has %d, observed has %d.", 
                 ncol(predicted), ncol(observed)))
  }
  
  # Base RSA metrics
  cormat     <- cor(predicted, observed)
  cors       <- diag(cormat)
  mean_cor   <- mean(cors)
  
  # Calculate mean of off-diagonal correlations
  n <- nrow(cormat)
  off_diag_cors <- (sum(cormat) - sum(cors)) / (n*n - n)  # Sum all elements minus diagonal, divided by number of off-diagonal elements
  
  # New metric: mean diagonal correlation minus mean off-diagonal correlation
  # This measures the spatial specificity of the model
  spatial_specificity <- mean_cor - off_diag_cors
  
  voxel_cor  <- cor(as.vector(predicted), as.vector(observed))
  mse        <- mean((predicted - observed)^2)
  rss        <- sum((observed - predicted)^2)
  tss        <- sum((observed - mean(observed))^2)
  r_squared  <- 1 - (rss / tss)
  
  perm_results <- NULL
  if (nperm > 0) {
    perm_results <- .perm_test_feature_rsa(
      observed = observed,
      predicted = predicted,
      nperm = nperm,
      save_distributions = save_distributions,
      mean_cor = mean_cor,
      spatial_specificity = spatial_specificity,
      voxel_cor = voxel_cor,
      mse = mse,
      r_squared = r_squared,
      cors = cors
    )
  }
  
  list(
    correlations        = cors,
    mean_correlation    = mean_cor,
    off_diag_correlation= off_diag_cors,
    spatial_specificity = spatial_specificity,
    voxel_correlation   = voxel_cor,
    mse                 = mse,
    r_squared           = r_squared,
    permutation_results = perm_results
  )
}


#' @export
train_model.feature_rsa_model <- function(obj, train_dat, ytrain, indices, ...) {
  # "train_dat" = X (the brain data for region or searchlight)
  # "ytrain"    = F (the feature set for these training rows)
  X <- as.matrix(train_dat)  
  Fsub <- as.matrix(ytrain)
  
  result <- list()
 
  if (obj$method == "pls") {
    # Number of observations:
    N <- nrow(X)
    V <- ncol(X)
    P <- ncol(Fsub)
    
    # 1) Combine X and Fsub in a single data frame
    #    * Brain data columns become v1..vV
    #    * Feature columns become f1..fP
    train_df <- data.frame(
      X,      # This creates V columns
      Fsub    # This creates P columns
    )
    
    # 3) Fit PLS with formula interface
    #    scale=TRUE => let pls automatically scale each column
    #    center=TRUE => mean-center
    fit <- pls::plsr(
      formula = X ~ Fsub,
      data    = train_df,
      ncomp   = obj$max_comps,  # or pick # of comps
      scale   = TRUE,
      center  = TRUE
      # Possibly specify validation="CV" etc. if you want
    )
    
    # 4) Store the model and any info
    result$trained_model <- fit
    # Use the actual # of comps from fit if needed
    result$ncomp <- min(obj$max_comps, fit$ncomp)
    
  } else if (obj$method == "scca") {
    #
    # ---- SCCA Train ----
    #
    # Standardize X, F
    sx <- .standardize(X)
    sf <- .standardize(Fsub)
    
    # SCCA requires enough rows
    if (nrow(sx$X_sc) < 2 || nrow(sf$X_sc) < 2) {
      stop("Cannot perform SCCA: insufficient observations.")
    }
    
    scca_res <- tryCatch(
      whitening::scca(sx$X_sc, sf$X_sc, scale=FALSE),
      error = function(e) {
      stop(paste("SCCA error:", e$message))
      }
    )
    
    # The effective # of comps from SCCA
    effective_ncomp <- sum(abs(scca_res$lambda) > 1e-6)
    if (effective_ncomp < 1) {
      scca_res <- whitening::scca(sx$X_sc, sf$X_sc, scale=FALSE, lambda.cor=.001)
      effective_ncomp <- sum(abs(scca_res$lambda) > 1e-6)
    }
    
    ncomp <- min(effective_ncomp, obj$max_comps)
    if (ncomp < 1) {
      stop("No effective canonical components available.")
    }
    
    # Save into result
    result$trained_model <- scca_res
    result$scca_x_mean <- sx$mean
    result$scca_x_sd   <- sx$sd
    result$scca_f_mean <- sf$mean
    result$scca_f_sd   <- sf$sd
    result$ncomp       <- ncomp
    
  } else if (obj$method == "pca") {
    #
    # -- PCA with possible caching --
    #
    # 1) Make a cache key from 'indices'
    # (assuming 'indices' is a vector of row indices used for training)
    key <- .hash_row_indices(indices)
    
    # 2) Check if caching is enabled AND if we have a cached PCA
    if (isTRUE(obj$cache_pca) && !is.null(obj$pca_cache[[key]])) {
      # -- REUSE the cached PCA info
      pca_info <- obj$pca_cache[[key]]
      # standardize brain data X with the stored means/sds
    sx <- .standardize(X)
      # we do not re-run prcomp, but reuse
      pca_res <- pca_info$pca_res
      sf_mean <- pca_info$f_mean
      sf_sd   <- pca_info$f_sd
      
      # figure out k
      available_k <- ncol(pca_res$x)
      k <- min(obj$max_comps, available_k)
      if (k < 1) {
        stop("No principal components available (check data).")
      }
      PC_train_subset <- pca_res$x[, seq_len(k), drop=FALSE]
      
      # Regress X_sc on these PCs
      df_pcs <- as.data.frame(PC_train_subset)
      if (nrow(df_pcs) <= k) {
        stop("Insufficient data for PCA regression (more comps than rows).")
      }
      
      fit <- lm(sx$X_sc ~ ., data=df_pcs)
      coefs <- coef(fit)
      
      result$pcarot     = pca_res$rotation[, seq_len(k), drop=FALSE]
      result$pca_f_mean = sf_mean
      result$pca_f_sd   = sf_sd
      result$pca_coefs  = coefs
      result$pca_x_mean = sx$mean
      result$pca_x_sd   = sx$sd
      result$ncomp      = k
      
    } else {
      # -- NO CACHE HIT (or caching off), so we do normal PCA
      sx <- .standardize(X)
      sf <- .standardize(Fsub)
      
      pca_res <- prcomp(sf$X_sc, scale.=FALSE)
      
      available_k <- ncol(pca_res$x)
      k <- min(obj$max_comps, available_k)
      if (k < 1) {
        stop("No principal components available (check data).")
      }
      PC_train_subset <- pca_res$x[, seq_len(k), drop=FALSE]
      
      # Regress X_sc on these PCs
      df_pcs <- as.data.frame(PC_train_subset)
      if (nrow(df_pcs) <= k) {
        stop("Insufficient data for PCA regression (more comps than rows).")
      }
      fit <- lm(sx$X_sc ~ ., data=df_pcs)
      coefs <- coef(fit)
      
      result$pcarot     = pca_res$rotation[, seq_len(k), drop=FALSE]
      result$pca_f_mean = sf$mean
      result$pca_f_sd   = sf$sd
      result$pca_coefs  = coefs
      result$pca_x_mean = sx$mean
      result$pca_x_sd   = sx$sd
      result$ncomp      = k
      
      # 3) Save a minimal set of info to the cache if cache is on
      if (isTRUE(obj$cache_pca)) {
        obj$pca_cache[[key]] <- list(
          pca_res = pca_res,
          f_mean  = sf$mean,
          f_sd    = sf$sd
          # we do *not* store the brain data standardization, 
          # because that depends on each ROI. 
          # We only store the PCA decomposition for the features. 
        )
      }
    }
    
  } else {
    stop("Unknown method in feature_rsa_model.")
  }
  
  result
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
      error         = TRUE, 
      error_message = error_message
    ))
  }
  
  Xobs  <- as.data.frame(context$test)
  Ftest <- as.matrix(context$ytest)
  
  Xpred <- tryCatch({
    predict_model.feature_rsa_model(obj, result, Ftest)
  }, error=function(e) NULL)
  
  if (is.null(Xpred)) {
    return(tibble::tibble(
      observed      = list(NULL), 
      predicted     = list(NULL), 
      error         = TRUE, 
      error_message = "Prediction failed"
    ))
  }
  
  # Evaluate WITHOUT permutations at the fold level
  perf <- evaluate_model.feature_rsa_model(
    object = obj,
    predicted = Xpred,
    observed  = Xobs,
    nperm = 0  # no permutation here
  )
  
  # Summarize
  perf_mat <- matrix(
    c(perf$mean_correlation,
      perf$spatial_specificity,
      perf$voxel_correlation,
      perf$mse,
      perf$r_squared),
    nrow = 1,
    ncol = 5,
    dimnames = list(NULL, c("mean_correlation", "spatial_specificity", "voxel_correlation", "mse", "r_squared"))
  )
  
  tibble::tibble(
    observed    = list(Xobs),
    predicted   = list(Xpred),
    result      = list(perf$correlations),
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
  
  combined_observed  <- do.call(rbind, observed_list)
  combined_predicted <- do.call(rbind, predicted_list)
  
  # Now we do permutations (if nperm>0 in the model spec)
  perf <- evaluate_model.feature_rsa_model(
    object    = obj,
    predicted = combined_predicted,
    observed  = combined_observed,
    nperm     = obj$nperm,
    save_distributions = obj$save_distributions
  )
  
  # Collate results
  if (is.null(perf$permutation_results)) {
    perf_mat <- matrix(
      c(perf$mean_correlation,
        perf$spatial_specificity,
        perf$voxel_correlation,
        perf$mse,
        perf$r_squared),
      nrow = 1,
      ncol = 5,
      dimnames = list(NULL, c("mean_correlation", "spatial_specificity", "voxel_correlation", "mse", "r_squared"))
    )
  } else {
    perf_values <- c(
      perf$mean_correlation,
      perf$spatial_specificity,
      perf$voxel_correlation,
      perf$mse,
      perf$r_squared,
      perf$permutation_results$p_values,
      perf$permutation_results$z_scores
    )
    perf_names <- c(
      "mean_correlation",
      "spatial_specificity",
      "voxel_correlation",
      "mse",
      "r_squared",
      paste0("p_", names(perf$permutation_results$p_values)),
      paste0("z_", names(perf$permutation_results$z_scores))
    )
    
    perf_mat <- matrix(
      perf_values,
      nrow = 1,
      ncol = length(perf_values),
      dimnames = list(NULL, perf_names)
    )
  }
  
  tibble::tibble(
    result      = list(NULL),
    indices     = list(indices),
    performance = list(perf_mat),
    id          = id,
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
  
  combined_observed  <- do.call(rbind, observed_list)
  combined_predicted <- do.call(rbind, predicted_list)
  
  # Now we do permutations (if nperm>0 in the model spec)
  perf <- evaluate_model.feature_rsa_model(
    object    = obj,
    predicted = combined_predicted,
    observed  = combined_observed,
    nperm     = obj$nperm,
    save_distributions = obj$save_distributions
  )
  
  # Collate results
  if (is.null(perf$permutation_results)) {
    perf_mat <- matrix(
      c(perf$mean_correlation,
        perf$spatial_specificity,
        perf$voxel_correlation,
        perf$mse,
        perf$r_squared),
      nrow = 1,
      ncol = 5,
      dimnames = list(NULL, c("mean_correlation", "spatial_specificity", "voxel_correlation", "mse", "r_squared"))
    )
  } else {
    perf_values <- c(
      perf$mean_correlation,
      perf$spatial_specificity,
      perf$voxel_correlation,
      perf$mse,
      perf$r_squared,
      perf$permutation_results$p_values,
      perf$permutation_results$z_scores
    )
    perf_names <- c(
      "mean_correlation",
      "spatial_specificity",
      "voxel_correlation",
      "mse",
      "r_squared",
      paste0("p_", names(perf$permutation_results$p_values)),
      paste0("z_", names(perf$permutation_results$z_scores))
    )
    
    perf_mat <- matrix(
      perf_values,
      nrow = 1,
      ncol = length(perf_values),
      dimnames = list(NULL, perf_names)
    )
  }
  
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

#' Run searchlight analysis with a feature RSA model
#'
#' This method provides specialized handling for feature RSA models during searchlight analysis.
#'
#' @param model_spec A \code{feature_rsa_model} object
#' @param radius Numeric radius for the searchlight spheres
#' @param method Method for searchlight: "standard" or "randomized"
#' @param niter Number of iterations for randomized searchlight
#' @param ... Additional arguments passed to run_searchlight_base
#'
#' @export
run_searchlight.feature_rsa_model <- function(model_spec, radius = 8, 
                                             method = c("standard", "randomized"),
                                             niter = 4, ...) {
  method <- match.arg(method)
  
  # No hacky adjustments needed since the component limits are now part of the design
  # and properly passed through to the model
  
  # Just log what we're using
  if (model_spec$method == "pca" && !is.null(model_spec$max_pca_comps)) {
    futile.logger::flog.info("Running searchlight with feature RSA (PCA): using max %d components", 
                           model_spec$max_pca_comps)
  } else if (model_spec$method == "scca" && !is.null(model_spec$max_scca_comps)) {
    futile.logger::flog.info("Running searchlight with feature RSA (SCCA): using max %d components", 
                           model_spec$max_scca_comps)
  } else if (model_spec$method == "pls" && !is.null(model_spec$max_pls_comps)) {
    futile.logger::flog.info("Running searchlight with feature RSA (PLS): using max %d components", 
                           model_spec$max_pls_comps)
  }
  
  # Use appropriate combiner based on method
  if (method == "standard") {
    futile.logger::flog.info("Running standard searchlight with feature RSA model")
    run_searchlight_base(
      model_spec = model_spec,
      radius = radius,
      method = method,
      combiner = combine_rsa_standard,  # Use RSA-specific combiner
      ...
    )
  } else {
    futile.logger::flog.info("Running randomized searchlight with feature RSA model")
    run_searchlight_base(
      model_spec = model_spec,
      radius = radius,
      method = method,
      niter = niter,
      ...  # Use default combiner for randomized
    )
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
  
  # Display component limits
  cat(crayon::bold(crayon::blue("Component Limits:\n")))
  cat("  ", crayon::bold(crayon::blue("PCA max components: ")), 
      if(!is.null(x$max_pca_comps)) x$max_pca_comps else "Not set", "\n")
  cat("  ", crayon::bold(crayon::blue("SCCA max components: ")), 
      if(!is.null(x$max_scca_comps)) x$max_scca_comps else "Not set", "\n")
  cat("  ", crayon::bold(crayon::blue("PLS max components: ")), 
      if(!is.null(x$max_pls_comps)) x$max_pls_comps else "Not set", "\n")
  
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
  
  # Display component limit for the current method
  if (x$method == "pca") {
    comp_limit <- if(!is.null(x$max_pca_comps)) {
      x$max_pca_comps
    } else if (!is.null(x$design$max_pca_comps)) {
      x$design$max_pca_comps
    } else {
      "Default"
    }
    cat(crayon::bold(crayon::blue("PCA max components:     ")), comp_limit, "\n")
    
    # Indicate if PCA has been precomputed for efficiency
    if (!is.null(x$precomputed_pca)) {
      cat(crayon::bold(crayon::green("PCA Optimization:       ")), 
          "PCA precomputed for efficiency\n")
    }
  } else if (x$method == "scca") {
    comp_limit <- if(!is.null(x$max_scca_comps)) {
      x$max_scca_comps
    } else if (!is.null(x$design$max_scca_comps)) {
      x$design$max_scca_comps
    } else {
      "Default"
    }
    cat(crayon::bold(crayon::blue("SCCA max components:    ")), comp_limit, "\n")
  } else if (x$method == "pls") {
    comp_limit <- if(!is.null(x$max_pls_comps)) {
      x$max_pls_comps
    } else if (!is.null(x$design$max_pls_comps)) {
      x$design$max_pls_comps
    } else {
      "Default"
    }
    cat(crayon::bold(crayon::blue("PLS max components:     ")), comp_limit, "\n")
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




