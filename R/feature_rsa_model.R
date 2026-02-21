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
#'          feature space F by subsequent `feature_rsa_model` methods (PCA, PLS).
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
#' methods (like PCA, PLS in `feature_rsa_model`) can use, and it cannot
#' exceed the final feature dimensionality `k`.
#'
#' @examples
#' \donttest{
#'   S <- as.matrix(dist(matrix(rnorm(5*3), 5, 3)))
#'   labels <- factor(letters[1:5])
#'   des <- feature_rsa_design(S = S, labels = labels)
#' }
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
    ret$cv_labels <- labels
    ret$targets <- ret$F
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
    ret$cv_labels <- labels
    ret$targets <- ret$F
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
#'     \item{pls}{Partial Least Squares regression predicting X from F (via \code{pls::plsr}).}
#'     \item{pca}{Principal Component Regression predicting X from PCs of F (via \code{pls::pcr}).}
#'     \item{glmnet}{Elastic net regression predicting X from F using glmnet with multivariate Gaussian response.}
#'   }
#' @param crossval Optional cross-validation specification.
#' @param ncomp_selection Character string controlling how the number of components
#'   is chosen for \code{pls} and \code{pca} methods.  One of:
#'   \describe{
#'     \item{loo}{(Default) Fit with leave-one-out validation and select the
#'       fewest components within one standard error of the minimum RMSEP
#'       (\code{pls::selectNcomp}, method \code{"onesigma"}).}
#'     \item{pve}{Keep the fewest components whose cumulative explained
#'       variance reaches \code{pve_threshold} of the total explained by all
#'       fitted components.}
#'     \item{max}{Use all \code{max_comps} components (legacy behaviour).}
#'   }
#'   Ignored when \code{method = "glmnet"}.
#' @param pve_threshold Numeric in (0, 1].  When \code{ncomp_selection = "pve"},
#'   the proportion of total explained X-variance at which to stop adding
#'   components.  Default 0.9.
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
#' @param ... Additional arguments (currently unused). Passing deprecated
#'   arguments such as \code{cache_pca} now results in an error.
#'
#' @return A \code{feature_rsa_model} object (S3 class).
#'
#' @examples
#' \donttest{
#'   S <- as.matrix(dist(matrix(rnorm(5*3), 5, 3)))
#'   labels <- factor(letters[1:5])
#'   des <- feature_rsa_design(S = S, labels = labels)
#'   # mdl <- feature_rsa_model(dataset, des, method="pls")
#' }
#' @details
#' Feature RSA models analyze how well a feature matrix \code{F} (defined in the `design`)
#' relates to neural data \code{X}. The `max_comps` parameter, inherited from the `design` object,
#' sets an upper limit on the number of components fitted:
#'   - \strong{pls}: PLS regression via \code{pls::plsr}. Fits up to `max_comps` components;
#'     the actual number used for prediction is chosen by \code{ncomp_selection}.
#'   - \strong{pca}: Principal Component Regression via \code{pls::pcr}. Fits up to
#'     `max_comps` components; selection controlled by \code{ncomp_selection}.
#'   - \strong{glmnet}: Elastic net regression via \code{glmnet} with multivariate Gaussian
#'     response. Regularisation (lambda) can be auto-selected via \code{cv_glmnet=TRUE}.
#'
#' For \code{pls} and \code{pca}, the \code{ncomp_selection} argument determines how many
#' of the fitted components are actually used for prediction.  The default
#' (\code{"loo"}) fits the model with leave-one-out cross-validation and picks
#' the fewest components within one SE of the minimum RMSEP.
#'
#' **Performance Metrics** (computed by `evaluate_model` after cross-validation):
#'
#' *Condition-pattern metrics* (trial x trial correlation matrix):
#'   - `pattern_correlation`: Average correlation between the predicted and observed
#'     spatial patterns for corresponding trials (diagonal of the trial x trial
#'     correlation matrix computed across voxels).
#'   - `pattern_discrimination`: `pattern_correlation` minus the mean off-diagonal
#'     correlation.  Measures how much better the model predicts the correct trial's
#'     pattern compared to incorrect trials.
#'   - `pattern_rank_percentile`: For each trial, percentile rank of the correct
#'     pattern match among all candidates.  0.5 = chance, 1 = perfect.
#'
#' *Representational geometry*:
#'   - `rdm_correlation`: Spearman correlation between the upper triangles of the
#'     observed and predicted RDMs (defined as 1 - trial-by-trial correlation
#'     across voxels).  Captures similarity of representational geometry.
#'
#' *Global reconstruction metrics*:
#'   - `voxel_correlation`: Correlation of the flattened predicted and observed
#'     matrices (all trials x all voxels).
#'   - `mse`: Mean Squared Error.
#'   - `r_squared`: 1 - RSS/TSS.
#'
#' *Voxel encoding fidelity*:
#'   - `mean_voxelwise_temporal_cor`: Average per-voxel temporal correlation
#'     between predicted and observed time courses.
#'
#'   - `p_*`, `z_*`: If `nperm > 0`, permutation-based p-values and z-scores for
#'     the above metrics.
#'
#' The number of components actually used (`ncomp`) for the region/searchlight is
#' also included in the performance output.
#'
#' @export
feature_rsa_model <- function(dataset,
                               design,
                               method = c("pls", "pca", "glmnet"),
                               crossval = NULL,
                               ncomp_selection = c("loo", "max", "pve"),
                               pve_threshold = 0.9,
                               alpha = 0.5,
                               cv_glmnet = FALSE,
                               lambda = NULL,
                               nperm = 0,
                               permute_by = c("features", "observations"),
                               save_distributions = FALSE,
                               ...) {
  
  method <- match.arg(method)
  ncomp_selection <- match.arg(ncomp_selection)
  permute_by <- match.arg(permute_by)
  assertthat::assert_that(inherits(dataset, "mvpa_dataset"))
  assertthat::assert_that(inherits(design, "feature_rsa_design"))
  extra_args <- list(...)
  if ("cache_pca" %in% names(extra_args)) {
    stop("`cache_pca` is no longer supported. PCA caching was removed; use `ncomp_selection` controls for `method='pca'`.")
  }

  if (ncomp_selection == "pve") {
    assertthat::assert_that(is.numeric(pve_threshold) && pve_threshold > 0 && pve_threshold <= 1,
                           msg = "pve_threshold must be in (0, 1]")
  }
  
  # Additional validation for dataset dimensions
  mask_dims <- dim(dataset$mask)[1:3]
  total_voxels <- prod(mask_dims)
  active_voxels <- sum(dataset$mask > 0)
  
  if (total_voxels <= 1) {
    stop("Invalid dataset for feature_rsa_model: Only 1 voxel detected (dimensions ",
         paste(mask_dims, collapse="x"),
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

  # Component selection strategy (for PLS/PCA)
  model_spec$ncomp_selection <- ncomp_selection
  model_spec$pve_threshold <- pve_threshold

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


#' Onesigma component selection for multivariate responses
#'
#' \code{pls::selectNcomp} only supports univariate responses.
#' This helper averages CV-MSEP across all response columns and
#' applies the same onesigma rule: pick the fewest components whose
#' mean MSEP is within one SE of the minimum.
#' @noRd
.selectNcomp_mv <- function(model, method = "onesigma") {
  msep_obj  <- pls::MSEP(model, estimate = "CV")
  msep_vals <- msep_obj$val                        # [estimate, response, ncomp+1]

  nresp     <- dim(msep_vals)[2]
  ncomp_max <- dim(msep_vals)[3] - 1L              # first slot = intercept-only
  if (ncomp_max < 1L) return(1L)

  # Average MSEP across responses for each ncomp (drop intercept slot)
  avg_msep <- drop(apply(msep_vals[1, , -1, drop = FALSE], 3, mean))

  # Guard against all-NaN MSEP (degenerate model)
  if (all(is.na(avg_msep) | !is.finite(avg_msep))) return(NA_integer_)

  min_idx <- which.min(avg_msep)
  if (length(min_idx) == 0L) return(NA_integer_)
  min_val <- avg_msep[min_idx]

  if (method == "onesigma" && nresp > 1L) {
    # SE of mean MSEP across responses at the optimum
    per_resp <- msep_vals[1, , min_idx + 1L]
    se       <- sd(per_resp) / sqrt(nresp)
    thresh   <- min_val + se
    selected <- which(avg_msep <= thresh)[1]
    if (is.na(selected)) return(as.integer(min_idx))
    return(as.integer(selected))
  }

  as.integer(min_idx)
}


#' @noRd
.standardize <- function(X) {
  cm <- colMeans(X)
  csd <- apply(X, 2, sd)
  csd[csd == 0] <- 1
  X_sc <- scale(X, center=cm, scale=csd)
  list(X_sc = X_sc, mean = cm, sd = csd)
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
  if (is.null(fit$trained_model)) {
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
    if (method %in% c("pls", "pca")) {
      # PLS and PCA (via pls::pcr) share the same prediction path
      pls_model    <- fit$trained_model
      f_mean       <- fit$f_mean
      f_sd         <- fit$f_sd
      x_mean       <- fit$x_mean
      x_sd         <- fit$x_sd
      ncomp_to_use <- fit$ncomp

      if (is.null(pls_model) || is.null(f_mean) || is.null(f_sd) ||
          is.null(x_mean) || is.null(x_sd) || is.null(ncomp_to_use)) {
        stop(sprintf("predict_model (%s): Missing essential components in the fit object.", method))
      }
      if (ncomp_to_use < 1) {
        stop(sprintf("predict_model (%s): ncomp (%d) < 1.", method, ncomp_to_use))
      }

      expected_cols <- length(f_mean)
      if (ncol(F_new) != expected_cols) {
        stop(sprintf("predict_model (%s): Feature column mismatch. Expected %d, got %d.",
                     method, expected_cols, ncol(F_new)))
      }

      sf_test <- scale(F_new, center = f_mean, scale = f_sd)
      if (any(is.nan(sf_test))) {
        stop(sprintf("predict_model (%s): NaNs after standardization.", method))
      }

      preds_raw <- predict(pls_model, newdata = sf_test, ncomp = ncomp_to_use)
      preds_sc  <- drop(preds_raw)

      if (!is.matrix(preds_sc)) {
        preds_sc <- matrix(preds_sc, nrow = nrow(F_new), ncol = length(x_mean))
      }

      preds <- sweep(sweep(preds_sc, 2, x_sd, "*"), 2, x_mean, "+")
      return(preds)

    } else if (method == "glmnet") {
      return(.predict_glmnet(fit, F_new))
    } else {
      stop(paste("Unknown method in predict_model.feature_rsa_model:", method))
    }
  }, error = function(e) {
      error_msg <- sprintf("predict_model (%s): Prediction failed - %s", method, e$message)
      futile.logger::flog.error(error_msg)
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

  # Early diagnostic: warn if predictions are constant across trials
  if (nrow(predictions) > 1) {
    pred_sds <- apply(predictions, 2, stats::sd)
    if (all(pred_sds < 1e-12, na.rm = TRUE)) {
      futile.logger::flog.debug(
        "predict_model (%s): all %d voxel predictions are constant across %d trials (max sd=%.2e). Model has no predictive power for this ROI.",
        method, ncol(predictions), nrow(predictions),
        max(pred_sds, na.rm = TRUE))
    }
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
#' @param pattern_cor Scalar: the observed pattern correlation
#' @param pattern_discrim Scalar: the observed pattern discrimination
#' @param pattern_rank Scalar: the observed pattern rank percentile
#' @param rdm_cor Scalar: the observed RDM correlation
#' @param voxel_cor Scalar: the observed voxel correlation
#' @param mse Scalar: the observed MSE
#' @param r_squared Scalar: the observed R^2
#' @param mean_voxelwise_temporal_cor Scalar: the observed mean voxelwise temporal correlation
#' @param valid_col Integer vector of valid voxel columns
#' @return A list with p-values, z-scores, and optionally a list of permutations
.perm_test_feature_rsa <- function(observed,
                                   predicted,
                                   nperm,
                                   save_distributions,
                                   pattern_cor,
                                   pattern_discrim,
                                   pattern_rank,
                                   rdm_cor,
                                   voxel_cor,
                                   mse,
                                   r_squared,
                                   mean_voxelwise_temporal_cor,
                                   valid_col)
{
  message("Performing permutation tests with ", nperm, " permutations... (feature_rsa_model)")

  n_rows <- nrow(predicted)
  sd_thresh <- 1e-12
  tss <- sum((observed - mean(observed, na.rm = TRUE))^2, na.rm = TRUE)

  ## Metric names and observed values (order matters)
  metric_names <- c("pattern_correlation", "pattern_discrimination",
                    "pattern_rank_percentile", "rdm_correlation",
                    "voxel_correlation",
                    "mse", "r_squared", "mean_voxelwise_temporal_cor")
  obs_vals <- c(pattern_cor, pattern_discrim, pattern_rank, rdm_cor,
                voxel_cor, mse, r_squared, mean_voxelwise_temporal_cor)
  names(obs_vals) <- metric_names
  n_met <- length(metric_names)

  ## Accumulators
  count_better <- rep(0L, n_met); names(count_better) <- metric_names
  sum_perm     <- rep(0,  n_met); names(sum_perm) <- metric_names
  sum_sq_perm  <- rep(0,  n_met); names(sum_sq_perm) <- metric_names
  n_valid_perm <- rep(0L, n_met); names(n_valid_perm) <- metric_names

  if (save_distributions) {
    dist_mat <- matrix(NA_real_, nrow = nperm, ncol = n_met,
                       dimnames = list(NULL, metric_names))
  }

  for (i in seq_len(nperm)) {
    perm_idx  <- sample(n_rows)
    perm_pred <- predicted[perm_idx, , drop = FALSE]

    ## -- Condition-pattern metrics (trial x trial) --
    ppc <- ppd <- ppr <- prdm <- NA_real_
    obs_row_sd  <- apply(observed[, valid_col, drop = FALSE], 1, stats::sd)
    pred_row_sd <- apply(perm_pred[, valid_col, drop = FALSE], 1, stats::sd)
    vr <- which(obs_row_sd > sd_thresh & pred_row_sd > sd_thresh)
    if (length(vr) >= 2) {
      cm <- tryCatch(stats::cor(t(perm_pred[vr, valid_col, drop = FALSE]),
                                t(observed[vr, valid_col, drop = FALSE]),
                                use = "pairwise.complete.obs"),
                     error = function(e) NULL)
      if (!is.null(cm)) {
        dc  <- diag(cm)
        ppc <- mean(dc, na.rm = TRUE)
        nc  <- nrow(cm)
        if (nc > 1) {
          off_vals <- cm[row(cm) != col(cm)]
          off_vals <- off_vals[!is.na(off_vals)]
          od <- if (length(off_vals) > 0) mean(off_vals) else NA_real_
        } else {
          od <- NA_real_
        }
        ppd <- if (is.finite(ppc) && is.finite(od)) ppc - od else NA_real_
        rnk <- numeric(nc)
        for (j in seq_len(nc)) {
          rc <- cm[j, ]
          dn <- sum(!is.na(rc)) - 1
          rnk[j] <- if (dn > 0 && is.finite(rc[j])) {
            (sum(rc <= rc[j], na.rm = TRUE) - 1) / dn
          } else {
            NA_real_
          }
        }
        ppr <- mean(rnk, na.rm = TRUE)
      }
    }

    ## -- Representational geometry (RDM correlation) --
    if (length(vr) >= 3) {
      pc <- tryCatch(stats::cor(t(perm_pred[vr, valid_col, drop = FALSE]),
                                use = "pairwise.complete.obs"),
                     error = function(e) NULL)
      oc <- tryCatch(stats::cor(t(observed[vr, valid_col, drop = FALSE]),
                                use = "pairwise.complete.obs"),
                     error = function(e) NULL)
      if (!is.null(pc) && !is.null(oc)) {
        pd <- 1 - pc
        od <- 1 - oc
        pv <- pd[upper.tri(pd)]
        ov <- od[upper.tri(od)]
        if (length(pv) >= 2 && length(pv) == length(ov)) {
          prdm <- tryCatch(stats::cor(pv, ov, method = "spearman", use = "complete.obs"),
                           error = function(e) NA_real_)
        }
      }
    }

    ## -- Global reconstruction --
    pvc  <- tryCatch(stats::cor(as.vector(perm_pred[, valid_col, drop = FALSE]),
                                as.vector(observed[, valid_col, drop = FALSE]),
                                use = "pairwise.complete.obs"),
                     error = function(e) NA_real_)
    pmse <- mean((perm_pred - observed)^2, na.rm = TRUE)
    prsq <- if (tss > 0) 1 - sum((observed - perm_pred)^2, na.rm = TRUE) / tss else NA_real_

    ## -- Voxel encoding --
    pmvtc <- NA_real_
    if (n_rows > 1 && length(valid_col) > 0) {
      pmvtc <- mean(vapply(valid_col, function(j) {
        tryCatch(stats::cor(observed[, j], perm_pred[, j], use = "pairwise.complete.obs"), error = function(e) NA_real_)
      }, numeric(1)), na.rm = TRUE)
    }

    pvals <- c(ppc, ppd, ppr, prdm, pvc, pmse, prsq, pmvtc)

    ## Update accumulators
    for (m in seq_len(n_met)) {
      pv <- pvals[m]
      ov <- obs_vals[m]
      if (!is.na(pv)) {
        n_valid_perm[m] <- n_valid_perm[m] + 1L
        sum_perm[m]    <- sum_perm[m] + pv
        sum_sq_perm[m] <- sum_sq_perm[m] + pv^2
        if (!is.na(ov)) {
          better <- if (metric_names[m] == "mse") pv <= ov else pv >= ov
          if (better) count_better[m] <- count_better[m] + 1L
        }
      }
    }
    if (save_distributions) dist_mat[i, ] <- pvals
  }

  ## --- p-values and z-scores ---
  eps <- .Machine$double.eps
  n_eff <- n_valid_perm

  p_values <- vapply(seq_len(n_met), function(m) {
    if (n_eff[m] > 0) (count_better[m] + 1) / (n_eff[m] + 1) else NA_real_
  }, numeric(1))
  names(p_values) <- metric_names

  z_scores <- vapply(seq_len(n_met), function(m) {
    if (n_eff[m] > 0) {
      mn <- sum_perm[m] / n_eff[m]
      sd_p <- sqrt(max(0, sum_sq_perm[m] / n_eff[m] - mn^2))
      sd_use <- max(sd_p, eps)
      if (metric_names[m] == "mse") (mn - obs_vals[m]) / sd_use
      else (obs_vals[m] - mn) / sd_use
    } else NA_real_
  }, numeric(1))
  names(z_scores) <- metric_names

  out <- list(p_values = p_values, z_scores = z_scores)
  if (save_distributions) {
    out$permutation_distributions <- as.list(as.data.frame(dist_mat))
  }
  out
}



#' Evaluate model performance for feature RSA
#'
#' Computes condition-pattern metrics (trial x trial correlation matrix),
#' voxel-level encoding metrics, global reconstruction metrics (MSE, R-squared),
#' and optionally performs permutation tests.
#'
#' @param object The feature RSA model
#' @param predicted Matrix of predicted values (observations x voxels)
#' @param observed Matrix of observed values (observations x voxels)
#' @param nperm Number of permutations for statistical testing (default: 0)
#' @param save_distributions Logical indicating whether to save full permutation distributions
#' @param ... Additional arguments
#'
#' @return A list containing:
#'   \describe{
#'     \item{pattern_correlation}{Mean diagonal of the trial x trial correlation
#'       matrix -- how well the predicted spatial pattern for each trial matches
#'       the correct observed pattern.}
#'     \item{pattern_discrimination}{Diagonal minus off-diagonal of the trial x
#'       trial correlation matrix -- how much better the correct trial is matched
#'       than incorrect trials.}
#'     \item{pattern_rank_percentile}{For each trial, percentile rank of the
#'       correct pattern among all candidates.  0.5 = chance, 1 = perfect.}
#'     \item{voxel_correlation}{Correlation of the flattened predicted and
#'       observed matrices (global reconstruction quality).}
#'     \item{mse}{Mean squared error.}
#'     \item{r_squared}{1 - RSS/TSS.}
#'     \item{mean_voxelwise_temporal_cor}{Average per-voxel temporal correlation
#'       (encoding fidelity).}
#'     \item{permutation_results}{If \code{nperm > 0}, a list with p-values and
#'       z-scores for each metric.}
#'   }
#' @examples
#' \dontrun{
#'   # Internal S3 method called after cross-validation
#'   # perf <- evaluate_model(feature_rsa_model, newdata, observed)
#' }
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

  if (nrow(observed) != nrow(predicted)) {
    stop(sprintf("Mismatch in rows: predicted has %d, observed has %d.",
                 nrow(predicted), nrow(observed)))
  }

  if (ncol(observed) != ncol(predicted)) {
    stop(sprintf("Mismatch in columns: predicted has %d, observed has %d.",
                 ncol(predicted), ncol(observed)))
  }

  sd_thresh <- 1e-12
  n_obs  <- nrow(observed)
  n_vox  <- ncol(observed)

  ## -- helpers to avoid repeating variance checks --
  obs_sd  <- apply(observed, 2, stats::sd)
  pred_sd <- apply(predicted, 2, stats::sd)
  valid_col <- which(obs_sd > sd_thresh & pred_sd > sd_thresh)

  ## -- check rows (trials) have variance across voxels --
  obs_row_sd  <- apply(observed[, valid_col, drop = FALSE], 1, stats::sd)
  pred_row_sd <- apply(predicted[, valid_col, drop = FALSE], 1, stats::sd)
  valid_row   <- which(obs_row_sd > sd_thresh & pred_row_sd > sd_thresh)

  na_result <- list(
    pattern_correlation        = NA_real_,
    pattern_discrimination     = NA_real_,
    pattern_rank_percentile    = NA_real_,
    rdm_correlation            = NA_real_,
    voxel_correlation          = NA_real_,
    mse                        = mean((predicted - observed)^2, na.rm = TRUE),
    r_squared                  = NA_real_,
    mean_voxelwise_temporal_cor = NA_real_,
    permutation_results        = NULL
  )

  if (length(valid_col) == 0) {
    n_obs_ok  <- sum(obs_sd  > sd_thresh, na.rm = TRUE)
    n_pred_ok <- sum(pred_sd > sd_thresh, na.rm = TRUE)
    futile.logger::flog.warn(
      paste0("evaluate_model: No columns with finite variance; returning NA metrics. ",
             "dims=%d obs x %d vox, obs cols with var=%d (max sd=%.2e), ",
             "pred cols with var=%d (max sd=%.2e). ",
             if (n_pred_ok == 0 && n_obs_ok > 0)
               "Predictions are constant across trials -- the model likely has no predictive power for this ROI."
             else if (n_obs_ok == 0)
               "Observed data has no cross-trial variance -- ROI may be too small or data is constant."
             else
               "Both predicted and observed lack variance."),
      n_obs, n_vox, n_obs_ok,
      if (all(is.na(obs_sd))) NA_real_ else max(obs_sd, na.rm = TRUE),
      n_pred_ok,
      if (all(is.na(pred_sd))) NA_real_ else max(pred_sd, na.rm = TRUE))
    return(na_result)
  }

  ## ================================================================
  ## 1. Condition-pattern metrics  (trial x trial cormat)
  ## ================================================================
  ## cor(t(predicted), t(observed)):  n_obs x n_obs
  ##   entry (i,j) = cor(predicted[i, ], observed[j, ])  across voxels
  pattern_cor     <- NA_real_
  pattern_discrim <- NA_real_
  pattern_rank    <- NA_real_
  rdm_cor         <- NA_real_

  if (length(valid_row) >= 2) {
    pmat <- predicted[valid_row, valid_col, drop = FALSE]
    omat <- observed[valid_row,  valid_col, drop = FALSE]
    cormat_cond <- stats::cor(t(pmat), t(omat), use = "pairwise.complete.obs")

    diag_cors   <- diag(cormat_cond)
    pattern_cor <- mean(diag_cors, na.rm = TRUE)

    nc <- nrow(cormat_cond)
    if (nc > 1) {
      off_vals <- cormat_cond[row(cormat_cond) != col(cormat_cond)]
      off_vals <- off_vals[!is.na(off_vals)]
      off_diag <- if (length(off_vals) > 0) mean(off_vals) else NA_real_
    } else {
      off_diag <- NA_real_
    }
    pattern_discrim <- if (is.finite(pattern_cor) && is.finite(off_diag)) {
      pattern_cor - off_diag
    } else {
      NA_real_
    }

    # rank percentile per trial
    ranks <- numeric(nc)
    for (i in seq_len(nc)) {
      row_cors <- cormat_cond[i, ]
      denom <- sum(!is.na(row_cors)) - 1
      ranks[i] <- if (denom > 0 && is.finite(row_cors[i])) {
        (sum(row_cors <= row_cors[i], na.rm = TRUE) - 1) / denom
      } else {
        NA_real_
      }
    }
    pattern_rank <- mean(ranks, na.rm = TRUE)
  }

  ## ================================================================
  ## 1b. Representational geometry metric (RDM correlation)
  ## ================================================================
  ## Compute within-space RDMs (1 - trial-by-trial correlation) and correlate
  ## upper triangles between predicted and observed.
  if (length(valid_row) >= 3) {
    pmat <- predicted[valid_row, valid_col, drop = FALSE]
    omat <- observed[valid_row,  valid_col, drop = FALSE]
    pc <- tryCatch(stats::cor(t(pmat), use = "pairwise.complete.obs"), error = function(e) NULL)
    oc <- tryCatch(stats::cor(t(omat), use = "pairwise.complete.obs"), error = function(e) NULL)
    if (!is.null(pc) && !is.null(oc)) {
      prdm <- 1 - pc
      ordm <- 1 - oc
      pv <- prdm[upper.tri(prdm)]
      ov <- ordm[upper.tri(ordm)]
      if (length(pv) >= 2 && length(pv) == length(ov)) {
        rdm_cor <- tryCatch(stats::cor(pv, ov, method = "spearman", use = "complete.obs"),
                            error = function(e) NA_real_)
      }
    }
  }
  if (!is.finite(rdm_cor)) rdm_cor <- NA_real_

  ## ================================================================
  ## 2. Global reconstruction metrics
  ## ================================================================
  pred_vec <- as.vector(predicted[, valid_col, drop = FALSE])
  obs_vec  <- as.vector(observed[,  valid_col, drop = FALSE])
  voxel_cor <- if (stats::sd(pred_vec) > sd_thresh && stats::sd(obs_vec) > sd_thresh) {
    stats::cor(pred_vec, obs_vec, use = "pairwise.complete.obs")
  } else {
    NA_real_
  }

  mse <- mean((predicted - observed)^2, na.rm = TRUE)
  rss <- sum((observed - predicted)^2, na.rm = TRUE)
  tss <- sum((observed - mean(observed, na.rm = TRUE))^2, na.rm = TRUE)
  r_squared <- if (tss == 0) NA_real_ else 1 - (rss / tss)

  ## ================================================================
  ## 3. Voxel encoding fidelity
  ## ================================================================
  mean_voxelwise_temporal_cor <- NA_real_
  if (n_obs > 1 && length(valid_col) > 0) {
    vcors <- vapply(valid_col, function(j) {
      if (stats::sd(observed[, j]) > sd_thresh && stats::sd(predicted[, j]) > sd_thresh) {
        tryCatch(stats::cor(observed[, j], predicted[, j], use = "pairwise.complete.obs"), error = function(e) NA_real_)
      } else {
        NA_real_
      }
    }, numeric(1))
    mean_voxelwise_temporal_cor <- mean(vcors, na.rm = TRUE)
  }
  if (!is.finite(mean_voxelwise_temporal_cor)) mean_voxelwise_temporal_cor <- NA_real_

  ## ================================================================
  ## 4. Permutation testing
  ## ================================================================
  perm_results <- NULL
  if (nperm > 0) {
    perm_results <- .perm_test_feature_rsa(
      observed  = observed,
      predicted = predicted,
      nperm     = nperm,
      save_distributions = save_distributions,
      pattern_cor    = pattern_cor,
      pattern_discrim = pattern_discrim,
      pattern_rank   = pattern_rank,
      rdm_cor        = rdm_cor,
      voxel_cor      = voxel_cor,
      mse            = mse,
      r_squared      = r_squared,
      mean_voxelwise_temporal_cor = mean_voxelwise_temporal_cor,
      valid_col      = valid_col
    )
  }

  list(
    pattern_correlation         = pattern_cor,
    pattern_discrimination      = pattern_discrim,
    pattern_rank_percentile     = pattern_rank,
    rdm_correlation             = rdm_cor,
    voxel_correlation           = voxel_cor,
    mse                         = mse,
    r_squared                   = r_squared,
    mean_voxelwise_temporal_cor = mean_voxelwise_temporal_cor,
    permutation_results         = perm_results
  )
}


#' @rdname train_model
#' @param obj An object of class \code{feature_rsa_model}.
#' @param X Brain data (samples x voxels).
#' @param y Feature matrix used for RSA (samples x features).
#' @param indices Spatial indices associated with the training data.
#' @param ... Additional arguments.
#' @return A list containing RSA metrics and, if requested, permutation results.
#' @method train_model feature_rsa_model
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
  
  # ---- PLS / PCA (unified via pls package) ----
  if (obj$method %in% c("pls", "pca")) {
    fit_func <- if (obj$method == "pls") pls::plsr else pls::pcr
    method_label <- toupper(obj$method)
    ncomp_sel <- obj$ncomp_selection %||% "max"

    pls_res <- tryCatch({
      # Check for near-zero variance before standardization
      var_X <- apply(X, 2, var, na.rm = TRUE)
      var_F <- apply(Fsub, 2, var, na.rm = TRUE)
      if (any(var_X < 1e-10) || any(var_F < 1e-10)) {
        stop("Near zero variance detected in X or F before standardization.")
      }

      sx <- .standardize(X)
      sf <- .standardize(Fsub)

      if (any(sx$sd < 1e-10) || any(sf$sd < 1e-10)) {
        stop("Near zero variance detected after standardization.")
      }

      max_k_possible <- min(nrow(sf$X_sc) - 1, ncol(sf$X_sc))
      k <- min(obj$max_comps, max_k_possible)

      if (k < 1) {
        stop(sprintf("Number of components (%d) < 1 (max_comps: %d, max_possible: %d).",
                     k, obj$max_comps, max_k_possible))
      }

      validation <- if (ncomp_sel == "loo") "LOO" else "none"
      model <- fit_func(sx$X_sc ~ sf$X_sc, ncomp = k, scale = FALSE,
                        validation = validation)

      # --- Component selection (always >= 1) ---
      ncomp_use <- k
      if (ncomp_sel == "loo" && k > 1) {
        ncomp_use <- tryCatch({
          nc <- .selectNcomp_mv(model, method = "onesigma")
          if (is.na(nc) || nc < 1L) {
            futile.logger::flog.warn(
              "train_model (%s): selectNcomp returned %s; falling back to %d components.",
              method_label, as.character(nc), k)
            k
          } else {
            nc
          }
        }, error = function(e) {
          futile.logger::flog.warn(
            "train_model (%s): selectNcomp failed (%s); using all %d components.",
            method_label, e$message, k)
          k
        })
      } else if (ncomp_sel == "pve") {
        xvar <- pls::explvar(model)
        cum_ratio <- cumsum(xvar) / sum(xvar)
        idx <- which(cum_ratio >= obj$pve_threshold)[1]
        if (is.na(idx)) {
          futile.logger::flog.warn(
            "train_model (%s): no component reaches pve_threshold=%.2f (max cum ratio=%.4f); using all %d components.",
            method_label, obj$pve_threshold, max(cum_ratio, na.rm = TRUE), k)
          ncomp_use <- k
        } else {
          ncomp_use <- max(1L, idx)
        }
      }

      list(model = model, sx = sx, sf = sf, ncomp_use = ncomp_use)
    }, error = function(e) {
      error_msg <- sprintf("train_model (%s): Error during training - %s",
                           method_label, e$message)
      futile.logger::flog.error(error_msg)
      list(error = error_msg)
    })

    if (!is.null(pls_res$error)) {
      result$error <- pls_res$error
      return(result)
    }

    result$trained_model <- pls_res$model
    result$x_mean        <- pls_res$sx$mean
    result$x_sd          <- pls_res$sx$sd
    result$f_mean        <- pls_res$sf$mean
    result$f_sd          <- pls_res$sf$sd
    result$ncomp         <- pls_res$ncomp_use

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
        } else if (!is.null(lambda_to_use)) {
           # User supplied an explicit lambda
           if (length(final_fit$lambda) > 0) final_fit$lambda[1] else lambda_to_use
        } else if (!is.null(final_fit$lambda) && length(final_fit$lambda) > 1) {
           # No lambda specified and no CV: glmnet auto-generated a sequence
           # (descending from lambda_max). Picking lambda[1] = lambda_max would
           # shrink all coefficients to zero.  Use 1% of lambda_max as a
           # heuristic that provides light regularisation without collapsing
           # predictions to the mean.
           heuristic <- final_fit$lambda[1] * 0.01
           futile.logger::flog.info(
             "train_model (GLMNet): no lambda specified; using heuristic lambda = %.4e (1%% of lambda_max). Set lambda explicitly or use cv_glmnet=TRUE for optimal selection.",
             heuristic)
           heuristic
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


#' @rdname y_train-methods
#' @export
y_train.feature_rsa_model <- function(obj) {
  obj$design$targets  # Feature matrix for cross-validation data splitting
}

#' @rdname y_train-methods
#' @export
y_train.feature_rsa_design <- function(obj) {
  obj$targets  # Feature matrix for cross-validation data splitting
}

#' @export
#' @rdname format_result
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

  if (nrow(Ftest) != nrow(Xobs)) {
    stop(sprintf("Mismatch in rows: feature rows (Ftest=%d) must match observed trial rows (Xobs=%d).",
                 nrow(Ftest), nrow(Xobs)))
  }
  
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
    c(perf$pattern_correlation,
      perf$pattern_discrimination,
      perf$pattern_rank_percentile,
      perf$rdm_correlation,
      perf$voxel_correlation,
      perf$mse,
      perf$r_squared,
      perf$mean_voxelwise_temporal_cor,
      ncomp_used),
    nrow = 1,
    ncol = 9,
    dimnames = list(NULL, c("pattern_correlation", "pattern_discrimination",
                            "pattern_rank_percentile", "rdm_correlation",
                            "voxel_correlation", "mse", "r_squared",
                            "mean_voxelwise_temporal_cor", "ncomp"))
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

#' @rdname merge_results-methods
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
    perf$pattern_correlation,
    perf$pattern_discrimination,
    perf$pattern_rank_percentile,
    perf$rdm_correlation,
    perf$voxel_correlation,
    perf$mse,
    perf$r_squared,
    perf$mean_voxelwise_temporal_cor,
    ncomp_used
  )
  base_names <- c(
    "pattern_correlation", "pattern_discrimination", "pattern_rank_percentile",
    "rdm_correlation",
    "voxel_correlation", "mse", "r_squared",
    "mean_voxelwise_temporal_cor",
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
  
  tibble::tibble(
    result      = list(NULL),
    indices     = list(indices),
    performance = list(perf_mat),
    id          = id,
    error       = FALSE,
    error_message = "~"
  )
}

#' @rdname fit_roi
#' @method fit_roi feature_rsa_model
#' @export
fit_roi.feature_rsa_model <- function(model, roi_data, context, ...) {
  roi <- list(train_roi = roi_data$train_roi, test_roi = roi_data$test_roi)
  id <- context$id
  center_global_id <- context$center_global_id %||% NA

  result_tbl <- internal_crossval(model, roi, id,
                                  center_global_id = center_global_id)

  if (isTRUE(result_tbl$error[1])) {
    roi_result(
      metrics = NULL,
      indices = roi_data$indices,
      id = id,
      result = result_tbl$result[[1]],
      error = TRUE,
      error_message = result_tbl$error_message[1]
    )
  } else {
    roi_result(
      metrics = result_tbl$performance[[1]],
      indices = roi_data$indices,
      id = id,
      result = result_tbl$result[[1]]
    )
  }
}


#' @rdname output_schema
#' @method output_schema feature_rsa_model
#' @keywords internal
#' @export
output_schema.feature_rsa_model <- function(model) {
  # When nperm > 0, permutation metrics are added dynamically; fall back to combine_standard.
  if (!is.null(model$nperm) && model$nperm > 0) {
    return(NULL)
  }
  nms <- c(
    "pattern_correlation", "pattern_discrimination", "pattern_rank_percentile",
    "rdm_correlation", "voxel_correlation", "mse", "r_squared",
    "mean_voxelwise_temporal_cor", "ncomp"
  )
  setNames(rep("scalar", length(nms)), nms)
}


#' Summary Method for Feature RSA Model
#'
#' @param object The feature RSA model
#' @param ... Additional args
#' @return A list of summary statistics for the feature RSA model (printed as side effect).
#' @examples
#' \dontrun{
#'   mdl <- feature_rsa_model(dataset, des)
#'   summary(mdl)
#' }
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
#' @return Invisibly returns the input object \code{x} (called for side effects).
#' @examples
#' \donttest{
#'   S <- as.matrix(dist(matrix(rnorm(5*3), 5, 3)))
#'   labels <- factor(letters[1:5])
#'   des <- feature_rsa_design(S = S, labels = labels)
#'   print(des)
#' }
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
#' @return Invisibly returns the input object \code{x} (called for side effects).
#' @examples
#' \dontrun{
#'   mdl <- feature_rsa_model(dataset, des)
#'   print(mdl)
#' }
#' @export
print.feature_rsa_model <- function(x, ...) {
  # Create a border line for styling
  border <- crayon::bold(crayon::cyan(strrep("=", 50)))
  
  # Header
  cat(border, "\n")
  cat(crayon::bold(crayon::cyan("          Feature RSA Model           \n")))
  cat(border, "\n\n")
  
  # Display the method used (e.g., pls, pca, or glmnet)
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

  if (x$method %in% c("pls", "pca")) {
    sel <- if (!is.null(x$ncomp_selection)) x$ncomp_selection else "max"
    cat(crayon::bold(crayon::blue("Component selection:    ")), sel, "\n")
    if (sel == "pve" && !is.null(x$pve_threshold)) {
      cat(crayon::bold(crayon::blue("PVE threshold:          ")), x$pve_threshold, "\n")
    }
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
