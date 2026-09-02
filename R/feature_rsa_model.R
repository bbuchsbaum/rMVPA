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
    kept_vals <- vals[seq_len(k)]
    kept_scales <- sqrt(pmax(kept_vals, 0))
    F <- sweep(vecs[, seq_len(k), drop = FALSE], 2L, kept_scales, "*")
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
#'     \item{pls}{Partial Least Squares regression predicting X from F using
#'       the numerical algorithm configured by \code{pls::pls.options()}.}
#'     \item{pca}{Principal Component Regression predicting X from PCs of F
#'       using the algorithm configured by \code{pls::pls.options()} (SVD-PCR
#'       by default).}
#'     \item{ridge}{Multi-response ridge regression predicting X from F with
#'       an economy-SVD solver and a normalized penalty.}
#'     \item{glmnet}{Elastic net regression predicting X from F using glmnet with multivariate Gaussian response.}
#'   }
#' @param crossval Optional cross-validation specification.
#' @param ncomp_selection Character string controlling how the number of components
#'   is chosen for \code{pls} and \code{pca} methods.  One of:
#'   \describe{
#'     \item{loo}{(Default) Use leave-one-observation-out validation and select
#'       the fewest components within one standard error of the minimum
#'       segment-wise MSE. This preserves the historical selection rule while
#'       streaming held-out errors instead of retaining a validation cube.}
#'     \item{pve}{Keep the fewest components whose cumulative explained
#'       variance reaches \code{pve_threshold} of the total explained by all
#'       fitted components.}
#'     \item{blocked}{Use leave-one-block-out validation within each outer
#'       training fold and apply the same one-standard-error rule as
#'       \code{"loo"}. Each held-out block contributes one segment MSE.
#'       Requires \code{design$block_var} and at least two training blocks in
#'       every outer fold.}
#'     \item{max}{Use all \code{max_comps} components (legacy behaviour).}
#'   }
#'   Ignored when \code{method} is \code{"ridge"} or \code{"glmnet"}.
#' @param ncomp_objective Character string naming the component-count tuning
#'   estimand for \code{pls} and \code{pca}. \code{"mse"} remains the default
#'   for backward compatibility. \code{"pattern_discrimination"} maximizes the
#'   held-out correct-pattern correlation minus the mean correlation with
#'   incorrect patterns. \code{"pattern_rank_percentile"} maximizes held-out
#'   identification rank. Both pattern-relative objectives require
#'   \code{ncomp_selection = "blocked"} with at least two observations per
#'   inner validation block. Ignored for other methods.
#' @param ncomp_one_se Optional logical controlling the one-standard-error rule
#'   for component-count tuning. \code{NULL} uses \code{TRUE} for MSE-based LOO
#'   or blocked tuning and \code{FALSE} otherwise. Pattern-relative tuning
#'   defaults to the empirical optimum because a simpler component model is not
#'   automatically preferable within one standard error on these objectives.
#'   Ignored for other methods.
#' @param pve_threshold Numeric in (0, 1].  When \code{ncomp_selection = "pve"},
#'   the proportion of total explained X-variance at which to stop adding
#'   components.  Default 0.9.
#' @param lambda_selection Character string controlling ridge penalty selection.
#'   \code{"gcv"} (default) minimizes generalized cross-validation error from
#'   one full-fold SVD. \code{"loo"} uses exact analytic leave-one-observation-out
#'   PRESS errors. \code{"blocked"} uses
#'   leave-one-block-out errors, with centering and scaling re-estimated inside
#'   every inner split; it requires
#'   \code{design$block_var}. \code{"fixed"} requires one supplied \code{lambda}.
#'   Selection applies the rule configured by \code{lambda_one_se}. Ignored for
#'   other methods.
#' @param lambda_objective Character string naming the ridge tuning estimand.
#'   \code{"mse"} (default) preserves GCV/LOO/blocked prediction-error tuning.
#'   \code{"pattern_discrimination"} maximizes the held-out correct-pattern
#'   correlation minus the mean correlation with incorrect patterns.
#'   \code{"pattern_rank_percentile"} maximizes identification rank among
#'   observations withheld together. Both pattern-relative objectives require
#'   \code{lambda_selection = "blocked"} with at least two observations per
#'   inner validation block. Ignored for other methods.
#' @param lambda_one_se Optional logical controlling the one-standard-error
#'   rule for ridge tuning. \code{NULL} uses \code{TRUE} for MSE-based LOO or
#'   blocked tuning and \code{FALSE} otherwise. Pattern-relative tuning defaults
#'   to the empirical optimum because stronger shrinkage is not inherently the
#'   safer error for a relative-pattern objective.
#' @param alpha Numeric value between 0 and 1, only used when method="glmnet". Controls the elastic net
#'   mixing parameter: 1 for lasso (default), 0 for ridge, values in between for a mixture.
#'   Defaults to 0.5 (equal mix of ridge and lasso).
#' @param cv_glmnet Logical, if TRUE and method="glmnet", use cv.glmnet to automatically select the
#'   optimal lambda value via cross-validation. Defaults to FALSE.
#' @param lambda Optional numeric value or sequence of regularization values.
#'   For \code{method = "ridge"}, values must be finite and non-negative and
#'   correspond to the normalized objective
#'   \eqn{n^{-1} ||Y-XB||_F^2 + \lambda ||B||_F^2}; \code{NULL} uses a fixed,
#'   data-independent grid. For \code{method = "glmnet"} with
#'   \code{cv_glmnet = FALSE}, values must be positive; \code{NULL} lets glmnet
#'   construct its path.
#' @param nperm Integer, number of permutations to run for statistical testing of model performance
#'   metrics after merging cross-validation folds. Default 0 (no permutation testing).
#' @param permute_by DEPRECATED. Permutation is always done by shuffling rows of the predicted matrix.
#' @param save_distributions Logical, if TRUE and nperm > 0, save the full null distributions
#'   from the permutation test. Defaults to FALSE.
#' @param return_rdm_vectors Logical; if TRUE, retain each ROI's predicted
#'   lower-triangle RDM vector in the regional result's `fits` slot. Cross-fold
#'   pairs are stored as missing because the two observations were not withheld
#'   together; the retained diagnostics include the observation order and fold
#'   assignment. This is off by default because it can add substantial memory
#'   use for long time series or many ROIs.
#' @param return_predictions Logical; if TRUE, retain each ROI's out-of-fold
#'   predicted patterns (`Yhat`) together with the matching observed patterns,
#'   observation order, outer-fold id, and voxel indices. These are stored in
#'   the regional result's `fits` slot and extracted with
#'   \code{\link{feature_rsa_predictions}}. Off by default because
#'   `n_obs x n_voxels` matrices across many ROIs can be large.
#' @param max_retained_mb Allocation contract (MiB) for retained out-of-fold
#'   predicted and observed patterns. The estimate counts both matrices for a
#'   partition of the active mask (one copy of each voxel). A refusal is
#'   preferred to a silent out-of-memory failure.
#' @param prediction_overflow What to do when the estimate exceeds
#'   \code{max_retained_mb}: \code{"error"} refuses the request;
#'   \code{"none"} disables prediction retention and records a notice.
#' @param feature_standardize How the feature matrix \code{F} is standardized
#'   inside each training fold: \code{"scale"} centers every column and
#'   divides it by its training-fold standard deviation; \code{"center"} only
#'   subtracts the training-fold column means. \code{NULL} (the default)
#'   resolves to \code{"center"} for designs built from a similarity matrix
#'   \code{S}, whose feature matrix consists of PCA-like scores, and to
#'   \code{"scale"} otherwise. Use \code{"center"} whenever the column
#'   variances of \code{F} are meaningful, e.g. PCA scores or other
#'   pre-whitened inputs. See the Feature standardization section.
#' @param ... Additional arguments (currently unused). Passing deprecated
#'   arguments such as \code{cache_pca} now results in an error.
#'
#' @section Feature standardization:
#' Inside every training fold the neural responses \code{X} are always
#' z-scored column-wise, and by default the feature matrix \code{F} is too
#' (\code{feature_standardize = "scale"}). Held-out rows are transformed with
#' the training-fold means and scales, and predictions are returned on the
#' original response scale. This differs from the default of
#' \code{pls::pcr()} (\code{scale = NULL}), so \code{method = "pca"} is
#' correlation-PCR rather than covariance-PCR.
#'
#' Column scaling discards the variance profile of \code{F}. That is harmless
#' for raw feature spaces but degenerate for inputs whose column variances
#' carry the structure, such as PCA scores or any whitened matrix: after
#' scaling every direction has unit variance, PCA component ordering is set by
#' noise, and \code{F = U \%*\% S} gives exactly the same fit as
#' \code{F = U}. Use \code{feature_standardize = "center"} for such inputs; it
#' centers \code{F} per training fold without rescaling columns, and applies
#' to the PLS/PCA, ridge, and glmnet paths alike, including nested blocked
#' tuning. Designs built from a similarity matrix \code{S} store exactly such
#' scores (eigenvectors weighted by the square roots of the eigenvalues), so
#' they default to \code{"center"}.
#'
#' For \code{method = "pca"} under column scaling, the constructor inspects
#' the correlation spectrum of \code{F} and warns when it is flat (condition
#' number within two percent of one), the signature of columns that are
#' exactly orthogonal on these rows, such as PCA scores computed from them.
#' Scores computed on a superset of rows, or otherwise nearly whitened inputs,
#' are not detected, so use \code{"center"} for any PCA-score input whether or
#' not a warning appears. The diagnostic is stored in
#' \code{model$feature_spectrum}.
#'
#' @return A \code{feature_rsa_model} object (S3 class).
#'
#' @seealso \code{\link{feature_rsa_predictions}}, \code{\link{feature_rsa_rdm_vectors}}
#'
#' @examples
#' \donttest{
#'   set.seed(79)
#'   sample <- gen_sample_dataset(c(4, 4, 4), nobs = 24, blocks = 3)
#'   Fmat <- matrix(rnorm(24 * 6), 24, 6)
#'   des <- feature_rsa_design(
#'     F = Fmat,
#'     labels = paste0("t", seq_len(24)),
#'     max_comps = 3,
#'     block_var = sample$design$block_var
#'   )
#'   mdl <- feature_rsa_model(
#'     sample$dataset, des, method = "pca",
#'     ncomp_selection = "max",
#'     return_predictions = TRUE
#'   )
#'   region_mask <- neuroim2::NeuroVol(
#'     sample(1:2, length(sample$dataset$mask), replace = TRUE),
#'     neuroim2::space(sample$dataset$mask)
#'   )
#'   res <- run_regional(mdl, region_mask)
#'   preds <- feature_rsa_predictions(res)
#'   dim(preds$predicted[[1]])
#' }
#' @details
#' Feature RSA models analyze how well a feature matrix \code{F} (defined in the `design`)
#' relates to neural data \code{X}. The `max_comps` parameter, inherited from the `design` object,
#' sets an upper limit on the number of components fitted:
#'   - \strong{pls}: PLS regression using the configured \pkg{pls} numerical
#'     kernel. Fits up to `max_comps` components; the actual number used for
#'     prediction is chosen by \code{ncomp_selection}.
#'   - \strong{pca}: Principal Component Regression using the configured
#'     \pkg{pls} PCR kernel (SVD-PCR by default). Fits up to `max_comps`
#'     components; selection is controlled by \code{ncomp_selection}.
#'   - \strong{ridge}: Multi-response ridge regression using one economy SVD
#'     per training matrix. This is a distinct estimator, not an approximation
#'     to PLS or elastic net. Its penalty is selected by \code{lambda_selection}.
#'   - \strong{glmnet}: Elastic net regression via \code{glmnet} with multivariate Gaussian
#'     response. Regularisation (lambda) can be auto-selected via \code{cv_glmnet=TRUE}.
#'
#' For \code{pls} and \code{pca}, the \code{ncomp_selection} argument determines how many
#' of the fitted components are actually used for prediction.  The default
#' (\code{"loo"}) uses leave-one-observation-out validation and picks the
#' fewest components within one SE of the minimum segment-wise MSE. The
#' \code{ncomp_objective = "mse"} default is retained for backward
#' compatibility. With \code{ncomp_selection = "blocked"}, callers may instead
#' maximize \code{"pattern_discrimination"} or
#' \code{"pattern_rank_percentile"}; these relative objectives default to the
#' empirical optimum rather than the one-SE rule. \code{"blocked"} uses
#' leave-one-block-out segments within each outer training fold, and centering
#' and scaling are estimated again from each inner training split. It is
#' usually the more faithful and much less expensive validation unit when
#' observations are dependent within acquisition runs, sessions, or subjects.
#' It is not interchangeable with \code{"pve"} or \code{"max"}, which do not
#' estimate held-out performance for component selection.
#'
#' **Performance Metrics** (computed by `evaluate_model` after cross-validation):
#'
#' *Condition-pattern metrics* (trial x trial correlation matrix):
#'   - `pattern_correlation`: Average correlation between the predicted and observed
#'     spatial patterns for corresponding trials (diagonal of the trial x trial
#'     correlation matrix computed across voxels).
#'   - `pattern_discrimination`: `pattern_correlation` minus the mean
#'     off-diagonal correlation among candidates withheld in the same outer
#'     fold. Measures how much better the model predicts the correct trial's
#'     pattern than eligible incorrect trials.
#'   - `pattern_rank_percentile`: For each trial, percentile rank of the correct
#'     pattern match among candidates withheld in the same outer fold.
#'     0.5 = chance, 1 = perfect.
#'
#' *Representational geometry*:
#'   - `rdm_correlation`: Spearman correlation between jointly withheld pairs
#'     in the observed and predicted RDMs (defined as 1 - trial-by-trial
#'     correlation across voxels). Captures similarity of held-out
#'     representational geometry without comparing a prediction with a target
#'     that trained that prediction.
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
#' For PLS, PCR, and glmnet, the number of components or its historical proxy
#' (`ncomp`) is included in the performance output. Ridge instead reports the
#' median selected penalty (`median_lambda`), mean effective degrees of freedom
#' (`mean_effective_df`), non-intercept degrees of freedom, and fractions of
#' folds selected at either end of the lambda grid.
#'
#' **Out-of-fold predictions** (`return_predictions = TRUE`):
#' Each ROI retains the merged out-of-fold `Yhat` and `Y` matrices, the
#' observation order, and the outer-fold id used by the built-in scoring
#' rules. Extract them with \code{\link{feature_rsa_predictions}} for
#' post-hoc scoring (temporal windowing, whitened distances, custom
#' identification rules) without crossing fold boundaries. This is a
#' different payload from `return_rdm_vectors` and from the classification
#' `prediction_table`. Overlapping searchlights are refused because they
#' would retain multiple copies of each voxel.
#'
#' @export
feature_rsa_model <- function(dataset,
                               design,
                               method = c("pls", "pca", "glmnet", "ridge"),
                               crossval = NULL,
                               ncomp_selection = c("loo", "max", "pve", "blocked"),
                               pve_threshold = 0.9,
                               alpha = 0.5,
                               cv_glmnet = FALSE,
                               lambda = NULL,
                               nperm = 0,
                               permute_by = c("features", "observations"),
                               save_distributions = FALSE,
                               return_rdm_vectors = FALSE,
                               lambda_selection = c("gcv", "loo", "blocked", "fixed"),
                               lambda_objective = c("mse", "pattern_discrimination", "pattern_rank_percentile"),
                               lambda_one_se = NULL,
                               ncomp_objective = c("mse", "pattern_discrimination", "pattern_rank_percentile"),
                               ncomp_one_se = NULL,
                               return_predictions = FALSE,
                               max_retained_mb = 1024,
                               prediction_overflow = c("error", "none"),
                               feature_standardize = NULL,
                               ...) {
  
  method <- match.arg(method)
  ncomp_selection <- match.arg(ncomp_selection)
  ncomp_objective <- match.arg(ncomp_objective)
  lambda_selection <- match.arg(lambda_selection)
  lambda_objective <- match.arg(lambda_objective)
  permute_by <- match.arg(permute_by)
  assertthat::assert_that(inherits(dataset, "mvpa_dataset"))
  assertthat::assert_that(inherits(design, "feature_rsa_design"))
  # A design built from a similarity matrix S holds PCA-like scores, whose
  # column variances are the eigenvalues; scaling them away is exactly the
  # degenerate case, so such designs default to centering only.
  if (is.null(feature_standardize)) {
    feature_standardize <- if (!is.null(design$S)) "center" else "scale"
  }
  feature_standardize <- match.arg(feature_standardize, c("scale", "center"))
  extra_args <- list(...)
  if ("cache_pca" %in% names(extra_args)) {
    stop("`cache_pca` is no longer supported. PCA caching was removed; use `ncomp_selection` controls for `method='pca'`.")
  }

  component_method <- method %in% c("pls", "pca")
  if (component_method && ncomp_selection == "pve") {
    assertthat::assert_that(is.numeric(pve_threshold) && pve_threshold > 0 && pve_threshold <= 1,
                           msg = "pve_threshold must be in (0, 1]")
  }
  if (component_method && ncomp_selection == "blocked" &&
      is.null(design$block_var)) {
    stop("ncomp_selection='blocked' requires design$block_var.", call. = FALSE)
  }
  if (component_method) {
    if (!is.null(ncomp_one_se) &&
        (!is.logical(ncomp_one_se) || length(ncomp_one_se) != 1L ||
         is.na(ncomp_one_se))) {
      stop("ncomp_one_se must be NULL, TRUE, or FALSE.", call. = FALSE)
    }
    if (ncomp_objective != "mse" && ncomp_selection != "blocked") {
      stop(
        paste0(
          "ncomp_objective='", ncomp_objective, "' requires ",
          "ncomp_selection='blocked' so candidates are withheld together."
        ),
        call. = FALSE
      )
    }
    if (isTRUE(ncomp_one_se) &&
        !ncomp_selection %in% c("loo", "blocked")) {
      stop(
        "The one-SE rule requires ncomp_selection='loo' or 'blocked'.",
        call. = FALSE
      )
    }
    if (is.null(ncomp_one_se)) {
      ncomp_one_se <-
        ncomp_objective == "mse" &&
        ncomp_selection %in% c("loo", "blocked")
    }
  }
  if (method == "ridge") {
    ridge_lambdas <- .feature_rsa_ridge_lambdas(lambda)
    if (!is.null(lambda_one_se) &&
        (!is.logical(lambda_one_se) || length(lambda_one_se) != 1L ||
         is.na(lambda_one_se))) {
      stop("lambda_one_se must be NULL, TRUE, or FALSE.", call. = FALSE)
    }
    if (lambda_selection == "fixed" && length(ridge_lambdas) != 1L) {
      stop(
        "lambda_selection='fixed' requires exactly one finite, non-negative lambda.",
        call. = FALSE
      )
    }
    if (lambda_selection == "blocked" && is.null(design$block_var)) {
      stop("lambda_selection='blocked' requires design$block_var.", call. = FALSE)
    }
    if (lambda_objective != "mse" &&
        lambda_selection != "blocked") {
      stop(
        paste0(
          "lambda_objective='", lambda_objective, "' requires ",
          "lambda_selection='blocked' so candidates are withheld together."
        ),
        call. = FALSE
      )
    }
    if (isTRUE(lambda_one_se) &&
        !lambda_selection %in% c("loo", "blocked")) {
      stop(
        "The one-SE rule requires lambda_selection='loo' or 'blocked'.",
        call. = FALSE
      )
    }
    if (is.null(lambda_one_se)) {
      lambda_one_se <-
        lambda_objective == "mse" &&
        lambda_selection %in% c("loo", "blocked")
    }
  }
  assertthat::assert_that(
    is.logical(return_rdm_vectors) && length(return_rdm_vectors) == 1L && !is.na(return_rdm_vectors),
    msg = "return_rdm_vectors must be TRUE or FALSE"
  )
  assertthat::assert_that(
    is.logical(return_predictions) && length(return_predictions) == 1L && !is.na(return_predictions),
    msg = "return_predictions must be TRUE or FALSE"
  )
  prediction_overflow <- tryCatch(
    match.arg(prediction_overflow),
    error = function(...) {
      stop('prediction_overflow must be "error" or "none".', call. = FALSE)
    }
  )
  if (!is.numeric(max_retained_mb) || length(max_retained_mb) != 1L ||
      !is.finite(max_retained_mb) || max_retained_mb <= 0) {
    stop("max_retained_mb must be one positive finite value.", call. = FALSE)
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
  n_obs <- nrow(design$F)
  if (is.null(n_obs) || !is.finite(n_obs) || n_obs < 1L) {
    stop("feature_rsa_model: design must contain a feature matrix with at least one row.",
         call. = FALSE)
  }
  retention <- .feature_rsa_validate_prediction_retention(
    n_obs = n_obs,
    n_voxels = active_voxels,
    return_predictions = return_predictions,
    max_retained_mb = max_retained_mb,
    prediction_overflow = prediction_overflow
  )
  return_predictions <- isTRUE(retention$return_predictions)
  
  model_spec <- create_model_spec(
    "feature_rsa_model", 
    dataset = dataset,
    design  = design, 
    method  = method, 
    crossval= crossval, 
    compute_performance = TRUE,
    return_predictions = return_predictions,
    return_fits = isTRUE(return_rdm_vectors) || return_predictions
  )
  
  # Single "max_comps" in use
  model_spec$max_comps <- max_comps
  model_spec$feature_standardize <- feature_standardize

  # Column scaling discards the variance profile of F. Warn when that leaves
  # PCA/PCR with no dominant direction to order components by.
  spectrum <- NULL
  if (component_method && feature_standardize == "scale") {
    spectrum <- .feature_rsa_spectrum_diagnostic(design$F)
    if (method == "pca" && !is.null(spectrum) && isTRUE(spectrum$flat)) {
      warning(sprintf(paste0(
        "feature_rsa_model: after per-fold standardization the correlation ",
        "spectrum of F is flat (condition number %.3f). Its columns are ",
        "orthogonal with equal variance once scaled, as with PCA scores or ",
        "whitened features, so PCA/PCR component ordering is arbitrary and ",
        "ncomp selection is degenerate. If the column variances of F are ",
        "meaningful, use feature_standardize = \"center\"."
      ), spectrum$condition_number), call. = FALSE)
    }
  }
  model_spec["feature_spectrum"] <- list(spectrum)

  # Component selection strategy (for PLS/PCA)
  model_spec$ncomp_selection <- ncomp_selection
  model_spec$ncomp_objective <- ncomp_objective
  if (component_method) {
    model_spec$ncomp_one_se <- isTRUE(ncomp_one_se)
  }
  model_spec$pve_threshold <- pve_threshold

  # GLMNet parameters
  if (method == "glmnet") {
    model_spec$alpha <- alpha
    model_spec$cv_glmnet <- cv_glmnet
    model_spec$lambda <- lambda
  } else if (method == "ridge") {
    model_spec$lambda_selection <- lambda_selection
    model_spec$lambda_objective <- lambda_objective
    model_spec$lambda_one_se <- isTRUE(lambda_one_se)
    # Preserve an explicit NULL entry. Without exact list assignment, `$lambda`
    # would partially match `lambda_selection` on ridge specifications.
    model_spec["lambda"] <- list(lambda)
  }
  
  model_spec$nperm <- nperm
  model_spec$permute_by <- permute_by
  model_spec$save_distributions <- save_distributions
  model_spec$return_rdm_vectors <- isTRUE(return_rdm_vectors)
  model_spec$max_retained_mb <- max_retained_mb
  model_spec$prediction_overflow <- prediction_overflow
  model_spec$retention_estimated_bytes <- retention$retention_estimated_bytes
  model_spec$retention_notice <- retention$retention_notice
  
  model_spec
}


#' Bytes retained for out-of-fold predicted and observed pattern matrices
#' @noRd
.feature_rsa_prediction_retention_bytes <- function(n_obs,
                                                    n_voxels,
                                                    return_predictions) {
  if (!isTRUE(return_predictions)) {
    return(0)
  }
  8 * as.numeric(n_obs) * as.numeric(n_voxels) * 2
}


#' Enforce the out-of-fold prediction allocation contract
#' @noRd
.feature_rsa_validate_prediction_retention <- function(n_obs,
                                                       n_voxels,
                                                       return_predictions,
                                                       max_retained_mb,
                                                       prediction_overflow = "error") {
  bytes <- .feature_rsa_prediction_retention_bytes(
    n_obs, n_voxels, return_predictions
  )
  notice <- NULL
  if (isTRUE(return_predictions) && bytes > max_retained_mb * 1024^2) {
    if (identical(prediction_overflow, "none")) {
      notice <- paste0(
        "Requested prediction retention was disabled because estimated retained ",
        "storage exceeded ", max_retained_mb, " MiB."
      )
      return_predictions <- FALSE
      bytes <- 0
    } else {
      stop(
        "feature_rsa_model: requested out-of-fold predictions are estimated at ",
        sprintf("%.3f", bytes / 1024^2),
        " MiB, above max_retained_mb = ", max_retained_mb,
        ". Reduce return_predictions or raise the explicit limit.",
        call. = FALSE
      )
    }
  }
  list(
    return_predictions = isTRUE(return_predictions),
    retention_estimated_bytes = c(predictions = bytes, total = bytes),
    retention_notice = notice
  )
}


#' Count voxels that a regional mask will actually analyze
#' @noRd
.feature_rsa_region_voxel_count <- function(dataset, region_mask) {
  region <- .ensure_dense_vol(region_mask)
  region_vals <- as.vector(region)
  n_region <- sum(region_vals > 0, na.rm = TRUE)
  mask <- dataset$mask
  if (is.null(mask)) {
    return(as.integer(n_region))
  }
  mask_vals <- as.vector(.ensure_dense_vol(mask))
  if (length(mask_vals) == length(region_vals)) {
    return(as.integer(sum(region_vals > 0 & mask_vals > 0, na.rm = TRUE)))
  }
  as.integer(n_region)
}


#' Apply a stored prediction-retention contract to a concrete voxel count
#' @noRd
.feature_rsa_apply_prediction_retention <- function(model_spec, n_voxels) {
  if (!isTRUE(model_spec$return_predictions)) {
    return(model_spec)
  }
  max_mb <- model_spec$max_retained_mb
  if (is.null(max_mb) || !is.numeric(max_mb) || !is.finite(max_mb)) {
    return(model_spec)
  }
  n_obs <- nrow(model_spec$design$F)
  overflow <- model_spec$prediction_overflow %||% "error"
  retention <- .feature_rsa_validate_prediction_retention(
    n_obs = n_obs,
    n_voxels = n_voxels,
    return_predictions = TRUE,
    max_retained_mb = max_mb,
    prediction_overflow = overflow
  )
  model_spec$return_predictions <- retention$return_predictions
  model_spec$return_fits <- isTRUE(model_spec$return_rdm_vectors) ||
    isTRUE(retention$return_predictions)
  model_spec$retention_estimated_bytes <- retention$retention_estimated_bytes
  if (!is.null(retention$retention_notice)) {
    model_spec$retention_notice <- retention$retention_notice
  }
  model_spec
}


#' Onesigma component selection for multivariate responses
#'
#' \code{pls::selectNcomp} only supports univariate responses.
#' This helper reconstructs segment-wise CV MSE from the stored validation
#' predictions and applies the standard onesigma rule across CV segments:
#' pick the fewest components whose mean segment MSE is within one SE of the
#' minimum. This avoids the earlier response-wise approximation, which could
#' be overly conservative when the response matrix has many columns.
#'
#' It computes CV errors directly from \code{model$validation$pred} and
#' \code{model$validation$segments}, avoiding fragile internal namespace
#' lookups in \pkg{pls} during parallel worker execution.
#' @noRd
.selectNcomp_mv <- function(model, method = "onesigma") {
  pred <- model$validation$pred
  segments <- model$validation$segments
  y_true <- tryCatch(stats::model.response(stats::model.frame(model)), error = function(...) NULL)

  if (is.null(pred) || is.null(segments) || is.null(y_true)) {
    return(NA_integer_)
  }

  y_true <- as.matrix(y_true)

  pred_dims <- dim(pred)
  if (is.null(pred_dims)) {
    return(NA_integer_)
  }
  if (length(pred_dims) == 2L) {
    pred <- array(pred, dim = c(pred_dims[1], 1L, pred_dims[2]))
    pred_dims <- dim(pred)
  }
  if (length(pred_dims) != 3L) {
    return(NA_integer_)
  }

  nobj <- pred_dims[1]
  nresp <- pred_dims[2]
  ncomp_max <- pred_dims[3]
  if (nobj != nrow(y_true) || nresp != ncol(y_true) || ncomp_max < 1L) {
    return(NA_integer_)
  }

  seg_mse <- matrix(NA_real_, nrow = length(segments), ncol = ncomp_max)
  for (s in seq_along(segments)) {
    idx <- as.integer(segments[[s]])
    idx <- idx[is.finite(idx) & idx >= 1L & idx <= nobj]
    if (!length(idx)) {
      next
    }

    y_seg <- y_true[idx, , drop = FALSE]
    for (k in seq_len(ncomp_max)) {
      pred_seg <- pred[idx, , k, drop = FALSE]
      seg_mse[s, k] <- mean((y_seg - pred_seg[, , 1L])^2, na.rm = TRUE)
    }
  }

  mean_mse <- colMeans(seg_mse, na.rm = TRUE)
  if (all(is.na(mean_mse) | !is.finite(mean_mse))) {
    return(NA_integer_)
  }

  min_idx <- which.min(mean_mse)
  if (length(min_idx) == 0L) {
    return(NA_integer_)
  }
  if (method != "onesigma") {
    return(as.integer(min_idx))
  }

  mse_at_min <- seg_mse[, min_idx]
  n_finite <- sum(is.finite(mse_at_min))
  if (n_finite < 2L) {
    return(as.integer(min_idx))
  }

  se <- stats::sd(mse_at_min, na.rm = TRUE) / sqrt(n_finite)
  thresh <- mean_mse[min_idx] + se
  selected <- which(mean_mse <= thresh)[1]
  if (is.na(selected)) {
    return(as.integer(min_idx))
  }

  as.integer(selected)
}


#' @noRd
.feature_rsa_col_sds <- function(X) {
  X <- as.matrix(X)
  n <- nrow(X)
  p <- ncol(X)

  if (p == 0L) {
    return(numeric(0))
  }
  if (n < 2L) {
    return(rep(NA_real_, p))
  }
  # Do not use E[X^2] - E[X]^2 here. That one-pass identity loses all
  # significant digits for constant or nearly constant non-zero columns and
  # can turn an exact zero variance into an apparent sd around 1e-8.
  matrixStats::colSds(X, na.rm = FALSE)
}

#' @noRd
.feature_rsa_row_sds <- function(X) {
  X <- as.matrix(X)
  n <- nrow(X)
  p <- ncol(X)

  if (n == 0L) {
    return(numeric(0))
  }
  if (p < 2L) {
    return(rep(NA_real_, n))
  }
  matrixStats::rowSds(X, na.rm = FALSE)
}

#' Correlate matrix rows with a no-missing-data BLAS fast path
#'
#' Rows are variables and columns are paired observations, matching
#' `stats::cor(t(X), t(Y), use = "pairwise.complete.obs")`. Missing data use
#' that reference implementation directly. Otherwise centering and one
#' `tcrossprod()` avoid the much slower generic pairwise-complete path.
#' @noRd
.feature_rsa_row_cor <- function(X, Y = NULL) {
  X <- as.matrix(X)
  same_matrix <- is.null(Y)
  if (same_matrix) {
    Y <- X
  } else {
    Y <- as.matrix(Y)
  }

  if (ncol(X) != ncol(Y)) {
    stop("feature RSA row correlation: matrices must have the same columns.",
         call. = FALSE)
  }
  if (anyNA(X) || anyNA(Y)) {
    return(stats::cor(
      t(X),
      if (same_matrix) NULL else t(Y),
      use = "pairwise.complete.obs"
    ))
  }
  if (ncol(X) < 2L) {
    out <- matrix(NA_real_, nrow = nrow(X), ncol = nrow(Y))
    if (!is.null(rownames(X)) || !is.null(rownames(Y))) {
      dimnames(out) <- list(rownames(X), rownames(Y))
    }
    return(out)
  }

  X_centered <- sweep(X, 1L, rowMeans(X), "-")
  X_norm <- sqrt(rowSums(X_centered * X_centered))
  X_valid <- is.finite(X_norm) & X_norm > 0
  if (any(X_valid)) {
    X_centered[X_valid, ] <- sweep(
      X_centered[X_valid, , drop = FALSE],
      1L,
      X_norm[X_valid],
      "/"
    )
  }
  X_centered[!X_valid, ] <- 0
  if (same_matrix) {
    Y_centered <- X_centered
    Y_norm <- X_norm
    Y_valid <- X_valid
  } else {
    Y_centered <- sweep(Y, 1L, rowMeans(Y), "-")
    Y_norm <- sqrt(rowSums(Y_centered * Y_centered))
    Y_valid <- is.finite(Y_norm) & Y_norm > 0
    if (any(Y_valid)) {
      Y_centered[Y_valid, ] <- sweep(
        Y_centered[Y_valid, , drop = FALSE],
        1L,
        Y_norm[Y_valid],
        "/"
      )
    }
    Y_centered[!Y_valid, ] <- 0
  }

  out <- tcrossprod(X_centered, Y_centered)
  if (any(!X_valid)) out[!X_valid, ] <- NA_real_
  if (any(!Y_valid)) out[, !Y_valid] <- NA_real_
  if (!is.null(rownames(X)) || !is.null(rownames(Y))) {
    dimnames(out) <- list(rownames(X), rownames(Y))
  }
  out
}


#' Normalize rows for correlation kernels
#' @noRd
.feature_rsa_normalize_rows <- function(X) {
  X <- as.matrix(X)
  if (anyNA(X)) {
    return(NULL)
  }
  centered <- sweep(X, 1L, rowMeans(X), "-")
  norms <- sqrt(rowSums(centered * centered))
  valid <- is.finite(norms) & norms > 0
  if (any(valid)) {
    centered[valid, ] <- sweep(
      centered[valid, , drop = FALSE],
      1L,
      norms[valid],
      "/"
    )
  }
  centered[!valid, ] <- 0
  list(values = centered, valid = valid)
}


#' Linear-time pattern-discrimination score
#'
#' The mean of every off-diagonal cross-correlation can be recovered from
#' column sums of row-normalized matrices. This computes the same
#' diagonal-minus-off-diagonal estimand as the dense correlation matrix without
#' allocating an observations-by-observations object.
#' @noRd
.feature_rsa_pattern_discrimination_linear <- function(predicted, observed) {
  predicted <- as.matrix(predicted)
  observed <- as.matrix(observed)
  if (!identical(dim(predicted), dim(observed))) {
    stop(
      "feature RSA pattern discrimination requires matching matrix dimensions.",
      call. = FALSE
    )
  }
  if (nrow(predicted) < 2L ||
      any(!is.finite(predicted)) || any(!is.finite(observed))) {
    return(NA_real_)
  }

  pred_norm <- .feature_rsa_normalize_rows(predicted)
  obs_norm <- .feature_rsa_normalize_rows(observed)
  pred_valid <- pred_norm$valid
  obs_valid <- obs_norm$valid
  diagonal_valid <- pred_valid & obs_valid
  diagonal_count <- sum(diagonal_valid)
  off_count <- sum(pred_valid) * sum(obs_valid) - diagonal_count
  if (diagonal_count < 1L || off_count < 1L) {
    return(NA_real_)
  }

  diagonal_values <- rowSums(
    pred_norm$values[diagonal_valid, , drop = FALSE] *
      obs_norm$values[diagonal_valid, , drop = FALSE]
  )
  all_sum <- sum(
    colSums(pred_norm$values[pred_valid, , drop = FALSE]) *
      colSums(obs_norm$values[obs_valid, , drop = FALSE])
  )
  diagonal_sum <- sum(diagonal_values)
  diagonal_sum / diagonal_count -
    (all_sum - diagonal_sum) / off_count
}


#' Score a held-out pattern-relative objective
#' @noRd
.feature_rsa_pattern_objective <- function(predicted,
                                           observed,
                                           objective) {
  objective <- match.arg(
    objective,
    c("pattern_discrimination", "pattern_rank_percentile")
  )
  if (objective == "pattern_discrimination") {
    .feature_rsa_pattern_discrimination_linear(predicted, observed)
  } else {
    .feature_rsa_pattern_metrics_blockwise(predicted, observed)[["rank"]]
  }
}


#' Pattern metrics without retaining an observations-by-observations matrix
#' @noRd
.feature_rsa_resolve_fold_id <- function(fold_id, n) {
  if (is.null(fold_id)) {
    return(rep.int(1L, n))
  }
  if (!is.atomic(fold_id) || length(fold_id) != n) {
    stop(sprintf("feature RSA fold_id must have length %d.", n), call. = FALSE)
  }
  if (anyNA(fold_id)) {
    stop("feature RSA fold_id cannot contain missing values.", call. = FALSE)
  }
  match(fold_id, unique(fold_id))
}


#' Fold-preserving row permutation
#' @noRd
.feature_rsa_permutation_index <- function(fold_id) {
  fold_id <- .feature_rsa_resolve_fold_id(fold_id, length(fold_id))
  out <- seq_along(fold_id)
  for (rows in unname(split(out, fold_id))) {
    # sample(5) samples from 1:5, so singleton groups need an explicit guard.
    out[rows] <- if (length(rows) > 1L) sample(rows) else rows
  }
  out
}


#' Pattern metrics without retaining an observations-by-observations matrix
#' @noRd
.feature_rsa_pattern_metrics_blockwise <- function(predicted,
                                                   observed,
                                                   fold_id = NULL) {
  predicted <- as.matrix(predicted)
  observed <- as.matrix(observed)
  if (!identical(dim(predicted), dim(observed)) || nrow(predicted) < 2L) {
    return(c(correlation = NA_real_, discrimination = NA_real_, rank = NA_real_))
  }
  n <- nrow(predicted)
  fold_id <- .feature_rsa_resolve_fold_id(fold_id, n)
  groups <- unname(split(seq_len(n), fold_id))
  pred_norm <- .feature_rsa_normalize_rows(predicted)
  obs_norm <- .feature_rsa_normalize_rows(observed)
  if (is.null(pred_norm) || is.null(obs_norm)) {
    dense <- .feature_rsa_row_cor(predicted, observed)
    diag_vals <- ranks <- rep(NA_real_, n)
    off_sum <- 0
    off_count <- 0L
    for (group in groups) {
      for (row_index in group) {
        correct <- match(row_index, group)
        row_values <- dense[row_index, group]
        diagonal <- row_values[[correct]]
        diag_vals[[row_index]] <- diagonal
        denominator <- sum(is.finite(row_values)) - 1L
        if (denominator > 0L && is.finite(diagonal)) {
          ranks[[row_index]] <-
            (sum(row_values <= diagonal, na.rm = TRUE) - 1L) /
            denominator
        }
        finite <- is.finite(row_values)
        if (is.finite(diagonal)) finite[[correct]] <- FALSE
        off_sum <- off_sum + sum(row_values[finite])
        off_count <- off_count + sum(finite)
      }
    }
    pattern_cor <- if (any(is.finite(diag_vals))) {
      mean(diag_vals, na.rm = TRUE)
    } else {
      NA_real_
    }
    off_mean <- if (off_count > 0L) off_sum / off_count else NA_real_
    return(c(
      correlation = pattern_cor,
      discrimination = if (is.finite(pattern_cor) && is.finite(off_mean)) {
        pattern_cor - off_mean
      } else {
        NA_real_
      },
      rank = if (any(is.finite(ranks))) mean(ranks, na.rm = TRUE) else NA_real_
    ))
  }

  block_size <- getOption("rMVPA.feature_rsa_metric_block_rows", 128L)
  block_size <- max(1L, min(n, as.integer(block_size)))
  diag_vals <- rep(NA_real_, n)
  ranks <- rep(NA_real_, n)
  off_sum <- 0
  off_count <- 0L

  for (group in groups) {
    group_size <- length(group)
    for (start in seq.int(1L, group_size, by = block_size)) {
      stop_at <- min(group_size, start + block_size - 1L)
      local_rows <- seq.int(start, stop_at)
      rows <- group[local_rows]
      block <- tcrossprod(
        pred_norm$values[rows, , drop = FALSE],
        obs_norm$values[group, , drop = FALSE]
      )
      for (local_row in seq_along(rows)) {
        row_index <- rows[[local_row]]
        correct <- local_rows[[local_row]]
        row_values <- block[local_row, ]
        if (!pred_norm$valid[[row_index]]) row_values[] <- NA_real_
        invalid_observed <- !obs_norm$valid[group]
        if (any(invalid_observed)) row_values[invalid_observed] <- NA_real_
        diagonal <- row_values[[correct]]
        diag_vals[[row_index]] <- diagonal
        denominator <- sum(is.finite(row_values)) - 1L
        if (denominator > 0L && is.finite(diagonal)) {
          ranks[[row_index]] <-
            (sum(row_values <= diagonal, na.rm = TRUE) - 1L) /
            denominator
        }
        finite <- is.finite(row_values)
        if (is.finite(diagonal)) finite[[correct]] <- FALSE
        off_sum <- off_sum + sum(row_values[finite])
        off_count <- off_count + sum(finite)
      }
    }
  }

  pattern_cor <- if (any(is.finite(diag_vals))) {
    mean(diag_vals, na.rm = TRUE)
  } else {
    NA_real_
  }
  off_mean <- if (off_count > 0L) off_sum / off_count else NA_real_
  c(
    correlation = pattern_cor,
    discrimination = if (is.finite(pattern_cor) && is.finite(off_mean)) {
      pattern_cor - off_mean
    } else {
      NA_real_
    },
    rank = if (any(is.finite(ranks))) mean(ranks, na.rm = TRUE) else NA_real_
  )
}


#' Lower-triangle correlation-distance vector using bounded row blocks
#' @noRd
.feature_rsa_rdm_vector_blockwise <- function(X) {
  X <- as.matrix(X)
  n <- nrow(X)
  n_pairs <- n * (n - 1L) / 2L
  if (n_pairs < 1L) return(numeric(0))
  normalized <- .feature_rsa_normalize_rows(X)
  if (is.null(normalized)) {
    dense <- .feature_rsa_row_cor(X)
    return(as.numeric((1 - dense)[lower.tri(dense)]))
  }

  out <- numeric(n_pairs)
  block_size <- getOption("rMVPA.feature_rsa_metric_block_rows", 128L)
  block_size <- max(1L, min(n - 1L, as.integer(block_size)))
  offset <- 0L
  for (start in seq.int(1L, n - 1L, by = block_size)) {
    stop_at <- min(n - 1L, start + block_size - 1L)
    columns <- seq.int(start, stop_at)
    block <- tcrossprod(
      normalized$values[columns, , drop = FALSE],
      normalized$values
    )
    for (local_column in seq_along(columns)) {
      column <- columns[[local_column]]
      rows <- seq.int(column + 1L, n)
      count <- length(rows)
      values <- block[local_column, rows]
      if (!normalized$valid[[column]]) values[] <- NA_real_
      invalid_rows <- !normalized$valid[rows]
      if (any(invalid_rows)) values[invalid_rows] <- NA_real_
      positions <- seq.int(offset + 1L, offset + count)
      out[positions] <- 1 - values
      offset <- offset + count
    }
  }
  out
}


#' Fold-restricted RDM vectors
#'
#' Compact output concatenates only within-fold pairs. Full output preserves
#' the ordinary lower-triangle length and marks cross-fold pairs as missing.
#' @noRd
.feature_rsa_grouped_rdm_vectors <- function(X,
                                             fold_id = NULL,
                                             full_length = FALSE) {
  X <- as.matrix(X)
  n <- nrow(X)
  fold_id <- .feature_rsa_resolve_fold_id(fold_id, n)
  groups <- unname(split(seq_len(n), fold_id))
  n_pairs <- n * (n - 1L) / 2L
  if (n_pairs < 1L) return(numeric(0))

  if (isTRUE(full_length)) {
    out <- rep(NA_real_, n_pairs)
    for (group in groups) {
      if (length(group) < 2L) next
      out[.feature_rsa_pair_indices(group, n)] <-
        .feature_rsa_rdm_vector_blockwise(X[group, , drop = FALSE])
    }
    return(out)
  }

  values <- lapply(groups, function(group) {
    if (length(group) < 2L) return(numeric(0))
    .feature_rsa_rdm_vector_blockwise(X[group, , drop = FALSE])
  })
  unlist(values, use.names = FALSE)
}


#' Map arbitrary row pairs into a full lower-triangle vector
#' @noRd
.feature_rsa_pair_indices <- function(rows, n) {
  rows <- as.integer(rows)
  m <- length(rows)
  if (m < 2L) return(integer(0))
  out <- integer(m * (m - 1L) / 2L)
  offset <- 0L
  for (column in seq_len(m - 1L)) {
    left <- rows[[column]]
    right <- rows[seq.int(column + 1L, m)]
    lower <- pmax.int(left, right)
    upper <- pmin.int(left, right)
    count <- length(right)
    positions <- seq.int(offset + 1L, offset + count)
    out[positions] <-
      (upper - 1L) * n - ((upper - 1L) * upper) %/% 2L + (lower - upper)
    offset <- offset + count
  }
  out
}


#' Full-vector pair indices restricted to jointly held-out folds
#' @noRd
.feature_rsa_grouped_pair_indices <- function(rows, fold_id, n) {
  if (length(rows) != length(fold_id)) {
    stop("feature RSA grouped pair rows and fold_id must have equal length.",
         call. = FALSE)
  }
  fold_id <- .feature_rsa_resolve_fold_id(fold_id, length(rows))
  groups <- unname(split(seq_along(rows), fold_id))
  unlist(lapply(groups, function(group) {
    .feature_rsa_pair_indices(rows[group], n)
  }), use.names = FALSE)
}


#' Pattern metrics from indexed rows of a cached cross-correlation matrix
#' @noRd
.feature_rsa_pattern_metrics_from_cross <- function(cross_cor,
                                                     predicted_rows,
                                                     observed_rows,
                                                     fold_id = NULL) {
  n <- length(predicted_rows)
  if (n < 1L || length(observed_rows) != n) {
    return(c(correlation = NA_real_, discrimination = NA_real_, rank = NA_real_))
  }
  fold_id <- .feature_rsa_resolve_fold_id(fold_id, n)
  groups <- unname(split(seq_len(n), fold_id))
  diag_vals <- ranks <- rep(NA_real_, n)
  off_sum <- 0
  off_count <- 0L
  for (group in groups) {
    for (i in group) {
      correct <- match(i, group)
      row_values <- cross_cor[predicted_rows[[i]], observed_rows[group]]
      diagonal <- row_values[[correct]]
      diag_vals[[i]] <- diagonal
      denominator <- sum(is.finite(row_values)) - 1L
      if (denominator > 0L && is.finite(diagonal)) {
        ranks[[i]] <-
          (sum(row_values <= diagonal, na.rm = TRUE) - 1L) / denominator
      }
      finite <- is.finite(row_values)
      if (is.finite(diagonal)) finite[[correct]] <- FALSE
      off_sum <- off_sum + sum(row_values[finite])
      off_count <- off_count + sum(finite)
    }
  }
  pattern_cor <- if (any(is.finite(diag_vals))) {
    mean(diag_vals, na.rm = TRUE)
  } else {
    NA_real_
  }
  off_mean <- if (off_count > 0L) off_sum / off_count else NA_real_
  c(
    correlation = pattern_cor,
    discrimination = if (is.finite(pattern_cor) && is.finite(off_mean)) {
      pattern_cor - off_mean
    } else {
      NA_real_
    },
    rank = if (any(is.finite(ranks))) mean(ranks, na.rm = TRUE) else NA_real_
  )
}

#' @noRd
.feature_rsa_offdiag_mean <- function(mat, diag_vals = NULL) {
  if (!is.matrix(mat) || nrow(mat) < 2L) {
    return(NA_real_)
  }

  if (is.null(diag_vals)) {
    diag_vals <- diag(mat)
  }

  finite_total <- is.finite(mat)
  finite_diag <- is.finite(diag_vals)
  n_off <- sum(finite_total) - sum(finite_diag)

  if (n_off <= 0L) {
    return(NA_real_)
  }

  (sum(mat[finite_total], na.rm = TRUE) - sum(diag_vals[finite_diag], na.rm = TRUE)) / n_off
}

#' @noRd
.feature_rsa_diag_col_cor <- function(X, Y) {
  X <- as.matrix(X)
  Y <- as.matrix(Y)

  if (!identical(dim(X), dim(Y))) {
    stop(".feature_rsa_diag_col_cor: matrices must have identical dimensions.", call. = FALSE)
  }
  if (ncol(X) == 0L) {
    return(numeric(0))
  }
  if (nrow(X) < 2L) {
    return(rep(NA_real_, ncol(X)))
  }

  if (anyNA(X) || anyNA(Y)) {
    return(vapply(seq_len(ncol(X)), function(j) {
      tryCatch(
        stats::cor(X[, j], Y[, j], use = "pairwise.complete.obs"),
        error = function(e) NA_real_
      )
    }, numeric(1)))
  }

  x_mean <- colMeans(X)
  y_mean <- colMeans(Y)
  x_ctr <- sweep(X, 2L, x_mean, "-")
  y_ctr <- sweep(Y, 2L, y_mean, "-")
  denom <- sqrt(colSums(x_ctr * x_ctr) * colSums(y_ctr * y_ctr))

  out <- rep(NA_real_, ncol(X))
  ok <- is.finite(denom) & denom > 0
  if (any(ok)) {
    out[ok] <- colSums(x_ctr[, ok, drop = FALSE] * y_ctr[, ok, drop = FALSE]) / denom[ok]
  }
  out
}

#' @noRd
.feature_rsa_global_cor <- function(X, Y) {
  X <- as.matrix(X)
  Y <- as.matrix(Y)

  if (!identical(dim(X), dim(Y))) {
    stop(".feature_rsa_global_cor: matrices must have identical dimensions.", call. = FALSE)
  }
  if (length(X) < 2L) {
    return(NA_real_)
  }

  if (anyNA(X) || anyNA(Y)) {
    return(tryCatch(
      stats::cor(as.vector(X), as.vector(Y), use = "pairwise.complete.obs"),
      error = function(e) NA_real_
    ))
  }

  x_ctr <- X - mean(X)
  y_ctr <- Y - mean(Y)
  denom <- sqrt(sum(x_ctr * x_ctr) * sum(y_ctr * y_ctr))

  if (!is.finite(denom) || denom <= 0) {
    return(NA_real_)
  }

  sum(x_ctr * y_ctr) / denom
}

#' @noRd
.feature_rsa_rdm_vector_from_cor <- function(cor_mat, n_obs) {
  n_pairs <- n_obs * (n_obs - 1L) / 2L

  if (n_pairs < 1L) {
    return(numeric(0))
  }
  if (is.null(cor_mat) || !is.matrix(cor_mat) || !all(dim(cor_mat) == c(n_obs, n_obs))) {
    return(rep(NA_real_, n_pairs))
  }

  as.numeric((1 - cor_mat)[lower.tri(cor_mat)])
}

#' @noRd
.standardize <- function(X, method = c("scale", "center")) {
  method <- match.arg(method)
  cm <- colMeans(X)
  if (method == "scale") {
    csd <- .feature_rsa_col_sds(X)
    csd[csd == 0] <- 1
  } else {
    csd <- stats::setNames(rep(1, ncol(X)), colnames(X))
  }
  X_sc <- scale(X, center=cm, scale=csd)
  list(X_sc = X_sc, mean = cm, sd = csd)
}


#' Is the correlation spectrum of F flat once its columns are scaled?
#'
#' Eigenvalues of the column-correlation matrix of F, accumulated over row
#' chunks so no full-size copy is made. A condition number of one means every
#' direction carries the same variance after column scaling, so PCA/PCR
#' component ordering is arbitrary. Exactly orthogonal inputs such as PCA
#' scores computed on these rows, or orthogonal equal-variance designs, give
#' a condition number of one; noise and structured feature spaces do not, and
#' PCA scores computed on a superset of rows are not detected. Columns are
#' thinned deterministically beyond `max_cols`, which preserves exact
#' orthogonality; rows are never thinned because that would not.
#' @noRd
.feature_rsa_spectrum_diagnostic <- function(F,
                                             max_cols = 1000L,
                                             chunk_rows = 4096L,
                                             flat_condition = 1.02) {
  F <- as.matrix(F)
  n <- nrow(F)
  p <- ncol(F)
  if (is.null(n) || n < 3L || p < 2L) return(NULL)
  cols <- if (p > max_cols) {
    unique(as.integer(round(seq(1, p, length.out = max_cols))))
  } else seq_len(p)
  chunks <- split(seq_len(n), ceiling(seq_len(n) / as.integer(chunk_rows)))
  sums <- numeric(length(cols))
  for (idx in chunks) {
    block <- F[idx, cols, drop = FALSE]
    if (any(!is.finite(block))) return(NULL)
    sums <- sums + colSums(block)
  }
  means <- sums / n
  C <- matrix(0, length(cols), length(cols))
  for (idx in chunks) {
    C <- C + crossprod(sweep(F[idx, cols, drop = FALSE], 2L, means, "-"))
  }
  d <- diag(C)
  keep <- is.finite(d) & d > 0
  if (sum(keep) < 2L) return(NULL)
  C <- C[keep, keep, drop = FALSE]
  s <- 1 / sqrt(diag(C))
  R <- C * tcrossprod(s)
  lambda <- eigen(R, symmetric = TRUE, only.values = TRUE)$values
  lam_max <- lambda[[1L]]
  lam_min <- max(lambda[[length(lambda)]], 0)
  condition <- if (lam_min > 0) lam_max / lam_min else Inf
  list(
    largest_eigenvalue = lam_max,
    smallest_eigenvalue = lam_min,
    condition_number = condition,
    flat_condition = flat_condition,
    n_rows = n,
    n_cols = ncol(C),
    flat = is.finite(condition) && condition <= flat_condition
  )
}


#' Resolve the low-level pls fitting kernel used by feature RSA
#' @noRd
.feature_rsa_kernel_function <- function(method) {
  method <- match.arg(method, c("pls", "pca"))
  opts <- pls::pls.options()
  algorithm <- if (identical(method, "pca")) opts$pcralg else opts$plsralg

  fit_name <- switch(
    algorithm,
    kernelpls = "kernelpls.fit",
    widekernelpls = "widekernelpls.fit",
    simpls = "simpls.fit",
    oscorespls = "oscorespls.fit",
    nipalspls = "nipals.fit",
    svdpc = "svdpc.fit",
    NULL
  )
  if (is.null(fit_name)) {
    stop(
      sprintf(
        "feature RSA does not support pls algorithm '%s' in the matrix kernel path.",
        algorithm
      ),
      call. = FALSE
    )
  }

  list(
    name = algorithm,
    fit = getExportedValue("pls", fit_name)
  )
}


#' Recover predictor variance explained from a coefficient path
#'
#' Each successive PLS/PCR coefficient increment is rank one.  Its left
#' direction is proportional to the corresponding effective score direction,
#' so the usual pls::explvar quantity can be recovered without asking the pls
#' kernels to retain fitted and residual arrays.
#' @noRd
.feature_rsa_x_variance_from_coefficients <- function(predictors,
                                                       coefficients,
                                                       xmeans) {
  predictors <- as.matrix(predictors)
  coefficients <- as.array(coefficients)
  dims <- dim(coefficients)
  if (length(dims) != 3L || dims[1L] != ncol(predictors)) {
    return(NULL)
  }

  centered <- sweep(predictors, 2L, xmeans, "-")
  ncomp <- dims[3L]
  xvar <- numeric(ncomp)
  tolerance <- .Machine$double.eps * max(
    1,
    sum(vapply(
      seq_len(ncomp),
      function(component) {
        values <- coefficients[, , component, drop = FALSE]
        sum(values * values)
      },
      numeric(1L)
    ))
  )

  for (component in seq_len(ncomp)) {
    increment <- coefficients[, , component, drop = FALSE][, , 1L]
    if (component > 1L) {
      increment <- increment - coefficients[, , component - 1L, drop = FALSE][, , 1L]
    }
    increment_ss <- sum(increment * increment)
    if (!is.finite(increment_ss) || increment_ss <= tolerance) {
      return(NULL)
    }
    direction <- La.svd(increment, nu = 1L, nv = 0L)$u[, 1L]
    score <- as.numeric(centered %*% direction)
    score_ss <- sum(score * score)
    if (!is.finite(score_ss) || score_ss <= tolerance) {
      return(NULL)
    }
    loading <- as.numeric(crossprod(centered, score)) / score_ss
    xvar[[component]] <- sum(loading * loading) * score_ss
  }

  xtotvar <- sum(vapply(
    seq_len(ncol(centered)),
    function(column) sum(centered[, column] * centered[, column]),
    numeric(1L)
  ))
  if (any(!is.finite(xvar)) || !is.finite(xtotvar)) {
    return(NULL)
  }
  list(Xvar = xvar, Xtotvar = xtotvar)
}


#' Fit a compact matrix-level PLS or PCR kernel
#'
#' The public `pls` formula interface retains model frames, fitted-value arrays,
#' and (for validation) a large prediction cube. Feature RSA needs only the
#' coefficient path and centering constants for the final prediction. This
#' helper calls the same exported numerical kernels used by `pls::plsr()` and
#' `pls::pcr()` and immediately discards all other material.
#' @noRd
.feature_rsa_fit_kernel <- function(predictors,
                                    responses,
                                    ncomp,
                                    method,
                                    keep_explained_variance = FALSE,
                                    retain_ncomp = NULL) {
  predictors <- as.matrix(predictors)
  responses <- as.matrix(responses)
  predictor_names <- colnames(predictors)
  response_names <- colnames(responses)
  method <- match.arg(method, c("pls", "pca"))
  ncomp <- as.integer(ncomp)

  if (nrow(predictors) != nrow(responses)) {
    stop("feature RSA kernel: predictors and responses must have the same rows.",
         call. = FALSE)
  }
  if (nrow(predictors) < 2L || ncol(predictors) < 1L ||
      ncol(responses) < 1L || length(ncomp) != 1L ||
      is.na(ncomp) || ncomp < 1L) {
    stop("feature RSA kernel: invalid matrix dimensions or component count.",
         call. = FALSE)
  }

  kernel <- .feature_rsa_kernel_function(method)
  raw_fit <- kernel$fit(
    predictors,
    responses,
    ncomp = ncomp,
    center = TRUE,
    stripped = TRUE
  )

  coefficients <- raw_fit$coefficients
  if (length(dim(coefficients)) == 2L) {
    coefficients <- array(
      coefficients,
      dim = c(nrow(coefficients), ncol(coefficients), 1L)
    )
  }
  component_numbers <- seq_len(dim(coefficients)[3L])
  if (!is.null(retain_ncomp)) {
    retain_ncomp <- as.integer(retain_ncomp)
    if (length(retain_ncomp) != 1L || is.na(retain_ncomp) ||
        retain_ncomp < 1L || retain_ncomp > length(component_numbers)) {
      stop("feature RSA kernel: invalid retained component.", call. = FALSE)
    }
    coefficients <- coefficients[, , retain_ncomp, drop = FALSE]
    component_numbers <- retain_ncomp
  }

  out <- list(
    coefficients = coefficients,
    Xmeans = as.numeric(raw_fit$Xmeans),
    Ymeans = as.numeric(raw_fit$Ymeans),
    ncomp = as.integer(max(component_numbers)),
    component_numbers = as.integer(component_numbers),
    method = method,
    algorithm = kernel$name,
    predictor_names = predictor_names,
    response_names = response_names
  )
  if (isTRUE(keep_explained_variance)) {
    explained <- .feature_rsa_x_variance_from_coefficients(
      predictors,
      raw_fit$coefficients,
      raw_fit$Xmeans
    )
    if (is.null(explained)) {
      # Degenerate response directions are extremely unusual after the
      # variance checks in train_model(), but the full pls path remains a
      # semantics-preserving fallback for them.
      variance_fit <- kernel$fit(
        predictors,
        responses,
        ncomp = ncomp,
        center = TRUE,
        stripped = FALSE
      )
      explained <- list(
        Xvar = as.numeric(variance_fit$Xvar),
        Xtotvar = as.numeric(variance_fit$Xtotvar)
      )
    }
    out$Xvar <- as.numeric(explained$Xvar)
    out$Xtotvar <- as.numeric(explained$Xtotvar)
  }
  class(out) <- "feature_rsa_kernel_fit"
  out
}


#' Predict from a compact feature RSA kernel
#' @noRd
.feature_rsa_predict_kernel <- function(fit, newdata, ncomp = fit$ncomp) {
  if (!inherits(fit, "feature_rsa_kernel_fit")) {
    stop("feature RSA kernel prediction requires a feature_rsa_kernel_fit.",
         call. = FALSE)
  }
  newdata <- as.matrix(newdata)
  ncomp <- as.integer(ncomp)
  component_numbers <- fit$component_numbers %||%
    seq_len(dim(fit$coefficients)[3L])
  coefficient_index <- match(ncomp, component_numbers)
  if (length(ncomp) != 1L || is.na(ncomp) || ncomp < 1L ||
      is.na(coefficient_index)) {
    stop(sprintf(
      "feature RSA kernel prediction: component %s is not retained (available: %s).",
      as.character(ncomp), paste(component_numbers, collapse = ", ")
    ), call. = FALSE)
  }
  if (ncol(newdata) != length(fit$Xmeans)) {
    stop(sprintf(
      "feature RSA kernel prediction: expected %d columns, got %d.",
      length(fit$Xmeans), ncol(newdata)
    ), call. = FALSE)
  }

  centered <- sweep(newdata, 2L, fit$Xmeans, "-")
  coefficients <- fit$coefficients[, , coefficient_index, drop = FALSE][, , 1L]
  predicted <- centered %*% coefficients
  if (!is.matrix(predicted)) {
    predicted <- matrix(predicted, nrow = nrow(newdata))
  }
  predicted <- sweep(predicted, 2L, fit$Ymeans, "+")
  if (!is.null(fit$response_names)) {
    colnames(predicted) <- fit$response_names
  }
  predicted
}


#' Construct leave-one-block-out segments for an outer training fold
#' @noRd
.feature_rsa_block_segments <- function(block_var,
                                        n_rows,
                                        observation_indices = NULL) {
  n_rows <- as.integer(n_rows)
  if (length(n_rows) != 1L || is.na(n_rows) || n_rows < 2L) {
    stop("feature RSA blocked selection: invalid training row count.",
         call. = FALSE)
  }
  if (is.null(block_var)) {
    stop("feature RSA blocked selection requires design$block_var.",
         call. = FALSE)
  }

  if (is.null(observation_indices)) {
    if (length(block_var) != n_rows) {
      stop(
        paste0(
          "feature RSA blocked selection needs original observation indices ",
          "when training on a subset of design$block_var."
        ),
        call. = FALSE
      )
    }
    fold_blocks <- block_var
  } else {
    observation_indices <- as.integer(observation_indices)
    if (length(observation_indices) != n_rows || anyNA(observation_indices) ||
        any(observation_indices < 1L | observation_indices > length(block_var))) {
      stop(
        "feature RSA blocked selection received invalid observation indices.",
        call. = FALSE
      )
    }
    fold_blocks <- block_var[observation_indices]
  }

  if (anyNA(fold_blocks)) {
    stop("feature RSA blocked selection does not allow missing block labels.",
         call. = FALSE)
  }
  block_key <- as.character(fold_blocks)
  segments <- unname(split(seq_len(n_rows), block_key, drop = TRUE))
  segments <- Filter(length, segments)
  if (length(segments) < 2L) {
    stop(
      "feature RSA blocked selection requires at least two non-empty blocks in every training fold.",
      call. = FALSE
    )
  }
  lapply(segments, as.integer)
}


#' Compute segment-wise validation scores without retaining a prediction cube
#'
#' `pls:::mvrCv()` predicts every observation for every validation segment and
#' only then retains the held-out rows. Feature RSA's component selector needs
#' only one scalar score per segment and component. Computing those scalars
#' while streaming segments preserves the selection estimand while bounding
#' memory. Pattern-relative objectives are evaluated in the original neural
#' response space after inner-fold predictions are unstandardized.
#' @noRd
.feature_rsa_cv_segment_scores <- function(predictors,
                                           responses,
                                           ncomp,
                                           method,
                                           segments,
                                           fold_standardize = FALSE,
                                           objective = c(
                                             "mse",
                                             "pattern_discrimination",
                                             "pattern_rank_percentile"
                                           ),
                                           feature_standardize = c("scale", "center")) {
  predictors <- as.matrix(predictors)
  responses <- as.matrix(responses)
  method <- match.arg(method, c("pls", "pca"))
  objective <- match.arg(objective)
  feature_standardize <- match.arg(feature_standardize)
  ncomp <- as.integer(ncomp)
  nobj <- nrow(predictors)

  if (!is.logical(fold_standardize) || length(fold_standardize) != 1L ||
      is.na(fold_standardize)) {
    stop("feature RSA CV: fold_standardize must be TRUE or FALSE.",
         call. = FALSE)
  }

  if (nobj != nrow(responses)) {
    stop("feature RSA CV: predictors and responses must have the same rows.",
         call. = FALSE)
  }
  if (!is.list(segments) || length(segments) < 1L) {
    stop("feature RSA CV: segments must be a non-empty list.", call. = FALSE)
  }

  segments <- lapply(segments, function(idx) {
    idx <- unique(as.integer(idx))
    if (!length(idx) || anyNA(idx) || any(idx < 1L | idx > nobj)) {
      stop("feature RSA CV: segment indices are empty or out of bounds.",
           call. = FALSE)
    }
    idx
  })
  if (objective != "mse") {
    if (!isTRUE(fold_standardize)) {
      stop(
        "feature RSA pattern-relative component tuning requires nested fold standardization.",
        call. = FALSE
      )
    }
    if (any(lengths(segments) < 2L)) {
      stop(
        paste0(
          "feature RSA pattern-relative component tuning requires at least ",
          "two held-out observations in every segment."
        ),
        call. = FALSE
      )
    }
  }
  max_segment <- max(lengths(segments))
  # Match pls:::mvrCv(), which reserves one residual degree of freedom in
  # every training split.
  ncomp_cv <- min(ncomp, nobj - max_segment - 1L)
  if (!is.finite(ncomp_cv) || ncomp_cv < 1L) {
    stop("feature RSA CV: too few training rows for component selection.",
         call. = FALSE)
  }

  nresp <- ncol(responses)
  all_rows <- seq_len(nobj)
  segment_scores <- matrix(
    NA_real_,
    nrow = length(segments),
    ncol = ncomp_cv
  )

  for (s in seq_along(segments)) {
    test_idx <- segments[[s]]
    train_idx <- all_rows[-test_idx]
    train_predictors <- predictors[train_idx, , drop = FALSE]
    train_responses <- responses[train_idx, , drop = FALSE]
    if (isTRUE(fold_standardize)) {
      sf <- .standardize(train_predictors, feature_standardize)
      sx <- .standardize(train_responses)
      train_predictors <- sf$X_sc
      train_responses <- sx$X_sc
      test_predictors <- scale(
        predictors[test_idx, , drop = FALSE],
        center = sf$mean,
        scale = sf$sd
      )
      observed_test <- scale(
        responses[test_idx, , drop = FALSE],
        center = sx$mean,
        scale = sx$sd
      )
    } else {
      test_predictors <- predictors[test_idx, , drop = FALSE]
      observed_test <- responses[test_idx, , drop = FALSE]
    }
    fit <- .feature_rsa_fit_kernel(
      train_predictors,
      train_responses,
      ncomp = ncomp_cv,
      method = method
    )

    test_centered <- sweep(
      test_predictors,
      2L,
      fit$Xmeans,
      "-"
    )
    target_bytes <- getOption(
      "rMVPA.feature_rsa_component_target_bytes",
      32 * 1024^2
    )
    bytes_per_component <- max(1, nrow(test_centered) * nresp * 8)
    component_block_size <- if (is.numeric(target_bytes) &&
                                length(target_bytes) == 1L &&
                                is.finite(target_bytes) &&
                                target_bytes > 0) {
      max(
        1L,
        min(ncomp_cv, as.integer(floor(target_bytes / bytes_per_component)))
      )
    } else {
      ncomp_cv
    }
    component_starts <- seq.int(1L, ncomp_cv, by = component_block_size)
    for (component_start in component_starts) {
      component_end <- min(
        ncomp_cv,
        component_start + component_block_size - 1L
      )
      component_index <- seq.int(component_start, component_end)
      coefficient_path <- matrix(
        fit$coefficients[, , component_index, drop = FALSE],
        nrow = dim(fit$coefficients)[1L],
        ncol = nresp * length(component_index)
      )
      pred_path <- test_centered %*% coefficient_path
      pred_path <- sweep(
        pred_path,
        2L,
        rep(fit$Ymeans, times = length(component_index)),
        "+"
      )
      for (local_component in seq_along(component_index)) {
        cols <- seq.int(
          (local_component - 1L) * nresp + 1L,
          local_component * nresp
        )
        predicted_component <- pred_path[, cols, drop = FALSE]
        score <- if (objective == "mse") {
          mean((observed_test - predicted_component)^2, na.rm = TRUE)
        } else {
          predicted_raw <- sweep(
            sweep(predicted_component, 2L, sx$sd, "*"),
            2L,
            sx$mean,
            "+"
          )
          .feature_rsa_pattern_objective(
            predicted_raw,
            responses[test_idx, , drop = FALSE],
            objective
          )
        }
        if (objective != "mse" && !is.finite(score)) {
          stop(
            sprintf(
              paste0(
                "feature RSA %s component tuning produced an undefined ",
                "score for segment %d, component %d."
              ),
              objective,
              s,
              component_index[[local_component]]
            ),
            call. = FALSE
          )
        }
        segment_scores[s, component_index[[local_component]]] <- score
      }
    }
  }

  segment_scores
}


#' Backward-compatible component MSE helper
#' @noRd
.feature_rsa_cv_segment_mse <- function(predictors,
                                        responses,
                                        ncomp,
                                        method,
                                        segments,
                                        fold_standardize = FALSE,
                                        feature_standardize = "scale") {
  .feature_rsa_cv_segment_scores(
    predictors = predictors,
    responses = responses,
    ncomp = ncomp,
    method = method,
    segments = segments,
    fold_standardize = fold_standardize,
    objective = "mse",
    feature_standardize = feature_standardize
  )
}


#' Select components from segment-wise validation scores
#' @noRd
.feature_rsa_select_from_segment_scores <- function(segment_scores,
                                                    maximize = FALSE,
                                                    one_se = TRUE) {
  segment_scores <- as.matrix(segment_scores)
  if (!nrow(segment_scores) || !ncol(segment_scores)) {
    return(NA_integer_)
  }
  mean_scores <- colMeans(segment_scores, na.rm = TRUE)
  valid <- is.finite(mean_scores)
  if (!any(valid)) {
    return(NA_integer_)
  }

  candidate_indices <- which(valid)
  best_idx <- if (isTRUE(maximize)) {
    candidate_indices[[which.max(mean_scores[valid])]]
  } else {
    candidate_indices[[which.min(mean_scores[valid])]]
  }
  if (!length(best_idx)) {
    return(NA_integer_)
  }
  if (!isTRUE(one_se)) {
    return(as.integer(best_idx))
  }

  score_at_best <- segment_scores[, best_idx]
  n_finite <- sum(is.finite(score_at_best))
  if (n_finite < 2L) {
    return(as.integer(best_idx))
  }

  se <- stats::sd(score_at_best, na.rm = TRUE) / sqrt(n_finite)
  eligible <- if (isTRUE(maximize)) {
    which(valid & mean_scores >= mean_scores[[best_idx]] - se)
  } else {
    which(valid & mean_scores <= mean_scores[[best_idx]] + se)
  }
  selected <- eligible[[1L]]
  if (is.na(selected)) as.integer(best_idx) else as.integer(selected)
}


#' Backward-compatible component MSE selector
#' @noRd
.feature_rsa_select_from_segment_mse <- function(segment_mse,
                                                 method = "onesigma") {
  .feature_rsa_select_from_segment_scores(
    segment_mse,
    maximize = FALSE,
    one_se = identical(method, "onesigma")
  )
}


#' Validate or construct the normalized ridge penalty grid
#'
#' Feature RSA ridge uses the objective
#' `mean((Y - X %*% B)^2) + lambda * ||B||_F^2`, up to the response-count
#' constant.  The corresponding normal equations add `n * lambda` to the
#' feature cross-product, making a supplied lambda comparable across folds of
#' different sizes.
#' @noRd
.feature_rsa_ridge_lambdas <- function(lambda = NULL) {
  if (is.null(lambda)) {
    return(10 ^ seq(-6, 2, length.out = 49L))
  }
  if (!is.numeric(lambda) || !length(lambda) || anyNA(lambda) ||
      any(!is.finite(lambda)) || any(lambda < 0)) {
    stop("feature RSA ridge lambda values must be finite and non-negative.",
         call. = FALSE)
  }
  sort(unique(as.numeric(lambda)))
}


#' Economy-SVD decomposition for multi-response feature RSA ridge
#'
#' The inputs are already on the preprocessing scale selected by the caller.
#' This helper estimates an unpenalized intercept by centering, but deliberately
#' does not rescale columns: LOO keeps the outer-fold scale fixed, while blocked
#' tuning supplies matrices standardized inside each inner training split.
#' @noRd
.feature_rsa_ridge_decomposition <- function(predictors, responses) {
  predictors <- as.matrix(predictors)
  responses <- as.matrix(responses)
  predictor_names <- colnames(predictors)
  response_names <- colnames(responses)

  if (nrow(predictors) != nrow(responses)) {
    stop("feature RSA ridge: predictors and responses must have the same rows.",
         call. = FALSE)
  }
  if (nrow(predictors) < 2L || ncol(predictors) < 1L ||
      ncol(responses) < 1L) {
    stop("feature RSA ridge: invalid matrix dimensions.", call. = FALSE)
  }
  if (any(!is.finite(predictors)) || any(!is.finite(responses))) {
    stop("feature RSA ridge requires finite predictors and responses.",
         call. = FALSE)
  }

  predictor_means <- colMeans(predictors)
  response_means <- colMeans(responses)
  predictors_centered <- sweep(predictors, 2L, predictor_means, "-")
  responses_centered <- sweep(responses, 2L, response_means, "-")
  spectral <- svd(
    predictors_centered,
    nu = min(dim(predictors_centered)),
    nv = min(dim(predictors_centered))
  )
  tolerance <- .brs_spectral_tolerance(
    spectral$d,
    max(dim(predictors_centered)),
    tolerance = NULL
  )

  out <- list(
    U = spectral$u,
    d = spectral$d,
    V = spectral$v,
    UY = crossprod(spectral$u, responses_centered),
    responses_centered = responses_centered,
    response_ss = sum(responses_centered ^ 2),
    predictor_means = predictor_means,
    response_means = response_means,
    predictor_names = predictor_names,
    response_names = response_names,
    n = nrow(predictors),
    p = ncol(predictors),
    n_response = ncol(responses),
    tolerance = tolerance,
    rank = sum(spectral$d > tolerance)
  )
  class(out) <- c("feature_rsa_ridge_decomposition", "list")
  out
}


.feature_rsa_ridge_shrinkage <- function(decomposition, lambda) {
  if (!inherits(decomposition, "feature_rsa_ridge_decomposition")) {
    stop("feature RSA ridge requires a spectral decomposition.", call. = FALSE)
  }
  lambda <- .feature_rsa_ridge_lambdas(lambda)
  if (length(lambda) != 1L) {
    stop("feature RSA ridge fit requires exactly one lambda.", call. = FALSE)
  }
  d <- decomposition$d
  if (lambda == 0) {
    keep <- d > decomposition$tolerance
    coefficient <- numeric(length(d))
    coefficient[keep] <- 1 / d[keep]
    leverage <- as.numeric(keep)
  } else {
    penalty <- decomposition$n * lambda
    coefficient <- d / (d ^ 2 + penalty)
    leverage <- d ^ 2 / (d ^ 2 + penalty)
  }
  list(
    lambda = lambda,
    coefficient = coefficient,
    leverage = leverage,
    effective_df = 1 + sum(leverage)
  )
}


#' Fit a compact multi-response ridge kernel
#' @noRd
.feature_rsa_ridge_fit_kernel <- function(predictors,
                                          responses,
                                          lambda,
                                          decomposition = NULL) {
  predictors <- as.matrix(predictors)
  responses <- as.matrix(responses)
  if (is.null(decomposition)) {
    decomposition <- .feature_rsa_ridge_decomposition(predictors, responses)
  } else if (!inherits(decomposition, "feature_rsa_ridge_decomposition") ||
             decomposition$n != nrow(predictors) ||
             decomposition$p != ncol(predictors) ||
             decomposition$n_response != ncol(responses)) {
    stop("feature RSA ridge decomposition does not match the supplied matrices.",
         call. = FALSE)
  }

  shrinkage <- .feature_rsa_ridge_shrinkage(decomposition, lambda)
  coefficients <- decomposition$V %*% sweep(
    decomposition$UY,
    1L,
    shrinkage$coefficient,
    "*"
  )
  dimnames(coefficients) <- list(
    decomposition$predictor_names,
    decomposition$response_names
  )

  out <- list(
    coefficients = coefficients,
    Xmeans = decomposition$predictor_means,
    Ymeans = decomposition$response_means,
    lambda = shrinkage$lambda,
    effective_df = shrinkage$effective_df,
    rank = decomposition$rank,
    predictor_names = decomposition$predictor_names,
    response_names = decomposition$response_names
  )
  class(out) <- c("feature_rsa_ridge_fit", "list")
  out
}


#' Predict from a compact feature RSA ridge kernel
#' @noRd
.feature_rsa_ridge_predict_kernel <- function(fit, newdata) {
  if (!inherits(fit, "feature_rsa_ridge_fit")) {
    stop("feature RSA ridge prediction requires a feature_rsa_ridge_fit.",
         call. = FALSE)
  }
  newdata <- as.matrix(newdata)
  if (ncol(newdata) != length(fit$Xmeans)) {
    stop(sprintf(
      "feature RSA ridge prediction: expected %d columns, got %d.",
      length(fit$Xmeans), ncol(newdata)
    ), call. = FALSE)
  }
  if (any(!is.finite(newdata))) {
    stop("feature RSA ridge prediction requires finite newdata.",
         call. = FALSE)
  }
  predicted <- sweep(newdata, 2L, fit$Xmeans, "-") %*%
    fit$coefficients
  predicted <- sweep(predicted, 2L, fit$Ymeans, "+")
  if (!is.null(fit$response_names)) {
    colnames(predicted) <- fit$response_names
  }
  predicted
}


#' Generalized cross-validation scores for normalized ridge penalties
#' @noRd
.feature_rsa_ridge_gcv_scores <- function(decomposition, lambdas) {
  lambdas <- .feature_rsa_ridge_lambdas(lambdas)
  response_projection_ss <- rowSums(decomposition$UY ^ 2)
  vapply(lambdas, function(lambda) {
    shrinkage <- .feature_rsa_ridge_shrinkage(decomposition, lambda)
    h <- shrinkage$leverage
    rss <- decomposition$response_ss -
      sum((2 * h - h ^ 2) * response_projection_ss)
    rss <- max(rss, 0)
    denominator <- 1 - shrinkage$effective_df / decomposition$n
    if (!is.finite(denominator) || denominator <= sqrt(.Machine$double.eps)) {
      return(Inf)
    }
    (rss / (decomposition$n * decomposition$n_response)) /
      denominator ^ 2
  }, numeric(1L))
}


#' Observation-wise analytic PRESS errors for normalized ridge penalties
#' @noRd
.feature_rsa_ridge_loo_mse <- function(decomposition, lambdas) {
  lambdas <- .feature_rsa_ridge_lambdas(lambdas)
  out <- matrix(
    Inf,
    nrow = decomposition$n,
    ncol = length(lambdas),
    dimnames = list(NULL, format(lambdas, scientific = TRUE))
  )
  squared_u <- decomposition$U ^ 2
  leave_one_out_scale <- decomposition$n / (decomposition$n - 1)

  for (j in seq_along(lambdas)) {
    # A normalized penalty is (n - 1) * lambda in every leave-one-out fit,
    # rather than n * lambda as in the full-data fit. Centering after deletion
    # gives S[-i] = S - n/(n-1) x_i x_i'. Applying Sherman-Morrison to that
    # downdate yields the exact refit residual below.
    d <- decomposition$d
    if (lambdas[[j]] == 0) {
      keep <- d > decomposition$tolerance
      shrinkage_leverage <- as.numeric(keep)
    } else {
      penalty <- (decomposition$n - 1) * lambdas[[j]]
      shrinkage_leverage <- d ^ 2 / (d ^ 2 + penalty)
    }
    fitted_centered <- decomposition$U %*% sweep(
      decomposition$UY,
      1L,
      shrinkage_leverage,
      "*"
    )
    residual <- decomposition$responses_centered - fitted_centered
    leverage <- leave_one_out_scale * rowSums(sweep(
      squared_u,
      2L,
      shrinkage_leverage,
      "*"
    ))
    denominator <- 1 - leverage
    valid <- is.finite(denominator) &
      denominator > sqrt(.Machine$double.eps)
    if (any(valid)) {
      loo_residual <- sweep(
        leave_one_out_scale * residual[valid, , drop = FALSE],
        1L,
        denominator[valid],
        "/"
      )
      out[valid, j] <- rowMeans(loo_residual ^ 2)
    }
  }
  out
}


#' Held-out MSE path from compact ridge spectral cross-products
#'
#' For blocked tuning, materializing one held-out observations-by-voxels
#' prediction matrix per lambda repeats a large matrix multiplication. The
#' squared-error identity below reduces the entire path to one feature cross-
#' product, one feature-response cross-product, and lambda-specific operations
#' on matrices bounded by the feature rank.
#' @noRd
.feature_rsa_ridge_prediction_mse_path <- function(decomposition,
                                                    newdata,
                                                    observed,
                                                    lambdas) {
  if (!inherits(decomposition, "feature_rsa_ridge_decomposition")) {
    stop("feature RSA ridge MSE path requires a spectral decomposition.",
         call. = FALSE)
  }
  newdata <- as.matrix(newdata)
  observed <- as.matrix(observed)
  lambdas <- .feature_rsa_ridge_lambdas(lambdas)
  if (nrow(newdata) != nrow(observed) ||
      ncol(newdata) != decomposition$p ||
      ncol(observed) != decomposition$n_response) {
    stop("feature RSA ridge MSE path received incompatible matrices.",
         call. = FALSE)
  }
  if (any(!is.finite(newdata)) || any(!is.finite(observed))) {
    stop("feature RSA ridge MSE path requires finite matrices.",
         call. = FALSE)
  }

  projected_features <- sweep(
    newdata,
    2L,
    decomposition$predictor_means,
    "-"
  ) %*% decomposition$V
  centered_observed <- sweep(
    observed,
    2L,
    decomposition$response_means,
    "-"
  )
  linear_term <- rowSums(
    crossprod(projected_features, centered_observed) * decomposition$UY
  )
  quadratic_kernel <- crossprod(projected_features) *
    tcrossprod(decomposition$UY)
  observed_ss <- sum(centered_observed ^ 2)
  denominator <- length(centered_observed)

  vapply(lambdas, function(lambda) {
    coefficient <- .feature_rsa_ridge_shrinkage(
      decomposition,
      lambda
    )$coefficient
    fitted_ss <- sum(quadratic_kernel * tcrossprod(coefficient))
    residual_ss <- observed_ss -
      2 * sum(coefficient * linear_term) + fitted_ss
    max(residual_ss, 0) / denominator
  }, numeric(1L))
}


#' Block-wise ridge errors with leakage-safe inner preprocessing
#' @noRd
.feature_rsa_ridge_block_mse <- function(features,
                                         responses,
                                         lambdas,
                                         segments,
                                         feature_standardize = c("scale", "center")) {
  features <- as.matrix(features)
  responses <- as.matrix(responses)
  lambdas <- .feature_rsa_ridge_lambdas(lambdas)
  feature_standardize <- match.arg(feature_standardize)
  n <- nrow(features)
  if (n != nrow(responses)) {
    stop("feature RSA ridge blocked CV requires matching row counts.",
         call. = FALSE)
  }
  if (!is.list(segments) || length(segments) < 2L) {
    stop("feature RSA ridge blocked CV requires at least two segments.",
         call. = FALSE)
  }
  segments <- lapply(segments, function(index) {
    index <- unique(as.integer(index))
    if (!length(index) || anyNA(index) || any(index < 1L | index > n)) {
      stop("feature RSA ridge blocked CV received invalid segment indices.",
           call. = FALSE)
    }
    index
  })

  out <- matrix(
    NA_real_,
    nrow = length(segments),
    ncol = length(lambdas),
    dimnames = list(NULL, format(lambdas, scientific = TRUE))
  )
  all_rows <- seq_len(n)
  for (s in seq_along(segments)) {
    test <- segments[[s]]
    train <- all_rows[-test]
    if (length(train) < 2L) {
      stop("feature RSA ridge blocked CV has fewer than two training rows.",
           call. = FALSE)
    }
    sf <- .standardize(features[train, , drop = FALSE], feature_standardize)
    sx <- .standardize(responses[train, , drop = FALSE])
    test_features <- scale(
      features[test, , drop = FALSE],
      center = sf$mean,
      scale = sf$sd
    )
    observed_test <- scale(
      responses[test, , drop = FALSE],
      center = sx$mean,
      scale = sx$sd
    )
    decomposition <- .feature_rsa_ridge_decomposition(sf$X_sc, sx$X_sc)
    out[s, ] <- .feature_rsa_ridge_prediction_mse_path(
      decomposition,
      test_features,
      observed_test,
      lambdas
    )
  }
  out
}


#' Block-wise held-out pattern-relative path for ridge selection
#'
#' Each score is computed only among observations withheld together, so no
#' candidate pattern used to train a prediction can enter its score. Unlike the
#' spectral MSE path, pattern-relative scoring requires materialized held-out
#' predictions for each candidate penalty. Predictions are returned to the raw
#' neural response scale used by the reported outer-fold metrics.
#' @noRd
.feature_rsa_ridge_block_pattern_scores <- function(features,
                                                    responses,
                                                    lambdas,
                                                    segments,
                                                    objective = c(
                                                      "pattern_discrimination",
                                                      "pattern_rank_percentile"
                                                    ),
                                                    feature_standardize = c("scale", "center")) {
  features <- as.matrix(features)
  responses <- as.matrix(responses)
  lambdas <- .feature_rsa_ridge_lambdas(lambdas)
  objective <- match.arg(objective)
  feature_standardize <- match.arg(feature_standardize)
  n <- nrow(features)
  if (n != nrow(responses)) {
    stop("feature RSA ridge blocked pattern tuning requires matching row counts.",
         call. = FALSE)
  }
  if (!is.list(segments) || length(segments) < 2L) {
    stop("feature RSA ridge blocked pattern tuning requires at least two segments.",
         call. = FALSE)
  }
  segments <- lapply(segments, function(index) {
    index <- unique(as.integer(index))
    if (!length(index) || anyNA(index) || any(index < 1L | index > n)) {
      stop("feature RSA ridge blocked pattern tuning received invalid segment indices.",
           call. = FALSE)
    }
    if (length(index) < 2L) {
      stop(
        paste0(
          "feature RSA ridge blocked pattern tuning requires at least two held-out ",
          "observations in every segment."
        ),
        call. = FALSE
      )
    }
    index
  })

  out <- matrix(
    NA_real_,
    nrow = length(segments),
    ncol = length(lambdas),
    dimnames = list(NULL, format(lambdas, scientific = TRUE))
  )
  all_rows <- seq_len(n)
  for (s in seq_along(segments)) {
    test <- segments[[s]]
    train <- all_rows[-test]
    if (length(train) < 2L) {
      stop("feature RSA ridge blocked pattern tuning has fewer than two training rows.",
           call. = FALSE)
    }
    sf <- .standardize(features[train, , drop = FALSE], feature_standardize)
    sx <- .standardize(responses[train, , drop = FALSE])
    test_features <- scale(
      features[test, , drop = FALSE],
      center = sf$mean,
      scale = sf$sd
    )
    decomposition <- .feature_rsa_ridge_decomposition(sf$X_sc, sx$X_sc)
    for (j in seq_along(lambdas)) {
      fit <- .feature_rsa_ridge_fit_kernel(
        sf$X_sc,
        sx$X_sc,
        lambdas[[j]],
        decomposition = decomposition
      )
      predicted_standardized <- .feature_rsa_ridge_predict_kernel(
        fit,
        test_features
      )
      predicted_test <- sweep(
        sweep(predicted_standardized, 2L, sx$sd, "*"),
        2L,
        sx$mean,
        "+"
      )
      score <- .feature_rsa_pattern_objective(
        predicted_test,
        responses[test, , drop = FALSE],
        objective
      )
      if (!is.finite(score)) {
        stop(
          sprintf(
            paste0(
              "feature RSA ridge %s tuning produced an undefined score ",
              "for segment %d, lambda %s."
            ),
            objective,
            s,
            format(lambdas[[j]], scientific = TRUE)
          ),
          call. = FALSE
        )
      }
      out[s, j] <- score
    }
  }
  out
}


#' Block-wise held-out pattern discrimination for ridge selection
#' @noRd
.feature_rsa_ridge_block_pattern_discrimination <- function(features,
                                                             responses,
                                                             lambdas,
                                                             segments,
                                                             feature_standardize = "scale") {
  .feature_rsa_ridge_block_pattern_scores(
    features,
    responses,
    lambdas,
    segments,
    objective = "pattern_discrimination",
    feature_standardize = feature_standardize
  )
}


#' Block-wise held-out identification rank for ridge selection
#' @noRd
.feature_rsa_ridge_block_pattern_rank <- function(features,
                                                  responses,
                                                  lambdas,
                                                  segments,
                                                  feature_standardize = "scale") {
  .feature_rsa_ridge_block_pattern_scores(
    features,
    responses,
    lambdas,
    segments,
    objective = "pattern_rank_percentile",
    feature_standardize = feature_standardize
  )
}


#' Select a normalized ridge penalty from scalar or segment-wise scores
#' @noRd
.feature_rsa_ridge_select_lambda <- function(scores,
                                             lambdas,
                                             one_se = FALSE) {
  lambdas <- .feature_rsa_ridge_lambdas(lambdas)
  if (is.matrix(scores)) {
    if (ncol(scores) != length(lambdas)) {
      stop("feature RSA ridge score columns must match lambda values.",
           call. = FALSE)
    }
    objective <- colMeans(scores, na.rm = TRUE)
  } else {
    scores <- as.numeric(scores)
    if (length(scores) != length(lambdas)) {
      stop("feature RSA ridge scores must match lambda values.",
           call. = FALSE)
    }
    objective <- scores
  }
  valid <- is.finite(objective)
  if (!any(valid)) {
    stop("feature RSA ridge tuning produced no finite candidate score.",
         call. = FALSE)
  }
  best_value <- min(objective[valid])
  tolerance <- sqrt(.Machine$double.eps) * max(1, abs(best_value))
  best <- which(valid & objective <= best_value + tolerance)
  threshold <- best_value + tolerance

  if (isTRUE(one_se) && is.matrix(scores)) {
    best_index <- best[[which.max(lambdas[best])]]
    values <- scores[, best_index]
    values <- values[is.finite(values)]
    standard_error <- if (length(values) >= 2L) {
      stats::sd(values) / sqrt(length(values))
    } else {
      0
    }
    threshold <- objective[[best_index]] + standard_error + tolerance
  }
  eligible <- which(valid & objective <= threshold)
  as.integer(eligible[[which.max(lambdas[eligible])]])
}





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
  # `lambda_used` is always numeric. For successful CV it is lambda.min from
  # cv.glmnet; the retained fit is cv_fit$glmnet.fit, so no second full-data
  # glmnet fit is needed merely to predict at that value.
  preds_mat <- drop(predict(
    model$trained_model,
    newx = Fsc,
    s = model$lambda_used
  ))


  
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

      preds_sc <- if (inherits(pls_model, "feature_rsa_kernel_fit")) {
        .feature_rsa_predict_kernel(
          pls_model,
          sf_test,
          ncomp = ncomp_to_use
        )
      } else {
        # Backward compatibility for fit objects saved by rMVPA versions that
        # retained a full pls::mvr object.
        drop(predict(pls_model, newdata = sf_test, ncomp = ncomp_to_use))
      }

      if (!is.matrix(preds_sc)) {
        preds_sc <- matrix(preds_sc, nrow = nrow(F_new), ncol = length(x_mean))
      }

      preds <- sweep(sweep(preds_sc, 2, x_sd, "*"), 2, x_mean, "+")
      return(preds)

    } else if (method == "ridge") {
      required <- c("f_mean", "f_sd", "x_mean", "x_sd")
      missing <- required[vapply(required, function(name) {
        is.null(fit[[name]])
      }, logical(1L))]
      if (length(missing)) {
        stop(sprintf(
          "predict_model (ridge): Missing essential components: %s.",
          paste(missing, collapse = ", ")
        ))
      }
      if (!inherits(fit$trained_model, "feature_rsa_ridge_fit")) {
        stop("predict_model (ridge): Invalid trained ridge model.")
      }
      if (ncol(F_new) != length(fit$f_mean)) {
        stop(sprintf(
          "predict_model (ridge): Feature column mismatch. Expected %d, got %d.",
          length(fit$f_mean), ncol(F_new)
        ))
      }

      features_scaled <- sweep(
        sweep(F_new, 2L, fit$f_mean, "-"),
        2L,
        fit$f_sd,
        "/"
      )
      predictions_scaled <- .feature_rsa_ridge_predict_kernel(
        fit$trained_model,
        features_scaled
      )
      return(sweep(
        sweep(predictions_scaled, 2L, fit$x_sd, "*"),
        2L,
        fit$x_mean,
        "+"
      ))
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
    pred_sds <- .feature_rsa_col_sds(predictions)
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
#' @param fold_id Vector assigning observations to outer test folds. Null
#'   permutations are restricted to rows in the same fold.
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
                                   valid_col,
                                   fold_id = NULL,
                                   predicted_rdm_cache = NULL,
                                   observed_rdm_cache = NULL)
{
  futile.logger::flog.info(
    "Performing permutation tests with %d permutations... (feature_rsa_model)",
    nperm
  )

  n_rows <- nrow(predicted)
  fold_id <- .feature_rsa_resolve_fold_id(fold_id, n_rows)
  sd_thresh <- 1e-12
  tss <- sum((observed - mean(observed, na.rm = TRUE))^2, na.rm = TRUE)
  observed_valid <- observed[, valid_col, drop = FALSE]
  predicted_valid <- predicted[, valid_col, drop = FALSE]
  obs_row_sd <- .feature_rsa_row_sds(observed)
  obs_row_ok <- obs_row_sd > sd_thresh
  pred_row_sd <- .feature_rsa_row_sds(predicted)
  pred_row_ok <- pred_row_sd > sd_thresh

  # Row permutation changes only matrix indexing, not pairwise correlation.
  # Keep one cross-correlation matrix for pattern matching and compact
  # lower-triangle vectors for the two within-space geometries.
  cross_trial_cor <- if (n_rows >= 2L) {
    tryCatch(
      .feature_rsa_row_cor(predicted, observed),
      error = function(e) NULL
    )
  } else NULL
  if (n_rows >= 2L && is.null(predicted_rdm_cache)) {
    predicted_rdm_cache <- tryCatch(
      .feature_rsa_grouped_rdm_vectors(
        predicted,
        fold_id = fold_id,
        full_length = TRUE
      ),
      error = function(e) NULL
    )
  }
  if (n_rows >= 2L && is.null(observed_rdm_cache)) {
    observed_rdm_cache <- tryCatch(
      .feature_rsa_grouped_rdm_vectors(
        observed,
        fold_id = fold_id,
        full_length = TRUE
      ),
      error = function(e) NULL
    )
  }

  no_na_global_fast_path <- n_rows > 1L &&
    !anyNA(observed) &&
    !anyNA(predicted)
  no_na_temporal_fast_path <- n_rows > 1L &&
    length(valid_col) > 0L &&
    !anyNA(observed_valid) &&
    !anyNA(predicted_valid)

  if (no_na_global_fast_path) {
    obs_centered_global <- observed - mean(observed)
    pred_centered_global <- predicted - mean(predicted)
    global_denom <- sqrt(
      sum(obs_centered_global * obs_centered_global) *
        sum(pred_centered_global * pred_centered_global)
    )
  } else {
    global_denom <- NA_real_
    obs_centered_global <- NULL
    pred_centered_global <- NULL
  }

  if (no_na_temporal_fast_path) {
    obs_centered_cols <- sweep(observed_valid, 2L, colMeans(observed_valid), "-")
    pred_centered_cols <- sweep(predicted_valid, 2L, colMeans(predicted_valid), "-")
    col_denom <- sqrt(
      colSums(obs_centered_cols * obs_centered_cols) *
      colSums(pred_centered_cols * pred_centered_cols)
    )
  } else {
    obs_centered_cols <- NULL
    pred_centered_cols <- NULL
    col_denom <- NULL
  }

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
    perm_idx  <- .feature_rsa_permutation_index(fold_id)
    perm_pred <- predicted[perm_idx, , drop = FALSE]
    perm_pred_valid <- perm_pred[, valid_col, drop = FALSE]

    ## -- Condition-pattern metrics (trial x trial) --
    ppc <- ppd <- ppr <- prdm <- NA_real_
    vr <- which(obs_row_ok & pred_row_ok[perm_idx])
    if (length(vr) >= 2) {
      if (!is.null(cross_trial_cor)) {
        pattern_metrics <- .feature_rsa_pattern_metrics_from_cross(
          cross_trial_cor,
          predicted_rows = perm_idx[vr],
          observed_rows = vr,
          fold_id = fold_id[vr]
        )
        ppc <- unname(pattern_metrics[["correlation"]])
        ppd <- unname(pattern_metrics[["discrimination"]])
        ppr <- unname(pattern_metrics[["rank"]])
      }
    }

    ## -- Representational geometry (RDM correlation) --
    if (length(vr) >= 2L) {
      if (!is.null(predicted_rdm_cache) && !is.null(observed_rdm_cache)) {
        predicted_pair_index <- .feature_rsa_grouped_pair_indices(
          perm_idx[vr],
          fold_id[vr],
          n_rows
        )
        observed_pair_index <- .feature_rsa_grouped_pair_indices(
          vr,
          fold_id[vr],
          n_rows
        )
        pv <- predicted_rdm_cache[predicted_pair_index]
        ov <- observed_rdm_cache[observed_pair_index]
        if (length(pv) >= 2 && length(pv) == length(ov)) {
          prdm <- tryCatch(stats::cor(pv, ov, method = "spearman", use = "complete.obs"),
                           error = function(e) NA_real_)
        }
      }
    }

    ## -- Global reconstruction --
    pvc <- if (no_na_global_fast_path &&
               is.finite(global_denom) && global_denom > 0) {
      sum(obs_centered_global * pred_centered_global[perm_idx, , drop = FALSE]) / global_denom
    } else {
      tryCatch(stats::cor(as.vector(perm_pred),
                          as.vector(observed),
                          use = "pairwise.complete.obs"),
               error = function(e) NA_real_)
    }
    pmse <- mean((perm_pred - observed)^2, na.rm = TRUE)
    prsq <- if (tss > 0) 1 - sum((observed - perm_pred)^2, na.rm = TRUE) / tss else NA_real_

    ## -- Voxel encoding --
    pmvtc <- NA_real_
    if (n_rows > 1 && length(valid_col) > 0) {
      if (no_na_temporal_fast_path) {
        vcors <- rep(NA_real_, length(valid_col))
        ok <- is.finite(col_denom) & col_denom > 0
        if (any(ok)) {
          vcors[ok] <- colSums(
            obs_centered_cols[, ok, drop = FALSE] *
              pred_centered_cols[perm_idx, ok, drop = FALSE]
          ) / col_denom[ok]
        }
        pmvtc <- mean(vcors, na.rm = TRUE)
      } else {
        pmvtc <- mean(vapply(seq_len(length(valid_col)), function(j) {
          tryCatch(
            stats::cor(observed_valid[, j], perm_pred_valid[, j], use = "pairwise.complete.obs"),
            error = function(e) NA_real_
          )
        }, numeric(1)), na.rm = TRUE)
      }
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
#' @param compute_rdm_vectors Logical; when TRUE, also return compact predicted
#'   and observed RDM vectors for reuse by downstream code.
#' @param fold_id Optional vector assigning each observation to its outer
#'   cross-validation test fold. Pattern discrimination, identification rank,
#'   and RDM correlation then use only candidates or pairs withheld together.
#'   Cross-validated workflows supply this automatically.
#' @param ... Additional arguments
#'
#' @return A list containing:
#'   \describe{
#'     \item{pattern_correlation}{Mean diagonal of the trial x trial correlation
#'       matrix -- how well the predicted spatial pattern for each trial matches
#'       the correct observed pattern.}
#'     \item{pattern_discrimination}{Diagonal minus off-diagonal of the trial x
#'       trial correlation matrix, restricted to candidates in the same outer
#'       test fold -- how much better the correct trial is matched than eligible
#'       incorrect trials.}
#'     \item{pattern_rank_percentile}{For each trial, percentile rank of the
#'       correct pattern among candidates in the same outer test fold.
#'       0.5 = chance, 1 = perfect.}
#'     \item{rdm_correlation}{Spearman correlation between predicted and
#'       observed correlation-distance pairs whose observations were withheld
#'       together.}
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
                                             compute_rdm_vectors = isTRUE(object$return_rdm_vectors),
                                             fold_id = NULL,
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

  fold_id <- .feature_rsa_resolve_fold_id(fold_id, nrow(observed))

  sd_thresh <- 1e-12
  n_obs  <- nrow(observed)
  n_pairs <- n_obs * (n_obs - 1L) / 2L

  ## -- helpers to avoid repeating variance checks --
  obs_sd  <- .feature_rsa_col_sds(observed)
  pred_sd <- .feature_rsa_col_sds(predicted)
  valid_col <- which(obs_sd > sd_thresh & pred_sd > sd_thresh)
  observed_valid <- observed[, valid_col, drop = FALSE]
  predicted_valid <- predicted[, valid_col, drop = FALSE]

  ## Pattern and RDM metrics correlate across voxels, so a voxel that is
  ## constant over observations can still carry spatial information. Temporal
  ## variance filtering applies only to the per-voxel temporal metric.
  obs_row_sd  <- .feature_rsa_row_sds(observed)
  pred_row_sd <- .feature_rsa_row_sds(predicted)
  valid_row   <- which(obs_row_sd > sd_thresh & pred_row_sd > sd_thresh)

  ## ================================================================
  ## 1. Condition-pattern metrics  (trial x trial cormat)
  ## ================================================================
  ## cor(t(predicted), t(observed)):  n_obs x n_obs
  ##   entry (i,j) = cor(predicted[i, ], observed[j, ])  across voxels
  pattern_cor     <- NA_real_
  pattern_discrim <- NA_real_
  pattern_rank    <- NA_real_
  rdm_cor         <- NA_real_
  predicted_rdm_subset <- NULL
  observed_rdm_subset <- NULL

  if (length(valid_row) >= 2) {
    pmat <- predicted[valid_row, , drop = FALSE]
    omat <- observed[valid_row,  , drop = FALSE]
    pattern_metrics <- .feature_rsa_pattern_metrics_blockwise(
      pmat,
      omat,
      fold_id = fold_id[valid_row]
    )
    pattern_cor <- unname(pattern_metrics[["correlation"]])
    pattern_discrim <- unname(pattern_metrics[["discrimination"]])
    pattern_rank <- unname(pattern_metrics[["rank"]])
  }

  ## ================================================================
  ## 1b. Representational geometry metric (RDM correlation)
  ## ================================================================
  ## Compute within-space RDMs (1 - trial-by-trial correlation) and correlate
  ## lower triangles between predicted and observed. With cross-validation,
  ## only pairs withheld together enter the estimand.
  if (length(valid_row) >= 2) {
    pmat <- predicted[valid_row, , drop = FALSE]
    omat <- observed[valid_row,  , drop = FALSE]
    predicted_rdm_subset <- tryCatch(
      .feature_rsa_grouped_rdm_vectors(
        pmat,
        fold_id = fold_id[valid_row]
      ),
      error = function(e) NULL
    )
    observed_rdm_subset <- tryCatch(
      .feature_rsa_grouped_rdm_vectors(
        omat,
        fold_id = fold_id[valid_row]
      ),
      error = function(e) NULL
    )
    if (!is.null(predicted_rdm_subset) && !is.null(observed_rdm_subset)) {
      if (length(predicted_rdm_subset) >= 2 &&
          length(predicted_rdm_subset) == length(observed_rdm_subset)) {
        rdm_cor <- tryCatch(stats::cor(
                              predicted_rdm_subset,
                              observed_rdm_subset,
                              method = "spearman",
                              use = "complete.obs"),
                            error = function(e) NA_real_)
      }
    }
  }
  if (!is.finite(rdm_cor)) rdm_cor <- NA_real_

  ## ================================================================
  ## 2. Global reconstruction metrics
  ## ================================================================
  voxel_cor <- .feature_rsa_global_cor(predicted, observed)

  mse <- mean((predicted - observed)^2, na.rm = TRUE)
  rss <- sum((observed - predicted)^2, na.rm = TRUE)
  tss <- sum((observed - mean(observed, na.rm = TRUE))^2, na.rm = TRUE)
  r_squared <- if (tss == 0) NA_real_ else 1 - (rss / tss)

  ## ================================================================
  ## 3. Voxel encoding fidelity
  ## ================================================================
  mean_voxelwise_temporal_cor <- NA_real_
  if (n_obs > 1 && length(valid_col) > 0) {
    vcors <- .feature_rsa_diag_col_cor(observed_valid, predicted_valid)
    mean_voxelwise_temporal_cor <- mean(vcors, na.rm = TRUE)
  }
  if (!is.finite(mean_voxelwise_temporal_cor)) mean_voxelwise_temporal_cor <- NA_real_

  predicted_rdm_cache <- observed_rdm_cache <- NULL
  need_full_rdm <- isTRUE(compute_rdm_vectors) || nperm > 0L
  if (need_full_rdm) {
    if (n_pairs < 1L) {
      predicted_rdm_cache <- numeric(0)
      observed_rdm_cache <- numeric(0)
    } else {
      predicted_rdm_cache <- tryCatch(
        suppressWarnings(.feature_rsa_grouped_rdm_vectors(
          predicted,
          fold_id = fold_id,
          full_length = TRUE
        )),
        error = function(e) NULL
      )
      observed_rdm_cache <- tryCatch(
        suppressWarnings(.feature_rsa_grouped_rdm_vectors(
          observed,
          fold_id = fold_id,
          full_length = TRUE
        )),
        error = function(e) NULL
      )
      if (is.null(predicted_rdm_cache)) {
        predicted_rdm_cache <- rep(NA_real_, n_pairs)
      }
      if (is.null(observed_rdm_cache)) {
        observed_rdm_cache <- rep(NA_real_, n_pairs)
      }
    }
  }
  predicted_rdm_vec <- if (isTRUE(compute_rdm_vectors)) predicted_rdm_cache else NULL
  observed_rdm_vec <- if (isTRUE(compute_rdm_vectors)) observed_rdm_cache else NULL

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
      valid_col      = valid_col,
      fold_id        = fold_id,
      predicted_rdm_cache = predicted_rdm_cache,
      observed_rdm_cache = observed_rdm_cache
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
    predicted_rdm_vec           = predicted_rdm_vec,
    observed_rdm_vec            = observed_rdm_vec,
    permutation_results         = perm_results
  )
}

.feature_rsa_predicted_rdm_vector <- function(predicted) {
  predicted <- as.matrix(predicted)
  n_obs <- nrow(predicted)
  n_pairs <- n_obs * (n_obs - 1L) / 2L

  if (n_pairs < 1L) {
    return(numeric(0))
  }

  out <- tryCatch(
    suppressWarnings(.feature_rsa_rdm_vector_blockwise(predicted)),
    error = function(e) NULL
  )

  if (is.null(out) || length(out) != n_pairs) {
    return(rep(NA_real_, n_pairs))
  }
  out
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
  training_context <- list(...)
  observation_indices <- training_context$observation_indices
  
  result <- list(method = obj$method)
  feature_standardize <- match.arg(obj$feature_standardize %||% "scale",
                                   c("scale", "center"))
  result$feature_standardize <- feature_standardize
  
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
    require_package("pls", "for PLS/PCA regression in feature RSA")
    method_label <- toupper(obj$method)
    ncomp_sel <- obj$ncomp_selection %||% "max"
    ncomp_objective <- obj$ncomp_objective %||% "mse"
    ncomp_one_se <- if (!is.null(obj$ncomp_one_se)) {
      isTRUE(obj$ncomp_one_se)
    } else {
      ncomp_objective == "mse" && ncomp_sel %in% c("loo", "blocked")
    }
    if (ncomp_objective != "mse" && ncomp_sel != "blocked") {
      stop(
        paste0(
          "ncomp_objective='", ncomp_objective, "' requires ",
          "ncomp_selection='blocked'."
        ),
        call. = FALSE
      )
    }

    pls_res <- tryCatch({
      # Check for near-zero variance before standardization
      var_X <- apply(X, 2, var, na.rm = TRUE)
      var_F <- apply(Fsub, 2, var, na.rm = TRUE)
      if (any(var_X < 1e-10) || any(var_F < 1e-10)) {
        stop("Near zero variance detected in X or F before standardization.")
      }

      sx <- .standardize(X)
      sf <- .standardize(Fsub, feature_standardize)

      if (any(sx$sd < 1e-10) || any(sf$sd < 1e-10)) {
        stop("Near zero variance detected after standardization.")
      }

      max_k_possible <- min(nrow(sf$X_sc) - 1, ncol(sf$X_sc))
      k <- min(obj$max_comps, max_k_possible)

      if (k < 1) {
        stop(sprintf("Number of components (%d) < 1 (max_comps: %d, max_possible: %d).",
                     k, obj$max_comps, max_k_possible))
      }

      # --- Component selection (always >= 1) ---
      ncomp_use <- k
      if (ncomp_sel %in% c("loo", "blocked") && k > 1) {
        ncomp_use <- tryCatch({
          segments <- if (identical(ncomp_sel, "blocked")) {
            .feature_rsa_block_segments(
              obj$design$block_var,
              n_rows = nrow(sf$X_sc),
              observation_indices = observation_indices
            )
          } else {
            lapply(seq_len(nrow(sf$X_sc)), as.integer)
          }
          blocked_selection <- identical(ncomp_sel, "blocked")
          segment_scores <- .feature_rsa_cv_segment_scores(
            if (blocked_selection) Fsub else sf$X_sc,
            if (blocked_selection) X else sx$X_sc,
            ncomp = k,
            method = obj$method,
            segments = segments,
            fold_standardize = blocked_selection,
            objective = ncomp_objective,
            feature_standardize = feature_standardize
          )
          nc <- .feature_rsa_select_from_segment_scores(
            segment_scores,
            maximize = ncomp_objective != "mse",
            one_se = ncomp_one_se
          )
          if (is.na(nc) || nc < 1L) {
            futile.logger::flog.warn(
              "train_model (%s): selectNcomp returned %s; falling back to %d components.",
              method_label, as.character(nc), k)
            k
          } else {
            nc
          }
        }, error = function(e) {
          if (identical(ncomp_sel, "blocked")) {
            stop(e)
          }
          futile.logger::flog.warn(
            "train_model (%s): selectNcomp failed (%s); using all %d components.",
            method_label, e$message, k)
          k
        })
      }

      if (ncomp_sel == "pve") {
        model <- .feature_rsa_fit_kernel(
          sf$X_sc,
          sx$X_sc,
          ncomp = k,
          method = obj$method,
          keep_explained_variance = TRUE
        )
        xvar <- 100 * model$Xvar / model$Xtotvar
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
        # Explained-variance diagnostics are needed only for component
        # selection; do not retain them in every fold fit.
        model$Xvar <- NULL
        model$Xtotvar <- NULL
        model$coefficients <- model$coefficients[, , ncomp_use, drop = FALSE]
        model$component_numbers <- as.integer(ncomp_use)
        model$ncomp <- as.integer(ncomp_use)
      } else {
        # Component selection is already complete. Ask the numerical kernel for
        # no more components than are needed and retain only the selected
        # coefficient matrix for prediction.
        model <- .feature_rsa_fit_kernel(
          sf$X_sc,
          sx$X_sc,
          ncomp = ncomp_use,
          method = obj$method,
          retain_ncomp = ncomp_use
        )
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
    result$ncomp_selection <- ncomp_sel
    result$ncomp_objective <- ncomp_objective
    result$ncomp_one_se <- ncomp_one_se

  } else if (obj$method == "ridge") {
    ridge_result <- tryCatch({
      var_X <- apply(X, 2L, stats::var, na.rm = TRUE)
      var_F <- apply(Fsub, 2L, stats::var, na.rm = TRUE)
      if (any(!is.finite(var_X)) || any(!is.finite(var_F)) ||
          any(var_X < 1e-10) || any(var_F < 1e-10)) {
        stop("Near zero or non-finite variance detected in X or F before standardization.")
      }

      sx <- .standardize(X)
      sf <- .standardize(Fsub, feature_standardize)
      if (any(!is.finite(sx$X_sc)) || any(!is.finite(sf$X_sc)) ||
          any(sx$sd < 1e-10) || any(sf$sd < 1e-10)) {
        stop("Near zero or non-finite variance detected after standardization.")
      }

      lambdas <- .feature_rsa_ridge_lambdas(obj$lambda)
      decomposition <- .feature_rsa_ridge_decomposition(
        sf$X_sc,
        sx$X_sc
      )
      selector <- obj$lambda_selection %||% "gcv"
      objective <- obj$lambda_objective %||% "mse"
      lambda_one_se <- if (!is.null(obj$lambda_one_se)) {
        isTRUE(obj$lambda_one_se)
      } else {
        objective == "mse" && selector %in% c("loo", "blocked")
      }
      if (objective != "mse" && selector != "blocked") {
        stop(
          paste0(
            "lambda_objective='", objective, "' requires ",
            "lambda_selection='blocked'."
          )
        )
      }
      selected_index <- switch(
        selector,
        fixed = {
          if (length(lambdas) != 1L) {
            stop("lambda_selection='fixed' requires exactly one lambda.")
          }
          1L
        },
        gcv = .feature_rsa_ridge_select_lambda(
          .feature_rsa_ridge_gcv_scores(decomposition, lambdas),
          lambdas
        ),
        loo = .feature_rsa_ridge_select_lambda(
          .feature_rsa_ridge_loo_mse(decomposition, lambdas),
          lambdas,
          one_se = lambda_one_se
        ),
        blocked = {
          segments <- .feature_rsa_block_segments(
            obj$design$block_var,
            n_rows = nrow(Fsub),
            observation_indices = observation_indices
          )
          scores <- switch(
            objective,
            mse = .feature_rsa_ridge_block_mse(
              Fsub, X, lambdas, segments,
              feature_standardize = feature_standardize
            ),
            pattern_discrimination =
              -.feature_rsa_ridge_block_pattern_discrimination(
                Fsub, X, lambdas, segments,
                feature_standardize = feature_standardize
              ),
            pattern_rank_percentile =
              1 - .feature_rsa_ridge_block_pattern_rank(
                Fsub, X, lambdas, segments,
                feature_standardize = feature_standardize
              ),
            stop(sprintf("Unknown ridge lambda objective: %s", objective))
          )
          .feature_rsa_ridge_select_lambda(
            scores,
            lambdas,
            one_se = lambda_one_se
          )
        },
        stop(sprintf("Unknown ridge lambda selection method: %s", selector))
      )
      selected_lambda <- lambdas[[selected_index]]
      model <- .feature_rsa_ridge_fit_kernel(
        sf$X_sc,
        sx$X_sc,
        selected_lambda,
        decomposition = decomposition
      )

      list(
        trained_model = model,
        x_mean = sx$mean,
        x_sd = sx$sd,
        f_mean = sf$mean,
        f_sd = sf$sd,
        lambda = selected_lambda,
        effective_df = model$effective_df,
        lambda_selection = selector,
        lambda_objective = objective,
        lambda_one_se = lambda_one_se,
        lambda_at_min_boundary = if (selector == "fixed") {
          NA_real_
        } else {
          as.numeric(selected_index == 1L)
        },
        lambda_at_max_boundary = if (selector == "fixed") {
          NA_real_
        } else {
          as.numeric(selected_index == length(lambdas))
        },
        ncomp = NA_real_
      )
    }, error = function(e) {
      error_msg <- sprintf(
        "train_model (Ridge): Error during training - %s",
        e$message
      )
      futile.logger::flog.error(error_msg)
      list(error = error_msg)
    })

    if (!is.null(ridge_result$error)) {
      result$error <- ridge_result$error
      return(result)
    }
    result <- c(result, ridge_result)

  } else if (obj$method == "glmnet") {
    #
    # ---- GLMNet Train ----
    #
    glm_result <- tryCatch({
        # Standardize X and F
        sx <- .standardize(X)
        sf <- .standardize(Fsub, feature_standardize)
        
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
                require_package("glmnet", "for regularized regression in feature RSA")
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
        
        # cv.glmnet already fits a full-data path. Reuse it after successful
        # CV; otherwise fit the requested path once as the primary/fallback.
        final_fit <- if (run_cv && !is.null(cv_results)) {
          cv_results$glmnet.fit
        } else {
          glmnet::glmnet(
            x = sf$X_sc,
            y = sx$X_sc,
            family = "mgaussian",
            alpha = obj$alpha,
            lambda = lambda_to_use,
            standardize = FALSE,
            intercept = TRUE
          )
        }

        # Determine lambda used for prediction
        lambda_used_for_pred <- if (run_cv && !is.null(cv_results)) {
           lambda_to_use # lambda.min from successful CV
        } else if (!is.null(lambda_to_use)) {
           # User supplied lambda explicitly. If a sequence was supplied,
           # predict at the weakest penalty rather than the strongest one
           # to avoid the degenerate all-zero solution at lambda_max.
           if (length(lambda_to_use) == 1L) {
             lambda_to_use
           } else if (length(final_fit$lambda) > 0) {
             chosen_lambda <- final_fit$lambda[length(final_fit$lambda)]
             futile.logger::flog.warn(
               paste0(
                 "train_model (GLMNet): lambda supplied as a sequence of length %d ",
                 "with cv_glmnet=FALSE; using the smallest fitted lambda (%.4e) ",
                 "for prediction. Supply a scalar lambda or set cv_glmnet=TRUE ",
                 "for automatic selection."
               ),
               length(final_fit$lambda),
               chosen_lambda
             )
             chosen_lambda
           } else {
             min(lambda_to_use)
           }
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
  # Parallel worker specs retain one target matrix in the fold cache and drop
  # the duplicate design copy. Controller specs continue to use design targets.
  obj$design$targets %||% obj$.cv_fold_cache$y
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
      test_index    = list(NULL),
      result        = list(NULL),
      performance   = list(NULL),
      error         = TRUE, 
      error_message = error_message
    ))
  }
  
  Xobs  <- as.matrix(context$test)
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
      test_index    = list(NULL),
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
      test_index    = list(NULL),
      result        = list(NULL),
      performance   = list(NULL),
      error         = TRUE, 
      error_message = "Prediction returned NULL unexpectedly."
    ))
  }
  
  # Fold metrics are intentionally deferred. merge_results() evaluates the
  # complete out-of-fold prediction exactly once, so computing geometry here
  # would be redundant and the result would be discarded.
  retained_result <- if (identical(obj$method, "ridge")) {
    list(
      lambda = result$lambda,
      effective_df = result$effective_df,
      lambda_at_min_boundary = result$lambda_at_min_boundary,
      lambda_at_max_boundary = result$lambda_at_max_boundary,
      ncomp = NA_real_
    )
  } else {
    list(ncomp = result$ncomp)
  }
  
  tibble::tibble(
    observed    = list(Xobs),
    predicted   = list(Xpred),
    test_index  = list(context$test_ind),
    # Prediction is complete, so the trained model and its coefficient path no
    # longer need to survive until fold merging.
    result      = list(retained_result),
    performance = list(NULL),
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
  
  oof <- attr(result_set, "feature_rsa_oof", exact = TRUE)
  observed_list <- result_set$observed
  predicted_list <- result_set$predicted
  test_index_list <- result_set$test_index
  combined_fold_id <- NULL
  
  # Get the list of results from each fold (contains ncomp)
  fold_results_list <- result_set$result
  
  if (is.list(oof) && is.matrix(oof$observed) && is.matrix(oof$predicted)) {
    combined_observed <- oof$observed
    combined_predicted <- oof$predicted
    combined_test_index <- oof$test_index
    if (!is.null(oof$fold_id) &&
        length(oof$fold_id) == nrow(combined_observed)) {
      combined_fold_id <- oof$fold_id
    }
  } else {
    combined_observed <- do.call(rbind, observed_list)
    combined_predicted <- do.call(rbind, predicted_list)
    combined_test_index <- NULL
    combined_fold_id <- unlist(lapply(seq_along(observed_list), function(i) {
      rep.int(i, nrow(as.matrix(observed_list[[i]])))
    }), use.names = FALSE)
    if (!is.null(test_index_list) && length(test_index_list) == length(observed_list) &&
        all(vapply(test_index_list, Negate(is.null), logical(1)))) {
      combined_test_index <- unlist(test_index_list, use.names = FALSE)
      ord <- order(combined_test_index)
      combined_test_index <- combined_test_index[ord]
      combined_observed <- combined_observed[ord, , drop = FALSE]
      combined_predicted <- combined_predicted[ord, , drop = FALSE]
      combined_fold_id <- combined_fold_id[ord]
    }
  }

  # Older preallocated result objects did not carry fold_id. Reconstruct it
  # from their per-fold test indices when possible rather than reverting to a
  # pooled candidate set.
  if (is.null(combined_fold_id) &&
      !is.null(test_index_list) &&
      length(test_index_list) == length(fold_results_list) &&
      all(vapply(test_index_list, Negate(is.null), logical(1)))) {
    candidate_index <- unlist(test_index_list, use.names = FALSE)
    candidate_fold <- rep.int(seq_along(test_index_list), lengths(test_index_list))
    if (length(candidate_index) == nrow(combined_observed)) {
      ord <- order(candidate_index)
      combined_fold_id <- candidate_fold[ord]
    }
  }
  if (is.null(combined_fold_id) ||
      length(combined_fold_id) != nrow(combined_observed)) {
    stop(
      "feature RSA cross-validation results do not identify outer test folds.",
      call. = FALSE
    )
  }
  
  scalar_from_fold <- function(fold, name) {
    value <- fold[[name]]
    if (is.numeric(value) && length(value) == 1L && is.finite(value)) {
      as.numeric(value)
    } else {
      NA_real_
    }
  }
  
  # Now we do permutations (if nperm>0 in the model spec)
  perf <- evaluate_model.feature_rsa_model(
    object    = obj,
    predicted = combined_predicted,
    observed  = combined_observed,
    nperm     = obj$nperm,
    save_distributions = obj$save_distributions,
    compute_rdm_vectors = isTRUE(obj$return_rdm_vectors),
    fold_id = combined_fold_id
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
    perf$mean_voxelwise_temporal_cor
  )
  base_names <- c(
    "pattern_correlation", "pattern_discrimination", "pattern_rank_percentile",
    "rdm_correlation",
    "voxel_correlation", "mse", "r_squared",
    "mean_voxelwise_temporal_cor"
  )
  if (identical(obj$method, "ridge")) {
    fold_lambdas <- vapply(
      fold_results_list,
      scalar_from_fold,
      numeric(1L),
      name = "lambda"
    )
    fold_effective_df <- vapply(
      fold_results_list,
      scalar_from_fold,
      numeric(1L),
      name = "effective_df"
    )
    finite_lambdas <- fold_lambdas[is.finite(fold_lambdas)]
    finite_effective_df <- fold_effective_df[is.finite(fold_effective_df)]
    fold_min_boundary <- vapply(
      fold_results_list,
      scalar_from_fold,
      numeric(1L),
      name = "lambda_at_min_boundary"
    )
    fold_max_boundary <- vapply(
      fold_results_list,
      scalar_from_fold,
      numeric(1L),
      name = "lambda_at_max_boundary"
    )
    finite_min_boundary <- fold_min_boundary[is.finite(fold_min_boundary)]
    finite_max_boundary <- fold_max_boundary[is.finite(fold_max_boundary)]
    diagnostics <- c(
      median_lambda = if (length(finite_lambdas)) {
        stats::median(finite_lambdas)
      } else {
        NA_real_
      },
      mean_effective_df = if (length(finite_effective_df)) {
        mean(finite_effective_df)
      } else {
        NA_real_
      },
      mean_nonintercept_df = if (length(finite_effective_df)) {
        mean(finite_effective_df) - 1
      } else {
        NA_real_
      },
      lambda_min_boundary_fraction = if (length(finite_min_boundary)) {
        mean(finite_min_boundary)
      } else {
        NA_real_
      },
      lambda_max_boundary_fraction = if (length(finite_max_boundary)) {
        mean(finite_max_boundary)
      } else {
        NA_real_
      }
    )
  } else {
    diagnostics <- c(
      ncomp = scalar_from_fold(fold_results_list[[1L]], "ncomp")
    )
  }
  base_metrics <- c(base_metrics, diagnostics)
  base_names <- c(base_names, names(diagnostics))

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

  predictor_obj <- NULL
  keep_rdm <- isTRUE(obj$return_rdm_vectors)
  keep_predictions <- isTRUE(obj$return_predictions)
  if (keep_predictions) {
    max_mb <- obj$max_retained_mb
    if (!is.null(max_mb) && is.numeric(max_mb) && is.finite(max_mb)) {
      overflow <- obj$prediction_overflow %||% "error"
      retention <- .feature_rsa_validate_prediction_retention(
        n_obs = nrow(combined_predicted),
        n_voxels = ncol(combined_predicted),
        return_predictions = TRUE,
        max_retained_mb = max_mb,
        prediction_overflow = overflow
      )
      keep_predictions <- isTRUE(retention$return_predictions)
    }
  }
  if (keep_rdm || keep_predictions) {
    # Store `predicted`/`observed` even when NULL so `$predicted` cannot
    # partial-match `predicted_rdm_vec`.
    predictor_obj <- list(
      predicted_rdm_vec = if (keep_rdm) perf$predicted_rdm_vec else NULL,
      observed_rdm_vec  = if (keep_rdm) perf$observed_rdm_vec else NULL,
      observation_index = combined_test_index,
      fold_id = combined_fold_id,
      n_obs = nrow(combined_predicted),
      predicted = if (keep_predictions) combined_predicted else NULL,
      observed = if (keep_predictions) combined_observed else NULL,
      voxel_index = if (keep_predictions) indices else NULL
    )
  }
  
  tibble::tibble(
    result      = list(list(predictor = predictor_obj)),
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
    perf <- result_tbl$performance[[1]]
    if (is.matrix(perf)) {
      perf <- setNames(as.numeric(perf[1, ]), colnames(perf))
    }
    roi_result(
      metrics = perf,
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
    "mean_voxelwise_temporal_cor",
    if (identical(model$method, "ridge")) {
      c(
        "median_lambda", "mean_effective_df", "mean_nonintercept_df",
        "lambda_min_boundary_fraction", "lambda_max_boundary_fraction"
      )
    } else {
      "ncomp"
    }
  )
  setNames(rep("scalar", length(nms)), nms)
}


#' Regional MVPA for `feature_rsa_model` Objects
#'
#' Out-of-fold predicted patterns are retained in `$fits` when
#' `feature_rsa_model(..., return_predictions = TRUE)`. Extract them with
#' \code{\link{feature_rsa_predictions}} rather than a classification-style
#' `prediction_table`.
#'
#' @rdname run_regional-methods
#' @export
run_regional.feature_rsa_model <- function(model_spec, region_mask,
                                           backend = c("default", "shard", "auto"),
                                           ...) {
  n_voxels <- .feature_rsa_region_voxel_count(model_spec$dataset, region_mask)
  model_spec <- .feature_rsa_apply_prediction_retention(model_spec, n_voxels)
  dots <- list(...)
  # Feature RSA predictions are observation-by-voxel matrices. The default
  # regional combiner expects scalar observed/predicted vectors, so skip it.
  dots$return_predictions <- FALSE
  backend <- match.arg(backend)
  do.call(run_regional_base, c(
    list(model_spec, region_mask, backend = backend),
    dots
  ))
}


#' @rdname run_searchlight
#' @export
run_searchlight.feature_rsa_model <- function(model_spec, radius,
                                              method = c("standard", "randomized", "resampled"),
                                              niter = 4,
                                              backend = c("default", "shard", "auto"),
                                              ...) {
  if (isTRUE(model_spec$return_predictions)) {
    stop(
      "feature_rsa_model(return_predictions = TRUE) is not supported for ",
      "run_searchlight(): overlapping spheres would retain multiple copies of ",
      "each voxel's out-of-fold predictions. Use run_regional() instead.",
      call. = FALSE
    )
  }
  run_searchlight.default(
    model_spec, radius, method = method, niter = niter, backend = backend, ...
  )
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
  cat(crayon::bold(crayon::green("Feature standardization: ")),
      x$feature_standardize %||% "scale", "\n")
  cat(
    crayon::bold(crayon::green("Return predictions:     ")),
    if (isTRUE(x$return_predictions)) "Yes" else "No",
    "\n"
  )
  if (!is.null(x$retention_notice)) {
    cat(crayon::bold(crayon::yellow("Retention notice:       ")),
        x$retention_notice, "\n")
  }
  
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
    objective <- x$ncomp_objective %||% "mse"
    one_se <- if (!is.null(x$ncomp_one_se)) {
      isTRUE(x$ncomp_one_se)
    } else {
      objective == "mse" && sel %in% c("loo", "blocked")
    }
    cat(crayon::bold(crayon::blue("Component selection:    ")), sel, "\n")
    cat(crayon::bold(crayon::blue("Component objective:    ")), objective, "\n")
    if (sel %in% c("loo", "blocked")) {
      cat(
        crayon::bold(crayon::blue("One-SE rule:            ")),
        if (one_se) "Yes" else "No",
        "\n"
      )
    }
    if (sel == "pve" && !is.null(x$pve_threshold)) {
      cat(crayon::bold(crayon::blue("PVE threshold:          ")), x$pve_threshold, "\n")
    }
  } else if (x$method == "ridge") {
    selector <- x$lambda_selection %||% "gcv"
    objective <- x$lambda_objective %||% "mse"
    one_se <- if (!is.null(x$lambda_one_se)) {
      isTRUE(x$lambda_one_se)
    } else {
      objective == "mse" && selector %in% c("loo", "blocked")
    }
    cat(crayon::bold(crayon::blue("Lambda selection:       ")), selector, "\n")
    cat(
      crayon::bold(crayon::blue("Lambda objective:       ")),
      objective,
      "\n"
    )
    cat(
      crayon::bold(crayon::blue("One-SE rule:            ")),
      if (one_se) "Yes" else "No",
      "\n"
    )
    ridge_lambdas <- .feature_rsa_ridge_lambdas(x$lambda)
    lambda_str <- if (length(ridge_lambdas) > 3L) {
      paste0(
        format(ridge_lambdas[[1L]], scientific = TRUE),
        " ... ",
        format(ridge_lambdas[[length(ridge_lambdas)]], scientific = TRUE),
        " (", length(ridge_lambdas), " values)"
      )
    } else {
      paste(format(ridge_lambdas, scientific = TRUE), collapse = ", ")
    }
    cat(crayon::bold(crayon::blue("Normalized lambda:      ")), lambda_str, "\n")
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
