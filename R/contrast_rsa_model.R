#' Constructor for contrast_rsa_model
#'
#' Creates a contrast_rsa_model specification object, which encapsulates the necessary
#' parameters and design information for a Multi-Dimensional Signed Representational
#' Voxel Encoding (MS-ReVE) style analysis.
#'
#' @param dataset An object of class \code{mvpa_dataset}, containing the neuroimaging
#'   data and associated metadata.
#' @param design An object of class \code{msreve_design}, containing the underlying
#'   \code{mvpa_design} and the contrast matrix. Created by \code{\\link{msreve_design}}.
#' @param estimation_method Character string specifying the method to estimate
#'   cross-validated condition means (\code{Û}) or distances. Supported:
#'   \itemize{
#'     \item \code{"average"}: Simple mean of training samples per condition (for \code{Û}).
#'     \item \code{"L2_norm"}: Like \code{"average"}, but \code{Û} rows are L2-normalized.
#'     \item \code{"crossnobis"}: Computes unbiased squared Euclidean distances directly using the Crossnobis method. Results in a distance vector for RSA, not a \code{Û} matrix for \(G_empirical\) construction in the same way. \code{U_hat} for \eqn{\Delta} calculation is still computed using "average" method internally when this is selected. Requires `return_folds=TRUE` from `compute_crossvalidated_means_sl`.
#'   }
#'   Passed to \code{\link{compute_crossvalidated_means_sl}} (for "average", "L2_norm", or to get per-fold means for "crossnobis").
#' @param regression_type Character string specifying the method for the RSA regression
#'   (regressing empirical RDM/Second Moment Matrix onto contrast RDMs). Options align
#'   with \code{rsa_model}: \code{"pearson"}, \code{"spearman"}, \code{"lm"}, \code{"rfit"}.
#'   Additional options for ridge regression: \code{"ridge_hkb"} (Hoerl-Kennard-Baldwin lambda).
#'   Default is \code{"lm"}.
#' @param output_metric Character vector specifying one or more output metrics to compute.
#'   Multiple metrics can be requested simultaneously and will be returned in a named list.
#'   Duplicates are removed while preserving the order of first occurrence.
#'   Supported options:
#'   \itemize{
#'     \item \code{"beta_delta"}: The product of the RSA regression coefficient (beta_q)
#'           and the voxel's projection onto the contrast (delta_q,v). This is the
#'           primary signed contribution metric.
#'     \item \code{"beta_only"}: Only the RSA regression coefficient (beta_q).
#'     \item \code{"delta_only"}: Only the voxel's projection onto the contrast (delta_q,v).
#'     \item \code{"recon_score"}: Voxel-specific RDM reconstruction score (r_v), correlating the RDM implied by the voxel's loadings with the empirical RDM.
#'     \item \code{"beta_delta_norm"}: Similar to `beta_delta`, but uses the L2-normalized voxel contribution vector (delta_q,v). Requires `normalize_delta=TRUE` to be meaningful.
#'     \item \code{"beta_delta_reliable"}: Reliability-weighted contributions, \eqn{\rho_{q,v} \beta_q \Delta_{q,v}}, where \eqn{\rho_{q,v}} reflects fold-wise stability.
#'     \item \code{"composite"}: Sum of beta-weighted, L2-normalized voxel contributions (Σ_q β_q ~Δ_q,v). Represents the net projection onto the positive diagonal of the contrast space. Interpretation requires caution if contrasts are not orthonormal.
#'   }
#'   Default is \code{c("beta_delta")}.
#' @param check_collinearity Logical, whether to check for collinearity among contrast RDMs
#'   when using \code{regression_type = "lm"}. Default is \code{FALSE}.
#' @param normalize_delta Logical. If TRUE, the voxel contribution vector (delta_q,v)
#'   is normalized to unit L2 norm before being potentially multiplied by beta_q
#'   for the "beta_delta" output metric. Default is \code{FALSE}.
#' @param allow_nonorth_composite Logical. If FALSE (default) the composite metric will return NA when the contrast matrix
#'   is not orthonormal, to avoid mis-interpretation. If TRUE, the composite score is returned regardless, but a warning
#'   is still emitted.
#' @param calc_reliability Logical. If TRUE, voxel-wise contribution reliability (ρ
#'   values) are estimated across cross-validation folds during training and can
#'   be incorporated into the output metrics. Default is \code{FALSE}.
#' @param whitening_matrix_W Optional V x V numeric matrix (voxels x voxels). Required if
#'   `estimation_method = "crossnobis"` and Mahalanobis (rather than Euclidean)
#'   distances are desired. This matrix (e.g., \eqn{\Sigma_{noise}^{-1/2}}) is passed to
#'   `compute_crossvalidated_means_sl` to whiten per-fold estimates before distance calculation.
#'   Default is `NULL` (Euclidean distances for Crossnobis).
#' @param ... Additional arguments passed to \code{create_model_spec}.
#'
#' @details
#' This model is designed for MS-ReVE style analyses where the goal is to understand
#' how different predefined contrasts contribute to the representational structure
#' observed in neural data, particularly at a fine-grained (e.g., voxel) level.
#' It involves:
#' 1. Estimating cross-validated condition mean patterns (Û) or distances (for Crossnobis).
#' 2. Constructing an empirical second-moment matrix (Ĝ) from Û or using the Crossnobis distances directly.
#' 3. Creating theoretical second-moment matrices (RDMs) from each contrast vector.
#' 4. Regressing the vectorized empirical RDM/distances onto the vectorized contrast RDMs to get β coefficients.
#' 5. Projecting voxel patterns (from Û) onto the contrast space to get Δ (delta) values.
#' 6. Combining β and Δ to form the final output metric (e.g., beta_delta).
#'
#' **Cross-Validation Compatibility:**
#' The `estimation_method` relies on `compute_crossvalidated_means_sl` which, in turn, requires
#' the cross-validation specification (`cv_spec` derived from `mvpa_design$crossval`)
#' to provide a deterministic, partition-based set of training indices for each fold.
#' Therefore, cross-validation schemes like `bootstrap_blocked_cross_validation` (which use
#' resampling with replacement for training folds) are **not suitable** for use with this model,
#' as they do not align with the assumptions of `compute_crossvalidated_means_sl`.
#' Schemes like `blocked_cross_validation`, `kfold_cross_validation`, and `custom_cross_validation`
#' (that define clear partitions) are appropriate.
#' For `twofold_blocked_cross_validation` and `sequential_blocked_cross_validation`,
#' their compatibility also depends on whether their `train_indices` methods can deterministically
#' define training sets for each fold iterated by `compute_crossvalidated_means_sl`.
#'
#' @return An object of class \code{contrast_rsa_model}, \code{mvpa_model_spec},
#'   \code{model_spec}, and \code{list}.
#'
#' @export
#' @importFrom assertthat assert_that
#' @seealso \code{\\link{msreve_design}}, \code{\\link{train_model.contrast_rsa_model}}, \code{\\link{run_searchlight}}
#' @examples
#' # --- Minimal Setup ---
#' # 1. Create dummy data and an mvpa_dataset
#'
#'   # Dummy data: 16 samples, 10 voxels, 4 conditions, 2 runs
#'   set.seed(123)
#'   n_samples <- 16
#'   n_voxels <- 10
#'   n_conditions <- 4 # condA, condB, condC, condD
#'   n_runs <- 2
#'
#'   dummy_sl_data <- matrix(rnorm(n_samples * n_voxels), n_samples, n_voxels)
#'   colnames(dummy_sl_data) <- paste0("V", 1:n_voxels)
#'
#'   dummy_mask <- neuroim2::NeuroVol(array(1, c(2,2,2)), neuroim2::NeuroSpace(c(2,2,2)))
#'
#'   condition_labels <- factor(rep(paste0("cond", LETTERS[1:n_conditions]), each = n_samples / n_conditions))
#'   run_labels <- factor(rep(1:n_runs, each = n_samples / n_runs))
#'
#'   # Create mvpa_dataset (without Y and block_var)
#'   mvpa_dat <- rMVPA::mvpa_dataset(
#'     train_data = dummy_sl_data,
#'     mask = dummy_mask
#'   )
#'
#'   # Create mvpa_design
#'   mvpa_des <- rMVPA::mvpa_design(
#'     train_design = data.frame(condition = condition_labels, run = run_labels),
#'     y_train = ~condition,
#'     block_var = ~run
#'   )
#'
#'   K <- mvpa_des$ncond # Use mvpa_des here
#'
#'   C_mat <- matrix(0, nrow = K, ncol = 2)
#'   rownames(C_mat) <- levels(mvpa_des$Y) # Use mvpa_des here
#'   C_mat["condA", 1] <- 1; C_mat["condB", 1] <- 1
#'   C_mat["condC", 1] <- -1; C_mat["condD", 1] <- -1
#'   C_mat["condA", 2] <- 1; C_mat["condB", 2] <- -1
#'   C_mat <- scale(C_mat, center = TRUE, scale = FALSE)
#'   colnames(C_mat) <- c("AB_vs_CD", "A_vs_B")
#'
#'   msreve_des <- rMVPA::msreve_design(
#'     mvpa_design = mvpa_des, # Use mvpa_des here
#'     contrast_matrix = C_mat
#'   )
#'
#'   # --- Example 1: Basic contrast_rsa_model ---
#'   model_basic <- contrast_rsa_model(
#'     dataset = mvpa_dat,
#'     design = msreve_des
#'   )
#'   print(model_basic)
#'
#'   # --- Example 1b: Requesting multiple metrics ---
#'   model_multi_metric <- contrast_rsa_model(
#'     dataset = mvpa_dat,
#'     design = msreve_des,
#'     output_metric = c("beta_delta", "recon_score", "beta_only")
#'   )
#'   print(model_multi_metric)
#'
#'   # --- Example 2: Using L2_norm for U_hat and normalize_delta ---
#'   model_l2_norm_delta <- contrast_rsa_model(
#'     dataset = mvpa_dat,
#'     design = msreve_des,
#'     estimation_method = "L2_norm",
#'     normalize_delta = TRUE,
#'     output_metric = "beta_delta_norm"
#'   )
#'   print(model_l2_norm_delta)
#'
#'   # --- Example 3: Ridge Regression (HKB) ---
#'   model_ridge <- contrast_rsa_model(
#'     dataset = mvpa_dat,
#'     design = msreve_des,
#'     regression_type = "ridge_hkb"
#'   )
#'   print(model_ridge)
#'
#'   # --- Example 4: Reconstruction Score Output ---
#'   model_recon <- contrast_rsa_model(
#'     dataset = mvpa_dat,
#'     design = msreve_des,
#'     output_metric = "recon_score"
#'   )
#'   print(model_recon)
#'
#'   # --- Example 5: Composite Score Output ---
#'   C_mat_ortho <- rMVPA::orthogonalize_contrasts(C_mat)
#'   msreve_des_ortho <- rMVPA::msreve_design(
#'       mvpa_design = mvpa_des, # Use mvpa_des here
#'       contrast_matrix = C_mat_ortho
#'   )
#'   print(paste("Is contrast matrix orthonormal:", attr(msreve_des_ortho, "is_orthonormal")))
#'
#'   model_composite <- contrast_rsa_model(
#'     dataset = mvpa_dat,
#'     design = msreve_des_ortho,
#'     output_metric = "composite",
#'     normalize_delta = TRUE 
#'   )
#'   print(model_composite)
#'
#'   # --- Example 6: Crossnobis estimation_method ---
#'   # This only shows setting the method. Actual training would require passing
#'   # a pre-computed whitening_matrix_W to compute_crossvalidated_means_sl,
#'   # which is called by train_model.contrast_rsa_model.
#'   model_crossnobis <- contrast_rsa_model(
#'       dataset = mvpa_dat,
#'       design = msreve_des,
#'       estimation_method = "crossnobis"
#'   )
#'   print(model_crossnobis)
contrast_rsa_model <- function(dataset,
                               design,
                               estimation_method = "average",
                               regression_type = "lm",
                               output_metric = c("beta_delta"),
                               check_collinearity = FALSE,
                               normalize_delta = FALSE,
                               allow_nonorth_composite = FALSE,
                               calc_reliability = FALSE,
                               whitening_matrix_W = NULL,
                               ...) {

  # --- Input Checks ---
  assert_that(inherits(dataset, "mvpa_dataset"))
  assert_that(inherits(design, "msreve_design"))

  estimation_method <- match.arg(estimation_method, c("average", "L2_norm", "crossnobis"))
  regression_type   <- match.arg(regression_type, c("pearson", "spearman", "lm", "rfit", "ridge_hkb"))
  
  # Validate output_metric
  allowed_metrics <- c("beta_delta", "beta_only", "delta_only", "recon_score", "beta_delta_norm", "beta_delta_reliable", "composite")
  if (!is.character(output_metric) || !all(output_metric %in% allowed_metrics)) {
    rlang::abort(paste0("`output_metric` must be a character vector containing only allowed metrics: ",
                       paste(allowed_metrics, collapse=", ")))
  }
  if (length(output_metric) == 0) {
    rlang::abort("`output_metric` cannot be an empty vector.")
  }
  # Ensure no duplicates, maintain user order for the first occurrences
  output_metric <- unique(output_metric)

  assert_that(is.logical(check_collinearity), length(check_collinearity) == 1)
  assert_that(is.logical(normalize_delta), length(normalize_delta) == 1)
  assert_that(is.logical(allow_nonorth_composite), length(allow_nonorth_composite)==1)
  assert_that(is.logical(calc_reliability), length(calc_reliability) == 1)
  
  if (!is.null(whitening_matrix_W)) {
    assert_that(is.matrix(whitening_matrix_W), is.numeric(whitening_matrix_W),
                msg = "`whitening_matrix_W` if provided, must be a numeric matrix.")
    # Further dimension checks will occur in compute_crossvalidated_means_sl against sl_data
  }

  # Check for collinearity when using linear regression and
  # `check_collinearity = TRUE`. Predictor matrix X is built from
  # vectorized contrast RDMs and the QR rank is tested; an abort is
  # triggered if the rank is less than the number of contrasts.
  if (regression_type == "lm" && check_collinearity) {
      # warning("Collinearity check requested but not yet implemented for contrast_rsa_model.")
      C <- design$contrast_matrix
      if (is.null(C) || !is.matrix(C) || ncol(C) < 2) {
          warning("Cannot check collinearity: Contrast matrix is missing, not a matrix, or has fewer than 2 columns.")
      } else {
          Q <- ncol(C)
          K <- nrow(C)
          contrast_names <- colnames(C)
          if (is.null(contrast_names)) contrast_names <- paste0("Contrast", 1:Q)
          
          # Create predictor matrix Xmat from contrast RDMs (lower triangles)
          # This mimics the logic in train_model.contrast_rsa_model but without data dependency
          # Note: This doesn't account for potential block exclusions, but checks the raw contrasts
          Xmat_list <- lapply(1:Q, function(q) {
              contrast_rdm <- C[, q, drop = FALSE] %*% t(C[, q, drop = FALSE])
              pred_vec <- contrast_rdm[lower.tri(contrast_rdm)]
              pred_vec
          })
          
          Xmat <- tryCatch({
              do.call(cbind, Xmat_list)
          }, error = function(e) {
              warning(paste("Could not bind contrast predictors into matrix for collinearity check:", e$message))
              NULL # Signal failure
          })
          
          if (!is.null(Xmat)) {
              # Check rank vs number of columns
              qr_Xmat <- qr(Xmat)
              if (qr_Xmat$rank < ncol(Xmat)) {
                   # Identify dependent columns
                   pivot_indices <- qr_Xmat$pivot
                   original_cols <- contrast_names
                   dependent_indices <- pivot_indices[(qr_Xmat$rank + 1):length(pivot_indices)]
                   dependent_names <- original_cols[dependent_indices]
                   
                   rlang::abort(paste0(
                      "Collinearity detected among contrast RDMs (when vectorized for regression): ",
                      "Rank (", qr_Xmat$rank, ") is less than the number of contrasts (", ncol(Xmat), ").",
                      "\nLinearly dependent contrast(s): ", paste(dependent_names, collapse=", "), ".",
                      "\nConsider orthogonalizing contrasts using `orthogonalize_contrasts()` or simplifying the model."
                   ))
              }
          }
      }
  }

  # --- Create Model Specification ---
  # Use the internal constructor helper
  create_model_spec(
    "contrast_rsa_model", # Class name
    dataset = dataset,
    design = design,
    estimation_method = estimation_method,
    regression_type = regression_type,
    output_metric = output_metric,
    check_collinearity = check_collinearity,
    normalize_delta = normalize_delta,
    allow_nonorth_composite = allow_nonorth_composite,
    calc_reliability = calc_reliability,
    whitening_matrix_W = whitening_matrix_W, # Store W in the spec
    ...
  )
}


#' @export
#' @method print contrast_rsa_model
#' @importFrom crayon bold cyan yellow white magenta italic blue
print.contrast_rsa_model <- function(x, ...) {
  header_style  <- bold$cyan
  section_style <- yellow
  info_style    <- white
  param_style   <- magenta
  class_style   <- italic$blue

  cat("\n", header_style("█▀▀ Contrast RSA Model Specification ▀▀█"), "\n\n")

  cat(section_style("├─ Dataset"), "\n")
  cat(info_style("│  └─ Class: "), class_style(class(x$dataset)[1]), "\n")

  cat(section_style("├─ Design (`msreve_design`)"), "\n")
  # Use ifelse for potentially NULL name
  design_name <- ifelse(is.null(x$design$name), "-", x$design$name)
  cat(info_style("│  ├─ Name: "), info_style(design_name), "\n")
  
  # Safely get dimensions of contrast_matrix
  if (!is.null(x$design$contrast_matrix) && is.matrix(x$design$contrast_matrix)) {
    dims_cm <- dim(x$design$contrast_matrix)
    cat(info_style("│  ├─ Contrast Matrix Dims: "), param_style(paste0(dims_cm[1], " Cond × ", dims_cm[2], " Contrasts")), "\n")
  } else {
    cat(info_style("│  ├─ Contrast Matrix Dims: "), param_style("Not a valid matrix or NULL"), "\n")
  }
  
  cat(info_style("│  └─ MVPA Design Class: "), class_style(class(x$design$mvpa_design)[1]), "\n")

  cat(section_style("└─ Analysis Parameters"), "\n")
  cat(info_style("   ├─ Estimation Method (Û/d): "), param_style(x$estimation_method), "\n")
  if (identical(x$estimation_method, "crossnobis")) {
    whiten_status <- if (!is.null(x$whitening_matrix_W)) "Provided (for Mahalanobis)" else "None (for Euclidean)"
    cat(info_style("   │  └─ Whitening Matrix (W): "), param_style(whiten_status), "\n")
  }
  cat(info_style("   ├─ Regression Type (β):   "), param_style(x$regression_type), "\n")
  cat(info_style("   ├─ Calc Reliability:     "), param_style(ifelse(x$calc_reliability, "TRUE", "FALSE")), "\n")
  cat(info_style("   └─ Output Metric(s):      "), param_style(paste(x$output_metric, collapse = ", ")), "\n")

  cat("\n")
  invisible(x)
}


#' Train method for contrast_rsa_model
#'
#' This function implements the core logic for the MS-ReVE analysis within a single
#' searchlight or region.
#'
#' @param obj An object of class \code{contrast_rsa_model}.
#' @param sl_data The data matrix for the current searchlight (samples x voxels).
#' @param sl_info A list containing information about the current searchlight,
#'   including \code{center_local_id}, the column index of the center voxel within \code{sl_data}.
#' @param cv_spec The cross-validation specification object (e.g., from \code{\\link{blocked_cross_validation}}).
#' @param ... Additional arguments (currently ignored).
#'
#' @return A named list where each element corresponds to a requested `output_metric` from the `obj$output_metric` vector.
#'   Each element is:
#'   \itemize{
#'     \item For metrics like "beta_delta", "beta_only", "delta_only": A Q-length named vector where
#'           values are indexed by contrast and names match the contrast_matrix column names (Q = number of contrasts)
#'     \item For metrics like "recon_score", "composite": A single numeric value
#'   }
#'   The list will have an attribute "na_reason" if any metric calculation failed, which can be used for diagnostics.
#'   
#'   For example, if `obj$output_metric = c("beta_delta", "recon_score")`, the returned list will have 
#'   two elements: `$beta_delta` (a Q-length vector) and `$recon_score` (a single value).
#'
#' @examples
#' # This example shows the structure of the returned list but doesn't actually run the function
#' # For a multi-metric model: output_metric = c("beta_delta", "recon_score", "beta_only")
#' @rdname train_model
#' @param obj An object of class \code{contrast_rsa_model}.
#' @param sl_data The data matrix for the current searchlight (samples x voxels).
#' @param sl_info A list containing information about the current searchlight, including \code{center_local_id}.
#' @param cv_spec The cross-validation specification.
#' @param ... Additional arguments (currently ignored).
#' @return A named list where each element corresponds to a requested metric from \code{obj$output_metric}.
#' @importFrom stats cor dist lm coef sd terms
#' @importFrom rlang abort
#' @method train_model contrast_rsa_model
#' @keywords internal
#' @export
train_model.contrast_rsa_model <- function(obj, sl_data, sl_info, cv_spec, ...) {

  # --- Input Checks & Setup ---
  if (is.null(sl_info) || is.null(sl_info$center_local_id) || is.na(sl_info$center_local_id)) {
    stop(paste0("`train_model.contrast_rsa_model` requires a valid center voxel ID (`sl_info$center_local_id`). ",
                "This model is typically used in searchlight analyses. For regional analyses, ",
                "a different aggregation approach might be needed."))
  }
  
  center_idx <- sl_info$center_local_id
  # This check should now be more robust due to the abort above if center_local_id is NA
  if (center_idx < 1 || center_idx > ncol(sl_data)) {
      rlang::abort(paste0("`sl_info$center_local_id` (", center_idx, ") is out of bounds for `sl_data` columns (", ncol(sl_data), "). Global ID: ", sl_info$center_global_id))
  }

  if (is.data.frame(sl_data)) {
    sl_data <- as.matrix(sl_data)
  } else if (inherits(sl_data, "Matrix")) {
    sl_data <- as.matrix(sl_data)
  } else if (inherits(sl_data, "neuroim2::ROIVec")) {
    sl_data <- neuroim2::values(sl_data)
  }
  if (is.matrix(sl_data) && !is.numeric(sl_data)) {
    storage.mode(sl_data) <- "double"
  }
  if (!is.matrix(sl_data) || !is.numeric(sl_data)) {
    rlang::abort("`sl_data` must be coercible to a numeric matrix (samples x features).")
  }

  if (missing(cv_spec) || is.null(cv_spec)) {
    cv_spec <- if (!is.null(obj$cv_spec)) obj$cv_spec else obj$crossval
  }
  if (is.null(cv_spec)) {
    rlang::abort("`cv_spec` must be supplied either as an argument or stored on the model specification.")
  }
  
  mvpa_des <- obj$design$mvpa_design
  C <- obj$design$contrast_matrix
  Q <- ncol(C)
  contrast_names <- colnames(C)
  if (is.null(contrast_names)) contrast_names <- paste0("Contrast", 1:Q)

  # Initialize results list
  results_list <- list()
  # Centralized NA reason tracking
  current_na_reason <- NULL
  set_na_reason <- function(reason) {
    if (is.null(current_na_reason)) current_na_reason <<- reason
  }

  # --- Step 1: Compute Cross-Validated Means (Û_sl) or Distances --- 
  # And also get U_hat for Delta calculation if using crossnobis for distances
  U_hat_for_delta_calc <- NULL
  dvec_sl              <- NULL

  if (obj$estimation_method == "crossnobis") {
    cv_outputs <- compute_crossvalidated_means_sl(
      sl_data,
      mvpa_des,
      cv_spec,
      estimation_method = "crossnobis", # Pass explicitly to ensure correct handling in cv_means
      whitening_matrix_W = obj$whitening_matrix_W,
      return_folds = TRUE
    )
    U_hat_for_delta_calc <- cv_outputs$mean_estimate # This is the overall mean estimate
    U_folds_data         <- cv_outputs$fold_estimates  # This is K x V x M for distance calc
    
    P_voxels <- ncol(sl_data) # Or ncol(U_hat_for_delta_calc)
    # Ensure crossnobis_helpers.R is sourced or function is available
    dvec_sl <- compute_crossnobis_distances_sl(U_folds_data, P_voxels)

  } else { # "average" or "L2_norm"
    cv_outputs <- compute_crossvalidated_means_sl(
      sl_data,
      mvpa_des,
      cv_spec,
      obj$estimation_method, # "average" or "L2_norm"
      whitening_matrix_W = NULL, # W is not used for these methods in compute_cv_means
      return_folds = obj$calc_reliability
    )
    if (is.list(cv_outputs) && obj$calc_reliability) {
      U_hat_for_delta_calc <- cv_outputs$mean_estimate
      U_folds_data <- cv_outputs$fold_estimates
    } else {
      U_hat_for_delta_calc <- cv_outputs
      U_folds_data <- NULL
    }
    # --- Step 2: Compute Empirical Second Moment Matrix (Ĝ_sl) ---
    # This step is only for "average" and "L2_norm" which produce U_hat directly for G_hat
    G_hat_sl <- U_hat_for_delta_calc %*% t(U_hat_for_delta_calc)
    dvec_sl <- G_hat_sl[lower.tri(G_hat_sl)]
  }
  
  # Assign K from the U_hat used for delta, which should always be valid
  K <- nrow(U_hat_for_delta_calc)
  
  # Check for NAs early in U_hat_for_delta_calc (used for Delta projections)
  if (anyNA(U_hat_for_delta_calc)) {
      warning("NA values present in cross-validated means (U_hat_for_delta_calc) used for Delta projections. Cannot proceed. Returning NAs for all metrics.")
      set_na_reason("NA in U_hat_for_delta_calc")
      # Populate all requested metrics with NAs
      for (metric_name in obj$output_metric) {
          num_outputs_metric <- if (metric_name %in% c("recon_score", "composite")) 1 else Q
          output_names_metric <- if (num_outputs_metric == 1) metric_name else contrast_names
          results_list[[metric_name]] <- setNames(rep(NA_real_, num_outputs_metric), output_names_metric)
      }
      attr(results_list, "na_reason") <- current_na_reason
      return(results_list)
  }
  # For dvec_sl, NA check for RSA regression happens later (complete.cases)

  # --- Step 3: Prepare RSA Regression --- (was Step 2 for G_hat_sl)
  include_vec <- NULL
  # Use condition_block_list from msreve_design if available
  if (!is.null(obj$design$condition_block_list) && isFALSE(mvpa_des$keep_intra_run)) {
      condition_blocks <- obj$design$condition_block_list
      K_conditions_in_list <- length(condition_blocks)
      
      if (K_conditions_in_list == K) { # K is from nrow(U_hat_for_delta_calc)
          # Ensure names in condition_block_list match the order of U_hat_for_delta_calc rows (conditions)
          # U_hat_for_delta_calc rows are named by unique_conditions from compute_cv_means_sl
          ordered_condition_names_uhat <- rownames(U_hat_for_delta_calc)
          if (all(names(condition_blocks) %in% ordered_condition_names_uhat) && 
              all(ordered_condition_names_uhat %in% names(condition_blocks))) {
              
              condition_blocks_ordered <- condition_blocks[ordered_condition_names_uhat]

              is_intra_block_pair <- matrix(FALSE, nrow = K, ncol = K)
              for (r in 1:(K-1)) {
                  for (c_col in (r+1):K) {
                      if (length(intersect(condition_blocks_ordered[[r]], condition_blocks_ordered[[c_col]])) > 0) {
                          is_intra_block_pair[r, c_col] <- TRUE
                          is_intra_block_pair[c_col, r] <- TRUE # Symmetric
                      }
                  }
              }
              include_vec <- !is_intra_block_pair[lower.tri(is_intra_block_pair)]
          } else {
              warning("Names/order in condition_block_list do not match conditions in U_hat_for_delta_calc. Cannot reliably exclude intra-block RDM cells.")
          }
      } else {
          warning(paste0("Length of condition_block_list (", K_conditions_in_list, 
                         ") does not match number of conditions in U_hat_for_delta_calc (", K, "). Cannot exclude intra-block RDM cells."))
      }
  } else if (isFALSE(mvpa_des$keep_intra_run) && is.null(obj$design$condition_block_list)){
      # Original warning if keep_intra_run is FALSE but we couldn't make the list
      warning("keep_intra_run is FALSE, but condition-to-block mapping (condition_block_list) was not available in the design. Intra-block RDM cells cannot be excluded.")
  }
  
  # Apply include_vec filtering if it was generated
  if (!is.null(include_vec)) {
      # Ensure dvec_sl has the expected number of elements before filtering
      # Expected length is K*(K-1)/2 if K is num conditions from U_hat_for_delta_calc
      expected_len_dvec <- K*(K-1)/2
      if (length(dvec_sl) != expected_len_dvec) {
          rlang::abort(paste0("Internal error: dvec_sl length (", length(dvec_sl), 
                               ") does not match expected K(K-1)/2 = ", expected_len_dvec, 
                               " before include_vec filtering."))
      }
      if (length(include_vec) != expected_len_dvec) {
          rlang::abort(paste0("Internal error: include_vec length (", length(include_vec), 
                               ") does not match expected K(K-1)/2 = ", expected_len_dvec, "."))
      }
      dvec_sl <- dvec_sl[include_vec]
  }

  # Create predictor matrix X_sl from contrast RDMs
  X_sl_list <- lapply(1:Q, function(q) {
      contrast_rdm <- C[, q, drop = FALSE] %*% t(C[, q, drop = FALSE])
      pred_vec <- contrast_rdm[lower.tri(contrast_rdm)]
      if (!is.null(include_vec)) {
          pred_vec <- pred_vec[include_vec]
      }
      pred_vec
  })
  names(X_sl_list) <- contrast_names
  
  # Add nuisance RDMs if present
  nuis_list <- list()
  if (!is.null(obj$design$nuisance_rdms)) {
    for (nm in names(obj$design$nuisance_rdms)) {
      M <- obj$design$nuisance_rdms[[nm]]
      # Extract lower triangle
      v <- if (inherits(M, "dist")) {
        as.vector(M)
      } else {
        M[lower.tri(M)]
      }
      # Apply include_vec mask if present
      if (!is.null(include_vec)) {
        v <- v[include_vec]
      }
      nuis_list[[paste0("nuis_", nm)]] <- v
    }
  }
  
  # Combine contrast and nuisance predictors
  X_sl_list <- c(X_sl_list, nuis_list)

  # Combine predictors into a matrix
  Xmat <- tryCatch({
    do.call(cbind, X_sl_list)
  }, error = function(e) {
    warning(paste("Could not bind contrast predictors into matrix:", e$message))
    NULL # Signal failure
  })
  
  if (is.null(Xmat)) {
      results_list <- list()
      for (metric_name in obj$output_metric) {
          results_list[[metric_name]] <- setNames(rep(NA_real_, Q), contrast_names)
      }
      attr(results_list, "na_reason") <- "Failed to create predictor matrix"
      return(results_list)
  }

  # Handle potential NAs in dvec_sl or predictors before regression
  # cbind will recycle dvec_sl if Xmat has more rows and dvec_sl is shorter - this should not happen
  # if include_vec was applied consistently to dvec_sl and to the construction of Xmat rows.
  if (!is.null(Xmat) && length(dvec_sl) != nrow(Xmat)) {
      rlang::abort(paste0("Dimension mismatch before complete.cases: length(dvec_sl) = ", length(dvec_sl), 
                         ", nrow(Xmat) = ", nrow(Xmat), ". This can happen if include_vec was not applied consistently."))
  }
  
  valid_idx <- complete.cases(cbind(dvec_sl, Xmat))
  n_valid <- sum(valid_idx)

  # Abort early if too few valid data points for regression
  min_samples_needed <- Q + 2 # Minimal number for lm
  if (n_valid < min_samples_needed) {
      warning(paste0("Insufficient valid data points (", n_valid, ") for RSA regression after handling NAs/exclusions. Need at least ", min_samples_needed, ". Returning NAs for all metrics."))
      set_na_reason("Insufficient valid RSA pairs")
      for (metric_name in obj$output_metric) {
          num_outputs_metric <- if (metric_name %in% c("recon_score", "composite")) 1 else Q
          output_names_metric <- if (num_outputs_metric == 1) metric_name else contrast_names
          results_list[[metric_name]] <- setNames(rep(NA_real_, num_outputs_metric), output_names_metric)
      }
      attr(results_list, "na_reason") <- current_na_reason
      return(results_list)
  }
  # Suggestion from audit: warn if low statistical power
  if (n_valid < 10 * Q) {
      warning(paste0("Low number of valid data points (", n_valid,") relative to number of contrasts (", Q, "). Regression results may be unstable."))
  }
  
  dvec_sl_valid <- dvec_sl[valid_idx]
  Xmat_valid <- Xmat[valid_idx, , drop = FALSE]
  colnames(Xmat_valid) <- contrast_names # Ensure names are preserved
  
  # Additional check for dimensional consistency after all filtering and NA removal
  if (nrow(Xmat_valid) != length(dvec_sl_valid)) {
      rlang::abort(paste0("Internal dimension mismatch after NA removal: nrow(Xmat_valid) = ", 
                         nrow(Xmat_valid), ", length(dvec_sl_valid) = ", length(dvec_sl_valid), "."))
  }
  
  # --- Optional collinearity check on the actual regression design ---
  if (isTRUE(obj$check_collinearity)) {
      qr_res <- qr(Xmat_valid)
      if (qr_res$rank < ncol(Xmat_valid)) {
          rlang::abort(paste0("Design matrix rank deficiency inside train_model: rank = ",
                              qr_res$rank, " < ", ncol(Xmat_valid),
                              ". Consider orthogonalizing contrasts or remove dependent ones."))
      }
  }

  # --- Step 4: Run Regression to get β_sl ---
  beta_sl <- NULL
  beta_all <- NULL  # Store all coefficients including nuisance
  n_obs <- nrow(Xmat_valid)
  p_preds <- ncol(Xmat_valid)
  Q_contrasts <- Q # Q from outer scope, number of original contrasts
  all_pred_names <- colnames(Xmat_valid)  # Includes both contrast and nuisance names
  contrast_idx <- seq_len(Q_contrasts)  # Indices for contrast predictors

  # Define SVD-based HKB ridge regression function (from audit)
  ridge_hkb_svd <- function(X, y, p_eff) {
    # X: n_obs x p_preds matrix
    # y: n_obs vector
    # p_eff: effective number of predictors for HKB lambda (usually p_preds)
    n_eff <- nrow(X)
    if (n_eff <= p_eff) {
        warning(paste0("HKB (SVD): n_obs (", n_eff, ") <= p_preds (", p_eff, "). May produce unstable lambda. Consider fixed lambda or more data."))
        # Fallback to a tiny lambda or return NA? For now, proceed but result might be mostly OLS-like if d is small.
    }
    
    # Handle case where X has zero columns after filtering (e.g. all NA predictors)
    if (ncol(X) == 0) {
        warning("HKB (SVD): Predictor matrix X has zero columns. Returning NAs.")
        return(setNames(rep(NA_real_, p_eff), colnames(X))) # or use original contrast_names if available
    }

    sv <- tryCatch(svd(X, nu = 0, nv = ncol(X)), # memory-efficient SVD (no U returned)
                   error = function(e) {
                     warning(paste0("SVD failed in ridge_hkb_svd: ", e$message))
                     return(NULL)
                   }
                   )
    
    if (is.null(sv)) {
        return(setNames(rep(NA_real_, p_eff), colnames(X)))
    }

    d_singular_values <- sv$d
    # Filter out very small singular values to avoid issues with d_singular_values^-1
    # and to define rank for sigma2 df calculation
    rank_X <- sum(d_singular_values > (max(d_singular_values) * .Machine$double.eps * 10)) # Practical rank
    if (rank_X == 0) {
        warning("HKB (SVD): Matrix X appears to be rank zero. Returning NAs.")
        return(setNames(rep(NA_real_, p_eff), colnames(X)))
    }

    d_inv <- ifelse(d_singular_values > (max(d_singular_values) * .Machine$double.eps * 10), 1/d_singular_values, 0)
    
    # Compute U^T y without explicit U:  uy = D^{-1} V^T (X^T y)
    Xt_y <- crossprod(X, y)                          # p x 1
    vt_Xty <- crossprod(sv$v, Xt_y)                 # p x 1 (V^T X^T y)
    uy_full <- (1 / d_singular_values) * vt_Xty     # D^{-1} V^T X^T y
    uy <- uy_full[1:rank_X, , drop=FALSE]           # first rank components
    beta_ols_for_sigma2_v_basis <- uy * d_inv[1:rank_X]
    # Convert to full beta_ols vector for residual calculation
    beta_ols_vec <- sv$v[, 1:rank_X, drop=FALSE] %*% beta_ols_for_sigma2_v_basis
    
    df_residual_sigma2 <- n_eff - rank_X # Use rank for df in sigma2
    if (df_residual_sigma2 <= 0) {
        warning(paste0("HKB (SVD): df for sigma2 (n_obs - rank(X) = ", df_residual_sigma2, ") is not positive. Using small default lambda."))
        lambda <- 1e-4 # Small default lambda as fallback
    } else {
        sigma2 <- sum((y - X %*% beta_ols_vec)^2) / df_residual_sigma2
        sum_beta_ols_scaled_sq <- sum((uy * d_inv[1:rank_X])^2) # beta_ols^T beta_ols in V-basis
        if (sum_beta_ols_scaled_sq < 1e-10) {
             warning("HKB (SVD): Sum of (scaled) OLS singular value components is near zero. Using small default lambda.")
             lambda <- 1e-4
        } else {
             lambda <- (p_eff * sigma2) / sum_beta_ols_scaled_sq # p_eff could be rank_X or original p_preds
             if (lambda < 0) {
                  warning("HKB (SVD): Negative lambda calculated, using absolute value. Check data quality.")
                  lambda <- abs(lambda)
             }
        }
    }
    
    # Ridge coefficients: beta = V * diag(d / (d^2 + lambda)) * (U^T y)
    d_eigen_sq <- d_singular_values^2
    coeff_v_basis <- (d_singular_values / (d_eigen_sq + lambda)) * uy_full  # uy_full = D^{-1} V^T X^T y
    beta_ridge <- sv$v %*% coeff_v_basis
    
    return(drop(beta_ridge))
  }
  
  if (obj$regression_type == "ridge_hkb") {
      beta_all_values <- ridge_hkb_svd(Xmat_valid, dvec_sl_valid, p_preds)
      if (length(beta_all_values) == p_preds) {
        # Extract only contrast coefficients
        beta_sl <- setNames(beta_all_values[contrast_idx], contrast_names)
      } else {
        warning("HKB (SVD) ridge regression did not return the expected number of coefficients. Returning NAs.")
        beta_sl <- setNames(rep(NA_real_, Q_contrasts), contrast_names)
        attr(beta_sl, "na_reason") <- "Ridge HKB (SVD) coefficient mismatch"
      }
  } else if (obj$regression_type == "lm") {
    if (n_obs <= p_preds) {
        warning(paste0("Number of observations (", n_obs, ") is not greater than number of predictors (", 
                       p_preds, ") for OLS. Regression may fail or be unstable. Returning NAs."))
        beta_sl <- setNames(rep(NA_real_, Q_contrasts), contrast_names)
        attr(beta_sl, "na_reason") <- "OLS n <= p"
    } else {
        xtx <- crossprod(Xmat_valid)
        xty <- crossprod(Xmat_valid, dvec_sl_valid)
        beta_ols_solve <- tryCatch({
            solve(xtx, xty)
        }, error = function(e) {
            warning(paste0("OLS regression failed (e.g. singular XTX): ", e$message, ". Returning NAs."))
            NULL
        })
        if (is.null(beta_ols_solve)) {
            beta_sl <- setNames(rep(NA_real_, Q_contrasts), contrast_names)
            attr(beta_sl, "na_reason") <- "OLS failed (e.g. singular XTX)"
        } else {
            # Extract only contrast coefficients
            beta_sl <- setNames(as.vector(beta_ols_solve)[contrast_idx], contrast_names)
        }
    }
  } else if (obj$regression_type %in% c("pearson", "spearman")) {
      # Centering for correlation types already handled for Xmat_valid by scale()
      # dvec_sl_valid needs to be centered for pearson, ranked for spearman
      dvec_sl_for_cor <- dvec_sl_valid
      Xmat_for_cor <- Xmat_valid
      if (obj$regression_type == "pearson") {
          dvec_sl_for_cor <- dvec_sl_valid - mean(dvec_sl_valid)
      } else if (obj$regression_type == "spearman") {
          # Rank–transform both outcome and predictors (no prior centering)
          dvec_sl_for_cor <- rank(dvec_sl_valid, ties.method = "average")
          Xmat_for_cor    <- apply(Xmat_valid, 2, rank, ties.method = "average")
      }
      
      temp_obj_for_helpers <- list(
          design = list(model_mat = as.list(data.frame(Xmat_for_cor))),
          distmethod = obj$regression_type 
      )
      beta_all <- tryCatch({
          run_cor(dvec_sl_for_cor, temp_obj_for_helpers)
      }, error = function(e) {
          warning(paste("Correlation-based RSA failed (", obj$regression_type, "):", e$message, ". Returning NAs."))
          setNames(rep(NA_real_, p_preds), all_pred_names)
      })
      # Extract only contrast coefficients
      if (!is.null(beta_all) && length(beta_all) >= Q_contrasts) {
          beta_sl <- setNames(beta_all[contrast_idx], contrast_names)
      } else {
          beta_sl <- setNames(rep(NA_real_, Q_contrasts), contrast_names)
      }
  } else if (obj$regression_type == "rfit") {
      temp_obj_for_helpers <- list(
          design = list(model_mat = Xmat_valid),
          distmethod = obj$regression_type
      )
      beta_all <- tryCatch({
          run_rfit(dvec_sl_valid, temp_obj_for_helpers)
      }, error = function(e) {
          warning(paste("Rfit RSA failed:", e$message, ". Returning NAs."))
          setNames(rep(NA_real_, p_preds), all_pred_names)
      })
      # Extract only contrast coefficients
      if (!is.null(beta_all) && length(beta_all) >= Q_contrasts) {
          beta_sl <- setNames(beta_all[contrast_idx], contrast_names)
      } else {
          beta_sl <- setNames(rep(NA_real_, Q_contrasts), contrast_names)
      }
  } else {
      rlang::abort(paste("Unsupported regression_type in train_model:", obj$regression_type))
  }

  # If regression failed, beta_sl might already be NA vector or NULL
  if (is.null(beta_sl) || (is.numeric(beta_sl) && length(beta_sl) != Q_contrasts)) {
      if (is.null(attr(beta_sl, "na_reason"))) {
          # Only set generic reason if a specific one wasn't set by the failing block
          # attr(beta_sl, "na_reason") <- "Regression resulted in NULL or incorrect length output" # This was for beta_sl itself
          set_na_reason("Regression resulted in NULL or incorrect length output for beta_sl")
      } else {
          set_na_reason(attr(beta_sl, "na_reason")) # Propagate specific reason from beta_sl
      }
      beta_sl <- setNames(rep(NA_real_, Q_contrasts), contrast_names) # Ensure it's a NA vector of correct form
  } else if (any(is.na(beta_sl)) && is.null(attr(beta_sl, "na_reason"))) {
      # attr(beta_sl, "na_reason") <- "Regression produced NA values"
      set_na_reason("Regression for beta_sl produced NA values")
  }
  
  # --- Step 5: Compute Voxel Projections (Δ_sl) ---
  # Ensure row order consistency between U_hat_for_delta_calc and C
  idx_match <- match(rownames(U_hat_for_delta_calc), rownames(C))
  if (anyNA(idx_match)) {
    rlang::abort(paste0("Row names of contrast matrix C must match those of U_hat_for_delta_calc (from ", obj$estimation_method, " path)."))
  }
  C_ord <- C[idx_match, , drop = FALSE]
  Delta_sl <- t(U_hat_for_delta_calc) %*% C_ord # V_sl x Q matrix

  # --- Step 6: Extract Center Voxel Projection (Δ_{v_c,sl}) ---
  delta_vc_sl <- Delta_sl[center_idx, , drop = TRUE] # Q-dimensional vector
  # Handle potential NAs from U_hat_for_delta_calc (though checked earlier) or C matrix
  if (anyNA(delta_vc_sl) && is.null(current_na_reason)) { # Only set if no prior reason
       set_na_reason("NA in delta_vc_sl projection")
  }

  # --- Step 6b: Optionally normalize delta_vc_sl for relevant metrics ---
  # This normalization is relevant for "beta_delta_norm" and "composite"
  # For "beta_delta", unnormalized delta is used unless obj$normalize_delta is TRUE for that specific metric (see below)
  
  delta_vc_sl_normalized <- delta_vc_sl # Keep original for some metrics
  if (any(c("beta_delta_norm", "composite") %in% obj$output_metric) || 
      (obj$normalize_delta && "beta_delta" %in% obj$output_metric) ) {
    if (anyNA(delta_vc_sl)) {
        warning("Cannot normalize delta_vc_sl due to NA values. Will use unnormalized version if possible or propagate NAs.")
        # NAs will propagate naturally
    } else {
        norm_delta_vc_val <- sqrt(sum(delta_vc_sl^2))
        if (norm_delta_vc_val < 1e-10) {
            warning("delta_vc_sl has near-zero L2 norm; using unnormalized or zero vector for normalized metrics.")
            # If norm is zero, normalized version is undefined or could be zeros.
            # For consistency, if not normalizing, use original. If normalizing, it might become NAs or zeros.
            # Let's ensure it doesn't produce NaNs if we proceed with division by zero norm.
            # If we intend to set to zero:
            # delta_vc_sl_normalized <- rep(0, length(delta_vc_sl))
            # For now, rely on NA propagation if delta_vc_sl was NA, or if beta_sl is NA.
            # If delta_vc_sl is all zeros, then delta_vc_sl_normalized will be all zeros if norm_delta_vc_val is treated as 1.
            # Let's make it so that if norm is too small, normalized delta is also all zeros, avoiding NaN.
            if (norm_delta_vc_val < 1e-10) {
                 delta_vc_sl_normalized <- rep(0, length(delta_vc_sl))
            } else {
                 delta_vc_sl_normalized <- delta_vc_sl / norm_delta_vc_val
            }
        } else {
            delta_vc_sl_normalized <- delta_vc_sl / norm_delta_vc_val
        }
    }
  }

  # --- Step 6c: Calculate Reliability Weights (rho) if requested ---
  rho_vc_sl <- rep(1, Q)
  if (isTRUE(obj$calc_reliability)) {
    S_total <- get_nfolds(cv_spec)
    if (exists("U_folds_data", inherits = FALSE) && !is.null(U_folds_data)) {
      fold_array <- U_folds_data
    } else {
      cv_tmp <- compute_crossvalidated_means_sl(
        sl_data,
        mvpa_des,
        cv_spec,
        obj$estimation_method,
        whitening_matrix_W = if (obj$estimation_method == "crossnobis") obj$whitening_matrix_W else NULL,
        return_folds = TRUE
      )
      fold_array <- cv_tmp$fold_estimates
    }

    if (is.array(fold_array) && length(dim(fold_array)) == 3) {
      S_eff <- dim(fold_array)[3]
      mean_delta <- rep(0, Q)
      M2_delta <- rep(0, Q)
      valid_folds <- 0
      for (s_idx in seq_len(S_eff)) {
        U_fold <- fold_array[,,s_idx]
        if (!is.matrix(U_fold)) next
        idx_match_fold <- match(rownames(U_fold), rownames(C))
        if (anyNA(idx_match_fold)) next
        C_fold <- C[idx_match_fold, , drop = FALSE]
        Delta_fold_sl <- t(U_fold) %*% C_fold
        delta_fold_center <- Delta_fold_sl[center_idx, , drop = TRUE]
        if (anyNA(delta_fold_center)) next
        valid_folds <- valid_folds + 1
        delta_diff <- delta_fold_center - mean_delta
        mean_delta <- mean_delta + delta_diff/valid_folds
        M2_delta <- M2_delta + delta_diff*(delta_fold_center - mean_delta)
      }

      if (valid_folds > 1) {
        var_delta <- M2_delta/(valid_folds - 1)
        sigma2_noise_param <- (valid_folds - 1) * var_delta
        denom <- var_delta + sigma2_noise_param
        rho_vc_sl <- ifelse(denom < 1e-10, 1, sigma2_noise_param/denom)
        rho_vc_sl[is.na(rho_vc_sl)] <- 0
      } else if (valid_folds == 1) {
        rho_vc_sl[M2_delta == 0] <- 1
        rho_vc_sl[M2_delta != 0] <- 0
      } else {
        rho_vc_sl <- rep(0, Q)
      }
    }
  }

  # --- Step 7: Calculate Final Metrics ---
  
  # Pre-check: if beta_sl or delta_vc_sl has NA due to critical failure,
  # then metrics depending on them will be NA.
  # The current_na_reason should capture this.

  for (metric_name in obj$output_metric) {
    # Skip calculation if a critical NA reason is already set and the metric depends on beta/delta
    if (!is.null(current_na_reason) && metric_name != "recon_score") { # recon_score has its own NA logic mostly
        num_outputs_metric <- if (metric_name %in% c("recon_score", "composite")) 1 else Q
        output_names_metric <- if (num_outputs_metric == 1) metric_name else contrast_names
        results_list[[metric_name]] <- setNames(rep(NA_real_, num_outputs_metric), output_names_metric)
        next
    }
    
    metric_value <- NULL
    
    if (metric_name == "recon_score") {
        r_v <- NA_real_
        na_reason_rv <- NULL
        if (anyNA(beta_sl)) { # recon_score needs valid beta_sl
            warning("Cannot calculate recon_score due to NA in beta_sl. Returning NA.")
            na_reason_rv <- attr(beta_sl, "na_reason") %||% "NA in beta_sl for recon_score"
        } else if (anyNA(delta_vc_sl)) { # and valid delta_vc_sl (unnormalized version)
             warning("Cannot calculate recon_score due to NA in delta_vc_sl (unnormalized). Returning NA.")
             na_reason_rv <- "NA in delta_vc_sl for recon_score"
        } else {
            # Build G_hat_v = C_ord %*% diag(beta_sl * delta_vc_sl) %*% t(C_ord)
            # Use unnormalized delta_vc_sl for recon_score's G_hat_v
            diag_vals <- beta_sl * delta_vc_sl 
            if (length(diag_vals) == 1 && Q_contrasts == 1) {
                loading_diag <- matrix(diag_vals, 1, 1)
            } else {
                loading_diag <- diag(diag_vals)
            }
            G_hat_v <- C_ord %*% loading_diag %*% t(C_ord)
            vec_G_hat_v <- G_hat_v[lower.tri(G_hat_v)]
            
            G_hat_sl_for_recon <- U_hat_for_delta_calc %*% t(U_hat_for_delta_calc)
            vec_G_hat_sl_empirical <- G_hat_sl_for_recon[lower.tri(G_hat_sl_for_recon)]
            
            if (!is.null(include_vec)) {
                if (length(vec_G_hat_v) == length(include_vec) && length(vec_G_hat_sl_empirical) == length(include_vec)) {
                    vec_G_hat_v <- vec_G_hat_v[include_vec]
                    vec_G_hat_sl_empirical <- vec_G_hat_sl_empirical[include_vec]
                } else {
                    rlang::abort("Internal error: Length mismatch for include_vec filtering in recon_score.")
                }
            }
            
            if (stats::var(vec_G_hat_v, na.rm=TRUE) < 1e-10 || stats::var(vec_G_hat_sl_empirical, na.rm=TRUE) < 1e-10) {
               warning("Cannot calculate recon_score: Zero variance in one or both RDM vectors.")
               na_reason_rv <- "Zero variance in RDM vectors for recon_score"
            } else {
                r_v <- stats::cor(vec_G_hat_v, vec_G_hat_sl_empirical, use = "pairwise.complete.obs")
                if (is.na(r_v)) {
                    na_reason_rv <- "Correlation for recon_score resulted in NA"
                }
            }
        }
        metric_value <- setNames(r_v, "recon_score")
        if (!is.null(na_reason_rv) && is.null(current_na_reason)) set_na_reason(na_reason_rv)

    } else if (metric_name == "beta_delta") {
        # Use unnormalized delta_vc_sl unless obj$normalize_delta is globally true
        current_delta <- if(obj$normalize_delta) delta_vc_sl_normalized else delta_vc_sl
    metric_value <- beta_sl * current_delta
    names(metric_value) <- contrast_names

  } else if (metric_name == "beta_delta_norm") {
        if (!obj$normalize_delta) {
            warning("Output metric is 'beta_delta_norm' but 'normalize_delta' is FALSE. Result for 'beta_delta_norm' will be equivalent to 'beta_delta' using unnormalized delta.")
        }
    metric_value <- beta_sl * delta_vc_sl_normalized # Always use normalized delta here
    names(metric_value) <- contrast_names

  } else if (metric_name == "beta_delta_reliable") {
        current_delta <- if(obj$normalize_delta) delta_vc_sl_normalized else delta_vc_sl
        metric_value <- beta_sl * current_delta * rho_vc_sl
        names(metric_value) <- contrast_names

    } else if (metric_name == "beta_only") {
        metric_value <- beta_sl
        names(metric_value) <- contrast_names

    } else if (metric_name == "delta_only") {
        # This should refer to the raw, unnormalized delta unless normalize_delta is globally true
        current_delta <- if(obj$normalize_delta) delta_vc_sl_normalized else delta_vc_sl
        metric_value <- current_delta
        names(metric_value) <- contrast_names
        
    } else if (metric_name == "composite") {
        is_orthonormal <- attr(obj$design, "is_orthonormal")
        if (is.null(is_orthonormal)) {
            warning("Could not determine if contrast matrix is orthonormal (attribute missing). Assuming not orthonormal for composite metric interpretation.")
            is_orthonormal <- FALSE 
        }
        if (!is_orthonormal) {
            warning("Composite metric (net-pull) calculated, but contrast matrix is not orthonormal.")
            if (!isTRUE(obj$allow_nonorth_composite)) {
                metric_value <- setNames(NA_real_, "composite")
                if (is.null(current_na_reason)) set_na_reason("Non-orthonormal contrast matrix for composite")
                results_list[[metric_name]] <- metric_value
                next # Skip to next metric
            }
        }
        
        # Composite always uses L2-normalized delta contributions, regardless of global obj$normalize_delta
        weighted_contributions <- beta_sl * delta_vc_sl_normalized 
        composite_score <- sum(weighted_contributions, na.rm = TRUE) # sum over Q contrasts
        
        if (all(is.na(weighted_contributions))) {
            composite_score <- NA_real_
            if (is.null(current_na_reason)) set_na_reason("All beta*norm_delta components were NA for composite")
        }
        metric_value <- setNames(composite_score, "composite")
    }
    
    results_list[[metric_name]] <- metric_value
  } # End loop over obj$output_metric

  # --- Step 8: Return Result ---
  if (!is.null(current_na_reason)) {
    attr(results_list, "na_reason") <- current_na_reason
    # Ensure all metrics in results_list are NA if a global NA reason was set mid-way
    # and some metrics were computed before the NA reason was established
    # This logic might be complex if some metrics could be valid while others are not.
    # For now, if current_na_reason is set, it implies a fairly global issue.
    # The individual metric assignment handles NA propagation.
  }
  
  # Add S3 class to the results list
  class(results_list) <- c("contrast_rsa_model_results", "list")
  
  results_list
}



#' @export
merge_results.contrast_rsa_model <- function(obj, result_set, indices, id, ...) {
  if (any(result_set$error)) {
    emessage <- result_set$error_message[which(result_set$error)[1]]
    return(tibble::tibble(
      result = list(NULL),
      indices = list(indices),
      performance = list(NULL),
      id = id,
      error = TRUE,
      error_message = emessage
    ))
  }

  metrics <- result_set$result[[1]]
  if (is.null(metrics)) {
    warning(sprintf('merge_results.contrast_rsa_model: Received NULL metrics for ROI %s.', id))
    return(tibble::tibble(
      result = list(NULL),
      indices = list(indices),
      performance = list(NULL),
      id = id,
      error = TRUE,
      error_message = 'train_model returned NULL metrics'
    ))
  }

  tibble::tibble(
    result = list(metrics),
    indices = list(indices),
    performance = list(metrics),
    id = id,
    error = FALSE,
    error_message = "~"
  )
}

#' Run Searchlight Analysis for Contrast RSA Model
#'
#' This is the S3 method for running a searchlight analysis specifically for a
#' \code{contrast_rsa_model} object. It performs the Multi-Dimensional Signed
#' Representational Voxel Encoding (MS-ReVE) style analysis across the brain
#' volume or surface.
#'
#' @param model_spec A \code{contrast_rsa_model} object created by
#'   \code{\\link{contrast_rsa_model}}.
#' @param radius The radius of the searchlight sphere (in mm for volume data, or
#'   geodesic distance/vertex count for surface data - interpretation depends
#'   on the \code{distance_metric} used in \code{neuroim2::searchlight} or
#'   \code{neurosurf::SurfaceSearchlight}).
#' @param method The type of searchlight procedure. Currently, only \code{"standard"}
#'   is fully supported and recommended for contrast RSA. The "standard" method
#'   iterates through all voxels/vertices in the mask, treating each as the center
#'   and calculating its specific contribution metric (e.g., beta_delta). The results
#'   are combined directly into Q output maps using \code{combine_msreve_standard}.
#'   \cr
#'   Using \code{"randomized"} is **not recommended** for this model. While the code
#'   will run, the default combiner (\code{combine_msreve_standard}) will produce sparse
#'   maps with values only at the random center locations. A proper interpretation of
#'   randomized searchlight (averaging results based on voxel coverage) would require a
#'   different, model-specific combiner (e.g., \code{combine_msreve_randomized}), which
#'   is not currently implemented due to conceptual challenges in averaging the center-voxel
#'   specific MS-ReVE metric.
#' @param niter The number of iterations if \code{method = "randomized"} is used (but see note above).
#' @param ... Additional arguments passed down to the underlying searchlight
#'   machinery (e.g., \code{mvpa_iterate}).
#'
#' @return A \code{searchlight_result} object (specifically with class
#'   \code{c("msreve_searchlight_result", "searchlight_result", "list")}), containing:
#'   \item{results}{A named list where each element is a \code{SparseNeuroVec} or
#'     \code{NeuroSurfaceVector} representing the map for one contrast.}
#'   \item{metrics}{A character vector of the contrast names.}
#'   \item{n_voxels}{Total number of voxels/vertices in the original mask space.}
#'   \item{active_voxels}{Number of voxels/vertices for which results were computed.}
#'
#' @examples
#' # Assuming 'spec' is a valid contrast_rsa_model object
#' # Standard (recommended) method:
#' # results <- run_searchlight(spec, radius = 8, method = "standard")
#' # plot(results$results[[1]]) # Plot the map for the first contrast
#'
#' @examples
#' # Assuming contrast_rsa_model examples have run and objects like
#' # 'mvpa_dat', 'msreve_des', 'model_basic', 'model_recon' are available.
#' # This requires the setup from contrast_rsa_model examples.
#' 
#' if (requireNamespace("neuroim2", quietly = TRUE) && 
#'     requireNamespace("rMVPA", quietly = TRUE) &&
#'     exists("model_basic") && inherits(model_basic, "contrast_rsa_model") &&
#'     exists("model_recon") && inherits(model_recon, "contrast_rsa_model")) {
#'
#'   # --- Example 1: Run searchlight with basic model ---
#'   # Use a very small radius for quick example run.
#'   # Actual searchlight analyses would use a more appropriate radius (e.g., 3-4 voxels).
#'   # With dummy data, results won't be meaningful; focus is on execution.
#'
#'   message("Running searchlight example 1 (basic model, radius=1)... May take a moment.")
#'   sl_results_basic <- tryCatch({
#'     run_searchlight(model_basic, radius = 1, method = "standard")
#'   }, error = function(e) {
#'     message("Searchlight (basic model) example failed: ", e$message)
#'     NULL
#'   })
#'   if (!is.null(sl_results_basic)) {
#'     print(sl_results_basic)
#'   }
#'
#'   # --- Example 2: Run searchlight with recon_score output ---
#'   message("Running searchlight example 2 (recon_score model, radius=1)... May take a moment.")
#'   sl_results_recon <- tryCatch({
#'     run_searchlight(model_recon, radius = 1, method = "standard")
#'   }, error = function(e) {
#'     message("Searchlight (recon_score model) example failed: ", e$message)
#'     NULL
#'   })
#'   if (!is.null(sl_results_recon)) {
#'     print(sl_results_recon)
#'   }
#'
#'   # Note on Crossnobis with searchlight:
#'   # To run a searchlight with 'estimation_method = "crossnobis"' from 'model_crossnobis',
#'   # the 'whitening_matrix_W' needs to be passed through the searchlight machinery
#'   # to 'compute_crossvalidated_means_sl'. This typically involves passing it via
#'   # the `...` argument of `run_searchlight` and ensuring `mvpa_iterate` and
#'   # `train_model` propagate it. This advanced usage is not shown here as it
#'   # requires modification to the general `mvpa_iterate` or a custom processor.
#'
#' } else {
#'  message("Skipping run_searchlight.contrast_rsa_model example execution here.")
#'  message("It can be time-consuming and depends on prior setup.")
#' }
#'
#' @export
#' @importFrom futile.logger flog.info
run_searchlight.contrast_rsa_model <- function(model_spec,
                                                 radius = NULL,
                                                 method = c("standard", "randomized"),
                                                 niter = NULL, # niter only relevant for randomized
                                                 ...) {
  method <- match.arg(method)

  # Currently, only standard method is fully recommended with the specific combiner
  if (method == "randomized") {
      # Warning strengthened based on discussion
      warning("Using 'randomized' method with contrast_rsa_model. This is NOT the recommended approach. The default combiner will produce sparse maps with values only at random centers, not averaged maps reflecting voxel coverage. Interpretation requires caution. See documentation for details.")
      if (is.null(niter) || niter < 1) {
          stop("'niter' must be provided and >= 1 for the 'randomized' method.")
      }
  }

  # The core processing function is train_model.contrast_rsa_model
  # The combiner function MUST be combine_msreve_standard for this model type
  the_combiner <- combine_msreve_standard 

  # Call the appropriate backend searchlight function
  if (method == "standard") {
    futile.logger::flog.info("Running standard MS-ReVE/Contrast RSA searchlight (radius = %s)", radius)
    # Pass the specific combiner for contrast RSA
    rMVPA:::do_standard(model_spec, radius, combiner = the_combiner, ...)
  } else { # method == "randomized"
    futile.logger::flog.info("Running randomized MS-ReVE/Contrast RSA searchlight (radius = %s, niter = %s)", radius, niter)
    # Pass the specific combiner for contrast RSA - Note: Applicability might need review
    rMVPA:::do_randomized(model_spec, radius, niter = niter, combiner = the_combiner, ...)
  }
} 