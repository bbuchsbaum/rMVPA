% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/contrast_rsa_model.R
\name{contrast_rsa_model}
\alias{contrast_rsa_model}
\title{Constructor for contrast_rsa_model}
\usage{
contrast_rsa_model(
  dataset,
  design,
  estimation_method = "average",
  regression_type = "lm",
  output_metric = c("beta_delta"),
  check_collinearity = FALSE,
  normalize_delta = FALSE,
  allow_nonorth_composite = FALSE,
  calc_reliability = FALSE,
  whitening_matrix_W = NULL,
  ...
)
}
\arguments{
\item{dataset}{An object of class \code{mvpa_dataset}, containing the neuroimaging
data and associated metadata.}

\item{design}{An object of class \code{msreve_design}, containing the underlying
\code{mvpa_design} and the contrast matrix. Created by \code{\\link{msreve_design}}.}

\item{estimation_method}{Character string specifying the method to estimate
cross-validated condition means (\code{Û}) or distances. Supported:
\itemize{
  \item \code{"average"}: Simple mean of training samples per condition (for \code{Û}).
  \item \code{"L2_norm"}: Like \code{"average"}, but \code{Û} rows are L2-normalized.
  \item \code{"crossnobis"}: Computes unbiased squared Euclidean distances directly using the Crossnobis method. Results in a distance vector for RSA, not a \code{Û} matrix for \(G_empirical\) construction in the same way. \code{U_hat} for \eqn{\Delta} calculation is still computed using "average" method internally when this is selected. Requires `return_folds=TRUE` from `compute_crossvalidated_means_sl`.
}
Passed to \code{\link{compute_crossvalidated_means_sl}} (for "average", "L2_norm", or to get per-fold means for "crossnobis").}

\item{regression_type}{Character string specifying the method for the RSA regression
(regressing empirical RDM/Second Moment Matrix onto contrast RDMs). Options align
with \code{rsa_model}: \code{"pearson"}, \code{"spearman"}, \code{"lm"}, \code{"rfit"}.
Additional options for ridge regression: \code{"ridge_hkb"} (Hoerl-Kennard-Baldwin lambda).
Default is \code{"lm"}.}

\item{output_metric}{Character vector specifying one or more output metrics to compute.
Multiple metrics can be requested simultaneously and will be returned in a named list.
Duplicates are removed while preserving the order of first occurrence.
Supported options:
\itemize{
  \item \code{"beta_delta"}: The product of the RSA regression coefficient (beta_q)
        and the voxel's projection onto the contrast (delta_q,v). This is the
        primary signed contribution metric.
  \item \code{"beta_only"}: Only the RSA regression coefficient (beta_q).
  \item \code{"delta_only"}: Only the voxel's projection onto the contrast (delta_q,v).
  \item \code{"recon_score"}: Voxel-specific RDM reconstruction score (r_v), correlating the RDM implied by the voxel's loadings with the empirical RDM.
  \item \code{"beta_delta_norm"}: Similar to `beta_delta`, but uses the L2-normalized voxel contribution vector (delta_q,v). Requires `normalize_delta=TRUE` to be meaningful.
  \item \code{"beta_delta_reliable"}: Reliability-weighted contributions, \eqn{\rho_{q,v} \beta_q \Delta_{q,v}}, where \eqn{\rho_{q,v}} reflects fold-wise stability.
  \item \code{"composite"}: Sum of beta-weighted, L2-normalized voxel contributions (Σ_q β_q ~Δ_q,v). Represents the net projection onto the positive diagonal of the contrast space. Interpretation requires caution if contrasts are not orthonormal.
}
Default is \code{c("beta_delta")}.}

\item{check_collinearity}{Logical, whether to check for collinearity among contrast RDMs
when using \code{regression_type = "lm"}. Default is \code{FALSE}.}

\item{normalize_delta}{Logical. If TRUE, the voxel contribution vector (delta_q,v)
is normalized to unit L2 norm before being potentially multiplied by beta_q
for the "beta_delta" output metric. Default is \code{FALSE}.}

\item{allow_nonorth_composite}{Logical. If FALSE (default) the composite metric will return NA when the contrast matrix
is not orthonormal, to avoid mis-interpretation. If TRUE, the composite score is returned regardless, but a warning
is still emitted.}

\item{calc_reliability}{Logical. If TRUE, voxel-wise contribution reliability (ρ
values) are estimated across cross-validation folds during training and can
be incorporated into the output metrics. Default is \code{FALSE}.}

\item{whitening_matrix_W}{Optional V x V numeric matrix (voxels x voxels). Required if
`estimation_method = "crossnobis"` and Mahalanobis (rather than Euclidean)
distances are desired. This matrix (e.g., \eqn{\Sigma_{noise}^{-1/2}}) is passed to
`compute_crossvalidated_means_sl` to whiten per-fold estimates before distance calculation.
Default is `NULL` (Euclidean distances for Crossnobis).}

\item{...}{Additional arguments passed to \code{create_model_spec}.}
}
\value{
An object of class \code{contrast_rsa_model}, \code{mvpa_model_spec},
  \code{model_spec}, and \code{list}.
}
\description{
Creates a contrast_rsa_model specification object, which encapsulates the necessary
parameters and design information for a Multi-Dimensional Signed Representational
Voxel Encoding (MS-ReVE) style analysis.
}
\details{
This model is designed for MS-ReVE style analyses where the goal is to understand
how different predefined contrasts contribute to the representational structure
observed in neural data, particularly at a fine-grained (e.g., voxel) level.
It involves:
1. Estimating cross-validated condition mean patterns (Û) or distances (for Crossnobis).
2. Constructing an empirical second-moment matrix (Ĝ) from Û or using the Crossnobis distances directly.
3. Creating theoretical second-moment matrices (RDMs) from each contrast vector.
4. Regressing the vectorized empirical RDM/distances onto the vectorized contrast RDMs to get β coefficients.
5. Projecting voxel patterns (from Û) onto the contrast space to get Δ (delta) values.
6. Combining β and Δ to form the final output metric (e.g., beta_delta).

**Cross-Validation Compatibility:**
The `estimation_method` relies on `compute_crossvalidated_means_sl` which, in turn, requires
the cross-validation specification (`cv_spec` derived from `mvpa_design$crossval`)
to provide a deterministic, partition-based set of training indices for each fold.
Therefore, cross-validation schemes like `bootstrap_blocked_cross_validation` (which use
resampling with replacement for training folds) are **not suitable** for use with this model,
as they do not align with the assumptions of `compute_crossvalidated_means_sl`.
Schemes like `blocked_cross_validation`, `kfold_cross_validation`, and `custom_cross_validation`
(that define clear partitions) are appropriate.
For `twofold_blocked_cross_validation` and `sequential_blocked_cross_validation`,
their compatibility also depends on whether their `train_indices` methods can deterministically
define training sets for each fold iterated by `compute_crossvalidated_means_sl`.
}
\examples{
# --- Minimal Setup ---
# 1. Create dummy data and an mvpa_dataset

  # Dummy data: 16 samples, 10 voxels, 4 conditions, 2 runs
  set.seed(123)
  n_samples <- 16
  n_voxels <- 10
  n_conditions <- 4 # condA, condB, condC, condD
  n_runs <- 2

  dummy_sl_data <- matrix(rnorm(n_samples * n_voxels), n_samples, n_voxels)
  colnames(dummy_sl_data) <- paste0("V", 1:n_voxels)

  dummy_mask <- neuroim2::NeuroVol(array(1, c(2,2,2)), neuroim2::NeuroSpace(c(2,2,2)))

  condition_labels <- factor(rep(paste0("cond", LETTERS[1:n_conditions]), each = n_samples / n_conditions))
  run_labels <- factor(rep(1:n_runs, each = n_samples / n_runs))

  # Create mvpa_dataset (without Y and block_var)
  mvpa_dat <- rMVPA::mvpa_dataset(
    train_data = dummy_sl_data,
    mask = dummy_mask
  )

  # Create mvpa_design
  mvpa_des <- rMVPA::mvpa_design(
    train_design = data.frame(condition = condition_labels, run = run_labels),
    y_train = ~condition,
    block_var = ~run
  )

  K <- mvpa_des$ncond # Use mvpa_des here

  C_mat <- matrix(0, nrow = K, ncol = 2)
  rownames(C_mat) <- levels(mvpa_des$Y) # Use mvpa_des here
  C_mat["condA", 1] <- 1; C_mat["condB", 1] <- 1
  C_mat["condC", 1] <- -1; C_mat["condD", 1] <- -1
  C_mat["condA", 2] <- 1; C_mat["condB", 2] <- -1
  C_mat <- scale(C_mat, center = TRUE, scale = FALSE)
  colnames(C_mat) <- c("AB_vs_CD", "A_vs_B")

  msreve_des <- rMVPA::msreve_design(
    mvpa_design = mvpa_des, # Use mvpa_des here
    contrast_matrix = C_mat
  )

  # --- Example 1: Basic contrast_rsa_model ---
  model_basic <- contrast_rsa_model(
    dataset = mvpa_dat,
    design = msreve_des
  )
  print(model_basic)

  # --- Example 1b: Requesting multiple metrics ---
  model_multi_metric <- contrast_rsa_model(
    dataset = mvpa_dat,
    design = msreve_des,
    output_metric = c("beta_delta", "recon_score", "beta_only")
  )
  print(model_multi_metric)

  # --- Example 2: Using L2_norm for U_hat and normalize_delta ---
  model_l2_norm_delta <- contrast_rsa_model(
    dataset = mvpa_dat,
    design = msreve_des,
    estimation_method = "L2_norm",
    normalize_delta = TRUE,
    output_metric = "beta_delta_norm"
  )
  print(model_l2_norm_delta)

  # --- Example 3: Ridge Regression (HKB) ---
  model_ridge <- contrast_rsa_model(
    dataset = mvpa_dat,
    design = msreve_des,
    regression_type = "ridge_hkb"
  )
  print(model_ridge)

  # --- Example 4: Reconstruction Score Output ---
  model_recon <- contrast_rsa_model(
    dataset = mvpa_dat,
    design = msreve_des,
    output_metric = "recon_score"
  )
  print(model_recon)

  # --- Example 5: Composite Score Output ---
  C_mat_ortho <- rMVPA::orthogonalize_contrasts(C_mat)
  msreve_des_ortho <- rMVPA::msreve_design(
      mvpa_design = mvpa_des, # Use mvpa_des here
      contrast_matrix = C_mat_ortho
  )
  print(paste("Is contrast matrix orthonormal:", attr(msreve_des_ortho, "is_orthonormal")))

  model_composite <- contrast_rsa_model(
    dataset = mvpa_dat,
    design = msreve_des_ortho,
    output_metric = "composite",
    normalize_delta = TRUE 
  )
  print(model_composite)

  # --- Example 6: Crossnobis estimation_method ---
  # This only shows setting the method. Actual training would require passing
  # a pre-computed whitening_matrix_W to compute_crossvalidated_means_sl,
  # which is called by train_model.contrast_rsa_model.
  model_crossnobis <- contrast_rsa_model(
      dataset = mvpa_dat,
      design = msreve_des,
      estimation_method = "crossnobis"
  )
  print(model_crossnobis)
}
\seealso{
\code{\\link{msreve_design}}, \code{\\link{train_model.contrast_rsa_model}}, \code{\\link{run_searchlight}}
}
